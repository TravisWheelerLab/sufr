//! A _.sufr_ file opened for queries

use crate::{
    format::{Header, SeedMask, SuffixSortType},
    io::{as_bytes, as_bytes_mut, FileAccess, Int},
};
use anyhow::{anyhow, bail, Result};
use chrono::{DateTime, Local};
use format_num::NumberFormat;
use home::home_dir;
use log::info;
use rayon::prelude::*;
use std::{
    borrow::Cow,
    cmp::min,
    fs::{self, File},
    io::{Read, Write},
    mem,
    ops::Range,
    path::{Path, PathBuf},
    time::Instant,
};

/// Where [`crate::SuffixArray::open_with`] puts the text or the suffix array.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Storage {
    /// Read into memory when the file is opened.
    Memory,

    /// Read from the file as each query needs it.
    Disk,
}

#[derive(Debug)]
pub(crate) enum Text {
    Memory(Vec<u8>),
    Disk(FileAccess<u8>),
}

/// Where a search reads suffixes from.
#[derive(Debug)]
pub(crate) enum Suffixes<T> {
    /// The whole suffix array.
    Memory(Vec<T>),

    /// The first suffix of each group sharing a `max_query_len` prefix,
    /// with its rank in the whole suffix array.
    Subsampled { suffixes: Vec<T>, ranks: Vec<T> },

    /// The suffix array on disk.
    Disk,
}

/// A summary of a _.sufr_ file.
#[derive(Debug, PartialEq)]
pub struct Metadata {
    /// The file's name.
    pub(crate) filename: String,

    /// When the file was last modified.
    pub modified: DateTime<Local>,

    /// The file's size in bytes.
    pub file_size: usize,

    /// The serialization version.
    pub file_version: usize,

    /// Whether the text holds nucleotides.
    pub is_dna: bool,

    /// Whether nucleotides other than A, C, G, and T were allowed.
    pub allow_ambiguity: bool,

    /// Whether softmasked (lowercase) nucleotides were ignored.
    pub ignore_softmask: bool,

    /// The length of the text.
    pub text_len: usize,

    /// The number of suffixes in the suffix array.
    pub num_suffixes: usize,

    /// The number of sequences in the text.
    pub num_sequences: usize,

    /// The start of each sequence in the text.
    pub sequence_starts: Vec<usize>,

    /// The name of each sequence.
    pub sequence_names: Vec<String>,

    /// How the suffixes were sorted.
    pub sort_type: SuffixSortType,

    /// Whether the file has an LCP array.
    pub has_lcp: bool,
}

impl Metadata {
    /// The labelled rows that `sufr summarize` prints.
    pub fn summary(&self) -> Vec<(&'static str, String)> {
        let num_fmt = NumberFormat::new();
        let number = |n: usize| num_fmt.format(",.0", n as f64);
        let sort_row = match &self.sort_type {
            SuffixSortType::Mask(seed_mask) => ("Seed mask", seed_mask.mask.clone()),
            SuffixSortType::MaxQueryLen(max_query_len) => {
                ("Max query len", number(*max_query_len))
            }
        };
        let sequence_starts = self
            .sequence_starts
            .iter()
            .map(|v| v.to_string())
            .collect::<Vec<_>>()
            .join(", ");
        let sequence_names = self.sequence_names.join(", ");
        vec![
            ("Filename", self.filename.clone()),
            (
                "Modified",
                self.modified.format("%Y-%m-%d %H:%M").to_string(),
            ),
            ("File Size", format!("{} bytes", number(self.file_size))),
            ("File Version", self.file_version.to_string()),
            ("DNA", self.is_dna.to_string()),
            ("Allow Ambiguity", self.allow_ambiguity.to_string()),
            ("Ignore Softmask", self.ignore_softmask.to_string()),
            ("Text Length", number(self.text_len)),
            ("Len Suffixes", number(self.num_suffixes)),
            (
                "LCP array",
                if self.has_lcp { "yes" } else { "no" }.to_string(),
            ),
            sort_row,
            ("Num sequences", number(self.num_sequences)),
            (
                "Sequence starts",
                textwrap::wrap(&sequence_starts, 40).join("\n"),
            ),
            (
                "Sequence names",
                textwrap::wrap(&sequence_names, 40).join("\n"),
            ),
        ]
    }
}

/// A _.sufr_ file opened for queries, with its text and suffixes in
/// memory or on disk as chosen by [`TypedSuffixArray::open`].
#[derive(Debug)]
pub(crate) struct TypedSuffixArray<T: Int> {
    /// The _.sufr_ filename.
    pub(crate) filename: String,

    /// The file's header.
    pub(crate) header: Header,

    /// The suffix array on disk.
    pub(crate) suffix_array_file: FileAccess<T>,

    /// The LCP array on disk, if the file has one.
    pub(crate) lcp_file: Option<FileAccess<T>>,

    text: Text,

    suffixes: Suffixes<T>,

    /// The seed mask the suffixes were sorted by, if any.
    pub(crate) seed_mask: Option<SeedMask>,

    /// The number of query characters a search compares, or 0 for all.
    pub(crate) max_query_len: usize,
}

impl<T: Int> TypedSuffixArray<T> {
    /// Open a _.sufr_ file whose arrays hold `T` values, with its
    /// `text` and `suffixes` stored as given. `max_query_len` limits how
    /// many characters of a query are compared, and when it is shorter
    /// than the sort, suffixes in memory are a subsample.
    pub(crate) fn open(
        filename: &str,
        text: Storage,
        suffixes: Storage,
        max_query_len: Option<usize>,
    ) -> Result<TypedSuffixArray<T>> {
        let header = Header::read(filename, mem::size_of::<T>())?;
        let text_file: FileAccess<u8> =
            FileAccess::new(filename, header.text_pos as u64, header.text_len)?;
        let text = match text {
            Storage::Disk => Text::Disk(text_file),
            Storage::Memory => Text::Memory(text_file.get_range(0..header.text_len)?),
        };
        Self::new(filename, header, text, suffixes, max_query_len)
    }

    /// The index over `filename`, whose `header` and `text` are already
    /// read, with its suffixes stored as given.
    pub(crate) fn new(
        filename: &str,
        header: Header,
        text: Text,
        suffixes: Storage,
        max_query_len: Option<usize>,
    ) -> Result<TypedSuffixArray<T>> {
        let suffix_array_file = FileAccess::new(
            filename,
            header.suffix_array_pos as u64,
            header.num_suffixes,
        )?;
        let lcp_file = match header.lcp_pos {
            Some(lcp_pos) => Some(FileAccess::new(
                filename,
                lcp_pos as u64,
                header.num_suffixes,
            )?),
            None => None,
        };

        let (seed_mask, search_max_query_len) = match &header.sort_type {
            SuffixSortType::Mask(seed_mask) => {
                (Some(seed_mask.clone()), max_query_len.unwrap_or(0))
            }
            SuffixSortType::MaxQueryLen(built) => {
                (None, searched_len(*built, max_query_len))
            }
        };

        let mut index = TypedSuffixArray {
            filename: filename.to_string(),
            header,
            suffix_array_file,
            lcp_file,
            text,
            suffixes: Suffixes::Disk,
            seed_mask,
            max_query_len: search_max_query_len,
        };
        if suffixes == Storage::Memory {
            index.suffixes = index.load_suffixes(max_query_len)?;
        }
        Ok(index)
    }

    /// Read the whole suffix array into memory when `max_query_len`
    /// matches how it was sorted, and otherwise a subsample of it.
    fn load_suffixes(&self, max_query_len: Option<usize>) -> Result<Suffixes<T>> {
        let built_max_query_len = match &self.header.sort_type {
            SuffixSortType::MaxQueryLen(0) => self.header.text_len,
            SuffixSortType::MaxQueryLen(mql) => *mql,
            SuffixSortType::Mask(seed_mask) => seed_mask.weight,
        };

        let mut max_query_len = max_query_len.unwrap_or(0);
        if built_max_query_len > 0 {
            max_query_len = if max_query_len > 0 {
                min(max_query_len, built_max_query_len)
            } else {
                built_max_query_len
            };
        }

        if max_query_len == built_max_query_len {
            let now = Instant::now();
            let suffixes = self.suffix_array_file.iter().collect();
            info!("Read entire SA from disk in {:?}", now.elapsed());
            return Ok(Suffixes::Memory(suffixes));
        }

        info!("Loading suffix_array_mem using max_query_len {max_query_len}");
        let (suffixes, ranks) = self.subsample_with_cache(max_query_len)?;
        Ok(Suffixes::Subsampled { suffixes, ranks })
    }

    /// Read the subsample for `max_query_len` from its cache file in
    /// `~/.sufr`, or compute it and write the cache file.
    fn subsample_with_cache(&self, max_query_len: usize) -> Result<(Vec<T>, Vec<T>)> {
        let cache_dir = cache_dir()?;
        let basename = Path::new(&self.filename)
            .file_name()
            .unwrap()
            .to_string_lossy()
            .into_owned();
        let cache_path = cache_dir.join(format!("locate-{max_query_len}-{basename}"));

        if let Ok(cache_meta) = fs::metadata(&cache_path) {
            let source_meta = fs::metadata(&self.filename)?;
            if let (Ok(source_modified), Ok(cache_modified)) =
                (source_meta.modified(), cache_meta.modified())
            {
                if source_modified > cache_modified {
                    info!("Removing stale cache {}", cache_path.display());
                    fs::remove_file(&cache_path)?;
                }
            }
        }

        if cache_path.is_file() {
            return self.read_cache(&cache_path);
        }

        let now = Instant::now();
        let (suffixes, ranks) = self.subsample(max_query_len);
        info!(
            "Loaded compressed SA ({}/{}) in {:?}",
            suffixes.len(),
            self.header.num_suffixes,
            now.elapsed()
        );
        if !suffixes.is_empty() {
            write_cache(&cache_path, &suffixes, &ranks)?;
        }
        Ok((suffixes, ranks))
    }

    fn read_cache(&self, cache_path: &Path) -> Result<(Vec<T>, Vec<T>)> {
        let now = Instant::now();
        let mut file = File::open(cache_path)
            .map_err(|e| anyhow!("{}: {e}", cache_path.display()))?;

        // the cache file holds the subsample length (8
        // bytes, little-endian), the suffixes, then the ranks
        let mut buffer = [0; 8];
        file.read_exact(&mut buffer)?;
        let len = u64::from_le_bytes(buffer) as usize;

        let mut suffixes = vec![T::from_usize(0); len];
        file.read_exact(as_bytes_mut(&mut suffixes))?;
        let mut ranks = vec![T::from_usize(0); len];
        file.read_exact(as_bytes_mut(&mut ranks))?;

        info!(
            "Read compressed SA ({}/{}) from cache file {} in {:?}",
            suffixes.len(),
            self.header.num_suffixes,
            cache_path.display(),
            now.elapsed()
        );
        Ok((suffixes, ranks))
    }

    /// The first suffix of each group sharing a `max_query_len` prefix,
    /// with its rank, found from the LCP array when the file has one and
    /// otherwise from the text.
    pub(crate) fn subsample(&self, max_query_len: usize) -> (Vec<T>, Vec<T>) {
        let Some(lcp_file) = &self.lcp_file else {
            return self.subsample_by_text(max_query_len);
        };

        let max_query_len = T::from_usize(max_query_len);
        let mut suffixes: Vec<T> = Vec::with_capacity(self.header.num_suffixes);
        let mut ranks: Vec<T> = Vec::with_capacity(self.header.num_suffixes);
        for (i, (lcp, suffix)) in lcp_file
            .iter()
            .zip(self.suffix_array_file.iter())
            .enumerate()
        {
            if lcp < max_query_len {
                suffixes.push(suffix);
                ranks.push(T::from_usize(i));
            }
        }

        (suffixes, ranks)
    }

    /// The subsample, keeping a suffix when the text of its first
    /// `max_query_len` bytes, or care positions under a seed mask,
    /// differs from its predecessor's.
    fn subsample_by_text(&self, max_query_len: usize) -> (Vec<T>, Vec<T>) {
        let now = Instant::now();
        let text_len = self.header.text_len;
        let text_owned: Vec<u8>;
        let text: &[u8] = match &self.text {
            Text::Memory(text) => text,
            Text::Disk(text_file) => {
                // reading the whole text cannot be out of range, so
                // this fails only on an I/O error, as the iterator
                // over the suffix array file does
                text_owned = text_file
                    .get_range(0..text_len)
                    .expect("failed to read text");
                &text_owned
            }
        };

        let (suffixes, ranks) = match &self.seed_mask {
            Some(seed_mask) => {
                let positions =
                    &seed_mask.positions[..min(max_query_len, seed_mask.weight)];
                self.keep_differing(|a, b| {
                    !positions.iter().all(|&offset| {
                        matches!(
                            (text.get(a + offset), text.get(b + offset)),
                            (Some(x), Some(y)) if x == y
                        )
                    })
                })
            }
            // compared as slices, which lowers to memcmp
            None => self.keep_differing(|a, b| {
                let end_a = min(a + max_query_len, text_len);
                let end_b = min(b + max_query_len, text_len);
                // a suffix shorter than max_query_len always starts a
                // new group
                end_a - a < max_query_len
                    || end_b - b < max_query_len
                    || text[a..end_a] != text[b..end_b]
            }),
        };
        info!(
            "Subsampled SA ({}/{}) from text (no LCP array) in {:?}",
            suffixes.len(),
            self.header.num_suffixes,
            now.elapsed()
        );

        (suffixes, ranks)
    }

    /// Keep each suffix, with its rank, for which `differs` is true of
    /// the previous suffix and it, and always the first suffix.
    /// The suffix array is streamed from disk in blocks, so memory is the
    /// subsample plus one block.
    fn keep_differing(
        &self,
        differs: impl Fn(usize, usize) -> bool + Sync,
    ) -> (Vec<T>, Vec<T>) {
        // compare within a block in parallel, carrying
        // the last suffix of the previous block across
        let block_len = 1 << 22;
        let mut block: Vec<T> = Vec::with_capacity(block_len);
        let mut suffixes: Vec<T> = vec![];
        let mut ranks: Vec<T> = vec![];
        let mut prev: Option<usize> = None;
        let mut seen = 0usize;
        let mut suffix_array = self.suffix_array_file.iter();
        loop {
            block.clear();
            block.extend(suffix_array.by_ref().take(block_len));
            if block.is_empty() {
                break;
            }
            let keep: Vec<bool> = block
                .par_iter()
                .enumerate()
                .map(|(i, &suffix)| {
                    let before = if i > 0 {
                        Some(block[i - 1].to_usize())
                    } else {
                        prev
                    };
                    match before {
                        None => true,
                        Some(a) => differs(a, suffix.to_usize()),
                    }
                })
                .collect();
            for (i, (&suffix, &kept)) in block.iter().zip(&keep).enumerate() {
                if kept {
                    suffixes.push(suffix);
                    ranks.push(T::from_usize(seen + i));
                }
            }
            seen += block.len();
            prev = block.last().map(|v| v.to_usize());
        }

        (suffixes, ranks)
    }

    /// The text byte at `pos`.
    pub(crate) fn text_byte(&self, pos: usize) -> Result<Option<u8>> {
        match &self.text {
            Text::Memory(text) => Ok(text.get(pos).copied()),
            Text::Disk(text_file) => text_file.get(pos),
        }
    }

    /// The text bytes in `range`.
    pub(crate) fn text_range(&self, range: Range<usize>) -> Result<Vec<u8>> {
        match &self.text {
            Text::Memory(text) => Ok(text[range].to_vec()),
            Text::Disk(text_file) => text_file.get_range(range),
        }
    }

    /// The text bytes in `range`, borrowed when the text is in memory.
    pub(crate) fn text_slice(&self, range: Range<usize>) -> Result<Cow<'_, [u8]>> {
        match &self.text {
            Text::Memory(text) => Ok(Cow::Borrowed(&text[range])),
            Text::Disk(text_file) => Ok(Cow::Owned(text_file.get_range(range)?)),
        }
    }

    /// The number of suffixes a search runs over.
    pub(crate) fn num_searched(&self) -> usize {
        match &self.suffixes {
            Suffixes::Memory(suffixes) | Suffixes::Subsampled { suffixes, .. } => {
                suffixes.len()
            }
            Suffixes::Disk => self.header.num_suffixes,
        }
    }

    /// The suffix at position `pos` of the suffixes a search runs over.
    pub(crate) fn searched_suffix(&self, pos: usize) -> Result<Option<T>> {
        match &self.suffixes {
            Suffixes::Memory(suffixes) | Suffixes::Subsampled { suffixes, .. } => {
                Ok(suffixes.get(pos).copied())
            }
            Suffixes::Disk => self.suffix_array_file.get(pos),
        }
    }

    /// The ranks in the whole suffix array of the searched positions
    /// `first..=last`.
    pub(crate) fn searched_ranks(&self, first: usize, last: usize) -> Range<usize> {
        let Suffixes::Subsampled { ranks, .. } = &self.suffixes else {
            return first..last + 1;
        };

        // the hit ends where the group after the last matched
        // one begins
        let end = match ranks.get(last + 1) {
            Some(next) => next.to_usize(),
            None => self.header.num_suffixes,
        };
        ranks[first].to_usize()..end
    }

    /// The suffixes at `ranks` in the whole suffix array.
    pub(crate) fn suffixes_at(&self, ranks: Range<usize>) -> Result<Vec<T>> {
        match &self.suffixes {
            Suffixes::Memory(suffixes) => Ok(suffixes[ranks].to_vec()),
            Suffixes::Subsampled { .. } | Suffixes::Disk => {
                self.suffix_array_file.get_range(ranks)
            }
        }
    }

    /// The LCPs at `ranks`, read from the file.
    pub(crate) fn lcps_at(&self, ranks: Range<usize>) -> Result<Vec<T>> {
        match &self.lcp_file {
            Some(lcp_file) => lcp_file.get_range(ranks),
            None => bail!("{} has no LCP array", self.filename),
        }
    }

    /// The file's metadata.
    pub(crate) fn metadata(&self) -> Result<Metadata> {
        let fs_meta = fs::metadata(&self.filename)?;
        let modified: DateTime<Local> = DateTime::from(fs_meta.modified()?);

        Ok(Metadata {
            filename: self.filename.clone(),
            modified,
            file_size: fs_meta.len() as usize,
            file_version: self.header.version as usize,
            is_dna: self.header.is_dna,
            allow_ambiguity: self.header.allow_ambiguity,
            ignore_softmask: self.header.ignore_softmask,
            text_len: self.header.text_len,
            num_suffixes: self.header.num_suffixes,
            num_sequences: self.header.sequence_starts.len(),
            sequence_starts: self.header.sequence_starts.clone(),
            sequence_names: self.header.sequence_names.clone(),
            sort_type: self.header.sort_type.clone(),
            has_lcp: self.lcp_file.is_some(),
        })
    }
}

/// The number of query characters a search compares, or 0 for all, in a
/// file sorted by its first `built` characters.
fn searched_len(built: usize, max_query_len: Option<usize>) -> usize {
    match max_query_len {
        Some(0) | None => built,
        Some(mql) if built > 0 => min(built, mql),
        Some(mql) => mql,
    }
}

/// Find or create the `~/.sufr` directory.
fn cache_dir() -> Result<PathBuf> {
    let home = home_dir().expect("Failed to get home directory");
    let cache_dir = home.join(".sufr");
    if !cache_dir.is_dir() {
        fs::create_dir(&cache_dir)?;
    }
    Ok(cache_dir)
}

fn write_cache<T: Int>(cache_path: &Path, suffixes: &[T], ranks: &[T]) -> Result<()> {
    let now = Instant::now();
    let mut file = File::create(cache_path)
        .map_err(|e| anyhow!("{}: {e}", cache_path.display()))?;
    file.write_all(&suffixes.len().to_le_bytes())?;
    file.write_all(as_bytes(suffixes))?;
    file.write_all(as_bytes(ranks))?;
    info!(
        "Wrote to cache {} in {:?}",
        cache_path.display(),
        now.elapsed()
    );
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::{Storage, Suffixes, Text, TypedSuffixArray};
    use crate::{
        format::{SeedMask, SuffixSortType},
        Alphabet, SuffixArrayBuilder,
    };
    use anyhow::Result;
    use pretty_assertions::assert_eq;
    use tempfile::NamedTempFile;

    #[test]
    fn test_text_byte() -> Result<()> {
        for very_low_memory in [true, false] {
            let index = TypedSuffixArray::<u32>::open(
                "../data/inputs/1.sufr",
                if very_low_memory {
                    Storage::Disk
                } else {
                    Storage::Memory
                },
                Storage::Disk,
                None,
            )?;
            assert_eq!(index.text_byte(index.header.text_len - 1)?, Some(b'$'));
            assert_eq!(index.text_byte(index.header.text_len)?, None);
        }
        Ok(())
    }

    #[test]
    fn test_text_range() -> Result<()> {
        for very_low_memory in [true, false] {
            let index = TypedSuffixArray::<u32>::open(
                "../data/inputs/1.sufr",
                if very_low_memory {
                    Storage::Disk
                } else {
                    Storage::Memory
                },
                Storage::Disk,
                None,
            )?;
            assert_eq!(index.text_range(0..5)?, b"ACGTN");
        }
        Ok(())
    }

    // 1.sufr is fully sorted, so no max_query_len, or one of at
    // least its text length (11), reads the whole suffix array; a
    // shorter one would subsample and write a cache file in ~/.sufr
    #[test]
    fn test_open_modes() -> Result<()> {
        let file = "../data/inputs/1.sufr";
        let index =
            TypedSuffixArray::<u32>::open(file, Storage::Disk, Storage::Disk, Some(2))?;
        assert!(matches!(index.text, Text::Disk(_)));
        assert!(matches!(index.suffixes, Suffixes::Disk));

        let index = TypedSuffixArray::<u32>::open(
            file,
            Storage::Memory,
            Storage::Disk,
            Some(2),
        )?;
        assert!(matches!(index.text, Text::Memory(_)));
        assert!(matches!(index.suffixes, Suffixes::Disk));

        for max_query_len in [None, Some(0), Some(11), Some(20)] {
            let index = TypedSuffixArray::<u32>::open(
                file,
                Storage::Memory,
                Storage::Memory,
                max_query_len,
            )?;
            assert!(matches!(index.suffixes, Suffixes::Memory(_)));
        }

        let index =
            TypedSuffixArray::<u32>::open(file, Storage::Disk, Storage::Memory, None)?;
        assert!(matches!(index.text, Text::Disk(_)));
        assert!(matches!(index.suffixes, Suffixes::Memory(_)));
        Ok(())
    }

    /// A small deterministic generator (xorshift64*) so the test needs no
    /// random-number crate
    struct Rng(u64);

    impl Rng {
        fn next(&mut self) -> u64 {
            self.0 ^= self.0 >> 12;
            self.0 ^= self.0 << 25;
            self.0 ^= self.0 >> 27;
            self.0.wrapping_mul(0x2545F4914F6CDD1D)
        }

        fn below(&mut self, n: usize) -> usize {
            (self.next() % n as u64) as usize
        }
    }

    fn random_text(rng: &mut Rng) -> Vec<u8> {
        let mut text = vec![];
        for rec in 0..3 {
            if rec > 0 {
                text.push(b'%');
            }
            let len = 300 + rng.below(400);
            text.extend((0..len).map(|_| b"ACGT"[rng.below(4)]));
            // A few Ns so some suffixes are skipped in DNA mode
            let at = rng.below(len - 5);
            let start = text.len() - len + at;
            for b in text[start..start + 4].iter_mut() {
                *b = b'N';
            }
        }
        text.push(b'$');
        text
    }

    /// Build with or without an LCP array and return the open file
    fn build(
        text: &[u8],
        seed_mask: Option<&str>,
        max_query_len: Option<usize>,
        write_lcp: bool,
    ) -> Result<(NamedTempFile, TypedSuffixArray<u32>)> {
        let outfile = NamedTempFile::new()?;
        let outpath = outfile.path().to_string_lossy().to_string();
        let mut builder =
            SuffixArrayBuilder::from_data(text.to_vec(), vec![(0, "s".to_string())])
                .with_alphabet(Alphabet::Dna {
                    allow_ambiguity: false,
                    ignore_softmask: false,
                })
                .with_num_partitions(3)
                .with_write_lcp(write_lcp);
        if let Some(n) = max_query_len {
            builder = builder.with_sort_type(SuffixSortType::MaxQueryLen(n));
        }
        if let Some(mask) = seed_mask {
            builder =
                builder.with_sort_type(SuffixSortType::Mask(SeedMask::new(mask)?));
        }
        builder.build_to_file(&outpath, true)?;
        let index: TypedSuffixArray<u32> =
            TypedSuffixArray::open(&outpath, Storage::Memory, Storage::Disk, None)?;
        Ok((outfile, index))
    }

    /// Independent definition: keep rank i when the first `mql` key symbols of
    /// SA[i] differ from those of SA[i-1] (a symbol past the text end differs
    /// from everything, including another past-the-end symbol).
    fn expected_subsample(
        text: &[u8],
        sa: &[u32],
        offsets: &[usize],
        mql: usize,
    ) -> (Vec<u32>, Vec<u32>) {
        let key = |pos: usize| -> Vec<Option<u8>> {
            offsets[..mql.min(offsets.len())]
                .iter()
                .map(|&off| text.get(pos + off).copied())
                .collect()
        };
        let mut sub_sa = vec![];
        let mut rank = vec![];
        for i in 0..sa.len() {
            let keep = i == 0 || {
                let a = key(sa[i - 1] as usize);
                let b = key(sa[i] as usize);
                !a.iter().zip(&b).all(|(x, y)| x.is_some() && x == y)
            };
            if keep {
                sub_sa.push(sa[i]);
                rank.push(i as u32);
            }
        }
        (sub_sa, rank)
    }

    #[test]
    fn subsample_without_lcp_matches_definition_and_lcp_version() -> Result<()> {
        let mut rng = Rng(11);
        let text = random_text(&mut rng);

        // Seed-mask index, query-time MQL shorter than the weight
        let mask = "11101101101111";
        let offsets: Vec<usize> = mask
            .bytes()
            .enumerate()
            .filter(|(_, b)| *b == b'1')
            .map(|(i, _)| i)
            .collect();
        for mql in [1, 2, 4, 7, 10] {
            let (_f1, with_lcp) = build(&text, Some(mask), None, true)?;
            let (_f2, without_lcp) = build(&text, Some(mask), None, false)?;
            assert!(with_lcp.lcp_file.is_some());
            assert!(without_lcp.lcp_file.is_none());
            let sa: Vec<u32> = without_lcp.suffix_array_file.iter().collect();
            let expected = expected_subsample(&text, &sa, &offsets, mql);
            let from_text = without_lcp.subsample(mql);
            assert_eq!(from_text, expected, "mask, mql={mql}, vs definition");
            let from_lcp = with_lcp.subsample(mql);
            assert_eq!(from_text, from_lcp, "mask, mql={mql}, vs LCP array");
        }

        // Max-query-len index, query-time MQL shorter than the build-time one
        let built_mql = 12;
        let offsets: Vec<usize> = (0..built_mql).collect();
        for mql in [1, 3, 6, 11] {
            let (_f1, with_lcp) = build(&text, None, Some(built_mql), true)?;
            let (_f2, without_lcp) = build(&text, None, Some(built_mql), false)?;
            let sa: Vec<u32> = without_lcp.suffix_array_file.iter().collect();
            let expected = expected_subsample(&text, &sa, &offsets, mql);
            let from_text = without_lcp.subsample(mql);
            assert_eq!(from_text, expected, "mql index, mql={mql}, vs definition");
            let from_lcp = with_lcp.subsample(mql);
            assert_eq!(from_text, from_lcp, "mql index, mql={mql}, vs LCP array");
        }
        Ok(())
    }
}
