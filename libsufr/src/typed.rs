//! A _.sufr_ file opened for queries

use crate::{
    format::{read_text_length, Header, SeedMask, SuffixSortType},
    io::{as_bytes, as_bytes_mut, Int},
    search::{Context, Extraction, Location},
    source::{DiskSource, MemorySource, Source},
    OpenOptions,
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
    fmt::Debug,
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

/// A _.sufr_ file opened for queries, with its text in `T` and its
/// suffix array in `S`.
#[derive(Debug)]
pub struct TypedSuffixArray<W: Int, T: Source<u8>, S: Source<W>> {
    /// The _.sufr_ filename.
    pub(crate) filename: String,

    /// The file's header.
    pub(crate) header: Header,

    text: T,

    suffixes: S,

    /// The LCP array, if the file has one.
    //
    // always on disk: no search reads it, and no mode loads it
    pub(crate) lcps: Option<DiskSource<W>>,

    /// The suffixes a search runs over instead of the whole array, when
    /// `max_query_len` is shorter than the sort.
    subsample: Option<Subsample<W>>,

    /// The seed mask the suffixes were sorted by, if any.
    pub(crate) seed_mask: Option<SeedMask>,

    /// The number of query characters a search compares, or 0 for all.
    pub(crate) max_query_len: usize,
}

/// The first suffix of each group sharing a `max_query_len` prefix,
/// with its rank in the whole suffix array.
#[derive(Debug, PartialEq)]
pub(crate) struct Subsample<W> {
    suffixes: Vec<W>,
    ranks: Vec<W>,
}

impl<W: Int> Subsample<W> {
    #[cfg(test)]
    pub(crate) fn pair(self) -> (Vec<W>, Vec<W>) {
        (self.suffixes, self.ranks)
    }
}

impl<W: Int, T: Source<u8>, S: Source<W>> TypedSuffixArray<W, T, S> {
    /// Open a _.sufr_ file whose arrays hold `W` values, with its text in
    /// `T` and its suffix array in `S`. `max_query_len` limits how many
    /// characters of a query are compared. With `subsample`, the first
    /// suffix of each group sharing that many characters is held in
    /// memory and searched instead of the whole array.
    pub fn open(
        filename: &str,
        max_query_len: Option<usize>,
        subsample: bool,
    ) -> Result<Self> {
        let header = Header::read(filename, mem::size_of::<W>())?;
        let path = Path::new(filename);
        let text = T::open(path, header.text_pos as u64, header.text_len)?;
        let suffixes =
            S::open(path, header.suffix_array_pos as u64, header.num_suffixes)?;
        let lcps = match header.lcp_pos {
            Some(lcp_pos) => {
                Some(DiskSource::open(path, lcp_pos as u64, header.num_suffixes)?)
            }
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
            text,
            suffixes,
            lcps,
            subsample: None,
            seed_mask,
            max_query_len: search_max_query_len,
        };
        if subsample {
            let len = subsample_len(&index.header, max_query_len);
            info!("Loading suffix_array_mem using max_query_len {len}");
            index.subsample = Some(index.subsample_with_cache(len)?);
        }
        Ok(index)
    }

    /// The text byte at `pos`.
    pub(crate) fn text_byte(&self, pos: usize) -> Result<Option<u8>> {
        self.text.get(pos)
    }

    /// The text bytes in `range`, borrowed when the text is in memory.
    pub(crate) fn text_slice(&self, range: Range<usize>) -> Result<Cow<'_, [u8]>> {
        self.text.range(range)
    }

    /// The number of suffixes a search runs over.
    pub(crate) fn num_searched(&self) -> usize {
        match &self.subsample {
            Some(subsample) => subsample.suffixes.len(),
            None => self.suffixes.len(),
        }
    }

    /// The suffix at position `pos` of the suffixes a search runs over.
    pub(crate) fn searched_suffix(&self, pos: usize) -> Result<Option<W>> {
        match &self.subsample {
            Some(subsample) => Ok(subsample.suffixes.get(pos).copied()),
            None => self.suffixes.get(pos),
        }
    }

    /// The ranks in the whole suffix array of the searched positions
    /// `first..=last`.
    pub(crate) fn searched_ranks(&self, first: usize, last: usize) -> Range<usize> {
        let Some(subsample) = &self.subsample else {
            return first..last + 1;
        };

        // the hit ends where the group after the last matched
        // one begins
        let end = match subsample.ranks.get(last + 1) {
            Some(next) => next.to_usize(),
            None => self.header.num_suffixes,
        };
        subsample.ranks[first].to_usize()..end
    }

    /// The positions in the text of the suffixes at `ranks`.
    pub(crate) fn suffixes(&self, ranks: Range<usize>) -> Result<Vec<usize>> {
        Ok(self
            .suffixes
            .range(ranks)?
            .iter()
            .map(|&v| v.to_usize())
            .collect())
    }

    /// The LCP of each suffix at `ranks` with the suffix before it. An
    /// error when the file has no LCP array.
    pub(crate) fn lcps(&self, ranks: Range<usize>) -> Result<Vec<usize>> {
        match &self.lcps {
            Some(lcps) => {
                Ok(lcps.range(ranks)?.iter().map(|&v| v.to_usize()).collect())
            }
            None => bail!("{} has no LCP array", self.filename),
        }
    }

    /// The text bytes in `range`.
    pub(crate) fn text(&self, range: Range<usize>) -> Result<Vec<u8>> {
        Ok(self.text.range(range)?.into_owned())
    }

    /// Read the subsample for `max_query_len` from its cache file in
    /// `~/.sufr`, or compute it and write the cache file.
    fn subsample_with_cache(&self, max_query_len: usize) -> Result<Subsample<W>> {
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
        let subsample = self.subsample(max_query_len)?;
        info!(
            "Loaded compressed SA ({}/{}) in {:?}",
            subsample.suffixes.len(),
            self.header.num_suffixes,
            now.elapsed()
        );
        if !subsample.suffixes.is_empty() {
            write_cache(&cache_path, &subsample)?;
        }
        Ok(subsample)
    }

    fn read_cache(&self, cache_path: &Path) -> Result<Subsample<W>> {
        let now = Instant::now();
        let mut file = File::open(cache_path)
            .map_err(|e| anyhow!("{}: {e}", cache_path.display()))?;

        // the cache file holds the subsample length (8
        // bytes, little-endian), the suffixes, then the ranks
        let mut buffer = [0; 8];
        file.read_exact(&mut buffer)?;
        let len = u64::from_le_bytes(buffer) as usize;

        let mut suffixes = vec![W::from_usize(0); len];
        file.read_exact(as_bytes_mut(&mut suffixes))?;
        let mut ranks = vec![W::from_usize(0); len];
        file.read_exact(as_bytes_mut(&mut ranks))?;

        info!(
            "Read compressed SA ({}/{}) from cache file {} in {:?}",
            suffixes.len(),
            self.header.num_suffixes,
            cache_path.display(),
            now.elapsed()
        );
        Ok(Subsample { suffixes, ranks })
    }

    /// The first suffix of each group sharing a `max_query_len` prefix,
    /// with its rank, found from the LCP array when the file has one and
    /// otherwise from the text.
    pub(crate) fn subsample(&self, max_query_len: usize) -> Result<Subsample<W>> {
        let Some(lcps) = &self.lcps else {
            return self.subsample_by_text(max_query_len);
        };

        let max_query_len = W::from_usize(max_query_len);
        let num_suffixes = self.header.num_suffixes;
        let mut suffixes: Vec<W> = Vec::with_capacity(num_suffixes);
        let mut ranks: Vec<W> = Vec::with_capacity(num_suffixes);
        for start in (0..num_suffixes).step_by(BLOCK_LEN) {
            let block = start..min(start + BLOCK_LEN, num_suffixes);
            let block_lcps = lcps.range(block.clone())?;
            let block_suffixes = self.suffixes.range(block.clone())?;
            for (i, (&lcp, &suffix)) in
                block_lcps.iter().zip(block_suffixes.iter()).enumerate()
            {
                if lcp < max_query_len {
                    suffixes.push(suffix);
                    ranks.push(W::from_usize(start + i));
                }
            }
        }

        Ok(Subsample { suffixes, ranks })
    }

    /// The subsample, keeping a suffix when the text of its first
    /// `max_query_len` bytes, or care positions under a seed mask,
    /// differs from its predecessor's.
    fn subsample_by_text(&self, max_query_len: usize) -> Result<Subsample<W>> {
        let now = Instant::now();
        let text_len = self.header.text_len;
        let text = self.text.range(0..text_len)?;
        let text: &[u8] = &text;

        let subsample = match &self.seed_mask {
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
                })?
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
            })?,
        };
        info!(
            "Subsampled SA ({}/{}) from text (no LCP array) in {:?}",
            subsample.suffixes.len(),
            self.header.num_suffixes,
            now.elapsed()
        );

        Ok(subsample)
    }

    /// Keep each suffix, with its rank, for which `differs` is true of
    /// the previous suffix and it, and always the first suffix.
    /// The suffix array is read from disk in blocks, so memory is the
    /// subsample plus one block.
    fn keep_differing(
        &self,
        differs: impl Fn(usize, usize) -> bool + Sync,
    ) -> Result<Subsample<W>> {
        // compare within a block in parallel, carrying
        // the last suffix of the previous block across
        let num_suffixes = self.header.num_suffixes;
        let mut suffixes: Vec<W> = vec![];
        let mut ranks: Vec<W> = vec![];
        let mut prev: Option<usize> = None;
        for start in (0..num_suffixes).step_by(BLOCK_LEN) {
            let block = self
                .suffixes
                .range(start..min(start + BLOCK_LEN, num_suffixes))?;
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
                    ranks.push(W::from_usize(start + i));
                }
            }
            prev = block.last().map(|v| v.to_usize());
        }

        Ok(Subsample { suffixes, ranks })
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
            has_lcp: self.lcps.is_some(),
        })
    }
}

/// What [`crate::SuffixArray`] asks of a [`TypedSuffixArray`] at whichever int
/// size and storage its file was opened with.
pub(crate) trait Queryable: Debug + Send + Sync {
    fn header(&self) -> &Header;

    fn metadata(&self) -> Result<Metadata>;

    fn matching_ranks(&self, query: &[u8]) -> Result<Range<usize>>;

    fn suffixes(&self, ranks: Range<usize>) -> Result<Vec<usize>>;

    fn lcps(&self, ranks: Range<usize>) -> Result<Vec<usize>>;

    fn text(&self, range: Range<usize>) -> Result<Vec<u8>>;

    fn locate(&self, ranks: Range<usize>) -> Result<Vec<Location<'_>>>;

    fn extract(
        &self,
        ranks: Range<usize>,
        context: Context,
    ) -> Result<Vec<Extraction<'_>>>;
}

impl<W: Int, T: Source<u8>, S: Source<W>> Queryable for TypedSuffixArray<W, T, S> {
    fn header(&self) -> &Header {
        &self.header
    }

    fn metadata(&self) -> Result<Metadata> {
        TypedSuffixArray::metadata(self)
    }

    fn matching_ranks(&self, query: &[u8]) -> Result<Range<usize>> {
        TypedSuffixArray::matching_ranks(self, query)
    }

    fn suffixes(&self, ranks: Range<usize>) -> Result<Vec<usize>> {
        TypedSuffixArray::suffixes(self, ranks)
    }

    fn lcps(&self, ranks: Range<usize>) -> Result<Vec<usize>> {
        TypedSuffixArray::lcps(self, ranks)
    }

    fn text(&self, range: Range<usize>) -> Result<Vec<u8>> {
        TypedSuffixArray::text(self, range)
    }

    fn locate(&self, ranks: Range<usize>) -> Result<Vec<Location<'_>>> {
        TypedSuffixArray::locate(self, ranks)
    }

    fn extract(
        &self,
        ranks: Range<usize>,
        context: Context,
    ) -> Result<Vec<Extraction<'_>>> {
        TypedSuffixArray::extract(self, ranks, context)
    }
}

/// Open `filename` as [`crate::SuffixArray::open_with`] asks, at the int
/// width its text length needs.
pub(crate) fn open(
    filename: &str,
    options: &OpenOptions,
) -> Result<Box<dyn Queryable>> {
    use Storage::{Disk as D, Memory as M};
    let cap = options.max_query_len;
    let text_len = read_text_length(filename)?;
    let wide = text_len as u64 >= u32::MAX as u64;

    // suffixes in memory under a cap shorter than the sort means
    // the subsample in memory and the whole array on disk
    let header = Header::read(filename, if wide { 8 } else { 4 })?;
    let subsample =
        options.suffixes == M && subsample_len(&header, cap) < header.sorted_len();
    let suffixes = if subsample { D } else { options.suffixes };

    let inner: Box<dyn Queryable> = match (wide, options.text, suffixes) {
        (false, D, D) => Box::new(TypedSuffixArray::<
            u32,
            DiskSource<u8>,
            DiskSource<u32>,
        >::open(filename, cap, subsample)?),
        (false, M, D) => Box::new(TypedSuffixArray::<
            u32,
            MemorySource<u8>,
            DiskSource<u32>,
        >::open(filename, cap, subsample)?),
        (false, D, M) => Box::new(TypedSuffixArray::<
            u32,
            DiskSource<u8>,
            MemorySource<u32>,
        >::open(filename, cap, subsample)?),
        (false, M, M) => Box::new(TypedSuffixArray::<
            u32,
            MemorySource<u8>,
            MemorySource<u32>,
        >::open(filename, cap, subsample)?),
        (true, D, D) => Box::new(TypedSuffixArray::<
            u64,
            DiskSource<u8>,
            DiskSource<u64>,
        >::open(filename, cap, subsample)?),
        (true, M, D) => Box::new(TypedSuffixArray::<
            u64,
            MemorySource<u8>,
            DiskSource<u64>,
        >::open(filename, cap, subsample)?),
        (true, D, M) => Box::new(TypedSuffixArray::<
            u64,
            DiskSource<u8>,
            MemorySource<u64>,
        >::open(filename, cap, subsample)?),
        (true, M, M) => Box::new(TypedSuffixArray::<
            u64,
            MemorySource<u8>,
            MemorySource<u64>,
        >::open(filename, cap, subsample)?),
    };
    Ok(inner)
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

fn write_cache<W: Int>(cache_path: &Path, subsample: &Subsample<W>) -> Result<()> {
    let now = Instant::now();
    let mut file = File::create(cache_path)
        .map_err(|e| anyhow!("{}: {e}", cache_path.display()))?;
    file.write_all(&subsample.suffixes.len().to_le_bytes())?;
    file.write_all(as_bytes(&subsample.suffixes))?;
    file.write_all(as_bytes(&subsample.ranks))?;
    info!(
        "Wrote to cache {} in {:?}",
        cache_path.display(),
        now.elapsed()
    );
    Ok(())
}

/// The number of suffixes read from the file at once while subsampling.
const BLOCK_LEN: usize = 1 << 22;

/// The prefix length a subsample groups by: `max_query_len`, clipped to
/// how the file was sorted.
pub(crate) fn subsample_len(header: &Header, max_query_len: Option<usize>) -> usize {
    let built = header.sorted_len();
    match max_query_len {
        Some(mql) if mql > 0 => min(mql, built),
        _ => built,
    }
}

#[cfg(test)]
mod tests {
    use super::TypedSuffixArray;
    use crate::{
        format::SeedMask,
        source::{DiskSource, MemorySource, Source},
        Alphabet, SuffixArrayBuilder,
    };
    use anyhow::Result;
    use pretty_assertions::assert_eq;
    use tempfile::NamedTempFile;

    type Index = TypedSuffixArray<u32, MemorySource<u8>, DiskSource<u32>>;

    #[test]
    fn test_text_byte() -> Result<()> {
        let file = "../data/inputs/1.sufr";
        let on_disk = TypedSuffixArray::<u32, DiskSource<u8>, DiskSource<u32>>::open(
            file, None, false,
        )?;
        let in_memory = Index::open(file, None, false)?;
        let text_len = on_disk.header.text_len;
        assert_eq!(on_disk.text_byte(text_len - 1)?, Some(b'$'));
        assert_eq!(on_disk.text_byte(text_len)?, None);
        assert_eq!(in_memory.text_byte(text_len - 1)?, Some(b'$'));
        assert_eq!(in_memory.text_byte(text_len)?, None);
        Ok(())
    }

    #[test]
    fn test_text_range() -> Result<()> {
        let file = "../data/inputs/1.sufr";
        let on_disk = TypedSuffixArray::<u32, DiskSource<u8>, DiskSource<u32>>::open(
            file, None, false,
        )?;
        let in_memory = Index::open(file, None, false)?;
        assert_eq!(on_disk.text(0..5)?, b"ACGTN");
        assert_eq!(in_memory.text(0..5)?, b"ACGTN");
        Ok(())
    }

    // 1.sufr is fully sorted, so a max_query_len of at least its text
    // length (11) subsamples to the whole array; a shorter one would
    // write a cache file in ~/.sufr
    #[test]
    fn test_subsample_at_the_sort_length_is_the_whole_array() -> Result<()> {
        let file = "../data/inputs/1.sufr";
        for max_query_len in [None, Some(0), Some(11), Some(20)] {
            let index =
                TypedSuffixArray::<u32, MemorySource<u8>, MemorySource<u32>>::open(
                    file,
                    max_query_len,
                    true,
                )?;
            let subsample = index.subsample.as_ref().unwrap();
            assert_eq!(subsample.suffixes, index.suffixes.0);
            assert_eq!(subsample.ranks, (0..9).collect::<Vec<u32>>());
            assert_eq!(index.num_searched(), 9);
        }
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
        write_lcp: bool,
    ) -> Result<(NamedTempFile, Index)> {
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
        if let Some(mask) = seed_mask {
            builder = builder.with_seed_mask(SeedMask::new(mask)?);
        }
        builder.build_to_file(&outpath)?;
        let index = Index::open(&outpath, None, false)?;
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
            let (_f1, with_lcp) = build(&text, Some(mask), true)?;
            let (_f2, without_lcp) = build(&text, Some(mask), false)?;
            assert!(with_lcp.lcps.is_some());
            assert!(without_lcp.lcps.is_none());
            let sa = without_lcp
                .suffixes
                .range(0..without_lcp.header.num_suffixes)?;
            let expected = expected_subsample(&text, &sa, &offsets, mql);
            let from_text = without_lcp.subsample(mql)?.pair();
            assert_eq!(from_text, expected, "mask, mql={mql}, vs definition");
            let from_lcp = with_lcp.subsample(mql)?.pair();
            assert_eq!(from_text, from_lcp, "mask, mql={mql}, vs LCP array");
        }
        Ok(())
    }
}
