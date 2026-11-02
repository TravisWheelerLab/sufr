//! A _.sufr_ file opened for queries

use crate::{
    format::Header,
    io::{slice_int_to_slice_u8, vec_to_slice_u8, FileAccess},
    types::{Int, ListOptions, SuffixSortType, SufrMetadata},
};
use anyhow::{anyhow, bail, Result};
use chrono::{DateTime, Local};
use home::home_dir;
use log::info;
use rayon::prelude::*;
use std::{
    cmp::min,
    fs::{self, File},
    io::{self, Read, Write},
    mem,
    ops::Range,
    path::{Path, PathBuf},
    time::Instant,
};

// --------------------------------------------------
/// Where the text is read from.
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

// --------------------------------------------------
/// A _.sufr_ file opened for queries, with its text and suffixes in
/// memory or on disk as chosen by [`Index::open`].
#[derive(Debug)]
pub struct Index<T: Int> {
    /// The _.sufr_ filename.
    pub filename: String,

    /// The file's header.
    pub header: Header,

    /// The suffix array on disk.
    pub suffix_array_file: FileAccess<T>,

    /// The LCP array on disk, if the file has one.
    pub lcp_file: Option<FileAccess<T>>,

    text: Text,

    suffixes: Suffixes<T>,

    /// The query-time maximum query length given to `open`.
    pub(crate) max_query_len: Option<usize>,
}

impl<T: Int> Index<T> {
    /// Open a _.sufr_ file whose arrays hold `T` values.
    /// `low_memory` leaves the suffix array on disk, `very_low_memory`
    /// leaves the text on disk, and `max_query_len` limits how many
    /// characters of a query are compared.
    pub fn open(
        filename: &str,
        low_memory: bool,
        very_low_memory: bool,
        max_query_len: Option<usize>,
    ) -> Result<Index<T>> {
        let header = Header::read(filename, mem::size_of::<T>())?;
        let text_file: FileAccess<u8> =
            FileAccess::new(filename, header.text_pos as u64, header.text_len)?;
        let text = if very_low_memory {
            Text::Disk(text_file)
        } else {
            Text::Memory(text_file.get_range(0..header.text_len)?)
        };
        let suffix_array_file = FileAccess::new(
            filename,
            header.suffix_array_pos as u64,
            header.len_suffixes,
        )?;
        let lcp_file = match header.lcp_pos {
            Some(lcp_pos) => Some(FileAccess::new(
                filename,
                lcp_pos as u64,
                header.len_suffixes,
            )?),
            None => None,
        };

        let mut index = Index {
            filename: filename.to_string(),
            header,
            suffix_array_file,
            lcp_file,
            text,
            suffixes: Suffixes::Disk,
            max_query_len,
        };
        if !low_memory {
            index.suffixes = index.load_suffixes(max_query_len)?;
        }
        Ok(index)
    }

    // --------------------------------------------------
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

    // --------------------------------------------------
    /// Read the subsample for `max_query_len` from its cache file in
    /// `~/.sufr`, or compute it and write the cache file.
    fn subsample_with_cache(&self, max_query_len: usize) -> Result<(Vec<T>, Vec<T>)> {
        let sufr_dir = sufr_dir()?;
        let basename = Path::new(&self.filename)
            .file_name()
            .unwrap()
            .to_string_lossy()
            .into_owned();
        let cache_path = sufr_dir.join(format!("locate-{max_query_len}-{basename}"));

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
        let (suffixes, ranks) = self.subsample_suffix_array(max_query_len);
        info!(
            "Loaded compressed SA ({}/{}) in {:?}",
            suffixes.len(),
            self.header.len_suffixes,
            now.elapsed()
        );
        if !suffixes.is_empty() {
            write_cache(&cache_path, &suffixes, &ranks)?;
        }
        Ok((suffixes, ranks))
    }

    // --------------------------------------------------
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
        file.read_exact(slice_int_to_slice_u8(&mut suffixes))?;
        let mut ranks = vec![T::from_usize(0); len];
        file.read_exact(slice_int_to_slice_u8(&mut ranks))?;

        info!(
            "Read compressed SA ({}/{}) from cache file {} in {:?}",
            suffixes.len(),
            self.header.len_suffixes,
            cache_path.display(),
            now.elapsed()
        );
        Ok((suffixes, ranks))
    }

    // --------------------------------------------------
    /// Subsample the suffix array for a maximum query length, using the
    /// LCP array when the file has one, and otherwise the text.
    pub(crate) fn subsample_suffix_array(
        &self,
        max_query_len: usize,
    ) -> (Vec<T>, Vec<T>) {
        let Some(lcp_file) = &self.lcp_file else {
            return self.subsample_suffix_array_by_text(max_query_len);
        };

        let max_query_len = T::from_usize(max_query_len);
        let mut suffixes: Vec<T> = Vec::with_capacity(self.header.len_suffixes);
        let mut ranks: Vec<T> = Vec::with_capacity(self.header.len_suffixes);
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

    // --------------------------------------------------
    /// Subsample without an LCP array: keep a suffix when its first
    /// `max_query_len` key symbols differ from its predecessor's, comparing
    /// the text directly. In max-query-len mode the key is the first
    /// `max_query_len` bytes (compared as slices, which lowers to memcmp);
    /// in seed-mask mode it is the first `max_query_len` care positions.
    fn subsample_suffix_array_by_text(&self, max_query_len: usize) -> (Vec<T>, Vec<T>) {
        let now = Instant::now();
        let text_len = self.header.text_len;
        let text_owned: Vec<u8>;
        let text: &[u8] = match &self.text {
            Text::Memory(text) => text,
            Text::Disk(text_file) => {
                // Reading the whole text cannot be out of range, so this only
                // fails on an I/O error, like the iterator over the SA file
                text_owned = text_file
                    .get_range(0..text_len)
                    .expect("failed to read text");
                &text_owned
            }
        };

        let (suffixes, ranks) = match &self.header.sort_type {
            SuffixSortType::Mask(seed_mask) => {
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
            SuffixSortType::MaxQueryLen(_) => self.keep_differing(|a, b| {
                let end_a = min(a + max_query_len, text_len);
                let end_b = min(b + max_query_len, text_len);
                // A suffix shorter than max_query_len always starts a new group
                end_a - a < max_query_len
                    || end_b - b < max_query_len
                    || text[a..end_a] != text[b..end_b]
            }),
        };
        info!(
            "Subsampled SA ({}/{}) from text (no LCP array) in {:?}",
            suffixes.len(),
            self.header.len_suffixes,
            now.elapsed()
        );

        (suffixes, ranks)
    }

    // --------------------------------------------------
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

    // --------------------------------------------------
    /// The text byte at `pos`.
    pub fn text_byte(&self, pos: usize) -> Option<u8> {
        match &self.text {
            Text::Memory(text) => text.get(pos).copied(),
            Text::Disk(text_file) => text_file.get(pos),
        }
    }

    /// The text bytes in `range`.
    pub fn text_range(&self, range: Range<usize>) -> Result<Vec<u8>> {
        match &self.text {
            Text::Memory(text) => Ok(text[range].to_vec()),
            Text::Disk(text_file) => text_file.get_range(range),
        }
    }

    /// The number of suffixes a search runs over.
    pub(crate) fn num_searched(&self) -> usize {
        match &self.suffixes {
            Suffixes::Memory(suffixes) | Suffixes::Subsampled { suffixes, .. } => {
                suffixes.len()
            }
            Suffixes::Disk => self.header.len_suffixes,
        }
    }

    /// The suffix at position `pos` of the suffixes a search runs over.
    pub(crate) fn suffix(&self, pos: usize) -> Option<T> {
        match &self.suffixes {
            Suffixes::Memory(suffixes) | Suffixes::Subsampled { suffixes, .. } => {
                suffixes.get(pos).copied()
            }
            Suffixes::Disk => self.suffix_array_file.get(pos),
        }
    }

    /// The ranks in the whole suffix array of the searched positions
    /// `first..=last`.
    pub(crate) fn rank_range(&self, first: usize, last: usize) -> Range<usize> {
        let Suffixes::Subsampled { ranks, .. } = &self.suffixes else {
            return first..last + 1;
        };

        // TODO: the last group ends at the subsample's length,
        //       not the suffix array's; should it be len_suffixes?
        let end = if first != last {
            ranks[last].to_usize() + 1
        } else if first == ranks.len() - 1 {
            ranks.len()
        } else {
            ranks[first + 1].to_usize()
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

    // --------------------------------------------------
    /// Calculate the length of the longest common prefix for two suffixes,
    /// comparing at most `len` positions. Does not consult the on-disk LCP.
    pub fn find_lcp(&self, start1: usize, start2: usize, len: usize) -> usize {
        let end1 = min(start1 + len, len);
        let end2 = min(start2 + len, len);
        (start1..end1)
            .zip(start2..end2)
            .take_while(|(a, b)| self.text_byte(*a) == self.text_byte(*b))
            .count()
    }

    // --------------------------------------------------
    /// The text from `pos`, of at most `len` bytes when given, and
    /// otherwise to the end of the text.
    pub fn string_at(&self, pos: usize, len: Option<usize>) -> Result<String> {
        let text_len = self.header.text_len;
        let end = len.map_or(text_len, |n| min(pos + n, text_len));
        let bytes = self.text_range(pos..end)?;
        Ok(String::from_utf8(bytes)?)
    }

    // --------------------------------------------------
    /// Retrieve file metadata
    pub fn metadata(&self) -> Result<SufrMetadata> {
        let fs_meta = fs::metadata(&self.filename)?;
        let modified: DateTime<Local> = DateTime::from(fs_meta.modified()?);

        Ok(SufrMetadata {
            filename: self.filename.clone(),
            modified,
            file_size: fs_meta.len() as usize,
            file_version: self.header.version as usize,
            is_dna: self.header.is_dna,
            allow_ambiguity: self.header.allow_ambiguity,
            ignore_softmask: self.header.ignore_softmask,
            text_len: self.header.text_len,
            len_suffixes: self.header.len_suffixes,
            num_sequences: self.header.sequence_starts.len(),
            sequence_starts: self.header.sequence_starts.clone(),
            sequence_names: self.header.sequence_names.clone(),
            sort_type: self.header.sort_type.clone(),
            has_lcp: self.lcp_file.is_some(),
        })
    }

    // --------------------------------------------------
    /// Print suffixes, in rank order or at the given ranks, to a file or
    /// stdout.
    pub fn list(&self, args: ListOptions) -> Result<()> {
        if args.show_lcp && self.lcp_file.is_none() {
            bail!(
                "{} has no LCP array; rebuild with --write-lcp to show LCP values",
                self.filename
            );
        }

        let mut output: Box<dyn Write> = match &args.output {
            Some(filename) => Box::new(
                File::create(filename).map_err(|e| anyhow!("{filename}: {e}"))?,
            ),
            _ => Box::new(io::stdout()),
        };

        let number = args.number.unwrap_or(0);
        if args.ranks.is_empty() {
            for (rank, suffix) in self.suffix_array_file.iter().enumerate() {
                writeln!(output, "{}", self.list_row(&args, rank, suffix)?)?;
                if number > 0 && rank == number - 1 {
                    break;
                }
            }
        } else {
            for &rank in &args.ranks {
                match self.suffix_array_file.get(rank) {
                    Some(suffix) => {
                        writeln!(output, "{}", self.list_row(&args, rank, suffix)?)?
                    }
                    None => eprintln!("Invalid rank: {rank}"),
                }
            }
        }

        Ok(())
    }

    /// One line of `list` output.
    fn list_row(&self, args: &ListOptions, rank: usize, suffix: T) -> Result<String> {
        let width = self.header.text_len.to_string().len();
        let text_len = self.header.text_len;
        let suffix = suffix.to_usize();
        let end = min(suffix + args.len.unwrap_or(text_len), text_len);

        let mut row = String::new();
        if args.show_rank {
            row += &format!("{rank:width$} ");
        }
        if args.show_suffix {
            row += &format!("{suffix:width$} ");
        }
        if args.show_lcp {
            let lcp = match &self.lcp_file {
                Some(lcp_file) => lcp_file.get(rank).unwrap(),
                None => T::from_usize(0),
            };
            row += &format!("{lcp:width$} ");
        }
        row += &String::from_utf8(self.text_range(suffix..end)?)?;
        Ok(row)
    }
}

// --------------------------------------------------
/// Find or create the `~/.sufr` directory.
fn sufr_dir() -> Result<PathBuf> {
    let home = home_dir().expect("Failed to get home directory");
    let sufr_dir = home.join(".sufr");
    if !sufr_dir.is_dir() {
        fs::create_dir(&sufr_dir)?;
    }
    Ok(sufr_dir)
}

fn write_cache<T: Int>(cache_path: &Path, suffixes: &[T], ranks: &[T]) -> Result<()> {
    let now = Instant::now();
    let mut file = File::create(cache_path)
        .map_err(|e| anyhow!("{}: {e}", cache_path.display()))?;
    file.write_all(&suffixes.len().to_le_bytes())?;
    file.write_all(vec_to_slice_u8(suffixes))?;
    file.write_all(vec_to_slice_u8(ranks))?;
    info!(
        "Wrote to cache {} in {:?}",
        cache_path.display(),
        now.elapsed()
    );
    Ok(())
}

// --------------------------------------------------
#[cfg(test)]
mod tests {
    use super::{Index, Suffixes, Text};
    use crate::{
        types::{ListOptions, SeedMask},
        SuffixArrayBuilder,
    };
    use anyhow::Result;
    use pretty_assertions::assert_eq;
    use std::fs;
    use tempfile::NamedTempFile;

    #[test]
    fn test_text_byte() -> Result<()> {
        for very_low_memory in [true, false] {
            let index = Index::<u32>::open(
                "../data/inputs/1.sufr",
                true,
                very_low_memory,
                None,
            )?;
            assert_eq!(index.text_byte(index.header.text_len - 1), Some(b'$'));
            assert_eq!(index.text_byte(index.header.text_len), None);
        }
        Ok(())
    }

    #[test]
    fn test_text_range() -> Result<()> {
        for very_low_memory in [true, false] {
            let index = Index::<u32>::open(
                "../data/inputs/1.sufr",
                true,
                very_low_memory,
                None,
            )?;
            assert_eq!(index.text_range(0..5)?, b"ACGTN");
        }
        Ok(())
    }

    #[test]
    fn test_find_lcp() -> Result<()> {
        // the text is ACGTNNACGT$
        let index = Index::<u32>::open("../data/inputs/1.sufr", true, true, None)?;
        assert_eq!(index.find_lcp(0, 6, index.header.text_len - 1), 4);
        Ok(())
    }

    #[test]
    fn test_string_at() -> Result<()> {
        let index = Index::<u32>::open("../data/inputs/1.sufr", true, true, None)?;
        assert_eq!(index.string_at(0, None)?, "ACGTNNACGT$");
        assert_eq!(index.string_at(0, Some(5))?, "ACGTN");
        assert_eq!(index.string_at(8, Some(5))?, "GT$");
        Ok(())
    }

    #[test]
    fn test_list_text_in_memory() -> Result<()> {
        let index = Index::<u32>::open("../data/inputs/1.sufr", true, false, None)?;
        let outfile = NamedTempFile::new()?;
        let list_opts = ListOptions {
            ranks: vec![2, 4, 99],
            show_rank: true,
            show_suffix: true,
            show_lcp: true,
            len: Some(3),
            number: None,
            output: Some(outfile.path().to_string_lossy().to_string()),
        };
        index.list(list_opts)?;
        assert_eq!(
            fs::read_to_string(&outfile)?,
            " 2  0  4 ACG\n 4  1  3 CGT\n"
        );
        Ok(())
    }

    // 1.sufr is fully sorted, so no max_query_len, or one of at
    // least its text length (11), reads the whole suffix array; a
    // shorter one would subsample and write a cache file in ~/.sufr
    #[test]
    fn test_open_modes() -> Result<()> {
        let file = "../data/inputs/1.sufr";
        let index = Index::<u32>::open(file, true, true, Some(2))?;
        assert!(matches!(index.text, Text::Disk(_)));
        assert!(matches!(index.suffixes, Suffixes::Disk));

        let index = Index::<u32>::open(file, true, false, Some(2))?;
        assert!(matches!(index.text, Text::Memory(_)));
        assert!(matches!(index.suffixes, Suffixes::Disk));

        for max_query_len in [None, Some(0), Some(11), Some(20)] {
            let index = Index::<u32>::open(file, false, false, max_query_len)?;
            assert!(matches!(index.suffixes, Suffixes::Memory(_)));
        }

        let index = Index::<u32>::open(file, false, true, None)?;
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
    ) -> Result<(NamedTempFile, Index<u32>)> {
        let outfile = NamedTempFile::new()?;
        let outpath = outfile.path().to_string_lossy().to_string();
        let mut builder =
            SuffixArrayBuilder::from_data(text.to_vec(), vec![(0, "s".to_string())])
                .with_dna_options(true, false, false)
                .with_num_partitions(3)
                .with_write_lcp(write_lcp);
        if let Some(n) = max_query_len {
            builder = builder.with_max_query_len(n);
        }
        if let Some(mask) = seed_mask {
            builder = builder.with_seed_mask(SeedMask::new(mask)?);
        }
        builder.build_to_file(&outpath, false)?;
        let index: Index<u32> = Index::open(&outpath, true, false, None)?;
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
            let from_text = without_lcp.subsample_suffix_array(mql);
            assert_eq!(from_text, expected, "mask, mql={mql}, vs definition");
            let from_lcp = with_lcp.subsample_suffix_array(mql);
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
            let from_text = without_lcp.subsample_suffix_array(mql);
            assert_eq!(from_text, expected, "mql index, mql={mql}, vs definition");
            let from_lcp = with_lcp.subsample_suffix_array(mql);
            assert_eq!(from_text, from_lcp, "mql index, mql={mql}, vs LCP array");
        }
        Ok(())
    }
}
