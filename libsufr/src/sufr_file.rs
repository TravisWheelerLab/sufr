//! Read and query on-disk suffix/LCP arrays
//!
//! This is a low-level interface that you should use if you want
//! fine-grained control over whether Sufr uses 32-bit or 64-bit
//! integers. Most likely, you should use [libsufr::suffix_array](super::suffix_array).
use crate::{
    file_access::FileAccess,
    sufr_search::{SufrSearch, SufrSearchArgs},
    types::{Int, Load, OpenOptions, SeedMask, SuffixSortType, SufrMetadata},
    util::{slice_int_to_slice_u8, vec_to_slice_u8},
};
use anyhow::{anyhow, bail, Result};
use chrono::{DateTime, Local};
use home::home_dir;
use log::info;
use rayon::prelude::*;
use std::{
    cmp::min,
    fs::{self, File},
    io::{Read, Seek, Write},
    ops::Range,
    path::{Path, PathBuf},
    time::Instant,
};

// --------------------------------------------------
/// Struct used to read a serialized _.sufr_ file representing
/// the suffix and LCP arrays for a given text.
/// Provides methods to query suffix array.
#[derive(Debug)]
pub struct SufrFile<T: Int> {
    /// The _.sufr_ filename
    pub filename: String,

    /// The serialization version.
    pub version: u8,

    /// Whether or not the text represents nucleotides.
    pub is_dna: bool,

    /// When `is_dna` is `true`, whether or not nucleotides other than
    /// A, C, G, or T were ignored.
    pub allow_ambiguity: bool,

    /// When `is_dna` is `true`, whether or not softmasked/lowercase
    /// nucleotides were ignored.
    pub ignore_softmask: bool,

    /// The byte position in the file where the text begins.
    pub text_pos: usize,

    /// The byte position in the file where the suffix array begins.
    pub suffix_array_pos: usize,

    /// The byte position in the file where the LCP array begins,
    /// or 0 when the file has no LCP array.
    pub lcp_pos: usize,

    /// Whether or not the file has an on-disk LCP array
    pub has_lcp: bool,

    /// How the suffixes were sorted (fully, max query length, spaced seeds)
    pub sort_type: SuffixSortType,

    /// The length of the text. Needed for parameterization of the struct
    /// as `u32` or `u64`. Cf. [util::read_text_length](super::util::read_text_length)
    pub text_len: T,

    /// The length of the suffix array.
    pub len_suffixes: T,

    /// The number of sequences in `text`.
    pub num_sequences: T,

    /// The start positions of the sequences in `text`
    pub sequence_starts: Vec<T>,

    /// The names of the sequences
    pub sequence_names: Vec<String>,

    /// The original text that was indexed.
    pub text: Vec<u8>,

    /// File access wrapper to the `text`
    pub text_file: FileAccess<u8>,

    /// File access wrapper to the suffix array
    pub suffix_array_file: FileAccess<T>,

    /// File access wrapper to the LCP array
    pub lcp_file: FileAccess<T>,

    /// In-memory access to the suffix array
    suffix_array_mem: Vec<T>,

    /// The maximum query length of the suffix array currently in
    /// `suffix_array_mem`
    suffix_array_mem_mql: Option<usize>,

    /// In-memory access to the suffix array ranks
    suffix_array_rank_mem: Vec<T>,

    /// The maximum query length every search uses
    max_query_len: Option<usize>,
}

// --------------------------------------------------
impl<T: Int> SufrFile<T> {
    /// Open a serialized Sufr file for searching
    ///
    /// Args:
    /// * `filename`: the _.sufr_ file
    /// * `opts`: where to read the text and suffix array from,
    ///   and the maximum query length
    ///
    /// You can use [util::read_text_length](super::util::read_text_length)
    /// to get the length of the text to parameterize the `SufrFile` with
    /// the correct integer size.
    ///
    /// ```
    /// use anyhow::Result;
    /// use libsufr::{sufr_file::SufrFile, types::OpenOptions, util::read_text_length};
    ///
    /// fn main() -> Result<()> {
    ///     let filename = "../data/inputs/1.sufr";
    ///     let text_len = read_text_length(filename)? as u64;
    ///     if text_len < u32::MAX as u64 {
    ///         let sufr_file: SufrFile<u32> = SufrFile::open(&filename, &OpenOptions::default())?;
    ///     } else {
    ///         let sufr_file: SufrFile<u64> = SufrFile::open(&filename, &OpenOptions::default())?;
    ///     }
    ///     Ok(())
    /// }
    /// ```
    pub fn open(filename: &str, opts: &OpenOptions) -> Result<SufrFile<T>> {
        let mut file = File::open(filename).map_err(|e| anyhow!("{filename}: {e}"))?;

        // Meta
        let mut buffer = [0u8; 4];
        file.read_exact(&mut buffer)?;
        let version = buffer[0];
        let is_dna = buffer[1] == 1;
        let allow_ambiguity = buffer[2] == 1;
        let ignore_softmask = buffer[3] == 1;

        // Length of text
        let mut buffer = [0; 8];
        file.read_exact(&mut buffer)?;
        let text_len = u64::from_le_bytes(buffer);
        if std::mem::size_of::<usize>() == 4 && text_len > usize::MAX as u64 {
            bail!("text_len exceeds native usize");
        }
        let text_len = text_len as usize;

        // Position of text
        let mut buffer = [0; 8];
        file.read_exact(&mut buffer)?;
        let text_pos = u64::from_le_bytes(buffer) as usize;

        // Position of suffix array
        let mut buffer = [0; 8];
        file.read_exact(&mut buffer)?;
        let suffix_array_pos = u64::from_le_bytes(buffer) as usize;

        // Position of LCP array
        let mut buffer = [0; 8];
        file.read_exact(&mut buffer)?;
        let lcp_pos = u64::from_le_bytes(buffer) as usize;

        // Number of suffixes
        let mut buffer = [0; 8];
        file.read_exact(&mut buffer)?;
        let len_suffixes = u64::from_le_bytes(buffer) as usize;

        // Max query length
        let mut buffer = [0; 8];
        file.read_exact(&mut buffer)?;
        let max_query_len = T::from_usize(u64::from_le_bytes(buffer) as usize);

        // Number of sequences
        let mut buffer = [0; 8];
        file.read_exact(&mut buffer)?;
        let num_sequences = T::from_usize(u64::from_le_bytes(buffer) as usize);

        // Sequence starts
        let mut sequence_starts = vec![T::from_usize(0); num_sequences.to_usize()];
        file.read_exact(slice_int_to_slice_u8(&mut sequence_starts))?;

        // Seed mask len
        let mut buffer = [0; 8];
        file.read_exact(&mut buffer)?;
        let seed_mask_len = u64::from_le_bytes(buffer) as usize;

        // Seed mask: stored as u8 for the 1s/0s
        let seed_mask: Vec<u8> = if seed_mask_len > 0 {
            let mut buffer = vec![0; seed_mask_len];
            file.read_exact(&mut buffer)?;
            buffer
        } else {
            vec![]
        };

        // Text
        let text = match opts.text {
            Load::Disk => {
                file.seek_relative(text_len as i64)?;
                vec![]
            }
            Load::Memory => {
                let mut tmp = vec![0; text_len];
                file.read_exact(&mut tmp)?;
                tmp
            }
        };

        // Text file access
        let text_file: FileAccess<u8> =
            FileAccess::new(filename, text_pos as u64, text_len)?;

        // Suffix Array
        let suffix_array_file: FileAccess<T> =
            FileAccess::new(filename, suffix_array_pos as u64, len_suffixes)?;
        file.seek_relative(suffix_array_file.size as i64)?;

        // LCP (lcp_pos == 0 means the file was written without an LCP array)
        let has_lcp = lcp_pos > 0;
        let lcp_file: FileAccess<T> = if has_lcp {
            FileAccess::new(filename, lcp_pos as u64, len_suffixes)?
        } else {
            FileAccess::new(filename, suffix_array_pos as u64, 0)?
        };
        file.seek_relative(lcp_file.size as i64)?;

        // Sequence names are variable in length so they are at the end
        let mut buffer = vec![];
        file.read_to_end(&mut buffer)?;
        let sequence_names: Vec<String> = bincode::deserialize(&buffer)?;

        let sort_type = if seed_mask.is_empty() {
            SuffixSortType::MaxQueryLen(max_query_len.to_usize())
        } else {
            let seed_mask = SeedMask::from_bytes(&seed_mask)?;
            SuffixSortType::Mask(seed_mask)
        };

        let mut sufr_file = SufrFile {
            filename: filename.to_string(),
            version,
            is_dna,
            allow_ambiguity,
            ignore_softmask,
            text_pos,
            suffix_array_pos,
            lcp_pos,
            has_lcp,
            text_len: T::from_usize(text_len),
            len_suffixes: T::from_usize(len_suffixes),
            sort_type,
            num_sequences,
            sequence_starts,
            sequence_names,
            text,
            suffix_array_file,
            lcp_file,
            text_file,
            suffix_array_mem: vec![],
            suffix_array_mem_mql: None,
            suffix_array_rank_mem: vec![],
            max_query_len: opts.max_query_len.filter(|&len| len > 0),
        };

        if opts.suffix_array == Load::Memory {
            sufr_file.set_suffix_array_mem(sufr_file.max_query_len)?;
        }

        Ok(sufr_file)
    }

    // --------------------------------------------------
    /// Get the suffix at a position
    ///
    /// Args:
    /// * `pos`: suffix position
    ///
    /// ```
    /// use anyhow::Result;
    /// use libsufr::{sufr_file::SufrFile, types::OpenOptions};
    ///
    /// fn main() -> Result<()> {
    ///     let mut sufr = SufrFile::<u32>::open("../data/inputs/1.sufr", &OpenOptions::default())?;
    ///
    ///     // The last character will be the sentinel `$`
    ///     assert_eq!(
    ///         Some('$'),
    ///         sufr.get_text(sufr.text_len as usize - 1).map(|v| v as char)
    ///     );
    ///     Ok(())
    /// }
    ///
    /// ```
    pub fn get_text(&mut self, pos: usize) -> Option<u8> {
        if self.text.is_empty() {
            self.text_file.get(pos)
        } else {
            self.text.get(pos).copied()
        }
    }

    // --------------------------------------------------
    /// Get a slice of text starting at a given position
    ///
    /// Args:
    /// * `pos`: a `Range` of start/stop suffix positions
    ///
    /// Given the text "ACGTNNACGT":
    ///
    /// ```
    /// use anyhow::Result;
    /// use libsufr::{sufr_file::SufrFile, types::OpenOptions};
    ///
    /// fn main() -> Result<()> {
    ///     let mut sufr = SufrFile::<u32>::open("../data/inputs/1.sufr", &OpenOptions::default())?;
    ///
    ///     let bytes: Vec<u8> = vec![65, 67, 71, 84, 78];
    ///     assert_eq!(bytes, sufr.get_text_range(0..5)?);
    ///
    ///     // Or convert to a string
    ///     assert_eq!(
    ///         "ACGTN",
    ///         String::from_utf8(sufr.get_text_range(0..5)?.to_vec())?
    ///     );
    ///
    ///     Ok(())
    /// }
    /// ```
    pub fn get_text_range(&self, pos: Range<usize>) -> Result<Vec<u8>> {
        if self.text.is_empty() {
            self.text_file.get_range(pos)
        } else {
            Ok(self.text[pos].to_vec())
        }
    }

    // --------------------------------------------------
    /// Calculate the length of the longest common prefix for two suffixes.
    /// Does not consult the on-disk LCP.
    ///
    /// Args:
    /// * `start1`: start position of first suffix
    /// * `start2`: start position of second suffix
    /// * `len`: end of the text
    ///
    /// Given the text "ACGTNNACGT", the LCP of suffixes 0 and 6 is 4:
    ///
    /// ```
    /// use anyhow::Result;
    /// use libsufr::{types::{Int, OpenOptions}, sufr_file::SufrFile};
    ///
    /// fn main() -> Result<()> {
    ///     let mut sufr = SufrFile::<u32>::open("../data/inputs/1.sufr", &OpenOptions::default())?;
    ///     let lcp = sufr.find_lcp(0, 6, sufr.text_len.to_usize() - 1);
    ///     assert_eq!(4, lcp);
    ///     Ok(())
    /// }
    /// ```
    ///
    pub fn find_lcp(&mut self, start1: usize, start2: usize, len: usize) -> usize {
        let end1 = min(start1 + len, len);
        let end2 = min(start2 + len, len);
        (start1..end1)
            .zip(start2..end2)
            .take_while(|(a, b)| self.get_text(*a) == self.get_text(*b))
            .count()
    }

    // --------------------------------------------------
    /// Return a suffix at a given position
    ///
    /// Args:
    /// * `pos`: suffix position
    /// * `len`: optional prefix length; without this, the entire suffix
    ///   will be returned, which could be the length of the `text`.
    ///
    /// Given a text of "ACGTNNACGT":
    ///
    /// ```
    /// use anyhow::Result;
    /// use libsufr::{sufr_file::SufrFile, types::OpenOptions};
    ///
    /// fn main() -> Result<()> {
    ///     let mut sufr = SufrFile::<u32>::open("../data/inputs/1.sufr", &OpenOptions::default())?;
    ///
    ///     // This will print the entire input text
    ///     assert_eq!("ACGTNNACGT$", sufr.string_at(0, None)?);
    ///
    ///     // This will print the first 5 characters
    ///     assert_eq!("ACGTN", sufr.string_at(0, Some(5))?);
    ///     Ok(())
    /// }
    /// ```
    pub fn string_at(&self, pos: usize, len: Option<usize>) -> Result<String> {
        let text_len = self.text_len.to_usize();
        let end = len.map_or(text_len, |n| {
            let end = pos + n;
            if end > text_len {
                text_len
            } else {
                end
            }
        });
        let bytes = self.get_text_range(pos..end)?;
        Ok(String::from_utf8(bytes.to_vec())?)
    }

    // --------------------------------------------------
    /// Find/create a hidden "~/.sufr" directory
    fn get_sufr_dir(&self) -> Result<PathBuf> {
        let home = home_dir().expect("Failed to get home directory");
        let sufr_dir = home.join(".sufr");
        if !sufr_dir.is_dir() {
            fs::create_dir(&sufr_dir)?;
        }
        Ok(sufr_dir)
    }

    // --------------------------------------------------
    /// Subsample a suffix array using a maximum query length.
    /// Uses the LCP array when the file has one, and otherwise compares
    /// each suffix with its predecessor in the text.
    ///
    /// Args:
    /// * `max_query_len`: prefix length
    pub(crate) fn subsample_suffix_array(
        &mut self,
        max_query_len: usize,
    ) -> (Vec<T>, Vec<T>) {
        if !self.has_lcp {
            return self.subsample_suffix_array_by_text(max_query_len);
        }

        let max_query_len = T::from_usize(max_query_len);

        let max_len = self.len_suffixes.to_usize();
        let mut suffix_array: Vec<T> = Vec::with_capacity(max_len);
        let mut rank: Vec<T> = Vec::with_capacity(max_len);

        for (i, (lcp, suffix)) in self
            .lcp_file
            .iter()
            .zip(self.suffix_array_file.iter())
            .enumerate()
        {
            if lcp < max_query_len {
                suffix_array.push(suffix);
                rank.push(T::from_usize(i));
            }
        }

        (suffix_array, rank)
    }

    // --------------------------------------------------
    /// Subsample without an LCP array: keep a suffix when its first
    /// `max_query_len` key symbols differ from its predecessor's, comparing
    /// the text directly. In max-query-len mode the key is the first
    /// `max_query_len` bytes (compared as slices, which lowers to memcmp);
    /// in seed-mask mode it is the first `max_query_len` care positions.
    /// The suffix array is streamed from disk in blocks, so memory is the
    /// text plus the subsample plus one block.
    ///
    /// Args:
    /// * `max_query_len`: prefix length
    fn subsample_suffix_array_by_text(
        &mut self,
        max_query_len: usize,
    ) -> (Vec<T>, Vec<T>) {
        let now = Instant::now();
        let text_len = self.text_len.to_usize();
        let text_owned: Vec<u8>;
        let text: &[u8] = if self.text.is_empty() {
            // Reading the whole text cannot be out of range, so this only
            // fails on an I/O error, like the iterator over the SA file
            text_owned = self
                .text_file
                .get_range(0..text_len)
                .expect("failed to read text");
            &text_owned
        } else {
            &self.text
        };

        let differs: Box<dyn Fn(usize, usize) -> bool + Send + Sync> = match &self
            .sort_type
        {
            SuffixSortType::Mask(seed_mask) => {
                let positions =
                    &seed_mask.positions[..min(max_query_len, seed_mask.weight)];
                Box::new(move |a: usize, b: usize| {
                    !positions.iter().all(|&offset| {
                        matches!(
                            (text.get(a + offset), text.get(b + offset)),
                            (Some(x), Some(y)) if x == y
                        )
                    })
                })
            }
            SuffixSortType::MaxQueryLen(_) => Box::new(move |a: usize, b: usize| {
                let end_a = min(a + max_query_len, text_len);
                let end_b = min(b + max_query_len, text_len);
                // A suffix shorter than max_query_len always starts a new group
                end_a - a < max_query_len
                    || end_b - b < max_query_len
                    || text[a..end_a] != text[b..end_b]
            }),
        };

        // Stream the suffix array in blocks; compare within a block in
        // parallel, carrying the last suffix of the previous block across.
        let block_len = 1 << 22;
        let mut block: Vec<T> = Vec::with_capacity(block_len);
        let mut suffix_array: Vec<T> = vec![];
        let mut rank: Vec<T> = vec![];
        let mut prev: Option<usize> = None;
        let mut seen = 0usize;
        let mut suffixes = self.suffix_array_file.iter();
        loop {
            block.clear();
            block.extend(suffixes.by_ref().take(block_len));
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
                    suffix_array.push(suffix);
                    rank.push(T::from_usize(seen + i));
                }
            }
            seen += block.len();
            prev = block.last().map(|v| v.to_usize());
        }
        info!(
            "Subsampled SA ({}/{seen}) from text (no LCP array) in {:?}",
            suffix_array.len(),
            now.elapsed()
        );

        (suffix_array, rank)
    }

    // --------------------------------------------------
    /// Retrieve file metadata
    ///
    /// ```
    /// use anyhow::Result;
    /// use libsufr::{sufr_file::SufrFile, types::{OpenOptions, SuffixSortType}};
    ///
    /// fn main() -> Result<()> {
    ///     let sufr = SufrFile::<u32>::open("../data/inputs/1.sufr", &OpenOptions::default())?;
    ///     let meta = sufr.metadata()?;
    ///     assert_eq!(meta.filename, "../data/inputs/1.sufr".to_string());
    ///     assert_eq!(meta.file_size, 172);
    ///     assert_eq!(meta.file_version, 6);
    ///     assert_eq!(meta.is_dna, true);
    ///     assert_eq!(meta.allow_ambiguity, false);
    ///     assert_eq!(meta.ignore_softmask, false);
    ///     assert_eq!(meta.text_len, 11);
    ///     assert_eq!(meta.len_suffixes, 9);
    ///     assert_eq!(meta.num_sequences, 1);
    ///     assert_eq!(meta.sequence_starts, vec![0]);
    ///     assert_eq!(meta.sequence_names, vec!["1".to_string()]);
    ///     assert_eq!(meta.sort_type, SuffixSortType::MaxQueryLen(0));
    ///     assert_eq!(meta.has_lcp, true);
    ///
    ///     Ok(())
    /// }
    /// ```
    pub fn metadata(&self) -> Result<SufrMetadata> {
        let fs_meta = fs::metadata(&self.filename)?;
        let modified: DateTime<Local> = DateTime::from(fs_meta.modified()?);

        Ok(SufrMetadata {
            filename: self.filename.clone(),
            modified,
            file_size: fs_meta.len().to_usize(),
            file_version: self.version as usize,
            is_dna: self.is_dna,
            allow_ambiguity: self.allow_ambiguity,
            ignore_softmask: self.ignore_softmask,
            text_len: self.text_len.to_usize(),
            len_suffixes: self.len_suffixes.to_usize(),
            num_sequences: self.num_sequences.to_usize(),
            sequence_starts: self
                .sequence_starts
                .iter()
                .map(|v| v.to_usize())
                .collect::<Vec<_>>(),
            sequence_names: self.sequence_names.clone(),
            sort_type: self.sort_type.clone(),
            has_lcp: self.has_lcp,
        })
    }

    // --------------------------------------------------
    /// Read a suffix array into memory
    ///
    /// Args:
    /// * `max_query_len`: prefix length
    fn set_suffix_array_mem(&mut self, max_query_len: Option<usize>) -> Result<()> {
        let mut max_query_len = max_query_len.unwrap_or(0);

        // If ".sufr" file was built with a nonzero max_query_len or seed mask
        // Then this is the value we must use
        let built_max_query_len = match &self.sort_type {
            SuffixSortType::MaxQueryLen(mql) => {
                if mql > &0 {
                    *mql
                } else {
                    self.text_len.to_usize()
                }
            }
            SuffixSortType::Mask(seed_mask) => seed_mask.weight,
        };

        if built_max_query_len > 0 {
            if max_query_len > 0 {
                max_query_len = min(max_query_len, built_max_query_len);
            } else {
                max_query_len = built_max_query_len;
            }
        }

        // Do nothing if we've already loaded the correct SA/MQL
        if !self.suffix_array_mem.is_empty()
            && self.suffix_array_mem_mql == Some(max_query_len)
        {
            info!("Using existing suffix_array_mem");
            return Ok(());
        } else if max_query_len == built_max_query_len {
            // The requested MQL matches how the SA was built
            // Stuff entire SA into memory
            let now = Instant::now();
            self.suffix_array_mem = self.suffix_array_file.iter().collect();
            info!("Read entire SA from disk in {:?}", now.elapsed());

            // There will be no ranks
            self.suffix_array_rank_mem = vec![];
        } else {
            info!("Loading suffix_array_mem using max_query_len {max_query_len}");

            let sufr_dir = &self.get_sufr_dir()?;
            let basename = Path::new(&self.filename)
                .file_name()
                .unwrap()
                .to_string_lossy()
                .into_owned();
            let cache_path =
                sufr_dir.join(format!("locate-{max_query_len}-{basename}"));

            // Check for stale cache
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
                let now = Instant::now();
                let mut file = File::open(&cache_path)
                    .map_err(|e| anyhow!("{}: {e}", cache_path.display()))?;

                let mut buffer = [0; 8];
                file.read_exact(&mut buffer)?;
                let suffix_array_len = u64::from_le_bytes(buffer) as usize;

                self.suffix_array_mem = {
                    let mut buffer = vec![T::from_usize(0); suffix_array_len];
                    file.read_exact(slice_int_to_slice_u8(&mut buffer))?;
                    buffer
                };

                self.suffix_array_rank_mem = {
                    let mut buffer = vec![T::from_usize(0); suffix_array_len];
                    file.read_exact(slice_int_to_slice_u8(&mut buffer))?;
                    buffer
                };

                info!(
                    "Read compressed SA ({}/{}) from cache file {} in {:?}",
                    self.suffix_array_mem.len(),
                    self.len_suffixes,
                    cache_path.display(),
                    now.elapsed()
                );
            } else {
                let now = Instant::now();
                let (sub_sa, sub_rank) = &self.subsample_suffix_array(max_query_len);
                self.suffix_array_mem = sub_sa.to_vec();
                self.suffix_array_rank_mem = sub_rank.to_vec();

                info!(
                    "Loaded compressed SA ({}/{}) in {:?}",
                    sub_sa.len(),
                    self.len_suffixes,
                    now.elapsed()
                );

                // Write cache file
                if !self.suffix_array_mem.is_empty() {
                    let now = Instant::now();
                    let mut file = File::create(&cache_path)
                        .map_err(|e| anyhow!("{}: {e}", cache_path.display()))?;
                    file.write_all(&self.suffix_array_mem.len().to_le_bytes())?;
                    file.write_all(vec_to_slice_u8(&self.suffix_array_mem))?;
                    file.write_all(vec_to_slice_u8(&self.suffix_array_rank_mem))?;

                    info!(
                        "Wrote to cache {} in {:?}",
                        cache_path.display(),
                        now.elapsed()
                    );
                }
            }
        }
        self.suffix_array_mem_mql = Some(max_query_len);

        Ok(())
    }

    // --------------------------------------------------
    /// Find the ranks of the suffixes matching a query,
    /// an empty range when none match
    ///
    /// Given the text "ACGTNNACGT":
    ///
    /// ```
    /// use anyhow::Result;
    /// use libsufr::{sufr_file::SufrFile, types::OpenOptions};
    ///
    /// fn main() -> Result<()> {
    ///     let sufr = SufrFile::<u32>::open("../data/inputs/1.sufr", &OpenOptions::default())?;
    ///     assert_eq!(sufr.search(b"ACG")?, 1..3);
    ///     assert!(sufr.search(b"GGC")?.is_empty());
    ///     Ok(())
    /// }
    /// ```
    pub fn search(&self, query: &[u8]) -> Result<Range<usize>> {
        let search = SufrSearch::new(SufrSearchArgs {
            text: &self.text,
            text_len: self.text_len.to_usize(),
            text_file: &self.text_file,
            file: &self.suffix_array_file,
            suffix_array: &self.suffix_array_mem,
            rank: &self.suffix_array_rank_mem,
            len_suffixes: self.len_suffixes.to_usize(),
            sort_type: &self.sort_type,
            max_query_len: self.max_query_len,
        });
        Ok(search.search(query))
    }

    // --------------------------------------------------
    /// The positions in the text of the suffixes at a range of ranks
    pub fn suffixes(&self, ranks: Range<usize>) -> Result<Vec<usize>> {
        // a sampled suffix array is missing ranks, so only
        // a complete one can answer from memory
        let values = if !self.suffix_array_mem.is_empty()
            && self.suffix_array_rank_mem.is_empty()
        {
            self.suffix_array_mem[ranks].to_vec()
        } else {
            self.suffix_array_file.get_range(ranks)?
        };
        Ok(values.iter().map(|v| v.to_usize()).collect())
    }

    // --------------------------------------------------
    /// The LCP values at a range of ranks
    pub fn lcps(&self, ranks: Range<usize>) -> Result<Vec<usize>> {
        if !self.has_lcp {
            bail!(
                "{} has no LCP array; rebuild with --write-lcp to show LCP values",
                self.filename
            );
        }
        let values = self.lcp_file.get_range(ranks)?;
        Ok(values.iter().map(|v| v.to_usize()).collect())
    }
}

// --------------------------------------------------
#[cfg(test)]
mod test {
    use crate::{
        sufr_file::SufrFile,
        types::{Load, OpenOptions},
    };
    use anyhow::Result;

    // --------------------------------------------------
    #[test]
    fn test_file_access() -> Result<()> {
        let input_file = "../data/expected/abba.sufr";
        let sufr_file: SufrFile<u32> = SufrFile::open(
            input_file,
            &OpenOptions { text: Load::Memory, ..Default::default() },
        )?;
        let suf_by_rank = [
            14, //  0: #
            0,  //  1: AABABABABBABAB#
            12, //  2: AB#
            10, //  3: ABAB#
            1,  //  4: ABABABABBABAB#
            3,  //  5: ABABABBABAB#
            5,  //  6: ABABBABAB#
            7,  //  7: ABBABAB#
            13, //  8: B#
            11, //  9: BAB#
            9,  // 10: BABAB#
            2,  // 11: BABABABBABAB#
            4,  // 12: BABABBABAB#
            6,  // 13: BABBABAB#
            8,  // 14: BBABAB#
        ];

        for (rank, &suffix) in suf_by_rank.iter().enumerate() {
            let res = sufr_file.suffix_array_file.get(rank);
            assert!(res.is_some());
            assert_eq!(res.unwrap(), suffix);
        }

        let res = sufr_file.suffix_array_file.get_range(1..100);
        assert!(res.is_err());
        assert_eq!(
            res.as_ref().unwrap_err().to_string(),
            "Invalid range: 1..100"
        );

        let res = sufr_file.suffix_array_file.get_range(8..9);
        assert!(res.is_ok());
        assert_eq!(res.as_ref().unwrap(), &[13]);

        let res = sufr_file.suffix_array_file.get_range(8..13);
        assert!(res.is_ok());
        assert_eq!(res.as_ref().unwrap(), &[13, 11, 9, 2, 4]);

        let res = sufr_file.suffix_array_file.get_range(1..8);
        assert!(res.is_ok());
        assert_eq!(res.as_ref().unwrap(), &[0, 12, 10, 1, 3, 5, 7]);

        let all: Vec<_> = sufr_file.suffix_array_file.iter().collect();
        assert_eq!(all, &[14, 0, 12, 10, 1, 3, 5, 7, 13, 11, 9, 2, 4, 6, 8]);

        for (i, suffix) in sufr_file.suffix_array_file.iter().enumerate() {
            assert_eq!(suf_by_rank[i], suffix);
        }

        Ok(())
    }

    // --------------------------------------------------
    // The "compare" function is now deeply nested inside the SuffixSearch
    // which is created inside the "suffix_search" function and I'm lost
    // how to untangle and test this.
    //#[test]
    //fn test_compare() -> Result<()> {
    //    // 0  1  2  3  4  5  6  7  8  9 10 11 12 13 14
    //    // A  A  B  A  B  A  B  A  B  B  A  B  A  B  #
    //    let sufr_file: SufrFile<u32> = SufrFile::read("../data/inputs/abba.sufr")?;
    //
    //    // Compare to B to B with no skip
    //    let query = "B".as_bytes();
    //    let res = sufr_file.compare(query, 13, 0);
    //    assert_eq!(res.cmp, Ordering::Equal);
    //    assert_eq!(res.lcp, 1);
    //
    //    // Compare to B to B with skip = 1
    //    let query = "B".as_bytes();
    //    let res = sufr_file.compare(query, 13, 1);
    //    assert_eq!(res.cmp, Ordering::Equal);
    //    assert_eq!(res.lcp, 1);
    //
    //    // Compare to B to AB
    //    let query = "B".as_bytes();
    //    let res = sufr_file.compare(query, 12, 0);
    //    assert_eq!(res.cmp, Ordering::Greater);
    //    assert_eq!(res.lcp, 0);
    //
    //    // Compare to ABABA to ABBABAB#
    //    let query = "ABABA".as_bytes();
    //    let res = sufr_file.compare(query, 7, 2);
    //    assert_eq!(res.cmp, Ordering::Less);
    //    assert_eq!(res.lcp, 2);
    //
    //    // Compare to ABAB to ABABBABAB#
    //    let query = "ABABA".as_bytes();
    //    let res = sufr_file.compare(query, 5, 2);
    //    assert_eq!(res.cmp, Ordering::Less);
    //    assert_eq!(res.lcp, 4);
    //
    //    Ok(())
    //}
}
