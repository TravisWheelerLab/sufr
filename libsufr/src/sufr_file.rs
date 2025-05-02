//! Read and query on-disk suffix/LCP arrays
//!
//! This is a low-level interface that you should use if you want
//! fine-grained control over whether Sufr uses 32-bit or 64-bit
//! integers. Most likely, you should use [libsufr::suffix_array](super::suffix_array).
use crate::{
    file_access::FileAccess,
    sufr_search::{SufrSearch, SufrSearchArgs},
    types::{
        BisectOptions, BisectResult, CountOptions, CountResult, ExtractOptions, ExtractResult, ExtractSequence,
        FromUsize, Int, ListOptions, LocateOptions, LocatePosition, LocateResult,
        SearchOptions, SearchResult, SeedMask, SuffixSortType, SufrMetadata,
    },
    util::{slice_u8_to_vec, usize_to_bytes},
};
use anyhow::{anyhow, Result, bail};
use chrono::{DateTime, Local};
use home::home_dir;
use log::info;
use rayon::prelude::*;
use std::{
    cell::RefCell,
    cmp::min,
    fs::{self, File},
    io::{self, Read, Seek, Write},
    mem,
    ops::Range,
    path::{Path, PathBuf},
    slice,
    time::Instant,
};
use thread_local::ThreadLocal;

// --------------------------------------------------
/// Struct used to read a serialized _.sufr_ file representing
/// the suffix and LCP arrays for a given text.
/// Provides methods to query suffix array.
#[derive(Debug)]
pub struct SufrFile<T>
where
    T: Int + FromUsize<T> + Sized + Send + Sync + serde::ser::Serialize,
{
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

    /// Whether or not to query the suffix array in-memory or on disk.
    pub query_low_memory: bool,

    /// The byte position in the file where the text begins.
    pub text_pos: usize,

    /// The byte position in the file where the suffix array begins.
    pub suffix_array_pos: usize,

    /// The byte position in the file where the LCP array begins.
    pub lcp_pos: usize,

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
}

// --------------------------------------------------
impl<T> SufrFile<T>
where
    T: Int + FromUsize<T> + Sized + Send + Sync,
{
    /// Read serialized Sufr file
    ///
    /// Args:
    /// * `filename`: the _.sufr_ file
    /// * `low_memory`: When `false`, the `text` will be loaded into memory.
    ///   When `true`, the `text` will be read from disk as needed.
    ///
    /// You can use [util::read_text_length](super::util::read_text_length)
    /// to get the length of the text to parameterize the `SufrFile` with
    /// the correct integer size.
    ///
    /// ```
    /// use anyhow::Result;
    /// use libsufr::{sufr_file::SufrFile, util::read_text_length};
    ///
    /// fn main() -> Result<()> {
    ///     let filename = "../data/inputs/1.sufr";
    ///     let text_len = read_text_length(filename)? as u64;
    ///     if text_len < u32::MAX as u64 {
    ///         let sufr_file: SufrFile<u32> = SufrFile::read(&filename, true)?;
    ///     } else {
    ///         let sufr_file: SufrFile<u64> = SufrFile::read(&filename, true)?;
    ///     }
    ///     Ok(())
    /// }
    /// ```
    pub fn read(filename: &str, low_memory: bool) -> Result<SufrFile<T>> {
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
        let text_len = u64::from_ne_bytes(buffer);
        if std::mem::size_of::<usize>() == 4 && text_len > usize::MAX as u64 {
            bail!("text_len exceeds native usize");
        }
        let text_len = text_len as usize;

        // Position of text
        let mut buffer = [0; 8];
        file.read_exact(&mut buffer)?;
        let text_pos = u64::from_ne_bytes(buffer) as usize;

        // Position of suffix array
        let mut buffer = [0; 8];
        file.read_exact(&mut buffer)?;
        let suffix_array_pos = u64::from_ne_bytes(buffer) as usize;

        // Position of LCP array
        let mut buffer = [0; 8];
        file.read_exact(&mut buffer)?;
        let lcp_pos = u64::from_ne_bytes(buffer) as usize;

        // Number of suffixes
        let mut buffer = [0; 8];
        file.read_exact(&mut buffer)?;
        let len_suffixes = u64::from_ne_bytes(buffer) as usize;

        // Max query length
        let mut buffer = [0; 8];
        file.read_exact(&mut buffer)?;
        let max_query_len = T::from_usize(u64::from_ne_bytes(buffer) as usize);

        // Number of sequences
        let mut buffer = [0; 8];
        file.read_exact(&mut buffer)?;
        let num_sequences = T::from_usize(u64::from_ne_bytes(buffer) as usize);

        // Sequence starts
        let mut buffer = vec![0; num_sequences.to_usize() * mem::size_of::<T>()];
        file.read_exact(&mut buffer)?;
        let sequence_starts: Vec<T> =
            slice_u8_to_vec(&buffer, num_sequences.to_usize());

        // Seed mask len
        let mut buffer = [0; 8];
        file.read_exact(&mut buffer)?;
        let seed_mask_len = u64::from_ne_bytes(buffer) as usize;

        // Seed mask: stored as u8 for the 1s/0s
        let seed_mask: Vec<u8> = if seed_mask_len > 0 {
            let mut buffer = vec![0; seed_mask_len];
            file.read_exact(&mut buffer)?;
            buffer
        } else {
            vec![]
        };

        // Text
        let text = if low_memory {
            file.seek_relative(text_len as i64)?;
            vec![]
        } else {
            let mut tmp = vec![0; text_len];
            file.read_exact(&mut tmp)?;
            tmp
        };

        // Text file access
        let text_file: FileAccess<u8> =
            FileAccess::new(filename, text_pos as u64, text_len)?;

        // Suffix Array
        let suffix_array_file: FileAccess<T> =
            FileAccess::new(filename, suffix_array_pos as u64, len_suffixes)?;
        file.seek_relative(suffix_array_file.size as i64)?;

        // LCP
        let lcp_file: FileAccess<T> =
            FileAccess::new(filename, lcp_pos as u64, len_suffixes)?;
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

        Ok(SufrFile {
            filename: filename.to_string(),
            version,
            is_dna,
            allow_ambiguity,
            ignore_softmask,
            query_low_memory: true,
            text_pos,
            suffix_array_pos,
            lcp_pos,
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
        })
    }

    // --------------------------------------------------
    /// Get the suffix at a position
    ///
    /// Args:
    /// * `pos`: suffix position
    ///
    /// ```
    /// use anyhow::Result;
    /// use libsufr::sufr_file::SufrFile;
    ///
    /// fn main() -> Result<()> {
    ///     let mut sufr = SufrFile::<u32>::read("../data/inputs/1.sufr", true)?;
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
    /// use libsufr::sufr_file::SufrFile;
    ///
    /// fn main() -> Result<()> {
    ///     let mut sufr = SufrFile::<u32>::read("../data/inputs/1.sufr", true)?;
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
    pub fn get_text_range(&mut self, pos: Range<usize>) -> Result<Vec<u8>> {
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
    /// use libsufr::{types::Int, sufr_file::SufrFile};
    ///
    /// fn main() -> Result<()> {
    ///     let mut sufr = SufrFile::<u32>::read("../data/inputs/1.sufr", true)?;
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
    /// use libsufr::sufr_file::SufrFile;
    ///
    /// fn main() -> Result<()> {
    ///     let mut sufr = SufrFile::<u32>::read("../data/inputs/1.sufr", true)?;
    ///
    ///     // This will print the entire input text
    ///     assert_eq!("ACGTNNACGT$", sufr.string_at(0, None)?);
    ///
    ///     // This will print the first 5 characters
    ///     assert_eq!("ACGTN", sufr.string_at(0, Some(5))?);
    ///     Ok(())
    /// }
    /// ```
    pub fn string_at(&mut self, pos: usize, len: Option<usize>) -> Result<String> {
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
    /// Subsample a suffix array using a maximum query length
    ///
    /// Args:
    /// * `max_query_len`: prefix length
    pub(crate) fn subsample_suffix_array(
        &mut self,
        max_query_len: usize,
    ) -> (Vec<T>, Vec<T>) {
        let max_query_len = T::from_usize(max_query_len);

        // Ensure we start from the beginning of the SA/LCP files
        self.lcp_file.reset();
        self.suffix_array_file.reset();

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
    /// Retrieve file metadata
    ///
    /// ```
    /// use anyhow::Result;
    /// use libsufr::{sufr_file::SufrFile, types::SuffixSortType};
    ///
    /// fn main() -> Result<()> {
    ///     let sufr = SufrFile::<u32>::read("../data/inputs/1.sufr", true)?;
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
            self.suffix_array_file.reset();
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
                let suffix_array_len = u64::from_ne_bytes(buffer) as usize;

                self.suffix_array_mem = if suffix_array_len == 0 {
                    vec![]
                } else {
                    let mut buffer = vec![0; suffix_array_len * mem::size_of::<T>()];
                    file.read_exact(&mut buffer)?;
                    slice_u8_to_vec(&buffer, suffix_array_len)
                };

                let mut buffer = vec![];
                file.read_to_end(&mut buffer)?;
                self.suffix_array_rank_mem = if buffer.is_empty() {
                    vec![]
                } else {
                    unsafe {
                        std::slice::from_raw_parts(
                            buffer.as_ptr() as *const _,
                            suffix_array_len,
                        )
                        .to_vec()
                    }
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
                    let _ = file.write(&usize_to_bytes(self.suffix_array_mem.len()))?;
                    let bytes = unsafe {
                        slice::from_raw_parts(
                            self.suffix_array_mem.as_ptr() as *const u8,
                            self.suffix_array_mem.len() * std::mem::size_of::<T>(),
                        )
                    };
                    file.write_all(bytes)?;

                    if !self.suffix_array_rank_mem.is_empty() {
                        let bytes = unsafe {
                            slice::from_raw_parts(
                                self.suffix_array_rank_mem.as_ptr() as *const u8,
                                self.suffix_array_rank_mem.len()
                                    * std::mem::size_of::<usize>(),
                            )
                        };
                        file.write_all(bytes)?;
                    }

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
    /// Bisect the index range of occurences of queries.
    /// If the index range of a prefix is already known,
    /// or if it is desirable to avoid enumerating every match,
    /// this method can be used as a faster stand-in for `count`
    /// ```
    /// use anyhow::Result;
    /// use libsufr::{types::{BisectOptions, BisectResult}, sufr_file::SufrFile};
    /// 
    /// fn main() -> Result<()> {
    ///     let mut sufr = SufrFile::<u32>::read("../data/inputs/1.sufr", false)?;
    ///     // bisect without a prefix result, searching the whole suffix array:
    ///     let opts_without_prefix = BisectOptions {
    ///         queries: vec!["AC".to_string(), "CG".to_string()],
    ///         max_query_len: None,
    ///         low_memory: false,
    ///         prefix_result: None,
    ///     };
    ///     let result_without_prefix = sufr.bisect(opts_without_prefix)?;
    ///     // ... both queries appear in the suffix array
    ///     assert_eq!(
    ///         result_without_prefix,
    ///         vec![
    ///             BisectResult { query_num: 0, query: "AC".to_string(), count: 2, first_position: 1, last_position: 2 }, 
    ///             BisectResult { query_num: 1, query: "CG".to_string(), count: 2, first_position: 3, last_position: 4 }]
    ///     );
    ///     // bisect within the range of a prefix result:
    ///     let prefix_opts = BisectOptions {
    ///         queries: vec!["A".to_string()],
    ///         max_query_len: None,
    ///         low_memory: false,
    ///         prefix_result: None,
    ///     };
    ///     let prefix_result = sufr.bisect(prefix_opts)?[0].clone();
    ///     let opts_with_prefix = BisectOptions {
    ///         queries: vec!["AC".to_string(), "CG".to_string()],
    ///         max_query_len: None,
    ///         low_memory: false,
    ///         prefix_result: Some(prefix_result),
    ///     };
    ///     let result_with_prefix = sufr.bisect(opts_with_prefix)?;
    ///     // ... the query AC is found within the range of the prefix result for A, but CG is not.
    ///     assert_eq!(
    ///         result_with_prefix,
    ///         vec![
    ///             BisectResult { query_num: 0, query: "AC".to_string(), count: 2, first_position: 1, last_position: 2 }, 
    ///             BisectResult { query_num: 1, query: "CG".to_string(), count: 0, first_position: 0, last_position: 0 }]
    ///     );
    ///     Ok(())
    /// }
    /// ```
    pub fn bisect(&mut self, args: BisectOptions) -> Result<Vec<BisectResult>> {
        // Set memory mode
        self.query_low_memory = args.low_memory;
        if !self.query_low_memory {
            self.set_suffix_array_mem(args.max_query_len)?;
        }

        // Construct SufrSearch factory 
        let now = Instant::now();
        let new_search = || -> Result<RefCell<SufrSearch<T>>> {
            let suffix_array_file: FileAccess<T> = FileAccess::new(
                &self.filename,
                self.suffix_array_pos as u64,
                self.len_suffixes.to_usize(),
            )?;
            let text_file: FileAccess<u8> = FileAccess::new(
                &self.filename,
                self.text_pos as u64,
                self.text_len.to_usize(),
            )?;
            let search_args = SufrSearchArgs {
                text: &self.text,
                text_len: self.text_len.to_usize(),
                text_file,
                file: suffix_array_file,
                suffix_array: &self.suffix_array_mem,
                rank: &self.suffix_array_rank_mem,
                len_suffixes: self.len_suffixes.to_usize(),
                sort_type: &self.sort_type,
                max_query_len: args.max_query_len,
            };
            Ok(RefCell::new(SufrSearch::new(search_args)))
        };

        // Retrieve the prefix result's index range. 
        // If no result was passed, deafult to the full range of the suffix array.
        let n = self.len_suffixes.to_usize() - 1;
        let search_range = match args.prefix_result {
            Some(result)    => (result.first_position, result.last_position),
            _               => (0, n),
        };
        
        // Bisect each query in its own thread
        let thread_local_search: ThreadLocal<RefCell<SufrSearch<T>>> =
            ThreadLocal::new();
        let mut res: Vec<_> = args
            .queries
            .clone()
            .into_par_iter()
            .enumerate()
            .flat_map(|(query_num, query)| -> Result<BisectResult> {
                let mut search = 
                    thread_local_search.get_or_try(new_search)?.borrow_mut();
                search.bisect(query_num, &query, search_range.0, search_range.1)
            })
            .collect(); 
        res.sort_by_key(|r| r.query_num);

        info!(
            "Bisection of {} queries finished in {:?}",
            args.queries.len(),
            now.elapsed()
        );
        
        Ok(res)
    }

    // --------------------------------------------------
    /// Count the occurrences of queries in a suffix array
    ///
    /// ```
    /// use anyhow::Result;
    /// use libsufr::{types::{CountOptions, CountResult}, sufr_file::SufrFile};
    ///
    /// fn main() -> Result<()> {
    ///     let mut sufr = SufrFile::<u32>::read("../data/inputs/1.sufr", false)?;
    ///     let opts = CountOptions {
    ///         queries: vec!["AC".to_string(), "AG".to_string(), "GT".to_string()],
    ///         max_query_len: None,
    ///         low_memory: true,
    ///     };
    ///     let res = sufr.count(opts)?;
    ///     let expected = vec![
    ///         CountResult {
    ///             query_num: 0,
    ///             query: "AC".to_string(),
    ///             count: 2,
    ///         },
    ///         CountResult {
    ///             query_num: 1,
    ///             query: "AG".to_string(),
    ///             count: 0,
    ///         },
    ///         CountResult {
    ///             query_num: 2,
    ///             query: "GT".to_string(),
    ///             count: 2,
    ///         },
    ///     ];
    ///     assert_eq!(expected, res);
    ///
    ///     Ok(())
    /// }
    /// ```
    pub fn count(&mut self, args: CountOptions) -> Result<Vec<CountResult>> {
        let search_args = SearchOptions {
            queries: args.queries,
            max_query_len: args.max_query_len,
            low_memory: args.low_memory,
            find_suffixes: false,
        };

        let counts: Vec<_> = self
            .suffix_search(&search_args)?
            .into_iter()
            .map(|res| CountResult {
                query_num: res.query_num,
                query: res.query.clone(),
                count: res.locations.map_or(0, |loc| loc.ranks.len()),
            })
            .collect();

        Ok(counts)
    }

    // --------------------------------------------------
    /// Search the suffix array.
    ///
    /// Args:
    /// * `args`: `SearchOptions`
    ///
    ///
    /// E.g., given a text of "ACGTNNACGT", a search for "ACG" and "GCC":
    ///
    /// ```
    /// use anyhow::Result;
    /// use libsufr::{
    ///     sufr_file::SufrFile,
    ///     types::{SearchOptions, SearchResult, SearchResultLocations}
    /// };
    ///
    /// fn main() -> Result<()> {
    ///     let mut sufr = SufrFile::<u32>::read("../data/inputs/1.sufr", false)?;
    ///     let opts = SearchOptions {
    ///         queries: vec!["ACG".to_string(), "GGC".to_string()],
    ///         max_query_len: None,
    ///         find_suffixes: true,
    ///         low_memory: true,
    ///     };
    ///     let expected: Vec<SearchResult<u32>> = vec![
    ///         SearchResult {
    ///             query_num: 0,
    ///             query: "ACG".to_string(),
    ///             locations: Some(
    ///                 SearchResultLocations {
    ///                     ranks: 1..3,
    ///                     suffixes: vec![
    ///                         6,
    ///                         0,
    ///                     ],
    ///                 },
    ///             ),
    ///         },
    ///         SearchResult {
    ///             query_num: 1,
    ///             query: "GGC".to_string(),
    ///             locations: None,
    ///         },
    ///     ];
    ///     assert_eq!(expected, sufr.suffix_search(&opts)?);
    ///     Ok(())
    /// }
    /// ```
    ///
    pub fn suffix_search(
        &mut self,
        args: &SearchOptions,
    ) -> Result<Vec<SearchResult<T>>> {
        self.query_low_memory = args.low_memory;

        if !self.query_low_memory {
            self.set_suffix_array_mem(args.max_query_len)?;
        }

        let now = Instant::now();
        let new_search = || -> Result<RefCell<SufrSearch<T>>> {
            let suffix_array_file: FileAccess<T> = FileAccess::new(
                &self.filename,
                self.suffix_array_pos as u64,
                self.len_suffixes.to_usize(),
            )?;
            let text_file: FileAccess<u8> = FileAccess::new(
                &self.filename,
                self.text_pos as u64,
                self.text_len.to_usize(),
            )?;
            let search_args = SufrSearchArgs {
                text: &self.text,
                text_len: self.text_len.to_usize(),
                text_file,
                file: suffix_array_file,
                suffix_array: &self.suffix_array_mem,
                rank: &self.suffix_array_rank_mem,
                len_suffixes: self.len_suffixes.to_usize(),
                sort_type: &self.sort_type,
                max_query_len: args.max_query_len,
            };
            Ok(RefCell::new(SufrSearch::new(search_args)))
        };

        let thread_local_search: ThreadLocal<RefCell<SufrSearch<T>>> =
            ThreadLocal::new();

        let mut res: Vec<_> = args
            .queries
            .clone()
            .into_par_iter()
            .enumerate()
            .flat_map(|(query_num, query)| -> Result<SearchResult<T>> {
                let mut search =
                    thread_local_search.get_or_try(new_search)?.borrow_mut();
                search.search(query_num, &query, args.find_suffixes)
            })
            .collect();
        res.sort_by_key(|r| r.query_num);

        info!(
            "Search of {} queries finished in {:?}",
            args.queries.len(),
            now.elapsed()
        );

        Ok(res)
    }

    // --------------------------------------------------
    /// Extract the suffixes for a given set of queries
    ///
    /// Args:
    /// * `args`: `ExtractOptions`
    ///
    /// Given the text "ACGTNNACGT":
    ///
    /// ```
    /// use anyhow::Result;
    /// use libsufr::{
    ///     sufr_file::SufrFile,
    ///     types::{ExtractOptions, ExtractResult, ExtractSequence}
    /// };
    ///
    /// fn main() -> Result<()> {
    ///     let mut sufr = SufrFile::<u32>::read("../data/inputs/1.sufr", true)?;
    ///     let extract_opts = ExtractOptions {
    ///         queries: vec!["ACG".to_string(), "GTC".to_string()],
    ///         max_query_len: None,
    ///         low_memory: true,
    ///         prefix_len: Some(1),
    ///         suffix_len: Some(3),
    ///     };
    ///
    ///     let expected = vec![
    ///         ExtractResult {
    ///             query_num: 0,
    ///             query: "ACG".to_string(),
    ///             sequences: vec![
    ///                 ExtractSequence {
    ///                     suffix: 6,
    ///                     rank: 1,
    ///                     sequence_name: "1".to_string(),
    ///                     sequence_start: 0,
    ///                     sequence_range: 5..9,
    ///                     suffix_offset: 1,
    ///                 },
    ///                 ExtractSequence {
    ///                     suffix: 0,
    ///                     rank: 2,
    ///                     sequence_name: "1".to_string(),
    ///                     sequence_start: 0,
    ///                     sequence_range: 0..3,
    ///                     suffix_offset: 0,
    ///                 },
    ///             ],
    ///         },
    ///         ExtractResult {
    ///             query_num: 1,
    ///             query: "GTC".to_string(),
    ///             sequences: vec![],
    ///         },
    ///     ];
    ///
    ///     assert_eq!(expected, sufr.extract(extract_opts)?);
    ///     Ok(())
    /// }
    /// ```
    ///
    pub fn extract(&mut self, args: ExtractOptions) -> Result<Vec<ExtractResult>> {
        let search_args = SearchOptions {
            queries: args.queries,
            max_query_len: args.max_query_len,
            low_memory: args.low_memory,
            find_suffixes: true,
        };
        let search_result = &self.suffix_search(&search_args)?;
        let seq_starts = self.sequence_starts.clone();
        let seq_names = self.sequence_names.clone();
        let text_len = self.text_len.to_usize();
        let now = Instant::now();

        // Augment the search with relative sequence positions
        let extract_result = search_result
            .into_par_iter()
            .map(|res| {
                let sequences: Vec<ExtractSequence> = match &res.locations {
                    Some(locs) => locs
                        .ranks
                        .clone()
                        .zip(locs.suffixes.clone())
                        .map(|(rank, suffix)| {
                            let i =
                                seq_starts.partition_point(|&val| val <= suffix) - 1;
                            let sequence_start = seq_starts[i].to_usize();
                            let seq_end = if i == seq_starts.len() - 1 {
                                text_len
                            } else {
                                seq_starts[i + 1].to_usize()
                            };
                            let suffix = suffix.to_usize();
                            let relative_suffix_start = suffix - sequence_start;
                            let context_start = relative_suffix_start
                                .saturating_sub(args.prefix_len.unwrap_or(0));
                            let context_end = min(
                                args.suffix_len
                                    .map_or(seq_end, |len| relative_suffix_start + len),
                                seq_end,
                            );
                            ExtractSequence {
                                rank,
                                suffix,
                                sequence_name: seq_names[i].clone(),
                                sequence_start,
                                sequence_range: (context_start..context_end),
                                suffix_offset: relative_suffix_start - context_start,
                            }
                        })
                        .collect(),
                    _ => vec![],
                };

                ExtractResult {
                    query_num: res.query_num,
                    query: res.query.clone(),
                    sequences,
                }
            })
            .collect();

        info!("Adding locate data finished in {:?}", now.elapsed());

        Ok(extract_result)
    }

    // --------------------------------------------------
    /// Print suffixes
    ///
    /// Args:
    /// * `args`: `SearchOptions`
    ///
    /// Given a text of "ACGTNNACGT":
    ///
    /// ```
    /// use anyhow::Result;
    /// use libsufr::{sufr_file::SufrFile, types::ListOptions};
    /// use std::fs;
    /// use tempfile::NamedTempFile;
    ///
    /// fn main() -> Result<()> {
    ///     let mut sufr = SufrFile::<u32>::read("../data/inputs/1.sufr", false)?;
    ///     let outfile = NamedTempFile::new()?;
    ///     let list_opts = ListOptions {
    ///         ranks: vec![],
    ///         show_rank: true,
    ///         show_suffix: true,
    ///         show_lcp: true,
    ///         len: None,
    ///         number: None,
    ///         output: Some(outfile.path().to_string_lossy().to_string()),
    ///     };
    ///     sufr.list(list_opts)?;
    ///
    ///     let expected = vec![
    ///         " 0 10  0 $",
    ///         " 1  6  0 ACGT$",
    ///         " 2  0  4 ACGTNNACGT$",
    ///         " 3  7  0 CGT$",
    ///         " 4  1  3 CGTNNACGT$",
    ///         " 5  8  0 GT$",
    ///         " 6  2  2 GTNNACGT$",
    ///         " 7  9  0 T$",
    ///         " 8  3  1 TNNACGT$",
    ///         "",
    ///     ]
    ///     .join("\n");
    ///     let output = fs::read_to_string(&outfile)?;
    ///     //let actual = String::from_utf8(output.to_vec())?;
    ///     assert_eq!(expected, output);
    ///     fs::remove_file(&outfile);
    ///     Ok(())
    /// }
    /// ```
    ///
    pub fn list(&mut self, args: ListOptions) -> Result<()> {
        let width = self.text_len.to_string().len();
        let text_len = self.text_len.to_usize();
        let suffix_len = args.len.unwrap_or(text_len);
        let mut output: Box<dyn Write> = match &args.output {
            Some(filename) => Box::new(
                File::create(filename).map_err(|e| anyhow!("{filename}: {e}"))?,
            ),
            _ => Box::new(io::stdout()),
        };

        let mut print = |rank: usize, suffix: usize, lcp: T| -> Result<()> {
            let end = if suffix + suffix_len > text_len {
                text_len
            } else {
                suffix + suffix_len
            };

            let rank_display = if args.show_rank {
                format!("{rank:width$} ")
            } else {
                "".to_string()
            };

            let suffix_display = if args.show_suffix {
                format!("{suffix:width$} ")
            } else {
                "".to_string()
            };
            let lcp_display = if args.show_lcp {
                format!("{:width$} ", lcp)
            } else {
                "".to_string()
            };

            writeln!(
                output,
                "{rank_display}{suffix_display}{lcp_display}{}",
                String::from_utf8(self.text_file.get_range(suffix..end)?)?
            )?;
            Ok(())
        };

        let number = args.number.unwrap_or(0);
        if args.ranks.is_empty() {
            self.suffix_array_file.reset();
            for (rank, suffix) in self.suffix_array_file.iter().enumerate() {
                print(rank, suffix.to_usize(), self.lcp_file.get(rank).unwrap())?;

                if number > 0 && rank == number - 1 {
                    break;
                }
            }
        } else {
            for rank in args.ranks {
                if let Some(suffix) = self.suffix_array_file.get(rank) {
                    print(rank, suffix.to_usize(), self.lcp_file.get(rank).unwrap())?;
                } else {
                    eprintln!("Invalid rank: {rank}");
                }
            }
        }

        Ok(())
    }

    // --------------------------------------------------
    /// Find the positions of queries in a suffix array.
    ///
    /// Args:
    /// * `args`: `SearchOptions`
    ///
    /// Given a text of "ACGTNNACGT":
    ///
    /// ```
    /// use anyhow::Result;
    /// use libsufr::{
    ///     sufr_file::SufrFile,
    ///     types::{LocateOptions, LocatePosition, LocateResult},
    /// };
    ///
    /// fn main() -> Result<()> {
    ///     let mut sufr = SufrFile::<u32>::read("../data/inputs/1.sufr", true)?;
    ///     let opts = LocateOptions {
    ///         queries: vec!["ACG".to_string(), "GGC".to_string()],
    ///         max_query_len: None,
    ///         low_memory: true,
    ///     };
    ///     let expected = vec![
    ///         LocateResult {
    ///             query_num: 0,
    ///             query: "ACG".to_string(),
    ///             positions: vec![
    ///                 LocatePosition {
    ///                     suffix: 6,
    ///                     rank: 1,
    ///                     sequence_name: "1".to_string(),
    ///                     sequence_position: 6,
    ///                 },
    ///                 LocatePosition {
    ///                     suffix: 0,
    ///                     rank: 2,
    ///                     sequence_name: "1".to_string(),
    ///                     sequence_position: 0,
    ///                 },
    ///             ],
    ///         },
    ///         LocateResult {
    ///             query_num: 1,
    ///             query: "GGC".to_string(),
    ///             positions: vec![],
    ///         },
    ///     ];
    ///
    ///     assert_eq!(expected, sufr.locate(opts)?);
    ///     Ok(())
    /// }
    /// ```
    ///
    pub fn locate(&mut self, args: LocateOptions) -> Result<Vec<LocateResult>> {
        let search_opts = SearchOptions {
            queries: args.queries,
            low_memory: args.low_memory,
            max_query_len: args.max_query_len,
            find_suffixes: true,
        };
        let search_result = &self.suffix_search(&search_opts)?;
        let seq_starts = self.sequence_starts.clone();
        let seq_names = self.sequence_names.clone();
        let mut locate_result: Vec<LocateResult> = vec![];
        let now = Instant::now();

        // Augment the search with relative sequence positions
        for res in search_result {
            let mut positions = vec![];
            if let Some(locs) = &res.locations {
                for (rank, suffix) in locs.ranks.clone().zip(locs.suffixes.clone()) {
                    let i = seq_starts.partition_point(|&val| val <= suffix) - 1;
                    positions.push(LocatePosition {
                        rank,
                        suffix: suffix.to_usize(),
                        sequence_name: seq_names[i].clone(),
                        sequence_position: (suffix - seq_starts[i]).to_usize(),
                    })
                }
            }
            locate_result.push(LocateResult {
                query_num: res.query_num,
                query: res.query.clone(),
                positions,
            });
        }

        info!("Adding locate data finished in {:?}", now.elapsed());

        Ok(locate_result)
    }
}

// --------------------------------------------------
#[cfg(test)]
mod test {
    use crate::{
        sufr_file::SufrFile,
        types::{
            ExtractOptions, ExtractResult, ExtractSequence, LocateOptions,
            LocatePosition, LocateResult, BisectOptions, BisectResult,
        },
    };
    use anyhow::Result;

    // --------------------------------------------------
    #[test]
    fn test_extract() -> Result<()> {
        // cargo run -- ex data/expected/1.sufr AC GT XX
        let mut sufr_file: SufrFile<u32> =
            SufrFile::read("../data/expected/1.sufr", false)?;

        let opts = ExtractOptions {
            queries: vec!["AC".to_string(), "GT".to_string(), "XX".to_string()],
            max_query_len: None,
            low_memory: true,
            prefix_len: Some(1),
            suffix_len: Some(3),
        };

        let expected = [
            ExtractResult {
                query_num: 0,
                query: "AC".to_string(),
                sequences: vec![
                    ExtractSequence {
                        suffix: 6,
                        rank: 1,
                        sequence_name: "1".to_string(),
                        sequence_start: 0,
                        sequence_range: 5..9,
                        suffix_offset: 1,
                    },
                    ExtractSequence {
                        suffix: 0,
                        rank: 2,
                        sequence_name: "1".to_string(),
                        sequence_start: 0,
                        sequence_range: 0..3,
                        suffix_offset: 0,
                    },
                ],
            },
            ExtractResult {
                query_num: 1,
                query: "GT".to_string(),
                sequences: vec![
                    ExtractSequence {
                        suffix: 8,
                        rank: 5,
                        sequence_name: "1".to_string(),
                        sequence_start: 0,
                        sequence_range: 7..11,
                        suffix_offset: 1,
                    },
                    ExtractSequence {
                        suffix: 2,
                        rank: 6,
                        sequence_name: "1".to_string(),
                        sequence_start: 0,
                        sequence_range: 1..5,
                        suffix_offset: 1,
                    },
                ],
            },
            ExtractResult {
                query_num: 2,
                query: "XX".to_string(),
                sequences: vec![],
            },
        ];

        let res = sufr_file.extract(opts);
        assert!(res.is_ok());
        assert_eq!(res.unwrap(), expected);

        Ok(())
    }

    // --------------------------------------------------
    #[test]
    fn test_locate() -> Result<()> {
        //  0  14: #
        //  1   0: AABABABABBABAB#
        //  2  12: AB#
        //  3  10: ABAB#
        //  4   1: ABABABABBABAB#
        //  5   3: ABABABBABAB#
        //  6   5: ABABBABAB#
        //  7   7: ABBABAB#
        //  8  13: B#
        //  9  11: BAB#
        // 10   9: BABAB#
        // 11   2: BABABABBABAB#
        // 12   4: BABABBABAB#
        // 13   6: BABBABAB#
        // 14   8: BBABAB#

        let mut sufr_file: SufrFile<u32> =
            SufrFile::read("../data/expected/abba.sufr", false)?;

        for low_memory in &[true, false] {
            let args = LocateOptions {
                queries: vec!["A".to_string()],
                max_query_len: None,
                low_memory: *low_memory,
            };
            let res = sufr_file.locate(args);
            assert!(res.is_ok());
            let res = res.unwrap();
            assert_eq!(res.len(), 1);

            assert_eq!(
                res,
                vec![LocateResult {
                    query_num: 0,
                    query: "A".to_string(),
                    positions: vec![
                        LocatePosition {
                            rank: 1,
                            suffix: 0,
                            sequence_name: "1".to_string(),
                            sequence_position: 0,
                        },
                        LocatePosition {
                            rank: 2,
                            suffix: 12,
                            sequence_name: "1".to_string(),
                            sequence_position: 12,
                        },
                        LocatePosition {
                            rank: 3,
                            suffix: 10,
                            sequence_name: "1".to_string(),
                            sequence_position: 10,
                        },
                        LocatePosition {
                            rank: 4,
                            suffix: 1,
                            sequence_name: "1".to_string(),
                            sequence_position: 1,
                        },
                        LocatePosition {
                            rank: 5,
                            suffix: 3,
                            sequence_name: "1".to_string(),
                            sequence_position: 3,
                        },
                        LocatePosition {
                            rank: 6,
                            suffix: 5,
                            sequence_name: "1".to_string(),
                            sequence_position: 5,
                        },
                        LocatePosition {
                            rank: 7,
                            suffix: 7,
                            sequence_name: "1".to_string(),
                            sequence_position: 7,
                        },
                    ]
                }]
            );
        }

        for low_memory in &[true, false] {
            let args = LocateOptions {
                queries: vec!["B".to_string()],
                max_query_len: None,
                low_memory: *low_memory,
            };
            let res = sufr_file.locate(args);
            assert!(res.is_ok());
            let res = res.unwrap();
            assert_eq!(res.len(), 1);

            assert_eq!(
                res,
                vec![LocateResult {
                    query_num: 0,
                    query: "B".to_string(),
                    positions: vec![
                        LocatePosition {
                            rank: 8,
                            suffix: 13,
                            sequence_name: "1".to_string(),
                            sequence_position: 13,
                        },
                        LocatePosition {
                            rank: 9,
                            suffix: 11,
                            sequence_name: "1".to_string(),
                            sequence_position: 11,
                        },
                        LocatePosition {
                            rank: 10,
                            suffix: 9,
                            sequence_name: "1".to_string(),
                            sequence_position: 9,
                        },
                        LocatePosition {
                            rank: 11,
                            suffix: 2,
                            sequence_name: "1".to_string(),
                            sequence_position: 2,
                        },
                        LocatePosition {
                            rank: 12,
                            suffix: 4,
                            sequence_name: "1".to_string(),
                            sequence_position: 4,
                        },
                        LocatePosition {
                            rank: 13,
                            suffix: 6,
                            sequence_name: "1".to_string(),
                            sequence_position: 6,
                        },
                        LocatePosition {
                            rank: 14,
                            suffix: 8,
                            sequence_name: "1".to_string(),
                            sequence_position: 8,
                        },
                    ]
                }]
            );
        }

        for low_memory in &[true, false] {
            let args = LocateOptions {
                queries: vec!["ABAB".to_string()],
                max_query_len: None,
                low_memory: *low_memory,
            };
            let res = sufr_file.locate(args);
            assert!(res.is_ok());
            let res = res.unwrap();
            assert_eq!(res.len(), 1);

            assert_eq!(
                res,
                vec![LocateResult {
                    query_num: 0,
                    query: "ABAB".to_string(),
                    positions: vec![
                        LocatePosition {
                            rank: 3,
                            suffix: 10,
                            sequence_name: "1".to_string(),
                            sequence_position: 10,
                        },
                        LocatePosition {
                            rank: 4,
                            suffix: 1,
                            sequence_name: "1".to_string(),
                            sequence_position: 1,
                        },
                        LocatePosition {
                            rank: 5,
                            suffix: 3,
                            sequence_name: "1".to_string(),
                            sequence_position: 3,
                        },
                        LocatePosition {
                            rank: 6,
                            suffix: 5,
                            sequence_name: "1".to_string(),
                            sequence_position: 5,
                        },
                    ]
                }]
            );
        }

        for low_memory in &[true, false] {
            let args = LocateOptions {
                queries: vec!["ABABB".to_string()],
                max_query_len: None,
                low_memory: *low_memory,
            };
            let res = sufr_file.locate(args);
            assert!(res.is_ok());
            let res = res.unwrap();
            assert_eq!(res.len(), 1);
            assert_eq!(
                res,
                vec![LocateResult {
                    query_num: 0,
                    query: "ABABB".to_string(),
                    positions: vec![LocatePosition {
                        rank: 6,
                        suffix: 5,
                        sequence_name: "1".to_string(),
                        sequence_position: 5,
                    },]
                }]
            );
        }

        for low_memory in &[true, false] {
            let args = LocateOptions {
                queries: vec!["BBBB".to_string()],
                max_query_len: None,
                low_memory: *low_memory,
            };
            let res = sufr_file.locate(args);
            assert!(res.is_ok());
            assert!(res.is_ok());
            let res = res.unwrap();
            assert_eq!(res.len(), 1);
            assert_eq!(
                res,
                vec![LocateResult {
                    query_num: 0,
                    query: "BBBB".to_string(),
                    positions: vec![],
                }]
            );
        }

        Ok(())
    }

    // --------------------------------------------------
    #[test]
    fn test_file_access() -> Result<()> {
        let input_file = "../data/expected/abba.sufr";
        let mut sufr_file: SufrFile<u32> = SufrFile::read(input_file, false)?;
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
    #[test]
    fn test_bisect() -> Result<()> {
        let mut sufr = SufrFile::<u32>::read("../data/inputs/3.sufr", false)?;
        // bisect "A"
        let prefix = vec!["A".to_string()];
        let prefix_result = sufr.bisect(BisectOptions {
            queries: prefix,
            max_query_len: None,
            low_memory: false,
            prefix_result: None
        })?[0].clone();
        // bisect "AA", "AC", "AG", "AT", "AN" within the range of "A".
        let queries = vec!["AA".to_string(), "AC".to_string(), "AG".to_string(), "AT".to_string(), "AN".to_string()];
        let queries_result = sufr.bisect(BisectOptions {
            queries: queries,
            max_query_len: None,
            low_memory: false,
            prefix_result: Some(prefix_result.clone()),
        })?;
        // because we queried all of the possible suffixes to "A", 
        // the count of "A" should be the sum of counts of queries. 
        assert_eq!(
            prefix_result.count,
            queries_result.iter().map(|res| res.count).sum(),
        );
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
