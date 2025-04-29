//! Create and query suffix arrays
use crate::{
    sufr_builder::SufrBuilder,
    sufr_file::SufrFile,
    types::{
        BisectOptions, BisectResult, CountOptions, CountResult, ExtractOptions, ExtractResult, ListOptions,
        LocateOptions, LocateResult, SufrBuilderArgs, SufrMetadata,
    },
};
use anyhow::Result;

// --------------------------------------------------
pub(crate) trait SuffixArrayTrait: Send + Sync {
    fn count(&mut self, args: CountOptions) -> Result<Vec<CountResult>>;
    fn extract(&mut self, args: ExtractOptions) -> Result<Vec<ExtractResult>>;
    fn list(&mut self, args: ListOptions) -> Result<()>;
    fn locate(&mut self, args: LocateOptions) -> Result<Vec<LocateResult>>;
    fn metadata(&self) -> Result<SufrMetadata>;
    fn string_at(&mut self, pos: usize, len: Option<usize>) -> Result<String>;
    fn bisect(&mut self, args: BisectOptions) -> Result<Vec<BisectResult>>; 
}

// --------------------------------------------------
pub(crate) struct SuffixArray32 {
    inner: SufrFile<u32>,
}

impl SuffixArrayTrait for SuffixArray32 {
    fn count(&mut self, args: CountOptions) -> Result<Vec<CountResult>> {
        self.inner.count(args)
    }

    fn extract(&mut self, args: ExtractOptions) -> Result<Vec<ExtractResult>> {
        self.inner.extract(args)
    }

    fn metadata(&self) -> Result<SufrMetadata> {
        self.inner.metadata()
    }

    fn list(&mut self, args: ListOptions) -> Result<()> {
        self.inner.list(args)
    }

    fn locate(&mut self, args: LocateOptions) -> Result<Vec<LocateResult>> {
        self.inner.locate(args)
    }

    fn string_at(&mut self, pos: usize, len: Option<usize>) -> Result<String> {
        self.inner.string_at(pos, len)
    }

    fn bisect(&mut self, args: BisectOptions) -> Result<Vec<BisectResult>> {
        self.inner.bisect(args)
    }
}

pub(crate) struct SuffixArray64 {
    inner: SufrFile<u64>,
}

// --------------------------------------------------
impl SuffixArrayTrait for SuffixArray64 {
    fn count(&mut self, args: CountOptions) -> Result<Vec<CountResult>> {
        self.inner.count(args)
    }

    fn extract(&mut self, args: ExtractOptions) -> Result<Vec<ExtractResult>> {
        self.inner.extract(args)
    }

    fn metadata(&self) -> Result<SufrMetadata> {
        self.inner.metadata()
    }

    fn list(&mut self, args: ListOptions) -> Result<()> {
        self.inner.list(args)
    }

    fn locate(&mut self, args: LocateOptions) -> Result<Vec<LocateResult>> {
        self.inner.locate(args)
    }

    fn string_at(&mut self, pos: usize, len: Option<usize>) -> Result<String> {
        self.inner.string_at(pos, len)
    }

    fn bisect(&mut self, args: BisectOptions) -> Result<Vec<BisectResult>> {
        self.inner.bisect(args)
    }
}

// --------------------------------------------------
/// Struct to create and read suffix arrays
pub struct SuffixArray {
    inner: Box<dyn SuffixArrayTrait>,
}

// --------------------------------------------------
impl SuffixArray {
    /// Create a new suffix array and return an in-memory struct.
    /// This is a combination of `write` and `read`.
    ///     
    /// ```
    /// use anyhow::Result;
    /// use libsufr::{
    ///     types::SufrBuilderArgs,
    ///     suffix_array::SuffixArray,
    ///     util::read_sequence_file,
    /// };
    /// use std::{fs, path::Path};
    ///
    /// fn main() -> Result<()> {
    ///     let sequence_delimiter = b'%';
    ///     let seq_data = read_sequence_file(Path::new("../data/inputs/3.fa"), sequence_delimiter)?;
    ///     let outfile = "3.sufr";
    ///     let builder_args = SufrBuilderArgs {
    ///         text: seq_data.seq,
    ///         path: Some(outfile.to_string()),
    ///         low_memory: true,
    ///         max_query_len: None,
    ///         is_dna: true,
    ///         allow_ambiguity: false,
    ///         ignore_softmask: true,
    ///         sequence_starts: seq_data.start_positions.into_iter().collect(),
    ///         sequence_names: seq_data.sequence_names,
    ///         num_partitions: 16,
    ///         seed_mask: None,
    ///         random_seed: 42,
    ///     };
    ///
    ///     let suffix_array = SuffixArray::new(builder_args)?;
    ///     let meta = suffix_array.metadata()?;
    ///     assert_eq!(meta.text_len, 113);
    ///     assert_eq!(meta.len_suffixes, 101);
    ///
    ///     fs::remove_file(&outfile)?;
    ///
    ///     Ok(())
    /// }
    /// ```
    ///
    pub fn new(args: SufrBuilderArgs) -> Result<SuffixArray> {
        let low_memory = args.low_memory;
        let path = Self::write(args)?;
        Self::read(&path, low_memory)
    }


    // --------------------------------------------------
    /// Bisect the index range of occurences of queries.
    /// If the index range of a prefix is already known,
    /// or if it is desirable to avoid enumerating every match,
    /// this method can be used as a faster stand-in for `count`
    /// ```
    /// use anyhow::Result;
    /// use libsufr::{types::{BisectOptions, BisectResult}, suffix_array::SuffixArray};
    ///
    /// fn main() -> Result<()> {
    ///     let mut suffix_array = SuffixArray::read("../data/inputs/1.sufr", false)?;
    ///     let opts = BisectOptions {
    ///         queries: vec!["A".to_string(), "AC".to_string(), "ACA".to_string(), "ACG".to_string(), "ACT".to_string(), "ACC".to_string()],
    ///         max_query_len: None,
    ///         low_memory: false,
    ///         prefix_result: None,
    ///     };
    ///     let res = suffix_array.bisect(opts)?;
    ///     let expected = vec![
    ///         BisectResult { 
    ///             query_num: 0, 
    ///             query: "A".to_string(),  
    ///             count: 2,  
    ///             first_position: 1,  
    ///             last_position: 2  
    ///         },  
    ///         BisectResult { 
    ///             query_num: 1, 
    ///             query: "AC".to_string(), 
    ///             count: 2, 
    ///             first_position: 1, 
    ///             last_position: 2 
    ///         }, 
    ///         BisectResult { 
    ///             query_num: 2, 
    ///             query: "ACA".to_string(), 
    ///             count: 0, 
    ///             first_position: 0, 
    ///             last_position: 0 
    ///         }, 
    ///         BisectResult { 
    ///             query_num: 3, 
    ///             query: "ACG".to_string(), 
    ///             count: 2, 
    ///             first_position: 1, 
    ///             last_position: 2 
    ///         }, 
    ///         BisectResult { 
    ///             query_num: 4, 
    ///             query: "ACT".to_string(), 
    ///             count: 0, 
    ///             first_position: 0, 
    ///             last_position: 0 
    ///         }, 
    ///         BisectResult {
    ///             query_num: 5, 
    ///             query: "ACC".to_string(), 
    ///             count: 0, 
    ///             first_position: 0, 
    ///             last_position: 0 
    ///         }
    ///     ];
    ///     assert_eq!(res, expected);
    ///     
    ///     let mut suffix_array = SuffixArray::read("../data/inputs/3.sufr", false)?;
    ///     let opts1 = BisectOptions {
    ///         queries: vec!["AC".to_string(), "ACA".to_string(), "ACG".to_string(), "ACT".to_string(), "ACC".to_string()],
    ///         max_query_len: None,
    ///         low_memory: false,
    ///         prefix_result: None,
    ///     };
    ///     let res1 = suffix_array.bisect(opts1)?;
    ///     let opts2 = BisectOptions {
    ///         queries: vec!["AC".to_string(), "ACA".to_string(), "ACG".to_string(), "ACT".to_string(), "ACC".to_string()],
    ///         max_query_len: None,
    ///         low_memory: false,
    ///         prefix_result: Some(res1[0].clone()),
    ///     };
    ///     let res2 = suffix_array.bisect(opts2)?;
    ///     assert_eq!(res1, res2);
    ///
    ///     Ok(())
    /// }
    /// ```
    pub fn bisect(&mut self, args: BisectOptions) -> Result<Vec<BisectResult>> {
        self.inner.bisect(args)
    }

    // --------------------------------------------------
    /// Count instances of queries
    ///
    /// ```
    /// use anyhow::Result;
    /// use libsufr::{
    ///     suffix_array::SuffixArray,
    ///     types::{CountOptions, CountResult}
    /// };
    ///
    /// fn main() -> Result<()> {
    ///     let mut suffix_array = SuffixArray::read("../data/inputs/1.sufr", true)?;
    ///     let count_args = CountOptions {
    ///         queries: vec!["AC".to_string(), "GG".to_string(), "CG".to_string()],
    ///         max_query_len: None,
    ///         low_memory: true
    ///     };
    ///     let res = suffix_array.count(count_args)?;
    ///     let expected = vec![
    ///         CountResult {
    ///             query_num: 0,
    ///             query: "AC".to_string(),
    ///             count: 2,
    ///         },
    ///         CountResult {
    ///             query_num: 1,
    ///             query: "GG".to_string(),
    ///             count: 0,
    ///         },
    ///         CountResult {
    ///             query_num: 2,
    ///             query: "CG".to_string(),
    ///             count: 2,
    ///         },
    ///     ];
    ///     assert_eq!(res, expected);
    ///
    ///     Ok(())
    /// }
    /// ```
    ///
    pub fn count(&mut self, args: CountOptions) -> Result<Vec<CountResult>> {
        self.inner.count(args)
    }

    // --------------------------------------------------
    /// Extract the suffixes matching given queries
    ///
    /// ```
    /// use anyhow::Result;
    /// use libsufr::{
    ///     suffix_array::SuffixArray,
    ///     types::{ExtractOptions, ExtractResult, ExtractSequence}
    /// };
    ///
    /// fn main() -> Result<()> {
    ///     let mut suffix_array = SuffixArray::read("../data/inputs/1.sufr", true)?;
    ///     let extract_args = ExtractOptions {
    ///         queries: vec!["CGT".to_string(), "GG".to_string()],
    ///         max_query_len: None,
    ///         low_memory: true,
    ///         prefix_len: Some(1),
    ///         suffix_len: None,
    ///     };
    ///     let expected = vec![
    ///         ExtractResult {
    ///             query_num: 0,
    ///             query: "CGT".to_string(),
    ///             sequences: vec![
    ///                 ExtractSequence {
    ///                     suffix: 7,
    ///                     rank: 3,
    ///                     sequence_name: "1".to_string(),
    ///                     sequence_start: 0,
    ///                     sequence_range: 6..11,
    ///                     suffix_offset: 1,
    ///                 },
    ///                 ExtractSequence {
    ///                     suffix: 1,
    ///                     rank: 4,
    ///                     sequence_name: "1".to_string(),
    ///                     sequence_start: 0,
    ///                     sequence_range: 0..11,
    ///                     suffix_offset: 1,
    ///                 },
    ///             ],
    ///         },
    ///         ExtractResult {
    ///             query_num: 1,
    ///             query: "GG".to_string(),
    ///             sequences: vec![],
    ///         },
    ///     ];
    ///     let res = suffix_array.extract(extract_args)?;
    ///     assert_eq!(res, expected);
    ///
    ///     Ok(())
    /// }
    /// ```
    ///
    pub fn extract(&mut self, args: ExtractOptions) -> Result<Vec<ExtractResult>> {
        self.inner.extract(args)
    }

    // --------------------------------------------------
    /// Extract the suffixes matching given queries
    ///
    /// ```
    /// use anyhow::Result;
    /// use libsufr::{suffix_array::SuffixArray, types::ListOptions};
    /// use std::fs;
    /// use tempfile::NamedTempFile;
    ///
    /// fn main() -> Result<()> {
    ///     let mut suffix_array = SuffixArray::read("../data/inputs/1.sufr", true)?;
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
    ///     suffix_array.list(list_opts)?;
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
    ///     assert_eq!(expected, output);
    ///     fs::remove_file(&outfile);
    ///     Ok(())
    ///
    /// }
    /// ```
    ///
    pub fn list(&mut self, args: ListOptions) -> Result<()> {
        self.inner.list(args)
    }

    // --------------------------------------------------
    /// Find the positions of queries in a suffix array
    ///
    /// ```
    /// use anyhow::Result;
    /// use libsufr::{
    ///     suffix_array::SuffixArray,
    ///     types::{LocateOptions, LocatePosition, LocateResult},
    /// };
    ///
    /// fn main() -> Result<()> {
    ///     let mut suffix_array = SuffixArray::read("../data/inputs/1.sufr", true)?;
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
    ///     assert_eq!(expected, suffix_array.locate(opts)?);
    ///     Ok(())
    /// }
    /// ```
    ///
    pub fn locate(&mut self, args: LocateOptions) -> Result<Vec<LocateResult>> {
        self.inner.locate(args)
    }

    // --------------------------------------------------
    /// Get suffix array metadata
    ///
    /// ```
    /// use anyhow::Result;
    /// use libsufr::{suffix_array::SuffixArray, types::SuffixSortType};
    ///
    /// fn main() -> Result<()> {
    ///     let mut suffix_array = SuffixArray::read("../data/inputs/1.sufr", true)?;
    ///     let meta = suffix_array.metadata()?;
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
    ///     Ok(())
    /// }
    /// ```
    ///
    pub fn metadata(&self) -> Result<SufrMetadata> {
        self.inner.metadata()
    }

    /// Read a _.sufr_ file
    ///
    /// Args:
    /// * `filename`: the _.sufr_ file!
    /// * `low_memory`: when `true`, leave text on disk; when `false`, read text into memory
    ///
    /// ```
    /// use anyhow::Result;
    /// use libsufr::suffix_array::SuffixArray;
    ///
    /// fn main() -> Result<()> {
    ///     let suffix_array = SuffixArray::read("../data/inputs/1.sufr", true)?;
    ///     Ok(())
    /// }
    /// ```
    ///
    pub fn read(filename: &str, low_memory: bool) -> Result<SuffixArray> {
        let text_len = crate::util::read_text_length(filename)? as u64;
        let sa: Box<dyn SuffixArrayTrait> = if text_len < u32::MAX as u64 {
            Box::new(SuffixArray32 {
                inner: SufrFile::read(filename, low_memory)?,
            })
        } else {
            Box::new(SuffixArray64 {
                inner: SufrFile::read(filename, low_memory)?,
            })
        };
        Ok(SuffixArray { inner: sa })
    }

    // --------------------------------------------------
    /// Retrieve a suffix
    ///
    /// ```
    /// use anyhow::Result;
    /// use libsufr::suffix_array::SuffixArray;
    ///
    /// fn main() -> Result<()> {
    ///     let mut suffix_array = SuffixArray::read("../data/inputs/1.sufr", true)?;
    ///     assert_eq!(suffix_array.string_at(0, None)?, "ACGTNNACGT$".to_string());
    ///     assert_eq!(suffix_array.string_at(6, Some(3))?, "ACG".to_string());
    ///     Ok(())
    /// }
    /// ```
    ///
    pub fn string_at(&mut self, pos: usize, len: Option<usize>) -> Result<String> {
        self.inner.string_at(pos, len)
    }

    // --------------------------------------------------
    /// Create a new suffix array and write to disk
    /// ```
    /// use anyhow::Result;
    /// use libsufr::{
    ///     types::SufrBuilderArgs,
    ///     suffix_array::SuffixArray,
    ///     util::read_sequence_file,
    /// };
    /// use std::{fs, path::Path};
    ///
    /// fn main() -> Result<()> {
    ///     let sequence_delimiter = b'%';
    ///     let seq_data = read_sequence_file(Path::new("../data/inputs/3.fa"), sequence_delimiter)?;
    ///     let outfile = "3.sufr";
    ///     let builder_args = SufrBuilderArgs {
    ///         text: seq_data.seq,
    ///         path: Some(outfile.to_string()),
    ///         low_memory: true,
    ///         max_query_len: None,
    ///         is_dna: true,
    ///         allow_ambiguity: false,
    ///         ignore_softmask: false,
    ///         sequence_starts: seq_data.start_positions.into_iter().collect(),
    ///         sequence_names: seq_data.sequence_names,
    ///         num_partitions: 16,
    ///         seed_mask: None,
    ///         random_seed: 42,
    ///     };
    ///
    ///     let outpath = SuffixArray::write(builder_args)?;
    ///     fs::remove_file(&outfile)?;
    ///
    ///     Ok(())
    /// }
    /// ```
    ///
    pub fn write(args: SufrBuilderArgs) -> Result<String> {
        let path = if (args.text.len() as u64) < u32::MAX as u64 {
            let builder: SufrBuilder<u32> = SufrBuilder::new(args)?;
            builder.path
        } else {
            let builder: SufrBuilder<u64> = SufrBuilder::new(args)?;
            builder.path
        };

        Ok(path)
    }
}
