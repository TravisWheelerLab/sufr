//! Create and query suffix arrays
use crate::{
    sufr_builder::SufrBuilder,
    sufr_file::SufrFile,
    types::{
        CountOptions, CountResult, ExtractOptions, ExtractResult, ListOptions,
        LocateOptions, LocateResult, SufrBuilderArgs, SufrMetadata,
    },
};
use anyhow::Result;

// --------------------------------------------------
pub(crate) trait SuffixArrayTrait {
    fn count(&mut self, args: CountOptions) -> Result<Vec<CountResult>>;
    fn extract(&mut self, args: ExtractOptions) -> Result<Vec<ExtractResult>>;
    fn list(&mut self, args: ListOptions) -> Result<()>;
    fn locate(&mut self, args: LocateOptions) -> Result<Vec<LocateResult>>;
    fn metadata(&self) -> Result<SufrMetadata>;
    fn string_at(&mut self, pos: usize, len: Option<usize>) -> Result<String>;
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
}

// --------------------------------------------------
/// Struct to create and read suffix arrays
pub struct SuffixArray {
    inner: Box<dyn SuffixArrayTrait>,
}

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
    /// use std::path::Path;
    ///
    /// fn main() -> Result<()> {
    ///     let sequence_delimiter = b'%';
    ///     let seq_data = read_sequence_file(Path::new("../data/inputs/3.fa"), sequence_delimiter)?;
    ///     let builder_args = SufrBuilderArgs {
    ///         text: seq_data.seq,
    ///         path: Some("3.sufr".to_string()),
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
    ///     Ok(())
    /// }
    /// ```
    ///
    pub fn new(args: SufrBuilderArgs) -> Result<SuffixArray> {
        let low_memory = args.low_memory;
        let path = Self::write(args)?;
        Self::read(&path, low_memory)
    }

    /// Create a new suffix array and write to disk
    /// ```
    /// use anyhow::Result;
    /// use libsufr::{
    ///     types::SufrBuilderArgs,
    ///     suffix_array::SuffixArray,
    ///     util::read_sequence_file,
    /// };
    /// use std::path::Path;
    ///
    /// fn main() -> Result<()> {
    ///     let sequence_delimiter = b'%';
    ///     let seq_data = read_sequence_file(Path::new("../data/inputs/3.fa"), sequence_delimiter)?;
    ///     let builder_args = SufrBuilderArgs {
    ///         text: seq_data.seq,
    ///         path: Some("3.sufr".to_string()),
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

    /// Count instances of queries
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

    pub fn extract(&mut self, args: ExtractOptions) -> Result<Vec<ExtractResult>> {
        self.inner.extract(args)
    }

    pub fn list(&mut self, args: ListOptions) -> Result<()> {
        self.inner.list(args)
    }

    pub fn locate(&mut self, args: LocateOptions) -> Result<Vec<LocateResult>> {
        self.inner.locate(args)
    }

    pub fn metadata(&self) -> Result<SufrMetadata> {
        self.inner.metadata()
    }

    pub fn string_at(&mut self, pos: usize, len: Option<usize>) -> Result<String> {
        self.inner.string_at(pos, len)
    }
}
