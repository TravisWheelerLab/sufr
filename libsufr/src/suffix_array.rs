use crate::{
    sufr_builder::{SufrBuilder, SufrBuilderArgs},
    sufr_file::SufrFile,
    types::{
        CountOptions, CountResult, ExtractOptions, ExtractResult, ListOptions,
        LocateResult, LocateOptions, SufrMetadata
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
pub struct SuffixArray {
    inner: Box<dyn SuffixArrayTrait>,
}

impl SuffixArray {
    pub fn new(args: SufrBuilderArgs) -> Result<SuffixArray> {
        let low_memory = args.low_memory;
        let path = Self::write(args)?;
        Self::read(&path, low_memory)
    }

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
