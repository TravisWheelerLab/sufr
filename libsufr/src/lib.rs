#![warn(missing_docs)]

//! # Parallel Construction of Suffix Arrays in Rust
//!
//! Inspired by the merge sort approach described in "Cache-friendly,
//! Parallel, and Samplesort-based Constructor for Suffix Arrays and LCP
//! Arrays." [^patro]
//!
//! [^patro]: <https://doi.org/10.4230/LIPIcs.WABI.2023.16>
//!
//! Use [SuffixArrayBuilder] to create a suffix array and [SuffixArray] to
//! read one. A file holds 32-bit or 64-bit values depending on the length
//! of its text, and [SuffixArray] reads either.
//!
//! ## Authors
//! * Ken Youens-Clark <kyclark@gmail.com>
//! * Jack Roddy <jroddy@arizona.edu>
//! * Travis Wheeler <twheeler@arizona.edu>

mod build;
mod format;
mod io;
mod lcp;
mod search;
mod typed;

pub use build::{Alphabet, SuffixArrayBuilder};
pub use format::{Header, SeedMask, SuffixSortType, OUTFILE_VERSION};
pub use search::{read_queries, Context, Extraction, Hits, Location};
pub use typed::{Metadata, Storage};

/// How [`SuffixArray::open_with`] opens a file. The default reads the text
/// and the suffix array from disk as needed and searches to the length
/// the file was sorted with.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct OpenOptions {
    /// Where the text is read from.
    pub text: Storage,

    /// Where the suffix array is read from.
    pub suffixes: Storage,

    /// The number of leading query bytes a suffix must match to be a hit.
    /// `None` uses the length the file was sorted to, and a longer value
    /// is lowered to it.
    pub max_query_len: Option<usize>,
}

impl Default for OpenOptions {
    fn default() -> Self {
        OpenOptions {
            text: Storage::Disk,
            suffixes: Storage::Disk,
            max_query_len: None,
        }
    }
}

use io::Int;
use std::ops::Range;
use typed::TypedSuffixArray;

use anyhow::{bail, Result};
use rayon::prelude::*;

/// A suffix array read from a _.sufr_ file.
#[derive(Debug)]
pub struct SuffixArray {
    inner: SuffixArrayKind,
}

/// A suffix array at whichever int size its file was written with.
#[derive(Debug)]
pub(crate) enum SuffixArrayKind {
    U32(TypedSuffixArray<u32>),
    U64(TypedSuffixArray<u64>),
}

impl SuffixArray {
    /// Open a _.sufr_ file with the default [`OpenOptions`]: the text and
    /// the suffix array stay on disk.
    pub fn open(filename: &str) -> Result<SuffixArray> {
        SuffixArray::open_with(filename, &OpenOptions::default())
    }

    /// Open a _.sufr_ file.
    pub fn open_with(filename: &str, options: &OpenOptions) -> Result<SuffixArray> {
        let OpenOptions {
            text,
            suffixes,
            max_query_len,
        } = *options;
        let text_len = format::read_text_length(filename)? as u64;
        let inner = if text_len < u32::MAX as u64 {
            SuffixArrayKind::U32(TypedSuffixArray::open(
                filename,
                text,
                suffixes,
                max_query_len,
            )?)
        } else {
            SuffixArrayKind::U64(TypedSuffixArray::open(
                filename,
                text,
                suffixes,
                max_query_len,
            )?)
        };
        Ok(SuffixArray { inner })
    }

    /// Find the suffixes matching `query`. An error is a failed read of
    /// the file.
    pub fn search(&self, query: &[u8]) -> Result<Hits<'_>> {
        let ranks = match &self.inner {
            SuffixArrayKind::U32(sa) => sa.matching_ranks(query)?,
            SuffixArrayKind::U64(sa) => sa.matching_ranks(query)?,
        };
        Ok(Hits::new(self, ranks))
    }

    /// Search each query in turn, stopping at the first failed read.
    pub fn search_all<Q>(
        &self,
        queries: impl IntoIterator<Item = Q>,
    ) -> Result<Vec<Hits<'_>>>
    where
        Q: AsRef<[u8]>,
    {
        queries
            .into_iter()
            .map(|q| self.search(q.as_ref()))
            .collect()
    }

    /// Search the queries on the current rayon pool, in the same order.
    pub fn par_search_all<Q>(
        &self,
        queries: impl IntoIterator<Item = Q>,
    ) -> Result<Vec<Hits<'_>>>
    where
        Q: AsRef<[u8]> + Sync,
    {
        // rayon parallelises slices, not arbitrary iterators
        let queries: Vec<Q> = queries.into_iter().collect();
        queries
            .par_iter()
            .map(|q| self.search(q.as_ref()))
            .collect()
    }

    /// The file's header.
    pub fn header(&self) -> &Header {
        match &self.inner {
            SuffixArrayKind::U32(sa) => &sa.header,
            SuffixArrayKind::U64(sa) => &sa.header,
        }
    }

    /// The file's metadata.
    pub fn metadata(&self) -> Result<Metadata> {
        match &self.inner {
            SuffixArrayKind::U32(sa) => sa.metadata(),
            SuffixArrayKind::U64(sa) => sa.metadata(),
        }
    }

    /// The positions in the text of the suffixes at `ranks`, in rank
    /// order.
    pub fn suffixes(&self, ranks: Range<usize>) -> Result<Vec<usize>> {
        match &self.inner {
            SuffixArrayKind::U32(sa) => Ok(to_usize(sa.suffixes_at(ranks)?)),
            SuffixArrayKind::U64(sa) => Ok(to_usize(sa.suffixes_at(ranks)?)),
        }
    }

    /// The LCP of each suffix at `ranks` with the suffix before it, in
    /// rank order. An error when the file has no LCP array.
    pub fn lcps(&self, ranks: Range<usize>) -> Result<Vec<usize>> {
        match &self.inner {
            SuffixArrayKind::U32(sa) => Ok(to_usize(sa.lcps_at(ranks)?)),
            SuffixArrayKind::U64(sa) => Ok(to_usize(sa.lcps_at(ranks)?)),
        }
    }

    /// The text in `range`, clipped to the end of the text.
    pub fn text(&self, range: Range<usize>) -> Result<Vec<u8>> {
        let end = range.end.min(self.header().text_len);
        if range.start > end {
            bail!("Invalid range: {range:?}");
        }
        match &self.inner {
            SuffixArrayKind::U32(sa) => sa.text_range(range.start..end),
            SuffixArrayKind::U64(sa) => sa.text_range(range.start..end),
        }
    }

    /// The suffixes at `ranks` within their sequences.
    pub(crate) fn locate(&self, ranks: Range<usize>) -> Result<Vec<Location<'_>>> {
        match &self.inner {
            SuffixArrayKind::U32(sa) => sa.locate(ranks),
            SuffixArrayKind::U64(sa) => sa.locate(ranks),
        }
    }

    /// The suffixes at `ranks` with the window of their sequences that
    /// `context` asks for.
    pub(crate) fn extract(
        &self,
        ranks: Range<usize>,
        context: Context,
    ) -> Result<Vec<Extraction<'_>>> {
        match &self.inner {
            SuffixArrayKind::U32(sa) => sa.extract(ranks, context),
            SuffixArrayKind::U64(sa) => sa.extract(ranks, context),
        }
    }
}

fn to_usize<T: Int>(values: Vec<T>) -> Vec<usize> {
    values.into_iter().map(|v| v.to_usize()).collect()
}

#[cfg(test)]
mod tests {
    use super::{
        Alphabet, Hits, SeedMask, Storage, SuffixArray, SuffixArrayBuilder,
        SuffixArrayKind, SuffixSortType, TypedSuffixArray, OUTFILE_VERSION,
    };
    use anyhow::Result;
    use std::path::Path;
    use tempfile::NamedTempFile;

    #[test]
    fn test_write_read_suffix_file_32() -> Result<()> {
        let seq_file = Path::new("../data/inputs/2.fa");
        let outfile = NamedTempFile::new()?;
        let outpath = outfile.path().to_string_lossy().to_string();

        SuffixArrayBuilder::from_sequence_file(seq_file, b'N')?
            .with_alphabet(Alphabet::Dna {
                allow_ambiguity: false,
                ignore_softmask: false,
            })
            .with_num_partitions(2)
            .with_write_lcp(true)
            .build_to_file(&outpath, true)?;
        assert!(outfile.path().exists());

        let index: TypedSuffixArray<u32> =
            TypedSuffixArray::open(&outpath, Storage::Memory, Storage::Disk, None)?;
        assert_eq!(index.header.version, OUTFILE_VERSION);
        assert!(index.header.is_dna);
        assert_eq!(index.header.text_len, 18);
        assert_eq!(index.header.sequence_starts, [0, 9]);
        assert_eq!(index.header.sequence_names, ["ABC", "DEF"]);
        assert_eq!(index.text_range(0..18)?, b"ACGTACGTNACGTACGT$");

        let file_sa: Vec<_> = index.suffix_array_file.iter().collect();
        let sorted_sa = [17, 13, 9, 0, 4, 14, 10, 1, 5, 15, 11, 2, 6, 16, 12, 3, 7];
        assert_eq!(file_sa, sorted_sa);

        let file_lcp: Vec<_> = index.lcp_file.as_ref().unwrap().iter().collect();
        let lcp = [0, 0, 4, 8, 4, 0, 3, 7, 3, 0, 2, 6, 2, 0, 1, 5, 1];
        assert_eq!(file_lcp, lcp);
        Ok(())
    }

    #[test]
    fn test_write_read_suffix_file_64() -> Result<()> {
        let seq_file = Path::new("../data/inputs/1.fa");
        let outfile = NamedTempFile::new()?;
        let outpath = outfile.path().to_string_lossy().to_string();

        SuffixArrayBuilder::from_sequence_file(seq_file, b'N')?
            .with_alphabet(Alphabet::Dna {
                allow_ambiguity: true,
                ignore_softmask: false,
            })
            .with_num_partitions(2)
            .with_force_64bit(true)
            .with_write_lcp(true)
            .build_to_file(&outpath, true)?;
        assert!(outfile.path().exists());

        let index: TypedSuffixArray<u64> =
            TypedSuffixArray::open(&outpath, Storage::Memory, Storage::Disk, None)?;
        assert_eq!(index.header.version, OUTFILE_VERSION);
        assert!(index.header.is_dna);
        assert_eq!(index.header.text_len, 11);
        assert_eq!(index.header.sequence_starts, [0]);
        assert_eq!(index.header.sequence_names, ["1"]);
        assert_eq!(index.text_range(0..11)?, b"ACGTNNACGT$");
        assert_eq!(index.header.num_suffixes, 11);

        let file_sa: Vec<_> = index.suffix_array_file.iter().collect();
        assert_eq!(file_sa, &[10, 6, 0, 7, 1, 8, 2, 5, 4, 9, 3]);
        let file_lcp: Vec<_> = index.lcp_file.as_ref().unwrap().iter().collect();
        assert_eq!(file_lcp, &[0, 0, 4, 0, 3, 0, 2, 0, 1, 0, 1]);
        Ok(())
    }

    #[test]
    fn test_subsample_suffix_array() -> Result<()> {
        let seq_file = Path::new("../data/inputs/smol.fa");
        let outfile = NamedTempFile::new()?;
        let outpath = outfile.path().to_string_lossy().to_string();

        SuffixArrayBuilder::from_sequence_file(seq_file, b'N')?
            .with_alphabet(Alphabet::Dna {
                allow_ambiguity: false,
                ignore_softmask: false,
            })
            .with_num_partitions(2)
            .with_write_lcp(true)
            .build_to_file(&outpath, true)?;
        assert!(outfile.path().exists());

        let index: TypedSuffixArray<u32> =
            TypedSuffixArray::open(&outpath, Storage::Memory, Storage::Disk, None)?;
        assert_eq!(index.header.num_suffixes, 364);
        let full_sa: Vec<_> = index.suffix_array_file.iter().collect();
        let full_lcp: Vec<_> = index.lcp_file.as_ref().unwrap().iter().collect();
        assert_eq!(full_sa.len(), 364);
        assert_eq!(full_lcp.len(), 364);

        let max_query_len = 1;
        let (sub_sa, sub_rank) = index.subsample(max_query_len);
        assert_eq!(sub_sa.len(), 5);
        assert_eq!(sub_rank.len(), 5);
        // $, A, C, G, T
        assert_eq!(sub_sa, vec![365, 364, 92, 224, 363]);
        assert_eq!(sub_rank, vec![0, 1, 94, 191, 284]);

        let max_query_len = 2;
        let (sub_sa, sub_rank) = index.subsample(max_query_len);
        assert_eq!(sub_sa.len(), 20);
        assert_eq!(sub_rank.len(), 20);
        // $, A$, AA, AC, AG, AN, AT, CA, CC, CG, CT, GA, GC, GG, GN, GT, TA, TC, TG, TT
        assert_eq!(
            sub_sa,
            vec![
                365, 364, 358, 91, 341, 255, 362, 92, 339, 233, 296, 224, 88, 129, 110,
                96, 363, 217, 223, 356
            ]
        );
        assert_eq!(
            sub_rank,
            vec![
                0, 1, 2, 38, 49, 70, 71, 94, 112, 143, 170, 191, 216, 252, 269, 270,
                284, 298, 315, 343
            ]
        );

        let max_query_len = 3;
        let (sub_sa, sub_rank) = index.subsample(max_query_len);
        assert_eq!(sub_sa.len(), 71);
        assert_eq!(sub_rank.len(), 71);

        let max_query_len = 5;
        let (sub_sa, sub_rank) = index.subsample(max_query_len);
        assert_eq!(sub_sa.len(), 293);
        assert_eq!(sub_rank.len(), 293);

        let max_query_len = 3;
        let (sub_sa, sub_rank) = index.subsample(max_query_len);
        assert_eq!(sub_sa.len(), 71);
        assert_eq!(sub_rank.len(), 71);

        Ok(())
    }

    #[test]
    fn test_spaced_seeds_1() -> Result<()> {
        let seq_file = Path::new("../data/inputs/mostlya1.fa");
        let outfile = NamedTempFile::new()?;
        let outpath = outfile.path().to_string_lossy().to_string();

        // 7 $
        // 6 A$
        // 5 AA$
        // 4 AAA$
        // 2 ACAAA$
        // 0 AAACAAA$
        // 1 AACAAA$
        // 3 CAAA$

        SuffixArrayBuilder::from_sequence_file(seq_file, b'N')?
            .with_alphabet(Alphabet::Dna {
                allow_ambiguity: false,
                ignore_softmask: false,
            })
            .with_num_partitions(1)
            .with_sort_type(SuffixSortType::Mask(SeedMask::new("101").unwrap()))
            .with_write_lcp(true)
            .build_to_file(&outpath, true)?;
        assert!(outfile.path().exists());

        let index: TypedSuffixArray<u32> =
            TypedSuffixArray::open(&outpath, Storage::Memory, Storage::Disk, None)?;
        assert_eq!(index.header.num_suffixes, 8);
        let suffix_array: Vec<_> = index.suffix_array_file.iter().collect();
        assert_eq!(suffix_array.len(), 8);
        assert_eq!(suffix_array, vec![7, 6, 5, 4, 2, 0, 1, 3]);

        Ok(())
    }

    #[test]
    fn test_spaced_seeds_2() -> Result<()> {
        let seq_file = Path::new("../data/inputs/mostlya2.fa");
        let outfile = NamedTempFile::new()?;
        let outpath = outfile.path().to_string_lossy().to_string();

        //  0 16 $
        //  1 13 AAC$
        //  2  5 AACAAACAAAC$
        //  3  1 AACAAACAAACAAAC$
        //  4  9 AACAAAC$
        //  5 12 AAAC$
        //  6  8 AAACAAAC$
        //  7  4 AAACAAACAAAC$
        //  8  0 AAACAAACAAACAAAC$
        //  9 14 AC$
        // 10 10 ACAAAC$
        // 11  6 ACAAACAAAC$
        // 12  2 ACAAACAAACAAAC$
        // 13 15 C$
        // 14 11 CAAAC$
        // 15  7 CAAACAAAC$
        // 16  3 CAAACAAACAAAC$

        SuffixArrayBuilder::from_sequence_file(seq_file, b'N')?
            .with_alphabet(Alphabet::Dna {
                allow_ambiguity: false,
                ignore_softmask: false,
            })
            .with_num_partitions(1)
            .with_sort_type(SuffixSortType::Mask(SeedMask::new("11011").unwrap()))
            .with_write_lcp(true)
            .build_to_file(&outpath, true)?;
        assert!(outfile.path().exists());

        let index = TypedSuffixArray::<u32>::open(
            &outpath,
            Storage::Memory,
            Storage::Disk,
            None,
        )?;
        assert_eq!(index.header.num_suffixes, 17);
        let suffix_array: Vec<_> = index.suffix_array_file.iter().collect();
        assert_eq!(suffix_array.len(), 17);
        assert_eq!(
            suffix_array,
            vec![16, 13, 9, 5, 1, 12, 8, 4, 0, 14, 10, 6, 2, 15, 11, 7, 3],
        );

        Ok(())
    }

    #[test]
    fn test_spaced_seeds_3() -> Result<()> {
        let seq_file = Path::new("../data/inputs/spaced_input.fa");
        let outfile = NamedTempFile::new()?;
        let outpath = outfile.path().to_string_lossy().to_string();

        SuffixArrayBuilder::from_sequence_file(seq_file, b'N')?
            .with_alphabet(Alphabet::Dna {
                allow_ambiguity: false,
                ignore_softmask: false,
            })
            .with_num_partitions(1)
            .with_sort_type(SuffixSortType::Mask(SeedMask::new("11000111").unwrap()))
            .with_write_lcp(true)
            .build_to_file(&outpath, true)?;
        assert!(outfile.path().exists());

        let index: TypedSuffixArray<u32> =
            TypedSuffixArray::open(&outpath, Storage::Memory, Storage::Disk, None)?;
        assert_eq!(index.header.num_suffixes, 43);
        let suffix_array: Vec<_> = index.suffix_array_file.iter().collect();
        assert_eq!(suffix_array.len(), 43);
        assert_eq!(
            suffix_array,
            vec![
                42, 18, 12, 0, 32, 29, 13, 23, 21, 6, 40, 1, 33, 19, 30, 10, 28, 9, 17,
                14, 4, 26, 39, 22, 25, 38, 24, 35, 7, 36, 15, 41, 5, 20, 31, 11, 27, 8,
                16, 3, 37, 34, 2
            ],
        );

        Ok(())
    }

    #[test]
    fn test_suffix_array_search() -> Result<()> {
        let suffix_array = SuffixArray::open("../data/inputs/1.sufr")?;
        let counts: Vec<_> = suffix_array
            .search_all(["AC", "GG", "CG"])?
            .iter()
            .map(Hits::count)
            .collect();
        assert_eq!(counts, [2, 0, 2]);
        let counts: Vec<_> = suffix_array
            .par_search_all(["AC", "GG", "CG"])?
            .iter()
            .map(Hits::count)
            .collect();
        assert_eq!(counts, [2, 0, 2]);
        assert_eq!(suffix_array.search(b"CG")?.ranks(), 3..5);
        Ok(())
    }

    // what: a compile-time check that SuffixArray is Send + Sync
    // why: nothing declares either; the compiler infers both from
    //      the fields, down through TypedSuffixArray, Text,
    //      Suffixes and FileAccess. par_search_all shares &self
    //      across the rayon pool, and a caller with their own
    //      threads is promised thread::scope with search or
    //      thread::spawn with the array moved in. A new field
    //      holding a Cell, RefCell or Rc would silently remove
    //      the trait, and the first error would come out of
    //      rayon's bounds on par_iter, or out of a downstream
    //      crate's build. This line fails instead, naming the
    //      field that lost it
    #[test]
    fn suffix_array_is_send_and_sync() {
        fn assert_send_sync<T: Send + Sync>() {}
        assert_send_sync::<SuffixArray>();
    }

    #[test]
    fn test_suffix_array_list() -> Result<()> {
        let suffix_array = SuffixArray::open("../data/inputs/1.sufr")?;
        assert_eq!(suffix_array.suffixes(0..9)?, [10, 6, 0, 7, 1, 8, 2, 9, 3]);
        assert_eq!(suffix_array.lcps(0..9)?, [0, 0, 4, 0, 3, 0, 2, 0, 1]);
        assert_eq!(suffix_array.suffixes(4..5)?, [1]);
        assert!(suffix_array.suffixes(9..10).is_err());
        assert!(suffix_array.lcps(9..10).is_err());
        Ok(())
    }

    #[test]
    fn test_suffix_array_metadata() -> Result<()> {
        let suffix_array = SuffixArray::open("../data/inputs/1.sufr")?;
        let meta = suffix_array.metadata()?;
        assert_eq!(meta.filename, "../data/inputs/1.sufr".to_string());
        assert_eq!(meta.file_size, 172);
        assert_eq!(meta.file_version, 6);
        assert!(meta.is_dna);
        assert!(!meta.allow_ambiguity);
        assert!(!meta.ignore_softmask);
        assert_eq!(meta.text_len, 11);
        assert_eq!(meta.num_suffixes, 9);
        assert_eq!(meta.num_sequences, 1);
        assert_eq!(meta.sequence_starts, vec![0]);
        assert_eq!(meta.sequence_names, vec!["1".to_string()]);
        assert_eq!(meta.sort_type, SuffixSortType::MaxQueryLen(0));
        assert!(meta.has_lcp);

        let summary = meta.summary();
        let labels: Vec<_> = summary.iter().map(|(label, _)| *label).collect();
        assert_eq!(
            labels,
            [
                "Filename",
                "Modified",
                "File Size",
                "File Version",
                "DNA",
                "Allow Ambiguity",
                "Ignore Softmask",
                "Text Length",
                "Len Suffixes",
                "LCP array",
                "Max query len",
                "Num sequences",
                "Sequence starts",
                "Sequence names",
            ]
        );
        assert_eq!(summary[2].1, "172 bytes");
        assert_eq!(summary[9].1, "yes");
        assert_eq!(summary[10].1, "0");
        Ok(())
    }

    #[test]
    fn test_suffix_array_text() -> Result<()> {
        let suffix_array = SuffixArray::open("../data/inputs/1.sufr")?;
        assert!(matches!(suffix_array.inner, SuffixArrayKind::U32(_)));
        assert_eq!(suffix_array.text(0..11)?, b"ACGTNNACGT$");
        assert_eq!(suffix_array.text(6..9)?, b"ACG");
        assert_eq!(suffix_array.text(8..20)?, b"GT$");
        assert!(suffix_array.text(12..13).is_err());
        Ok(())
    }
}
