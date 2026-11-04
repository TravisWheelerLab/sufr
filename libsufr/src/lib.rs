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
pub use search::{read_queries, Extraction, Hits, Location};
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

use anyhow::Result;
use log::info;
use rayon::prelude::*;
use std::time::Instant;

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

    /// Find the suffixes matching `query`.
    pub fn search(&self, query: &[u8]) -> Hits<'_> {
        let ranks = match &self.inner {
            SuffixArrayKind::U32(sa) => sa.matching_ranks(query),
            SuffixArrayKind::U64(sa) => sa.matching_ranks(query),
        };
        Hits::new(self, ranks)
    }

    /// Search for each of `queries` in parallel, returning their hits in
    /// the same order.
    pub fn search_all<Q: AsRef<[u8]> + Sync>(&self, queries: &[Q]) -> Vec<Hits<'_>> {
        let now = Instant::now();
        let hits = queries
            .par_iter()
            .map(|query| self.search(query.as_ref()))
            .collect();
        info!(
            "Search of {} queries finished in {:?}",
            queries.len(),
            now.elapsed()
        );
        hits
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

    /// The text from `pos`, of at most `len` bytes when given, and
    /// otherwise to the end of the text.
    pub fn string_at(&self, pos: usize, len: Option<usize>) -> Result<String> {
        match &self.inner {
            SuffixArrayKind::U32(sa) => sa.string_at(pos, len),
            SuffixArrayKind::U64(sa) => sa.string_at(pos, len),
        }
    }

    /// Every suffix in rank order, read from the file.
    pub fn suffixes(&self) -> Box<dyn Iterator<Item = usize> + '_> {
        match &self.inner {
            SuffixArrayKind::U32(sa) => {
                Box::new(sa.suffix_array_file.iter().map(|v| v.to_usize()))
            }
            SuffixArrayKind::U64(sa) => {
                Box::new(sa.suffix_array_file.iter().map(|v| v.to_usize()))
            }
        }
    }

    /// The suffix at `rank`, read from the file.
    pub fn suffix_at(&self, rank: usize) -> Option<usize> {
        match &self.inner {
            SuffixArrayKind::U32(sa) => {
                sa.suffix_array_file.get(rank).map(|v| v.to_usize())
            }
            SuffixArrayKind::U64(sa) => {
                sa.suffix_array_file.get(rank).map(|v| v.to_usize())
            }
        }
    }

    /// The LCP of the suffix at `rank` with the one before it, read from
    /// the file, or `None` when the file has no LCP array.
    pub fn lcp_at(&self, rank: usize) -> Option<usize> {
        match &self.inner {
            SuffixArrayKind::U32(sa) => {
                sa.lcp_file.as_ref()?.get(rank).map(|v| v.to_usize())
            }
            SuffixArrayKind::U64(sa) => {
                sa.lcp_file.as_ref()?.get(rank).map(|v| v.to_usize())
            }
        }
    }

    /// The positions of the suffixes at `ranks`, in rank order.
    pub(crate) fn suffixes_at(&self, ranks: Range<usize>) -> Result<Vec<usize>> {
        match &self.inner {
            SuffixArrayKind::U32(sa) => Ok(to_usize(sa.suffixes_at(ranks)?)),
            SuffixArrayKind::U64(sa) => Ok(to_usize(sa.suffixes_at(ranks)?)),
        }
    }

    /// The suffixes at `ranks` within their sequences.
    pub(crate) fn locate(&self, ranks: Range<usize>) -> Result<Vec<Location<'_>>> {
        match &self.inner {
            SuffixArrayKind::U32(sa) => sa.locate(ranks),
            SuffixArrayKind::U64(sa) => sa.locate(ranks),
        }
    }

    /// The suffixes at `ranks` with `before` and `after` bytes of context.
    pub(crate) fn extract(
        &self,
        ranks: Range<usize>,
        before: Option<usize>,
        after: Option<usize>,
    ) -> Result<Vec<Extraction<'_>>> {
        match &self.inner {
            SuffixArrayKind::U32(sa) => sa.extract(ranks, before, after),
            SuffixArrayKind::U64(sa) => sa.extract(ranks, before, after),
        }
    }
}

fn to_usize<T: Int>(values: Vec<T>) -> Vec<usize> {
    values.into_iter().map(|v| v.to_usize()).collect()
}

#[cfg(test)]
mod tests {
    use super::{
        Alphabet, SeedMask, Storage, SuffixArray, SuffixArrayBuilder, SuffixArrayKind,
        SuffixSortType, TypedSuffixArray, OUTFILE_VERSION,
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
            .search_all(&["AC", "GG", "CG"])
            .iter()
            .map(|hits| hits.count())
            .collect();
        assert_eq!(counts, [2, 0, 2]);
        assert_eq!(suffix_array.search(b"CG").ranks(), 3..5);
        Ok(())
    }

    #[test]
    fn test_suffix_array_list() -> Result<()> {
        let suffix_array = SuffixArray::open("../data/inputs/1.sufr")?;
        let suffixes: Vec<_> = suffix_array.suffixes().collect();
        assert_eq!(suffixes, [10, 6, 0, 7, 1, 8, 2, 9, 3]);
        let lcps: Vec<_> = (0..9).map(|rank| suffix_array.lcp_at(rank)).collect();
        let expected = [0, 0, 4, 0, 3, 0, 2, 0, 1].map(Some);
        assert_eq!(lcps, expected);
        assert_eq!(suffix_array.suffix_at(4), Some(1));
        assert_eq!(suffix_array.suffix_at(9), None);
        assert_eq!(suffix_array.lcp_at(9), None);
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
    fn test_suffix_array_string_at() -> Result<()> {
        let suffix_array = SuffixArray::open("../data/inputs/1.sufr")?;
        assert!(matches!(suffix_array.inner, SuffixArrayKind::U32(_)));
        assert_eq!(suffix_array.string_at(0, None)?, "ACGTNNACGT$".to_string());
        assert_eq!(suffix_array.string_at(6, Some(3))?, "ACG".to_string());
        Ok(())
    }
}
