//! Create and query suffix arrays
use crate::{
    sufr_builder::SufrBuilder,
    sufr_file::SufrFile,
    types::{
        Context, Extraction, Int, Location, OpenOptions, SufrBuilderArgs, SufrMetadata,
    },
};
use anyhow::Result;
use rayon::prelude::*;
use std::{cmp::min, ops::Range};

// --------------------------------------------------
pub(crate) trait SuffixArrayTrait: Send + Sync {
    fn search(&self, query: &[u8]) -> Result<Range<usize>>;
    fn suffixes(&self, ranks: Range<usize>) -> Result<Vec<usize>>;
    fn lcps(&self, ranks: Range<usize>) -> Result<Vec<usize>>;
    fn text(&self, range: Range<usize>) -> Result<Vec<u8>>;
    fn metadata(&self) -> Result<SufrMetadata>;
}

impl<T: Int> SuffixArrayTrait for SufrFile<T> {
    fn search(&self, query: &[u8]) -> Result<Range<usize>> {
        SufrFile::search(self, query)
    }

    fn suffixes(&self, ranks: Range<usize>) -> Result<Vec<usize>> {
        SufrFile::suffixes(self, ranks)
    }

    fn lcps(&self, ranks: Range<usize>) -> Result<Vec<usize>> {
        SufrFile::lcps(self, ranks)
    }

    fn text(&self, range: Range<usize>) -> Result<Vec<u8>> {
        self.get_text_range(range)
    }

    fn metadata(&self) -> Result<SufrMetadata> {
        SufrFile::metadata(self)
    }
}

// --------------------------------------------------
/// Struct to create and read suffix arrays
pub struct SuffixArray {
    inner: Box<dyn SuffixArrayTrait>,
    text_len: usize,
    sequence_starts: Vec<usize>,
    sequence_names: Vec<String>,
}

// --------------------------------------------------
impl SuffixArray {
    // --------------------------------------------------
    /// Create a new suffix array and return an in-memory struct.
    /// This is a combination of `write` and `open`.
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
    ///         write_lcp: true,
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
        let path = Self::write(args)?;
        Self::open(&path)
    }

    // --------------------------------------------------
    /// Open a _.sufr_ file with the text and suffix array on disk
    /// and the maximum query length it was built with
    ///
    /// ```
    /// use anyhow::Result;
    /// use libsufr::suffix_array::SuffixArray;
    ///
    /// fn main() -> Result<()> {
    ///     let suffix_array = SuffixArray::open("../data/inputs/1.sufr")?;
    ///     Ok(())
    /// }
    /// ```
    pub fn open(filename: &str) -> Result<SuffixArray> {
        Self::open_with(filename, &OpenOptions::default())
    }

    // --------------------------------------------------
    /// Open a _.sufr_ file
    ///
    /// ```
    /// use anyhow::Result;
    /// use libsufr::{suffix_array::SuffixArray, types::{Load, OpenOptions}};
    ///
    /// fn main() -> Result<()> {
    ///     let opts = OpenOptions {
    ///         text: Load::Memory,
    ///         suffix_array: Load::Memory,
    ///         max_query_len: Some(2),
    ///     };
    ///     let suffix_array = SuffixArray::open_with("../data/inputs/1.sufr", &opts)?;
    ///     assert_eq!(suffix_array.search(b"ACX")?.count(), 2);
    ///     Ok(())
    /// }
    /// ```
    pub fn open_with(filename: &str, opts: &OpenOptions) -> Result<SuffixArray> {
        let text_len = crate::util::read_text_length(filename)? as u64;
        let inner: Box<dyn SuffixArrayTrait> = if text_len < u32::MAX as u64 {
            Box::new(SufrFile::<u32>::open(filename, opts)?)
        } else {
            Box::new(SufrFile::<u64>::open(filename, opts)?)
        };
        let meta = inner.metadata()?;
        Ok(SuffixArray {
            inner,
            text_len: meta.text_len,
            sequence_starts: meta.sequence_starts,
            sequence_names: meta.sequence_names,
        })
    }

    // --------------------------------------------------
    /// Find the suffixes matching a query
    ///
    /// Given the text "ACGTNNACGT":
    ///
    /// ```
    /// use anyhow::Result;
    /// use libsufr::suffix_array::SuffixArray;
    ///
    /// fn main() -> Result<()> {
    ///     let suffix_array = SuffixArray::open("../data/inputs/1.sufr")?;
    ///     let hits = suffix_array.search(b"ACG")?;
    ///     assert_eq!(hits.count(), 2);
    ///     assert_eq!(hits.ranks(), 1..3);
    ///     assert_eq!(hits.suffixes()?, vec![6, 0]);
    ///     assert!(suffix_array.search(b"GGC")?.is_empty());
    ///     Ok(())
    /// }
    /// ```
    pub fn search(&self, query: &[u8]) -> Result<Hits<'_>> {
        Ok(Hits {
            suffix_array: self,
            ranks: self.inner.search(query)?,
        })
    }

    // --------------------------------------------------
    /// Search for many queries in parallel.
    /// The results are in the order of the queries.
    pub fn search_many<Q>(&self, queries: &[Q]) -> Vec<Result<Hits<'_>>>
    where
        Q: AsRef<[u8]> + Sync,
    {
        queries.par_iter().map(|q| self.search(q.as_ref())).collect()
    }

    // --------------------------------------------------
    /// The positions in the text of the suffixes at a range of ranks
    pub fn suffixes(&self, ranks: Range<usize>) -> Result<Vec<usize>> {
        self.inner.suffixes(ranks)
    }

    // --------------------------------------------------
    /// The LCP values at a range of ranks
    pub fn lcps(&self, ranks: Range<usize>) -> Result<Vec<usize>> {
        self.inner.lcps(ranks)
    }

    // --------------------------------------------------
    /// The text in a range of positions, clipped to the end of the text
    ///
    /// ```
    /// use anyhow::Result;
    /// use libsufr::suffix_array::SuffixArray;
    ///
    /// fn main() -> Result<()> {
    ///     let suffix_array = SuffixArray::open("../data/inputs/1.sufr")?;
    ///     assert_eq!(suffix_array.text(0..100)?, b"ACGTNNACGT$");
    ///     assert_eq!(suffix_array.text(6..9)?, b"ACG");
    ///     Ok(())
    /// }
    /// ```
    pub fn text(&self, range: Range<usize>) -> Result<Vec<u8>> {
        let end = min(range.end, self.text_len);
        self.inner.text(min(range.start, end)..end)
    }

    // --------------------------------------------------
    /// Read metadata from a _.sufr_ file
    ///
    /// ```
    /// use anyhow::Result;
    /// use libsufr::{suffix_array::SuffixArray, types::SuffixSortType};
    ///
    /// fn main() -> Result<()> {
    ///     let suffix_array = SuffixArray::open("../data/inputs/1.sufr")?;
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
    pub fn metadata(&self) -> Result<SufrMetadata> {
        self.inner.metadata()
    }

    // --------------------------------------------------
    /// The index of the sequence containing a text position
    fn sequence_index(&self, pos: usize) -> usize {
        self.sequence_starts.partition_point(|&start| start <= pos) - 1
    }

    // --------------------------------------------------
    /// The end of a sequence in the text
    fn sequence_end(&self, index: usize) -> usize {
        self.sequence_starts
            .get(index + 1)
            .copied()
            .unwrap_or(self.text_len)
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
    ///         write_lcp: true,
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

// --------------------------------------------------
/// The suffixes matching one query, as a range of ranks
/// in the suffix array they came from
#[derive(Clone)]
pub struct Hits<'a> {
    suffix_array: &'a SuffixArray,
    ranks: Range<usize>,
}

impl std::fmt::Debug for Hits<'_> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Hits").field("ranks", &self.ranks).finish()
    }
}

// --------------------------------------------------
impl<'a> Hits<'a> {
    /// The number of matching suffixes
    pub fn count(&self) -> usize {
        self.ranks.len()
    }

    /// Whether no suffixes matched
    pub fn is_empty(&self) -> bool {
        self.ranks.is_empty()
    }

    /// The ranks of the matching suffixes
    pub fn ranks(&self) -> Range<usize> {
        self.ranks.clone()
    }

    // --------------------------------------------------
    /// The positions in the text of the matching suffixes, in rank order
    pub fn suffixes(&self) -> Result<Vec<usize>> {
        if self.ranks.is_empty() {
            return Ok(vec![]);
        }
        self.suffix_array.suffixes(self.ranks())
    }

    // --------------------------------------------------
    /// The matching suffixes in the context of their sequences, in rank order
    ///
    /// Given the text "ACGTNNACGT":
    ///
    /// ```
    /// use anyhow::Result;
    /// use libsufr::suffix_array::SuffixArray;
    ///
    /// fn main() -> Result<()> {
    ///     let suffix_array = SuffixArray::open("../data/inputs/1.sufr")?;
    ///     let locs = suffix_array.search(b"ACG")?.locate()?;
    ///     let found: Vec<_> = locs
    ///         .iter()
    ///         .map(|loc| (loc.sequence_name, loc.sequence_position))
    ///         .collect();
    ///     assert_eq!(found, vec![("1", 6), ("1", 0)]);
    ///     Ok(())
    /// }
    /// ```
    pub fn locate(&self) -> Result<Vec<Location<'a>>> {
        let sa = self.suffix_array;
        Ok(self
            .ranks
            .clone()
            .zip(self.suffixes()?)
            .map(|(rank, suffix)| {
                let i = sa.sequence_index(suffix);
                let sequence_start = sa.sequence_starts[i];
                Location {
                    rank,
                    suffix,
                    sequence_name: &sa.sequence_names[i],
                    sequence_start,
                    sequence_position: suffix - sequence_start,
                }
            })
            .collect())
    }

    // --------------------------------------------------
    /// The matching suffixes with a window of their sequences, in rank order
    ///
    /// Given the text "ACGTNNACGT":
    ///
    /// ```
    /// use anyhow::Result;
    /// use libsufr::{suffix_array::SuffixArray, types::Context};
    ///
    /// fn main() -> Result<()> {
    ///     let suffix_array = SuffixArray::open("../data/inputs/1.sufr")?;
    ///     let context = Context { before: 1, after: Some(3) };
    ///     let found: Vec<_> = suffix_array
    ///         .search(b"ACG")?
    ///         .extract(context)?
    ///         .into_iter()
    ///         .map(|ex| (ex.sequence_range, ex.suffix_offset))
    ///         .collect();
    ///     assert_eq!(found, vec![(5..9, 1), (0..3, 0)]);
    ///     Ok(())
    /// }
    /// ```
    pub fn extract(&self, context: Context) -> Result<Vec<Extraction<'a>>> {
        let sa = self.suffix_array;
        Ok(self
            .locate()?
            .into_iter()
            .map(|location| {
                let i = sa.sequence_index(location.suffix);
                let sequence_end = sa.sequence_end(i);
                let start = location.sequence_position.saturating_sub(context.before);
                let end = context.after.map_or(sequence_end, |after| {
                    min(location.sequence_position + after, sequence_end)
                });
                Extraction {
                    suffix_offset: location.sequence_position - start,
                    sequence_range: start..end,
                    location,
                }
            })
            .collect())
    }
}

// --------------------------------------------------
#[cfg(test)]
mod test {
    use super::SuffixArray;
    use crate::types::{Context, Extraction, Load, Location, OpenOptions};
    use anyhow::Result;

    // --------------------------------------------------
    #[test]
    fn test_extract() -> Result<()> {
        // cargo run -- ex data/expected/1.sufr AC GT XX
        let opts = OpenOptions {
            text: Load::Memory,
            suffix_array: Load::Disk,
            max_query_len: None,
        };
        let sa = SuffixArray::open_with("../data/expected/1.sufr", &opts)?;
        let context = || Context {
            before: 1,
            after: Some(3),
        };

        let res = sa.search(b"AC")?.extract(context())?;
        assert_eq!(
            res,
            vec![
                Extraction {
                    location: Location {
                        rank: 1,
                        suffix: 6,
                        sequence_name: "1",
                        sequence_start: 0,
                        sequence_position: 6,
                    },
                    sequence_range: 5..9,
                    suffix_offset: 1,
                },
                Extraction {
                    location: Location {
                        rank: 2,
                        suffix: 0,
                        sequence_name: "1",
                        sequence_start: 0,
                        sequence_position: 0,
                    },
                    sequence_range: 0..3,
                    suffix_offset: 0,
                },
            ]
        );

        let res = sa.search(b"GT")?.extract(context())?;
        assert_eq!(
            res,
            vec![
                Extraction {
                    location: Location {
                        rank: 5,
                        suffix: 8,
                        sequence_name: "1",
                        sequence_start: 0,
                        sequence_position: 8,
                    },
                    sequence_range: 7..11,
                    suffix_offset: 1,
                },
                Extraction {
                    location: Location {
                        rank: 6,
                        suffix: 2,
                        sequence_name: "1",
                        sequence_start: 0,
                        sequence_position: 2,
                    },
                    sequence_range: 1..5,
                    suffix_offset: 1,
                },
            ]
        );

        let res = sa.search(b"XX")?.extract(context())?;
        assert_eq!(res, vec![]);

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

        let loc = |rank: usize, suffix: usize| Location {
            rank,
            suffix,
            sequence_name: "1",
            sequence_start: 0,
            sequence_position: suffix,
        };

        let cases: Vec<(&[u8], Vec<(usize, usize)>)> = vec![
            (
                b"A",
                vec![(1, 0), (2, 12), (3, 10), (4, 1), (5, 3), (6, 5), (7, 7)],
            ),
            (
                b"B",
                vec![(8, 13), (9, 11), (10, 9), (11, 2), (12, 4), (13, 6), (14, 8)],
            ),
            (b"ABAB", vec![(3, 10), (4, 1), (5, 3), (6, 5)]),
            (b"ABABB", vec![(6, 5)]),
            (b"BBBB", vec![]),
        ];

        for suffix_array in [Load::Disk, Load::Memory] {
            let opts = OpenOptions {
                text: Load::Memory,
                suffix_array,
                max_query_len: None,
            };
            let sa = SuffixArray::open_with("../data/expected/abba.sufr", &opts)?;
            for (query, expected) in &cases {
                let res = sa.search(query)?.locate()?;
                let expected: Vec<_> =
                    expected.iter().map(|&(r, s)| loc(r, s)).collect();
                assert_eq!(res, expected);
            }
        }

        Ok(())
    }
}
