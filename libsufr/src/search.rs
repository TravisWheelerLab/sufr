//! Search a suffix array

use crate::{
    io::Int,
    lcp::{lcp, mask_offset, masked_lcp},
    typed::TypedSuffixArray,
    SuffixArray,
};
use anyhow::Result;
use regex::Regex;
use std::{
    cmp::{min, Ordering},
    fs,
    ops::Range,
    path::Path,
};

/// The suffixes matching one query, as a range of ranks in the
/// suffix array.
#[derive(Debug)]
pub struct Hits<'a> {
    suffix_array: &'a SuffixArray,
    ranks: Range<usize>,
}

/// Where a hit falls in the sequences of the text.
#[derive(Debug, PartialEq)]
pub struct Location<'a> {
    /// The rank of the suffix in the suffix array.
    pub rank: usize,

    /// The position of the hit in the text.
    pub suffix: usize,

    /// The name of the sequence holding the hit.
    pub sequence_name: &'a str,

    /// The start of that sequence in the text.
    pub sequence_start: usize,

    /// The position of the hit in its sequence.
    pub position: usize,
}

/// How much of a hit's sequence [`Hits::extract`] takes around it.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Context {
    /// Bytes before the hit, fewer at the start of the sequence.
    pub before: usize,

    /// Bytes from the start of the hit, or to the end of the sequence
    /// when `None`; fewer at the end of the sequence.
    pub after: Option<usize>,
}

/// A hit with the window of its sequence around it.
#[derive(Debug, PartialEq)]
pub struct Extraction<'a> {
    /// The hit.
    pub location: Location<'a>,

    /// The window, relative to the start of the sequence. Its bytes are
    /// `text(sequence_start + window)`.
    pub sequence_range: Range<usize>,

    /// The position of the hit in the window.
    pub suffix_offset: usize,
}

impl<'a> Hits<'a> {
    pub(crate) fn new(suffix_array: &'a SuffixArray, ranks: Range<usize>) -> Self {
        Hits {
            suffix_array,
            ranks,
        }
    }

    /// The number of hits.
    pub fn count(&self) -> usize {
        self.ranks.len()
    }

    /// Whether nothing matched.
    pub fn is_empty(&self) -> bool {
        self.ranks.is_empty()
    }

    /// The ranks of the hits in the suffix array, empty when nothing
    /// matched.
    pub fn ranks(&self) -> Range<usize> {
        self.ranks.clone()
    }

    /// The positions of the hits in the text, in rank order.
    pub fn suffixes(&self) -> Result<Vec<usize>> {
        self.suffix_array.suffixes(self.ranks())
    }

    /// The hits within their sequences, sorted by sequence name and then
    /// by position.
    pub fn locate(&self) -> Result<Vec<Location<'a>>> {
        self.suffix_array.locate(self.ranks())
    }

    /// The hits in rank order, each with the window of its sequence that
    /// `context` asks for.
    pub fn extract(&self, context: Context) -> Result<Vec<Extraction<'a>>> {
        self.suffix_array.extract(self.ranks(), context)
    }
}

/// The queries, with each one that names an existing file replaced by
/// the whitespace-separated queries in that file.
pub fn read_queries(queries: &[String]) -> Result<Vec<String>> {
    let whitespace = Regex::new(r"\s+").unwrap();
    let mut ret = vec![];
    for query in queries {
        if Path::new(&query).exists() {
            let contents = fs::read_to_string(query)?;
            let mut vals: Vec<String> = whitespace
                .split(&contents)
                .filter(|v| !v.is_empty())
                .map(|v| v.to_string())
                .collect();
            ret.append(&mut vals);
        } else {
            ret.push(query.to_string());
        }
    }

    Ok(ret)
}

/// How a query compares to a suffix.
#[derive(Debug)]
struct Comparison {
    /// Whether the query sorts before, after, or with the suffix.
    cmp: Ordering,

    /// The LCP of the query and the suffix.
    lcp: usize,
}

impl<T: Int> TypedSuffixArray<T> {
    /// The ranks of the suffixes matching `query`, empty when none match.
    pub(crate) fn matching_ranks(&self, query: &[u8]) -> Result<Range<usize>> {
        let n = self.num_searched();
        let Some(start) = self.first_match(query, 0, n - 1, 0, 0)? else {
            return Ok(0..0);
        };
        let end = self
            .last_match(query, start, n - 1, n, 0, 0)?
            .unwrap_or(start);
        Ok(self.searched_ranks(start, end))
    }

    fn first_match(
        &self,
        qry: &[u8],
        low: usize,
        high: usize,
        left_lcp: usize,
        right_lcp: usize,
    ) -> Result<Option<usize>> {
        if high >= low {
            let mid = low + ((high - low) / 2);
            let Some(mid_val) = self.searched_suffix(mid)? else {
                return Ok(None);
            };
            let mid_val = mid_val.to_usize();
            let mid_cmp = self.compare(qry, mid_val, min(left_lcp, right_lcp))?;

            let mid_minus_one = if mid > 0 {
                match self.searched_suffix(mid - 1)? {
                    Some(suffix) => suffix.to_usize(),
                    None => return Ok(None),
                }
            } else {
                mid_val
            };

            if mid_cmp.cmp == Ordering::Equal
                && (mid == 0
                    || self.compare(qry, mid_minus_one, 0)?.cmp == Ordering::Greater)
            {
                Ok(Some(mid))
            } else if mid_cmp.cmp == Ordering::Greater {
                self.first_match(qry, mid + 1, high, mid_cmp.lcp, right_lcp)
            } else {
                // Ordering::Less
                self.first_match(qry, low, mid - 1, left_lcp, mid_cmp.lcp)
            }
        } else {
            Ok(None)
        }
    }

    fn last_match(
        &self,
        qry: &[u8],
        low: usize,
        high: usize,
        n: usize,
        left_lcp: usize,
        right_lcp: usize,
    ) -> Result<Option<usize>> {
        if high >= low {
            let mid = low + ((high - low) / 2);
            let Some(mid_val) = self.searched_suffix(mid)? else {
                return Ok(None);
            };
            let mid_val = mid_val.to_usize();
            let mid_cmp = self.compare(qry, mid_val, min(left_lcp, right_lcp))?;

            let mid_plus_one = if mid < n - 1 {
                match self.searched_suffix(mid + 1)? {
                    Some(suffix) => suffix.to_usize(),
                    None => return Ok(None),
                }
            } else {
                mid_val
            };

            if mid_cmp.cmp == Ordering::Equal
                && (mid == n - 1
                    || self.compare(qry, mid_plus_one, 0)?.cmp == Ordering::Less)
            {
                Ok(Some(mid))
            } else if mid_cmp.cmp == Ordering::Less {
                self.last_match(qry, low, mid - 1, n, left_lcp, mid_cmp.lcp)
            } else {
                self.last_match(qry, mid + 1, high, n, mid_cmp.lcp, right_lcp)
            }
        } else {
            Ok(None)
        }
    }

    /// How `query` compares to the suffix at `suffix_pos`, which share
    /// their first `skip` characters.
    fn compare(
        &self,
        query: &[u8],
        suffix_pos: usize,
        skip: usize,
    ) -> Result<Comparison> {
        let max_query_len = self.max_query_len;
        let lcp = self.query_lcp(query, suffix_pos, skip)?;
        let cmp = if (max_query_len > 0) && (lcp >= max_query_len) {
            Ordering::Equal
        } else {
            let full_offset = match &self.seed_mask {
                Some(seed_mask) => mask_offset(lcp, seed_mask),
                None => lcp,
            };
            match (
                query.get(full_offset),
                &self.text_byte(suffix_pos + full_offset)?,
            ) {
                (None, _) => Ordering::Equal,
                (Some(a), Some(b)) => a.cmp(b),
                _ => unreachable!(),
            }
        };

        Ok(Comparison { lcp, cmp })
    }

    /// The LCP of `query` and the suffix at `suffix_pos`, which share their
    /// first `skip` characters, counted in care positions under a seed mask.
    fn query_lcp(&self, query: &[u8], suffix_pos: usize, skip: usize) -> Result<usize> {
        let text_len = self.header.text_len;
        let max_query_len = self.max_query_len;
        Ok(match &self.seed_mask {
            None => {
                if (max_query_len > 0) && (skip >= max_query_len) {
                    skip
                } else {
                    let text_start = suffix_pos + skip;
                    let text_end = if max_query_len > 0 {
                        min(text_len, text_start + max_query_len)
                    } else {
                        min(text_len, text_start + query.len())
                    };
                    let text = self.text_slice(text_start..text_end)?;
                    skip + lcp(query.get(skip..).unwrap_or_default(), &text, None)
                }
            }
            Some(seed_mask) => {
                if skip >= seed_mask.weight
                    || (max_query_len > 0 && skip >= max_query_len)
                {
                    skip
                } else {
                    let end = if max_query_len > 0 {
                        min(max_query_len, seed_mask.weight)
                    } else {
                        seed_mask.weight
                    };
                    let positions = &seed_mask.positions[skip..end];
                    let last = positions[positions.len() - 1];
                    let text_end = min(text_len, suffix_pos + last + 1);
                    let text = self.text_slice(suffix_pos..text_end)?;
                    skip + masked_lcp(query, &text, positions)
                }
            }
        })
    }

    pub(crate) fn locate(&self, ranks: Range<usize>) -> Result<Vec<Location<'_>>> {
        let seq_starts = &self.header.sequence_starts;
        let suffixes = self.suffixes_at(ranks.clone())?;
        let mut locations: Vec<Location> = ranks
            .zip(suffixes)
            .map(|(rank, suffix)| {
                let suffix = suffix.to_usize();
                let i = seq_starts.partition_point(|&val| val <= suffix) - 1;
                Location {
                    rank,
                    suffix,
                    sequence_name: &self.header.sequence_names[i],
                    sequence_start: seq_starts[i],
                    position: suffix - seq_starts[i],
                }
            })
            .collect();
        locations.sort_by(|a, b| {
            a.sequence_name
                .cmp(b.sequence_name)
                .then(a.position.cmp(&b.position))
        });
        Ok(locations)
    }

    pub(crate) fn extract(
        &self,
        ranks: Range<usize>,
        context: Context,
    ) -> Result<Vec<Extraction<'_>>> {
        let suffixes = self.suffixes_at(ranks.clone())?;
        Ok(ranks
            .zip(suffixes)
            .map(|(rank, suffix)| self.extraction(rank, suffix.to_usize(), context))
            .collect())
    }

    fn extraction(
        &self,
        rank: usize,
        suffix: usize,
        context: Context,
    ) -> Extraction<'_> {
        let seq_starts = &self.header.sequence_starts;
        let text_len = self.header.text_len;
        let i = seq_starts.partition_point(|&val| val <= suffix) - 1;
        let sequence_start = seq_starts[i];

        // TODO: seq_end is a position in the text but is compared
        //       with positions relative to the sequence start, so
        //       context past the first sequence runs past its end
        let seq_end = if i == seq_starts.len() - 1 {
            text_len
        } else {
            seq_starts[i + 1]
        };
        let relative_suffix_start = suffix - sequence_start;
        let context_start = relative_suffix_start.saturating_sub(context.before);
        let context_end = min(
            context
                .after
                .map_or(seq_end, |len| relative_suffix_start + len),
            seq_end,
        );
        Extraction {
            location: Location {
                rank,
                suffix,
                sequence_name: &self.header.sequence_names[i],
                sequence_start,
                position: relative_suffix_start,
            },
            sequence_range: context_start..context_end,
            suffix_offset: relative_suffix_start - context_start,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::{read_queries, Context, Extraction, Location};
    use crate::{OpenOptions, Storage, SuffixArray};
    use anyhow::Result;
    use pretty_assertions::assert_eq;
    use std::io::Write;
    use tempfile::NamedTempFile;

    fn open(file: &str, text: Storage, suffixes: Storage) -> Result<SuffixArray> {
        SuffixArray::open_with(
            file,
            &OpenOptions {
                text,
                suffixes,
                max_query_len: None,
            },
        )
    }

    /// An extraction from a one-sequence file, whose sequence starts at 0.
    fn extraction(
        suffix: usize,
        rank: usize,
        sequence_range: std::ops::Range<usize>,
        suffix_offset: usize,
    ) -> Extraction<'static> {
        Extraction {
            location: Location {
                rank,
                suffix,
                sequence_name: "1",
                sequence_start: 0,
                position: suffix,
            },
            sequence_range,
            suffix_offset,
        }
    }

    /// The bytes of an extraction's window.
    fn window(sa: &SuffixArray, ex: &Extraction) -> Result<Vec<u8>> {
        let start = ex.location.sequence_start;
        sa.text(start + ex.sequence_range.start..start + ex.sequence_range.end)
    }

    const ONE_THREE: Context = Context {
        before: 1,
        after: Some(3),
    };

    #[test]
    fn test_count() -> Result<()> {
        let sa = open("../data/inputs/1.sufr", Storage::Memory, Storage::Disk)?;
        let counts = sa
            .search_all(&["AC", "AG", "GT"])
            .into_iter()
            .map(|hits| Ok(hits?.count()))
            .collect::<Result<Vec<_>>>()?;
        assert_eq!(counts, [2, 0, 2]);
        Ok(())
    }

    #[test]
    fn test_search() -> Result<()> {
        // the text is ACGTNNACGT$
        let sa = open("../data/inputs/1.sufr", Storage::Memory, Storage::Disk)?;
        let hits = sa.search(b"ACG")?;
        assert_eq!(hits.ranks(), 1..3);
        assert_eq!(hits.suffixes()?, [6, 0]);

        let hits = sa.search(b"GGC")?;
        assert_eq!(hits.ranks(), 0..0);
        assert!(hits.is_empty());
        assert_eq!(hits.suffixes()?, []);
        Ok(())
    }

    #[test]
    fn test_extract_text_on_disk() -> Result<()> {
        let sa = open("../data/inputs/1.sufr", Storage::Disk, Storage::Disk)?;
        let extractions = sa.search(b"ACG")?.extract(ONE_THREE)?;
        assert_eq!(
            extractions,
            [extraction(6, 1, 5..9, 1), extraction(0, 2, 0..3, 0)]
        );
        assert_eq!(window(&sa, &extractions[0])?, b"NACG");
        assert_eq!(window(&sa, &extractions[1])?, b"ACG");
        assert_eq!(sa.search(b"GTC")?.extract(ONE_THREE)?, []);
        Ok(())
    }

    #[test]
    fn test_extract() -> Result<()> {
        let sa = open("../data/expected/1.sufr", Storage::Memory, Storage::Disk)?;
        let extractions = sa
            .search_all(&["AC", "GT", "XX"])
            .into_iter()
            .map(|hits| hits?.extract(ONE_THREE))
            .collect::<Result<Vec<_>>>()?;
        assert_eq!(
            extractions,
            [
                vec![extraction(6, 1, 5..9, 1), extraction(0, 2, 0..3, 0)],
                vec![extraction(8, 5, 7..11, 1), extraction(2, 6, 1..5, 1)],
                vec![],
            ]
        );
        assert_eq!(window(&sa, &extractions[1][0])?, b"CGT$");
        assert_eq!(window(&sa, &extractions[1][1])?, b"CGTN");
        Ok(())
    }

    #[test]
    fn test_extract_without_after() -> Result<()> {
        let sa = open("../data/inputs/1.sufr", Storage::Memory, Storage::Memory)?;
        let extractions = sa.search(b"CGT")?.extract(Context {
            before: 1,
            after: None,
        })?;
        assert_eq!(
            extractions,
            [extraction(7, 3, 6..11, 1), extraction(1, 4, 0..11, 1)]
        );
        assert_eq!(window(&sa, &extractions[0])?, b"ACGT$");
        assert_eq!(window(&sa, &extractions[1])?, b"ACGTNNACGT$");
        Ok(())
    }

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

        // (query, [(rank, suffix)] sorted by position)
        let expected: [(&str, &[(usize, usize)]); 5] = [
            (
                "A",
                &[(1, 0), (4, 1), (5, 3), (6, 5), (7, 7), (3, 10), (2, 12)],
            ),
            (
                "B",
                &[
                    (11, 2),
                    (12, 4),
                    (13, 6),
                    (14, 8),
                    (10, 9),
                    (9, 11),
                    (8, 13),
                ],
            ),
            ("ABAB", &[(4, 1), (5, 3), (6, 5), (3, 10)]),
            ("ABABB", &[(6, 5)]),
            ("BBBB", &[]),
        ];
        for suffixes in [Storage::Disk, Storage::Memory] {
            let sa = open("../data/expected/abba.sufr", Storage::Memory, suffixes)?;
            for (query, hits) in expected {
                let expected: Vec<_> = hits
                    .iter()
                    .map(|&(rank, suffix)| Location {
                        rank,
                        suffix,
                        sequence_name: "1",
                        sequence_start: 0,
                        position: suffix,
                    })
                    .collect();
                assert_eq!(sa.search(query.as_bytes())?.locate()?, expected, "{query}");
            }
        }
        Ok(())
    }

    #[test]
    fn test_read_queries() -> Result<()> {
        let mut file = NamedTempFile::new()?;
        write!(file, "AC GT\n  ACG\tXX\n\nT\n")?;
        let path = file.path().to_string_lossy().to_string();
        let queries = ["A".to_string(), path, "C".to_string()];
        assert_eq!(
            read_queries(&queries)?,
            ["A", "AC", "GT", "ACG", "XX", "T", "C"]
        );
        Ok(())
    }
}
