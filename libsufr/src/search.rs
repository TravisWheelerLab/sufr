//! Search a suffix array

use crate::{
    index::Index,
    types::{
        Comparison, CountOptions, CountResult, ExtractOptions, ExtractResult,
        ExtractSequence, Int, LocateOptions, LocatePosition, LocateResult,
        SearchOptions, SearchResult, SearchResultLocations, SuffixSortType,
    },
    util::find_lcp_full_offset,
};
use anyhow::Result;
use log::info;
use rayon::prelude::*;
use std::{
    cmp::{min, Ordering},
    time::Instant,
};

impl<T: Int> Index<T> {
    /// Count the occurrences of queries in a suffix array
    pub fn count(&self, args: CountOptions) -> Result<Vec<CountResult>> {
        let search_args = SearchOptions {
            queries: args.queries,
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
    pub fn suffix_search(&self, args: &SearchOptions) -> Result<Vec<SearchResult<T>>> {
        let now = Instant::now();
        let mut res: Vec<_> = args
            .queries
            .clone()
            .into_par_iter()
            .enumerate()
            .flat_map(|(query_num, query)| -> Result<SearchResult<T>> {
                self.search(query_num, &query, args.find_suffixes)
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
    /// Find one query, and its suffixes when `find_suffixes`.
    fn search(
        &self,
        query_num: usize,
        query: &str,
        find_suffixes: bool,
    ) -> Result<SearchResult<T>> {
        let qry = query.as_bytes();
        let n = self.num_searched();
        let Some(start) = self.suffix_search_first(qry, 0, n - 1, 0, 0) else {
            return Ok(SearchResult {
                query_num,
                query: query.to_string(),
                locations: None,
            });
        };

        let end = self
            .suffix_search_last(qry, start, n - 1, n, 0, 0)
            .unwrap_or(start);
        let ranks = self.rank_range(start, end);

        // Getting suffixes may require going to disk
        // but "count" doesn't need them so this is an option
        let suffixes = if find_suffixes {
            self.suffixes_at(ranks.clone())?
        } else {
            vec![]
        };

        Ok(SearchResult {
            query_num,
            query: query.to_string(),
            locations: Some(SearchResultLocations { suffixes, ranks }),
        })
    }

    // --------------------------------------------------
    fn suffix_search_first(
        &self,
        qry: &[u8],
        low: usize,
        high: usize,
        left_lcp: usize,
        right_lcp: usize,
    ) -> Option<usize> {
        if high >= low {
            let mid = low + ((high - low) / 2);
            let mid_val = self.suffix(mid)?.to_usize();
            let mid_cmp = self.compare(qry, mid_val, min(left_lcp, right_lcp));

            let mid_minus_one = if mid > 0 {
                self.suffix(mid - 1)?.to_usize()
            } else {
                mid_val
            };

            if mid_cmp.cmp == Ordering::Equal
                && (mid == 0
                    || self.compare(qry, mid_minus_one, 0).cmp == Ordering::Greater)
            {
                Some(mid)
            } else if mid_cmp.cmp == Ordering::Greater {
                self.suffix_search_first(qry, mid + 1, high, mid_cmp.lcp, right_lcp)
            } else {
                // Ordering::Less
                self.suffix_search_first(qry, low, mid - 1, left_lcp, mid_cmp.lcp)
            }
        } else {
            None
        }
    }

    // --------------------------------------------------
    fn suffix_search_last(
        &self,
        qry: &[u8],
        low: usize,
        high: usize,
        n: usize,
        left_lcp: usize,
        right_lcp: usize,
    ) -> Option<usize> {
        if high >= low {
            let mid = low + ((high - low) / 2);
            let mid_val = self.suffix(mid)?.to_usize();
            let mid_cmp = self.compare(qry, mid_val, min(left_lcp, right_lcp));

            // Weird hack because I cannot embed this call in the "if"
            let mid_plus_one = if mid < n - 1 {
                self.suffix(mid + 1)?.to_usize()
            } else {
                mid_val
            };

            if mid_cmp.cmp == Ordering::Equal
                && (mid == n - 1
                    || self.compare(qry, mid_plus_one, 0).cmp == Ordering::Less)
            {
                Some(mid)
            } else if mid_cmp.cmp == Ordering::Less {
                self.suffix_search_last(qry, low, mid - 1, n, left_lcp, mid_cmp.lcp)
            } else {
                self.suffix_search_last(qry, mid + 1, high, n, mid_cmp.lcp, right_lcp)
            }
        } else {
            None
        }
    }

    // --------------------------------------------------
    /// Compare a query string to a suffix.
    /// When the suffix array was sorted using a seed mask, only
    /// compare the "care" positions; otherwise, compare all the positions
    /// up to an optional maximum query length
    fn compare(&self, query: &[u8], suffix_pos: usize, skip: usize) -> Comparison {
        let text_len = self.header.text_len;
        let sort_type = &self.header.sort_type;
        let (lcp, max_query_len) = match sort_type {
            SuffixSortType::MaxQueryLen(mql) => {
                // The "MaxQueryLen(mql)" refers to how the suffix array
                // was built, but there may be a runtime value in self.max_query_len.
                // If both are present, take the lowest.
                let max_query_len: usize = if mql > &0 && self.max_query_len.is_some() {
                    min(*mql, self.max_query_len.unwrap_or(0))
                } else if let Some(val) = self.max_query_len {
                    val
                } else {
                    *mql
                };

                let lcp = if (max_query_len > 0) && (skip >= max_query_len) {
                    // If we've already seen enough
                    skip
                } else {
                    let text_start = suffix_pos + skip;
                    let text_end = if max_query_len > 0 {
                        min(text_len, text_start + max_query_len)
                    } else {
                        min(text_len, text_start + query.len())
                    };

                    query
                        .iter()
                        .skip(skip)
                        .zip(self.text_range(text_start..text_end).unwrap())
                        .map_while(|(a, b)| (a == &b).then_some(a))
                        .count()
                        + skip
                };
                (lcp, max_query_len)
            }
            SuffixSortType::Mask(seed_mask) => {
                let max_query_len = self.max_query_len.unwrap_or(0);
                let lcp = if skip >= seed_mask.weight
                    || (max_query_len > 0 && skip >= max_query_len)
                {
                    skip
                } else {
                    let end = if max_query_len > 0 {
                        min(max_query_len, seed_mask.weight)
                    } else {
                        seed_mask.weight
                    };
                    let mask_pos = &seed_mask.positions[skip..end];
                    let query_len =
                        mask_pos.iter().filter(|&v| v < &query.len()).count();
                    let suffix_len = mask_pos
                        .iter()
                        .map(|offset| suffix_pos + offset) // add offset to start of suffix
                        .filter(|&v| v < text_len) // don't go off end of text
                        .count();
                    let len = min(query_len, suffix_len);

                    if len > 0 {
                        skip + seed_mask.positions[skip..skip + len]
                            .iter()
                            .take_while(|&offset| {
                                self.text_byte(suffix_pos + offset)
                                    .is_some_and(|val| query[*offset] == val)
                            })
                            .count()
                    } else {
                        skip
                    }
                };
                (lcp, max_query_len)
            }
        };

        let cmp = if (max_query_len > 0) && (lcp >= max_query_len) {
            // We've seen enough
            Ordering::Equal
        } else {
            // Get the next chars
            let full_offset = find_lcp_full_offset(lcp, sort_type);
            match (
                query.get(full_offset),
                &self.text_byte(suffix_pos + full_offset),
            ) {
                // Entire query matched
                (None, _) => Ordering::Equal,

                // Compare next char
                (Some(a), Some(b)) => a.cmp(b),

                // Panic at the disco
                _ => unreachable!(),
            }
        };

        Comparison { lcp, cmp }
    }

    // --------------------------------------------------
    /// Extract the suffixes for a given set of queries
    pub fn extract(&self, args: ExtractOptions) -> Result<Vec<ExtractResult>> {
        let search_args = SearchOptions {
            queries: args.queries,
            find_suffixes: true,
        };
        let search_result = &self.suffix_search(&search_args)?;
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
                            self.extract_sequence(
                                rank,
                                suffix.to_usize(),
                                args.prefix_len,
                                args.suffix_len,
                            )
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

    /// The sequence context of one hit, from `prefix_len` before the
    /// suffix to `suffix_len` after its start, within its sequence.
    fn extract_sequence(
        &self,
        rank: usize,
        suffix: usize,
        prefix_len: Option<usize>,
        suffix_len: Option<usize>,
    ) -> ExtractSequence {
        let seq_starts = &self.header.sequence_starts;
        let i = seq_starts.partition_point(|&val| val <= suffix) - 1;
        let sequence_start = seq_starts[i];
        let seq_end = if i == seq_starts.len() - 1 {
            self.header.text_len
        } else {
            seq_starts[i + 1]
        };
        let relative_suffix_start = suffix - sequence_start;
        let context_start =
            relative_suffix_start.saturating_sub(prefix_len.unwrap_or(0));
        let context_end = min(
            suffix_len.map_or(seq_end, |len| relative_suffix_start + len),
            seq_end,
        );
        ExtractSequence {
            rank,
            suffix,
            sequence_name: self.header.sequence_names[i].clone(),
            sequence_start,
            sequence_range: (context_start..context_end),
            suffix_offset: relative_suffix_start - context_start,
        }
    }

    // --------------------------------------------------
    /// Find the positions of queries in a suffix array.
    pub fn locate(&self, args: LocateOptions) -> Result<Vec<LocateResult>> {
        let search_opts = SearchOptions {
            queries: args.queries,
            find_suffixes: true,
        };
        let search_result = &self.suffix_search(&search_opts)?;
        let seq_starts = &self.header.sequence_starts;
        let seq_names = &self.header.sequence_names;
        let mut locate_result: Vec<LocateResult> = vec![];
        let now = Instant::now();

        // Augment the search with relative sequence positions
        for res in search_result {
            let mut positions = vec![];
            if let Some(locs) = &res.locations {
                for (rank, suffix) in locs.ranks.clone().zip(locs.suffixes.clone()) {
                    let suffix = suffix.to_usize();
                    let i = seq_starts.partition_point(|&val| val <= suffix) - 1;
                    positions.push(LocatePosition {
                        rank,
                        suffix,
                        sequence_name: seq_names[i].clone(),
                        sequence_position: suffix - seq_starts[i],
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
mod tests {
    use crate::{
        index::Index,
        types::{
            CountOptions, CountResult, ExtractOptions, ExtractResult, ExtractSequence,
            LocateOptions, LocatePosition, LocateResult, SearchOptions, SearchResult,
            SearchResultLocations,
        },
    };
    use anyhow::Result;

    #[test]
    fn test_count() -> Result<()> {
        let sufr = Index::<u32>::open("../data/inputs/1.sufr", true, false, None)?;
        let opts = CountOptions {
            queries: vec!["AC".to_string(), "AG".to_string(), "GT".to_string()],
        };
        let expected = vec![
            CountResult {
                query_num: 0,
                query: "AC".to_string(),
                count: 2,
            },
            CountResult {
                query_num: 1,
                query: "AG".to_string(),
                count: 0,
            },
            CountResult {
                query_num: 2,
                query: "GT".to_string(),
                count: 2,
            },
        ];
        assert_eq!(sufr.count(opts)?, expected);
        Ok(())
    }

    #[test]
    fn test_suffix_search() -> Result<()> {
        // the text is ACGTNNACGT$
        let sufr = Index::<u32>::open("../data/inputs/1.sufr", true, false, None)?;
        let opts = SearchOptions {
            queries: vec!["ACG".to_string(), "GGC".to_string()],
            find_suffixes: true,
        };
        let expected: Vec<SearchResult<u32>> = vec![
            SearchResult {
                query_num: 0,
                query: "ACG".to_string(),
                locations: Some(SearchResultLocations {
                    ranks: 1..3,
                    suffixes: vec![6, 0],
                }),
            },
            SearchResult {
                query_num: 1,
                query: "GGC".to_string(),
                locations: None,
            },
        ];
        assert_eq!(sufr.suffix_search(&opts)?, expected);
        Ok(())
    }

    #[test]
    fn test_extract_text_on_disk() -> Result<()> {
        let sufr = Index::<u32>::open("../data/inputs/1.sufr", true, true, None)?;
        let extract_opts = ExtractOptions {
            queries: vec!["ACG".to_string(), "GTC".to_string()],
            prefix_len: Some(1),
            suffix_len: Some(3),
        };
        let expected = vec![
            ExtractResult {
                query_num: 0,
                query: "ACG".to_string(),
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
                query: "GTC".to_string(),
                sequences: vec![],
            },
        ];
        assert_eq!(sufr.extract(extract_opts)?, expected);
        Ok(())
    }

    // --------------------------------------------------
    #[test]
    fn test_extract() -> Result<()> {
        // cargo run -- ex data/expected/1.sufr AC GT XX
        let sufr_file: Index<u32> =
            Index::open("../data/expected/1.sufr", true, false, None)?;

        let opts = ExtractOptions {
            queries: vec!["AC".to_string(), "GT".to_string(), "XX".to_string()],
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

        for low_memory in [true, false] {
            let sufr_file: Index<u32> =
                Index::open("../data/expected/abba.sufr", low_memory, false, None)?;
            let args = LocateOptions {
                queries: vec!["A".to_string()],
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

        for low_memory in [true, false] {
            let sufr_file: Index<u32> =
                Index::open("../data/expected/abba.sufr", low_memory, false, None)?;
            let args = LocateOptions {
                queries: vec!["B".to_string()],
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

        for low_memory in [true, false] {
            let sufr_file: Index<u32> =
                Index::open("../data/expected/abba.sufr", low_memory, false, None)?;
            let args = LocateOptions {
                queries: vec!["ABAB".to_string()],
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

        for low_memory in [true, false] {
            let sufr_file: Index<u32> =
                Index::open("../data/expected/abba.sufr", low_memory, false, None)?;
            let args = LocateOptions {
                queries: vec!["ABABB".to_string()],
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

        for low_memory in [true, false] {
            let sufr_file: Index<u32> =
                Index::open("../data/expected/abba.sufr", low_memory, false, None)?;
            let args = LocateOptions {
                queries: vec!["BBBB".to_string()],
            };
            let res = sufr_file.locate(args);
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
}
