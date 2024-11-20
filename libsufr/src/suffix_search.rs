use crate::{Comparison, FileAccess, FromUsize, Int};
use anyhow::Result;
use std::{
    cmp::{min, Ordering},
    ops::Range,
};

// --------------------------------------------------
#[derive(Debug, Clone)]
pub struct SearchOptions {
    pub queries: Vec<String>,
    pub max_query_len: Option<usize>,
    pub low_memory: bool,
    pub find_suffixes: bool,
}

// --------------------------------------------------
#[derive(Debug)]
pub struct SearchResult<T>
where
    T: Int + FromUsize<T> + Sized + Send + Sync + serde::ser::Serialize,
{
    pub query_num: usize,
    pub query: String,
    pub locations: Option<SearchResultLocations<T>>,
}

// --------------------------------------------------
#[derive(Debug)]
pub struct SearchResultLocations<T> {
    pub ranks: Range<usize>,
    pub suffixes: Vec<T>,
}

// --------------------------------------------------
#[derive(Debug)]
pub struct SuffixSearch<'a, T>
where
    T: Int + FromUsize<T> + Sized + Send + Sync + serde::ser::Serialize,
{
    text: &'a [u8],
    suffix_array_file: FileAccess<T>,
    suffix_array_mem: &'a [T],
    suffix_array_rank_mem: &'a [usize],
    query_low_memory: bool,
    num_suffixes: usize,
    //suffix_array_mem_mql: Option<usize>,
}

// --------------------------------------------------
impl<'a, T> SuffixSearch<'a, T>
where
    T: Int + FromUsize<T> + Sized + Send + Sync + serde::ser::Serialize,
{
    pub fn new(
        text: &'a [u8],
        file: FileAccess<T>,
        sa: &'a [T],
        rank: &'a [usize],
        query_low_memory: bool,
        num_suffixes: usize,
    ) -> SuffixSearch<'a, T> {
        SuffixSearch {
            text,
            suffix_array_file: file,
            suffix_array_mem: sa,
            suffix_array_rank_mem: rank,
            query_low_memory,
            num_suffixes: if sa.is_empty() { num_suffixes } else { sa.len() },
        }
    }

    // --------------------------------------------------
    pub fn search(
        &mut self,
        query_num: usize,
        query: &str,
        find_suffixes: bool
    ) -> Result<SearchResult<T>> {
        let qry = query.as_bytes();
        let n = self.num_suffixes;
        if let Some(start) = self.suffix_search_first(qry, 0, n - 1, 0, 0) {
            let end = self
                .suffix_search_last(qry, start, n - 1, n, 0, 0)
                .unwrap_or(start);

            // Rank is empty when we have the full SA in memory
            // AND when doing low-memory searches
            let ranks = if self.suffix_array_rank_mem.is_empty() {
                start..end + 1
            } else {
                // This is the case for the compressed/in-memory SA
                let start_rank = self.suffix_array_rank_mem[start];
                let end_rank = if start == end {
                    if start == self.suffix_array_rank_mem.len() - 1 {
                        // We're on the last rank, so go to end
                        self.num_suffixes
                    } else {
                        // Use the next LCP rank
                        self.suffix_array_rank_mem[start + 1]
                    }
                } else {
                    self.suffix_array_rank_mem[end] + 1
                };
                start_rank..end_rank
            };

            // This requires going to disk but "count" doesn't need them
            let suffixes = if find_suffixes {
                self.suffix_array_file.get_range(ranks.clone())?
            } else {
                vec![]
            };

            Ok(SearchResult {
                query_num,
                query: query.to_string(),
                locations: Some(SearchResultLocations { suffixes, ranks }),
            })
        } else {
            Ok(SearchResult {
                query_num,
                query: query.to_string(),
                locations: None,
            })
        }
    }

    // --------------------------------------------------
    fn suffix_search_first(
        &mut self,
        qry: &[u8],
        low: usize,
        high: usize,
        left_lcp: usize,
        right_lcp: usize,
    ) -> Option<usize> {
        if high >= low {
            let mid = low + ((high - low) / 2);
            let mid_val = self.get_suffix(mid)?.to_usize();
            let mid_cmp = self.compare(qry, mid_val, min(left_lcp, right_lcp));

            let mid_minus_one = if mid > 0 {
                self.get_suffix(mid - 1)?.to_usize()
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
        &mut self,
        qry: &[u8],
        low: usize,
        high: usize,
        n: usize,
        left_lcp: usize,
        right_lcp: usize,
    ) -> Option<usize> {
        if high >= low {
            let mid = low + ((high - low) / 2);
            let mid_val = self.get_suffix(mid)?.to_usize();
            let mid_cmp = self.compare(qry, mid_val, min(left_lcp, right_lcp));

            // Weird hack because I cannot embed this call in the "if"
            let mid_plus_one = if mid < n - 1 {
                self.get_suffix(mid + 1)?.to_usize()
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
    pub fn compare(&self, query: &[u8], suffix_pos: usize, skip: usize) -> Comparison {
        let lcp = query
            .iter()
            .zip(self.text.get(suffix_pos..).unwrap())
            .skip(skip)
            .map_while(|(a, b)| (a == b).then_some(a))
            .count()
            + skip;

        let cmp = match (query.get(lcp), self.text.get(suffix_pos + lcp)) {
            // Entire query matched
            (None, _) => Ordering::Equal,
            // Compare next char
            (Some(a), Some(b)) => a.cmp(b),
            // Panic at the disco
            _ => unreachable!(),
        };

        Comparison { lcp, cmp }
    }

    // --------------------------------------------------
    fn get_suffix(&mut self, pos: usize) -> Option<T> {
        if self.query_low_memory {
            self.suffix_array_file.get(pos)
        } else {
            self.suffix_array_mem.get(pos).copied()
        }
    }

    // --------------------------------------------------
    //pub fn locate(&mut self, args: Locate) -> Result<Vec<Result<LocateResult<T>>>> {
    //    self.query_low_memory = args.low_memory;
    //    let n = if self.query_low_memory {
    //        self.num_suffixes.to_usize()
    //    } else {
    //        let max_query_len =
    //            args.max_query_len.unwrap_or(self.max_query_len.to_usize());
    //        self.set_suffix_array_mem(max_query_len)?;
    //        self.suffix_array_mem.len()
    //    };
    //    let seq_starts = self.sequence_starts.clone();
    //    let seq_names = self.headers.clone();
    //
    //    let now = Instant::now();
    //    let res: Vec<_> = args
    //        .queries
    //        .into_iter()
    //        .map(|query| -> Result<LocateResult<T>> {
    //            let qry = query.as_bytes();
    //
    //            if let Some(start) = self.suffix_search_first(qry, 0, n - 1, 0, 0) {
    //                let end = self
    //                    .suffix_search_last(qry, start, n - 1, n, 0, 0)
    //                    .unwrap_or(start);
    //
    //                // Rank is empty when we have the full SA in memory
    //                // AND when doing low-memory searches
    //                let (suffixes, ranks) = if self.suffix_array_rank_mem.is_empty() {
    //                    let (start_rank, end_rank) = (start, end + 1);
    //                    // For low-memory, go to disk
    //                    let suffixes = if self.suffix_array_mem.is_empty() {
    //                        self.suffix_array_file.get_range(start_rank..end_rank)?
    //                    } else {
    //                        // Otherwise, get from memory
    //                        self.suffix_array_mem[start_rank..end_rank].to_vec()
    //                    };
    //                    (suffixes, start_rank..end_rank)
    //                } else {
    //                    // This is the case for the compressed/in-memory SA
    //                    let start_rank = self.suffix_array_rank_mem[start];
    //                    let end_rank = if start == end {
    //                        if start == self.suffix_array_rank_mem.len() - 1 {
    //                            // We're on the last rank, so go to end
    //                            self.num_suffixes.to_usize()
    //                        } else {
    //                            // Use the next LCP rank
    //                            self.suffix_array_rank_mem[start + 1]
    //                        }
    //                    } else {
    //                        self.suffix_array_rank_mem[end] + 1
    //                    };
    //
    //                    // I have to go to disk to get the actual suffixes
    //                    let suffixes =
    //                        self.suffix_array_file.get_range(start_rank..end_rank)?;
    //                    (suffixes, start_rank..end_rank)
    //                };
    //
    //                let positions: Vec<_> = suffixes
    //                    .iter()
    //                    .map(|&suffix| {
    //                        let i =
    //                            seq_starts.partition_point(|&val| val <= suffix) - 1;
    //                        LocateResultPosition {
    //                            suffix,
    //                            sequence_name: seq_names[i].clone(),
    //                            sequence_position: suffix - seq_starts[i],
    //                        }
    //                    })
    //                    .collect();
    //
    //                Ok(LocateResult {
    //                    query: query.to_string(),
    //                    positions,
    //                    ranks,
    //                })
    //            } else {
    //                Err(anyhow!("{query}"))
    //            }
    //        })
    //        .collect();
    //
    //    info!("Search finished in {:?}", now.elapsed());
    //    Ok(res)
    //}
}
