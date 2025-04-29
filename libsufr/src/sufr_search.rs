//! # Search a suffix array

use crate::{
    file_access::FileAccess,
    types::{
        Comparison, FromUsize, Int, BisectResult, SearchResult, SearchResultLocations, SuffixSortType,
    },
    util::find_lcp_full_offset,
};
use anyhow::Result;
use std::{
    cmp::{min, Ordering},
    ops::Range,
};

// --------------------------------------------------
/// Arguments to create a search
#[derive(Debug)]
pub struct SufrSearchArgs<'a, T>
where
    T: Int + FromUsize<T> + Sized + Send + Sync + serde::ser::Serialize,
{
    /// The text being searched (which may be empty in the case of very
    /// low memory)
    pub text: &'a [u8],

    /// The length of the text
    pub text_len: usize,

    /// A `FileAccess` to the text
    pub text_file: FileAccess<u8>,

    /// File?
    pub file: FileAccess<T>,

    /// In-memory suffix array
    pub suffix_array: &'a [T],

    /// In-memory rank array
    pub rank: &'a [T],

    /// The number of suffixes 
    pub len_suffixes: usize,

    /// How the suffixes were built?
    pub sort_type: &'a SuffixSortType,

    /// The maximum query length to use when querying
    pub max_query_len: Option<usize>,
}

// --------------------------------------------------
#[derive(Debug)]
pub struct SufrSearch<'a, T>
where
    T: Int + FromUsize<T> + Sized + Send + Sync + serde::ser::Serialize,
{
    text: &'a [u8],
    text_len: usize,
    text_file: FileAccess<u8>,
    suffix_array_file: FileAccess<T>,
    suffix_array_mem: &'a [T],
    suffix_array_rank_mem: &'a [T],
    len_suffixes: usize,
    sort_type: &'a SuffixSortType,
    max_query_len: Option<usize>,
}

// --------------------------------------------------
impl<'a, T> SufrSearch<'a, T>
where
    T: Int + FromUsize<T> + Sized + Send + Sync + serde::ser::Serialize,
{
    /// Create a `SufrSearch` struct
    ///
    /// Args:
    /// * `args`: a `SufrSearchArgs` struct
    pub fn new(args: SufrSearchArgs<'a, T>) -> SufrSearch<'a, T> {
        SufrSearch {
            text: args.text,
            text_len: args.text_len,
            text_file: args.text_file,
            suffix_array_file: args.file,
            suffix_array_mem: args.suffix_array,
            suffix_array_rank_mem: args.rank,
            len_suffixes: if args.suffix_array.is_empty() {
                args.len_suffixes
            } else {
                args.suffix_array.len()
            },
            sort_type: args.sort_type,
            max_query_len: args.max_query_len,
        }
    }

    // --------------------------------------------------
    /// Find a query string in a suffix array.
    /// Returns a `SearchResult`
    ///
    /// Args:
    /// * `query_num`: ordinal number of the query
    /// * `query`: a string to search for
    /// * `find_suffixes`: whether or not to return the suffixes locations
    pub fn search(
        &mut self,
        query_num: usize,
        query: &str,
        find_suffixes: bool,
    ) -> Result<SearchResult<T>> {
        let qry = query.as_bytes();
        let n = self.len_suffixes;
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
                        self.len_suffixes
                    } else {
                        // Use the next LCP rank
                        self.suffix_array_rank_mem[start + 1].to_usize()
                    }
                } else {
                    self.suffix_array_rank_mem[end].to_usize() + 1
                };
                start_rank.to_usize()..end_rank
            };

            // Getting suffixes may require going to disk
            // but "count" doesn't need them so this is an option
            let suffixes = if find_suffixes {
                // First condition is that there is no in-memory SA
                if self.suffix_array_mem.is_empty()
                    // Second condition is there is an in-memory SA
                    // BUT the ranks are also there meaning the SA
                    // is compressed, forcing us to go to disk
                    || !self.suffix_array_rank_mem.is_empty()
                {
                    self.suffix_array_file.get_range(ranks.clone())?
                } else {
                    // Otherwise, get from memory
                    self.suffix_array_mem[ranks.clone()].to_vec()
                }
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

    pub fn bisect(
        &mut self,
        query_num: usize,
        query: &str,
        low: usize,
        high: usize,
    ) -> Result<BisectResult> {
        let qry = query.as_bytes();
        if let Some(start) = self.suffix_search_first(qry, low, high, 0, 0) {
            // something was found
            let end = self
                .suffix_search_last(qry, start, high, high + 1, 0, 0)
                .unwrap_or(start);
            Ok(BisectResult {
                query_num: query_num,
                query: query.to_string(),
                count: end - start + 1,
                first_position: start,
                last_position: end,
            })
        } else {
            // nothing was found
            Ok(BisectResult {
                query_num: query_num,
                query: query.to_string(),
                count: 0,
                first_position: 0,
                last_position: 0,
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
    /// Compare a query string to a suffix.
    /// When the suffix array was sorted using a seed mask, only
    /// compare the "care" positions; otherwise, compare all the positions
    /// up to an optional maximum query length
    /// Returns a `Comparison`
    ///
    /// Args:
    /// * `query`: string to search for
    /// * `suffix_pos`: suffix position
    /// * `skip`: number of places to skip
    fn compare(&mut self, query: &[u8], suffix_pos: usize, skip: usize) -> Comparison {
        let (lcp, max_query_len) = match &self.sort_type {
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
                        min(self.text_len, text_start + max_query_len)
                    } else {
                        min(self.text_len, text_start + query.len())
                    };

                    query
                        .iter()
                        .skip(skip)
                        .zip(self.get_text_range(text_start..text_end).unwrap())
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
                        .filter(|&v| v < self.text_len) // don't go off end of text
                        .count();
                    let len = min(query_len, suffix_len);

                    if len > 0 {
                        skip + seed_mask.positions[skip..skip + len]
                            .iter()
                            .take_while(|&offset| {
                                self.get_text(suffix_pos + offset)
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
            let full_offset = find_lcp_full_offset(lcp, self.sort_type);
            match (
                query.get(full_offset),
                &self.get_text(suffix_pos + full_offset),
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
    fn get_text(&mut self, pos: usize) -> Option<u8> {
        if self.text.is_empty() {
            self.text_file.get(pos)
        } else {
            self.text.get(pos).copied()
        }
    }

    // --------------------------------------------------
    pub fn get_text_range(&mut self, pos: Range<usize>) -> Result<Vec<u8>> {
        // this is too expensive, copies loooooong stretches of text into vec
        if self.text.is_empty() {
            self.text_file.get_range(pos.clone())
        } else {
            Ok(self.text.get(pos.clone()).expect("text").to_vec())
        }
    }

    // --------------------------------------------------
    fn get_suffix(&mut self, pos: usize) -> Option<T> {
        if self.suffix_array_mem.is_empty() {
            self.suffix_array_file.get(pos)
        } else {
            self.suffix_array_mem.get(pos).copied()
        }
    }
}
