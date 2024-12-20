use crate::{
    file_access::FileAccess,
    types::{
        Comparison, FromUsize, Int, LowMemoryUsage, SearchResult,
        SearchResultLocations, SuffixSortType,
    },
    util::find_lcp_full_offset,
};
use anyhow::Result;
use std::{
    cmp::{min, Ordering},
    ops::Range,
};

// --------------------------------------------------
#[derive(Debug)]
pub struct SufrSearchArgs<'a, T>
where
    T: Int + FromUsize<T> + Sized + Send + Sync + serde::ser::Serialize,
{
    pub text: &'a [u8],
    pub text_len: usize,
    pub text_file: FileAccess<u8>,
    pub file: FileAccess<T>,
    pub suffix_array: &'a [T],
    pub rank: &'a [T],
    pub query_low_memory: Option<LowMemoryUsage>,
    pub len_suffixes: usize,
    pub sort_type: &'a SuffixSortType,
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
    query_low_memory: Option<LowMemoryUsage>,
    len_suffixes: usize,
    sort_type: &'a SuffixSortType,
    max_query_len: Option<usize>,
}

// --------------------------------------------------
impl<'a, T> SufrSearch<'a, T>
where
    T: Int + FromUsize<T> + Sized + Send + Sync + serde::ser::Serialize,
{
    pub fn new(args: SufrSearchArgs<'a, T>) -> SufrSearch<'a, T> {
        SufrSearch {
            text: args.text,
            text_len: args.text_len,
            text_file: args.text_file,
            suffix_array_file: args.file,
            suffix_array_mem: args.suffix_array,
            suffix_array_rank_mem: args.rank,
            query_low_memory: args.query_low_memory,
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
    pub fn compare(
        &mut self,
        query: &[u8],
        suffix_pos: usize,
        skip: usize,
    ) -> Comparison {
        //let end = min(suffix_pos + query.len(), self.text_len);
        //eprintln!(
        //    "\nskip {skip} query {:?} suffix {suffix_pos} = {:?} \
        //    max_query_len {:?}",
        //    String::from_utf8(query.to_vec()),
        //    String::from_utf8(
        //        self.get_text_range(suffix_pos..end)
        //            .expect("OK")
        //            .to_vec()
        //    ),
        //    self.max_query_len,
        //);
        //eprintln!("sort_type {:?} mem {:?}", self.sort_type, self.query_low_memory);

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
                                self.text
                                    .get(suffix_pos + offset)
                                    .map_or(false, |&val| query[*offset] == val)
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
            //println!(
            //    "full_offset {full_offset} next query {:?} next text {:?} = {res:?}",
            //    query.get(full_offset).map(|v| *v as char),
            //    self.text.get(suffix_pos + full_offset).map(|v| *v as char),
            //);
            //res
            match (
                query.get(full_offset),
                //self.text.get(suffix_pos + full_offset),
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
        match self.query_low_memory {
            Some(LowMemoryUsage::VeryLow) => self.text_file.get(pos),
            _ => self.text.get(pos).copied(),
        }
    }

    // --------------------------------------------------
    pub fn get_text_range(&mut self, pos: Range<usize>) -> Result<Vec<u8>> {
        match self.query_low_memory {
            // this is too expensive, copies loooooong stretches of text into vec
            Some(LowMemoryUsage::VeryLow) => self.text_file.get_range(pos),
            _ => Ok(self.text.get(pos).expect("foo").to_vec()),
        }
    }

    // --------------------------------------------------
    fn get_suffix(&mut self, pos: usize) -> Option<T> {
        if self.query_low_memory.is_some() {
            self.suffix_array_file.get(pos)
        } else {
            self.suffix_array_mem.get(pos).copied()
        }
    }
}
