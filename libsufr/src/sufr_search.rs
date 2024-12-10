use crate::{
    file_access::FileAccess,
    types::{Comparison, FromUsize, Int, SearchResult, SearchResultLocations},
    util::{find_lcp_full_offset, seed_mask_difference, seed_mask_positions},
};
use anyhow::Result;
use std::cmp::{min, Ordering};

// --------------------------------------------------
#[derive(Debug)]
pub struct SufrSearchArgs<'a, T>
where
    T: Int + FromUsize<T> + Sized + Send + Sync + serde::ser::Serialize,
{
    pub text: &'a [u8],
    pub file: FileAccess<T>,
    pub suffix_array: &'a [T],
    pub rank: &'a [usize],
    pub query_low_memory: bool,
    pub max_query_len: usize,
    pub num_suffixes: usize,
    pub seed_mask: Vec<u8>,
}

// --------------------------------------------------
#[derive(Debug)]
pub struct SufrSearch<'a, T>
where
    T: Int + FromUsize<T> + Sized + Send + Sync + serde::ser::Serialize,
{
    text: &'a [u8],
    suffix_array_file: FileAccess<T>,
    suffix_array_mem: &'a [T],
    suffix_array_rank_mem: &'a [usize],
    query_low_memory: bool,
    max_query_len: usize,
    num_suffixes: usize,
    //seed_mask: Vec<u8>,
    seed_mask_pos: Vec<usize>,
    seed_mask_diff: Vec<usize>,
}

// --------------------------------------------------
impl<'a, T> SufrSearch<'a, T>
where
    T: Int + FromUsize<T> + Sized + Send + Sync + serde::ser::Serialize,
{
    pub fn new(args: SufrSearchArgs<'a, T>) -> SufrSearch<'a, T> {
        let seed_mask_pos = seed_mask_positions(&args.seed_mask);
        let seed_mask_diff = seed_mask_difference(&seed_mask_pos);
        SufrSearch {
            text: args.text,
            suffix_array_file: args.file,
            suffix_array_mem: args.suffix_array,
            suffix_array_rank_mem: args.rank,
            query_low_memory: args.query_low_memory,
            max_query_len: args.max_query_len,
            num_suffixes: if args.suffix_array.is_empty() {
                args.num_suffixes
            } else {
                args.suffix_array.len()
            },
            //seed_mask: args.seed_mask.clone(),
            seed_mask_pos,
            seed_mask_diff,
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
    pub fn compare(&self, query: &[u8], suffix_pos: usize, skip: usize) -> Comparison {
        //println!(
        //    "\nskip {skip} query {:?} suffix {suffix_pos} = {:?} seed_mask_pos = {:?} mql {}",
        //    String::from_utf8(query.to_vec()),
        //    String::from_utf8(self.text.get(suffix_pos..suffix_pos + query.len()).expect("OK").to_vec()),
        //    self.seed_mask_pos,
        //    self.max_query_len
        //);

        let lcp = if (self.max_query_len > 0) && (skip >= self.max_query_len) {
            // If we've already seen enough
            skip
        } else if self.seed_mask_pos.is_empty() {
            let text_start = suffix_pos + skip;
            let text_end = if self.max_query_len > 0 {
                min(self.text.len(), text_start + self.max_query_len)
            } else {
                self.text.len()
            };

            query
                .iter()
                .skip(skip)
                .zip(self.text.get(text_start..text_end).unwrap())
                .map_while(|(a, b)| (a == b).then_some(a))
                .count()
                + skip
        } else {
            let mask_end = if self.max_query_len > 0 {
                self.max_query_len
            } else {
                self.seed_mask_pos.len()
            };
            let mask_pos = &self.seed_mask_pos[skip..mask_end];
            let query_len = mask_pos.iter().filter(|&v| v < &query.len()).count();
            let suffix_len = mask_pos
                .iter()
                .map(|offset| suffix_pos + offset) // add offset to start of suffix
                .filter(|&v| v < self.text.len()) // don't go off end of text
                .count();
            //let query_len = self
            //    .seed_mask_pos
            //    .iter()
            //    .skip(skip)
            //    .filter(|&v| v < &query.len())
            //    .count();
            //let suffix_len = self
            //    .seed_mask_pos
            //    .iter()
            //    .skip(skip)
            //    .map(|offset| suffix_pos + offset) // add offset to start of suffix
            //    .filter(|&v| v < self.text.len()) // don't go off end of text
            //    .count();
            let len = min(query_len, suffix_len);
            //println!("query_len {query_len} suffix_len {suffix_len} actual_len {len}");
            if len > 0 {
                skip + self.seed_mask_pos[skip..skip + len]
                    .iter()
                    .take_while(|&offset| {
                        query[*offset] == self.text[suffix_pos + offset]
                    })
                    .count()
            } else {
                skip
            }
        };
        let full_offset = find_lcp_full_offset(lcp, &self.seed_mask_pos);
        //println!(
        //    "lcp {lcp} full_offset {full_offset} query {:?} text {:?}\n",
        //    query.get(full_offset).map(|&v| v as char),
        //    self.text.get(suffix_pos + full_offset).map(|&v| v as char)
        //);

        let cmp = if (self.max_query_len > 0) && (lcp >= self.max_query_len) {
            // We've seen enough
            Ordering::Equal
        } else {
            // Get the next chars
            match (
                query.get(full_offset),
                self.text.get(suffix_pos + full_offset),
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
    fn get_suffix(&mut self, pos: usize) -> Option<T> {
        if self.query_low_memory {
            self.suffix_array_file.get(pos)
        } else {
            self.suffix_array_mem.get(pos).copied()
        }
    }
}
