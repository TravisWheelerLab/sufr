use crate::{
    file_access::FileAccess,
    sufr_search::{SufrSearch, SufrSearchArgs},
    types::{
        ExtractOptions, ExtractResult, ExtractSequence, FromUsize, Int, LocatePosition,
        LocateResult, LowMemoryUsage, SearchOptions, SearchResult, SeedMask,
        SuffixSortType,
    },
    util::{slice_u8_to_vec, usize_to_bytes},
};
use anyhow::{anyhow, Result};
use home::home_dir;
use log::info;
use rayon::prelude::*;
use std::{
    cell::RefCell,
    cmp::min,
    fs::{self, File},
    io::{Read, Seek, Write},
    mem,
    ops::Range,
    path::{Path, PathBuf},
    slice,
    time::Instant,
};
use thread_local::ThreadLocal;

// --------------------------------------------------
#[derive(Debug)]
pub struct SufrFile<T>
where
    T: Int + FromUsize<T> + Sized + Send + Sync + serde::ser::Serialize,
{
    pub filename: String,
    pub version: u8,
    pub is_dna: bool,
    pub allow_ambiguity: bool,
    pub ignore_softmask: bool,
    pub query_low_memory: Option<LowMemoryUsage>,
    pub text_pos: usize,
    pub suffix_array_pos: usize,
    pub lcp_pos: usize,
    pub sort_type: SuffixSortType,
    pub text_len: T,
    pub len_suffixes: T,
    pub num_sequences: T,
    pub sequence_starts: Vec<T>,
    pub headers: Vec<String>,
    pub text: Vec<u8>,
    pub suffix_array_mem: Vec<T>,
    pub suffix_array_mem_mql: Option<usize>,
    pub suffix_array_rank_mem: Vec<T>,
    pub text_file: FileAccess<u8>,
    pub suffix_array_file: FileAccess<T>,
    pub lcp_file: FileAccess<T>,
}

// --------------------------------------------------
impl<T> SufrFile<T>
where
    T: Int + FromUsize<T> + Sized + Send + Sync,
{
    // Read serialized ".sufr" file
    pub fn read(filename: &str, low_memory: bool) -> Result<SufrFile<T>> {
        let mut file = File::open(filename).map_err(|e| anyhow!("{filename}: {e}"))?;

        // Meta
        let mut buffer = [0u8; 4];
        file.read_exact(&mut buffer)?;
        let version = buffer[0];
        let is_dna = buffer[1] == 1;
        let allow_ambiguity = buffer[2] == 1;
        let ignore_softmask = buffer[3] == 1;

        // Length of text
        let mut buffer = [0; 8];
        file.read_exact(&mut buffer)?;
        let text_len = usize::from_ne_bytes(buffer);

        // Position of text
        let mut buffer = [0; 8];
        file.read_exact(&mut buffer)?;
        let text_pos = usize::from_ne_bytes(buffer);

        // Position of suffix array
        let mut buffer = [0; 8];
        file.read_exact(&mut buffer)?;
        let suffix_array_pos = usize::from_ne_bytes(buffer);

        // Position of LCP array
        let mut buffer = [0; 8];
        file.read_exact(&mut buffer)?;
        let lcp_pos = usize::from_ne_bytes(buffer);

        // Number of suffixes
        let mut buffer = [0; 8];
        file.read_exact(&mut buffer)?;
        let len_suffixes = usize::from_ne_bytes(buffer);

        // Max query length
        let mut buffer = [0; 8];
        file.read_exact(&mut buffer)?;
        let max_query_len = T::from_usize(usize::from_ne_bytes(buffer));

        // Number of sequences
        let mut buffer = [0; 8];
        file.read_exact(&mut buffer)?;
        let num_sequences = T::from_usize(usize::from_ne_bytes(buffer));

        // Sequence starts
        let mut buffer = vec![0; num_sequences.to_usize() * mem::size_of::<T>()];
        file.read_exact(&mut buffer)?;
        let sequence_starts: Vec<T> =
            slice_u8_to_vec(&buffer, num_sequences.to_usize());

        // Seed mask len
        let mut buffer = [0; 8];
        file.read_exact(&mut buffer)?;
        let seed_mask_len = usize::from_ne_bytes(buffer);

        // Seed mask: stored as u8 for the 1s/0s
        let seed_mask: Vec<u8> = if seed_mask_len > 0 {
            let mut buffer = vec![0; seed_mask_len];
            file.read_exact(&mut buffer)?;
            buffer
        } else {
            vec![]
        };

        // Text
        let text = if low_memory {
            file.seek_relative(text_len as i64)?;
            vec![]
        } else {
            let mut tmp = vec![0; text_len];
            file.read_exact(&mut tmp)?;
            tmp
        };

        // Text file access
        let text_file: FileAccess<u8> =
            FileAccess::new(filename, text_pos as u64, text_len)?;

        // Suffix Array
        let suffix_array_file: FileAccess<T> =
            FileAccess::new(filename, suffix_array_pos as u64, len_suffixes)?;
        file.seek_relative(suffix_array_file.size as i64)?;

        // LCP
        let lcp_file: FileAccess<T> =
            FileAccess::new(filename, lcp_pos as u64, len_suffixes)?;
        file.seek_relative(lcp_file.size as i64)?;

        // Headers are variable in length so they are at the end
        let mut buffer = vec![];
        file.read_to_end(&mut buffer)?;
        let headers: Vec<String> = bincode::deserialize(&buffer)?;

        let sort_type = if seed_mask.is_empty() {
            SuffixSortType::MaxQueryLen(max_query_len.to_usize())
        } else {
            let seed_mask = SeedMask::from_bytes(&seed_mask)?;
            SuffixSortType::Mask(seed_mask)
        };

        Ok(SufrFile {
            filename: filename.to_string(),
            version,
            is_dna,
            allow_ambiguity,
            ignore_softmask,
            query_low_memory: None,
            text_pos,
            suffix_array_pos,
            lcp_pos,
            text_len: T::from_usize(text_len),
            len_suffixes: T::from_usize(len_suffixes),
            sort_type,
            num_sequences,
            sequence_starts,
            headers,
            text,
            suffix_array_file,
            lcp_file,
            text_file,
            suffix_array_mem: vec![],
            suffix_array_mem_mql: None,
            suffix_array_rank_mem: vec![],
        })
    }

    // --------------------------------------------------
    pub fn get_text(&mut self, pos: usize) -> Option<u8> {
        match self.query_low_memory {
            Some(LowMemoryUsage::VeryLow) => self.text_file.get(pos),
            _ => self.text.get(pos).copied(),
        }
    }

    // --------------------------------------------------
    pub fn get_text_range(&mut self, pos: Range<usize>) -> Result<Vec<u8>> {
        match self.query_low_memory {
            Some(LowMemoryUsage::VeryLow) => self.text_file.get_range(pos),
            _ => Ok(self.text[pos].to_vec()),
        }
    }

    // --------------------------------------------------
    pub fn find_lcp(&self, start1: usize, start2: usize, len: usize) -> usize {
        let end1 = min(start1 + len, len);
        let end2 = min(start2 + len, len);
        unsafe {
            (start1..end1)
                .zip(start2..end2)
                .take_while(|(a, b)| {
                    self.text.get_unchecked(*a) == self.text.get_unchecked(*b)
                })
                .count()
        }
    }

    // --------------------------------------------------
    pub fn check(&mut self) -> Result<Vec<String>> {
        // TODO: Handle MQL/spaced seeds
        //let mut previous: Option<usize> = None;
        //let mut errors: Vec<String> = vec![];
        //let text_len = self.text_len.to_usize();
        //let len_suffixes = self.len_suffixes.to_usize();

        //for i in 0..len_suffixes {
        //    if i > 0 && i % 1_000_000 == 0 {
        //        info!("Checked {i}");
        //    }
        //    let cur_sa = self.suffix_array_file.get(i).expect("sa").to_usize();
        //    let cur_lcp = self.lcp_file.get(i).expect("lcp").to_usize();
        //
        //    if let Some(prev_sa) = previous {
        //        println!("Check prev {prev_sa} cur {cur_sa}");
        //        let check_lcp = self.find_lcp(cur_sa, prev_sa, text_len);
        //        if check_lcp != cur_lcp {
        //            errors.push(format!(
        //                "{cur_sa} (r. {i}): LCP {cur_lcp} should be {check_lcp}"
        //            ));
        //        }
        //
        //        let is_less = match (
        //            self.text.get(prev_sa + cur_lcp),
        //            self.text.get(cur_sa + cur_lcp),
        //        ) {
        //            (Some(a), Some(b)) => a < b,
        //            (None, Some(_)) => true,
        //            _ => false,
        //        };
        //
        //        if !is_less {
        //            errors.push(format!("{cur_sa} (r. {i}): greater than previous"));
        //        }
        //
        //        if !errors.is_empty() {
        //            dbg!(errors);
        //            panic!("blah");
        //        }
        //    }
        //    previous = Some(cur_sa);
        //}
        //Ok(errors)
        Ok(vec![])
    }

    // --------------------------------------------------
    pub fn string_at(&self, pos: usize, len: Option<usize>) -> String {
        let text_len = self.text_len.to_usize();
        let end = len.map_or(text_len, |n| {
            let end = pos + n;
            if end > text_len {
                text_len
            } else {
                end
            }
        });
        self.text
            .get(pos..end)
            .map(|v| String::from_utf8(v.to_vec()).unwrap())
            .unwrap()
    }

    // --------------------------------------------------
    fn get_sufr_dir(&self) -> Result<PathBuf> {
        let home = home_dir().expect("Failed to get home directory");
        let sufr_dir = home.join(".sufr");
        if !sufr_dir.is_dir() {
            fs::create_dir(&sufr_dir)?;
        }
        Ok(sufr_dir)
    }

    // --------------------------------------------------
    pub fn subsample_suffix_array(&mut self, max_query_len: usize) -> (Vec<T>, Vec<T>) {
        let max_query_len = T::from_usize(max_query_len);

        // Ensure we start from the beginning of the SA/LCP files
        self.lcp_file.reset();
        self.suffix_array_file.reset();

        let max_len = self.len_suffixes.to_usize();
        let mut suffix_array: Vec<T> = Vec::with_capacity(max_len);
        let mut rank: Vec<T> = Vec::with_capacity(max_len);

        for (i, (lcp, suffix)) in self
            .lcp_file
            .iter()
            .zip(self.suffix_array_file.iter())
            .enumerate()
        {
            if lcp < max_query_len {
                suffix_array.push(suffix);
                rank.push(T::from_usize(i));
                //rank.push(i);
            }
        }

        (suffix_array, rank)
    }

    // --------------------------------------------------
    pub fn set_suffix_array_mem(&mut self, max_query_len: Option<usize>) -> Result<()> {
        let mut max_query_len = max_query_len.unwrap_or(0);

        // If ".sufr" file was built with a nonzero max_query_len or seed mask
        // Then this is the value we must use
        let built_max_query_len = match &self.sort_type {
            SuffixSortType::MaxQueryLen(mql) => {
                if mql > &0 {
                    *mql
                } else {
                    self.text_len.to_usize()
                }
            }
            SuffixSortType::Mask(seed_mask) => seed_mask.weight,
        };

        if built_max_query_len > 0 {
            if max_query_len > 0 {
                max_query_len = min(max_query_len, built_max_query_len);
            } else {
                max_query_len = built_max_query_len;
            }
        }

        // The requested MQL matches how the SA was built
        if max_query_len == built_max_query_len {
            // Stuff entire SA into memory
            let now = Instant::now();
            self.suffix_array_file.reset();
            self.suffix_array_mem = self.suffix_array_file.iter().collect();
            info!("Read entire SA from disk in {:?}", now.elapsed());

            // There will be no ranks
            self.suffix_array_rank_mem = vec![];
        } else {
            // Do nothing if we've already loaded the correct SA/MQL
            if !self.suffix_array_mem.is_empty()
                && self
                    .suffix_array_mem_mql
                    .map_or(false, |cur_mql| cur_mql == max_query_len)
            {
                info!("Using existing suffix_array_mem");
                return Ok(());
            }

            info!("Loading suffix_array_mem using max_query_len {max_query_len}");

            let sufr_dir = &self.get_sufr_dir()?;
            let basename = Path::new(&self.filename)
                .file_name()
                .unwrap()
                .to_string_lossy()
                .into_owned();
            let cache_path =
                sufr_dir.join(format!("locate-{max_query_len}-{basename}"));

            // Check for stale cache
            if let Ok(cache_meta) = fs::metadata(&cache_path) {
                let source_meta = fs::metadata(&self.filename)?;
                if let (Ok(source_modified), Ok(cache_modified)) =
                    (source_meta.modified(), cache_meta.modified())
                {
                    if source_modified > cache_modified {
                        info!("Removing stale cache {}", cache_path.display());
                        fs::remove_file(&cache_path)?;
                    }
                }
            }

            if cache_path.is_file() {
                let now = Instant::now();
                let mut file = File::open(&cache_path)
                    .map_err(|e| anyhow!("{}: {e}", cache_path.display()))?;

                let mut buffer = [0; 8];
                file.read_exact(&mut buffer)?;
                let suffix_array_len = usize::from_ne_bytes(buffer);

                self.suffix_array_mem = if suffix_array_len == 0 {
                    vec![]
                } else {
                    let mut buffer = vec![0; suffix_array_len * mem::size_of::<T>()];
                    file.read_exact(&mut buffer)?;
                    slice_u8_to_vec(&buffer, suffix_array_len)
                };

                let mut buffer = vec![];
                file.read_to_end(&mut buffer)?;
                self.suffix_array_rank_mem = if buffer.is_empty() {
                    vec![]
                } else {
                    unsafe {
                        std::slice::from_raw_parts(
                            buffer.as_ptr() as *const _,
                            suffix_array_len,
                        )
                        .to_vec()
                    }
                };

                info!(
                    "Read compressed SA ({}/{}) from cache file {} in {:?}",
                    self.suffix_array_mem.len(),
                    self.len_suffixes,
                    cache_path.display(),
                    now.elapsed()
                );
            } else {
                let now = Instant::now();
                let (sub_sa, sub_rank) = &self.subsample_suffix_array(max_query_len);
                self.suffix_array_mem_mql = Some(max_query_len);
                self.suffix_array_mem = sub_sa.to_vec();
                self.suffix_array_rank_mem = sub_rank.to_vec();

                info!(
                    "Loaded compressed SA ({}/{}) in {:?}",
                    sub_sa.len(),
                    self.len_suffixes,
                    now.elapsed()
                );

                // Write cache file
                if !self.suffix_array_mem.is_empty() {
                    let now = Instant::now();
                    let mut file = File::create(&cache_path)
                        .map_err(|e| anyhow!("{}: {e}", cache_path.display()))?;
                    let _ = file.write(&usize_to_bytes(self.suffix_array_mem.len()))?;
                    let bytes = unsafe {
                        slice::from_raw_parts(
                            self.suffix_array_mem.as_ptr() as *const u8,
                            self.suffix_array_mem.len() * std::mem::size_of::<T>(),
                        )
                    };
                    file.write_all(bytes)?;

                    if !self.suffix_array_rank_mem.is_empty() {
                        let bytes = unsafe {
                            slice::from_raw_parts(
                                self.suffix_array_rank_mem.as_ptr() as *const u8,
                                self.suffix_array_rank_mem.len()
                                    * std::mem::size_of::<usize>(),
                            )
                        };
                        file.write_all(bytes)?;
                    }

                    info!(
                        "Wrote to cache {} in {:?}",
                        cache_path.display(),
                        now.elapsed()
                    );
                }
            }
        }

        Ok(())
    }

    // --------------------------------------------------
    pub fn suffix_search(
        &mut self,
        args: &SearchOptions,
    ) -> Result<Vec<SearchResult<T>>> {
        self.query_low_memory = args.low_memory.clone();

        if self.query_low_memory.is_none() {
            self.set_suffix_array_mem(args.max_query_len)?;
        }

        let now = Instant::now();
        let new_search = || -> Result<RefCell<SufrSearch<T>>> {
            let suffix_array_file: FileAccess<T> = FileAccess::new(
                &self.filename,
                self.suffix_array_pos as u64,
                self.len_suffixes.to_usize(),
            )?;
            let text_file: FileAccess<u8> = FileAccess::new(
                &self.filename,
                self.text_pos as u64,
                self.text_len.to_usize(),
            )?;
            let search_args = SufrSearchArgs {
                text: &self.text,
                text_len: self.text_len.to_usize(),
                text_file,
                file: suffix_array_file,
                suffix_array: &self.suffix_array_mem,
                rank: &self.suffix_array_rank_mem,
                query_low_memory: args.low_memory.clone(),
                len_suffixes: self.len_suffixes.to_usize(),
                sort_type: &self.sort_type,
                max_query_len: args.max_query_len,
            };
            Ok(RefCell::new(SufrSearch::new(search_args)))
        };

        let thread_local_search: ThreadLocal<RefCell<SufrSearch<T>>> =
            ThreadLocal::new();

        let mut res: Vec<_> = args
            .queries
            .clone()
            .into_par_iter()
            .enumerate()
            .flat_map(|(query_num, query)| -> Result<SearchResult<T>> {
                let mut search =
                    thread_local_search.get_or_try(new_search)?.borrow_mut();
                //dbg!(&search);
                search.search(query_num, &query, args.find_suffixes)
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
    pub fn extract(&mut self, args: ExtractOptions) -> Result<Vec<ExtractResult>> {
        let search_args = SearchOptions {
            queries: args.queries,
            max_query_len: args.max_query_len,
            low_memory: args.low_memory,
            find_suffixes: true,
        };
        let search_result = &self.suffix_search(&search_args)?;
        let seq_starts = self.sequence_starts.clone();
        let seq_names = self.headers.clone();
        let text_len = self.text_len.to_usize();
        let mut extract_result: Vec<ExtractResult> = vec![];
        let now = Instant::now();

        // Augment the search with relative sequence positions
        for res in search_result {
            let mut sequences = vec![];
            if let Some(locs) = &res.locations {
                for (rank, suffix) in locs.ranks.clone().zip(locs.suffixes.clone()) {
                    let i = seq_starts.partition_point(|&val| val <= suffix) - 1;
                    let sequence_start = seq_starts[i].to_usize();
                    let seq_end = if i == seq_starts.len() - 1 {
                        text_len
                    } else {
                        seq_starts[i + 1].to_usize()
                    };
                    let suffix = suffix.to_usize();
                    let relative_suffix_start = suffix - sequence_start;
                    let context_start = relative_suffix_start
                        .saturating_sub(args.prefix_len.unwrap_or(0));
                    let context_end = min(
                        args.suffix_len
                            .map_or(seq_end, |len| relative_suffix_start + len),
                        seq_end,
                    );
                    sequences.push(ExtractSequence {
                        rank,
                        suffix,
                        sequence_name: seq_names[i].clone(),
                        sequence_start,
                        sequence_range: (context_start..context_end),
                        suffix_offset: relative_suffix_start - context_start,
                    })
                }
            }

            extract_result.push(ExtractResult {
                query_num: res.query_num,
                query: res.query.clone(),
                sequences,
            });
        }

        info!("Adding locate data finished in {:?}", now.elapsed());

        Ok(extract_result)
    }

    // --------------------------------------------------
    pub fn locate(&mut self, args: SearchOptions) -> Result<Vec<LocateResult<T>>> {
        let search_result = &self.suffix_search(&args)?;
        let seq_starts = self.sequence_starts.clone();
        let seq_names = self.headers.clone();
        let mut locate_result: Vec<LocateResult<T>> = vec![];
        let now = Instant::now();

        // Augment the search with relative sequence positions
        for res in search_result {
            let mut positions = vec![];
            if let Some(locs) = &res.locations {
                for (rank, suffix) in locs.ranks.clone().zip(locs.suffixes.clone()) {
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
mod test {
    use crate::{
        sufr_file::SufrFile,
        types::{
            ExtractOptions, ExtractResult, ExtractSequence, LocatePosition,
            LocateResult, LowMemoryUsage, SearchOptions,
        },
    };
    use anyhow::Result;

    // --------------------------------------------------
    #[test]
    fn test_extract() -> Result<()> {
        // cargo run -- ex data/expected/1.sufr AC GT XX
        let mut sufr_file: SufrFile<u32> =
            SufrFile::read("../data/expected/1.sufr", false)?;

        let opts = ExtractOptions {
            queries: vec!["AC".to_string(), "GT".to_string(), "XX".to_string()],
            max_query_len: None,
            low_memory: None,
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

        let mut sufr_file: SufrFile<u32> =
            SufrFile::read("../data/expected/abba.sufr", false)?;

        for val in &[true, false] {
            let args = SearchOptions {
                queries: vec!["A".to_string()],
                max_query_len: None,
                low_memory: val.then_some(LowMemoryUsage::Low),
                find_suffixes: true,
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

        for val in &[true, false] {
            let args = SearchOptions {
                queries: vec!["B".to_string()],
                max_query_len: None,
                low_memory: val.then_some(LowMemoryUsage::Low),
                find_suffixes: true,
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

        for val in &[true, false] {
            let args = SearchOptions {
                queries: vec!["ABAB".to_string()],
                max_query_len: None,
                low_memory: val.then_some(LowMemoryUsage::Low),
                find_suffixes: true,
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

        for val in &[true, false] {
            let args = SearchOptions {
                queries: vec!["ABABB".to_string()],
                max_query_len: None,
                low_memory: val.then_some(LowMemoryUsage::Low),
                find_suffixes: true,
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

        for val in &[true, false] {
            let args = SearchOptions {
                queries: vec!["BBBB".to_string()],
                max_query_len: None,
                low_memory: val.then_some(LowMemoryUsage::Low),
                find_suffixes: true,
            };
            let res = sufr_file.locate(args);
            assert!(res.is_ok());
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

    // --------------------------------------------------
    #[test]
    fn test_file_access() -> Result<()> {
        let input_file = "../data/expected/abba.sufr";
        let mut sufr_file: SufrFile<u32> = SufrFile::read(input_file, false)?;
        let suf_by_rank = [
            14, //  0: #
            0,  //  1: AABABABABBABAB#
            12, //  2: AB#
            10, //  3: ABAB#
            1,  //  4: ABABABABBABAB#
            3,  //  5: ABABABBABAB#
            5,  //  6: ABABBABAB#
            7,  //  7: ABBABAB#
            13, //  8: B#
            11, //  9: BAB#
            9,  // 10: BABAB#
            2,  // 11: BABABABBABAB#
            4,  // 12: BABABBABAB#
            6,  // 13: BABBABAB#
            8,  // 14: BBABAB#
        ];

        for (rank, &suffix) in suf_by_rank.iter().enumerate() {
            let res = sufr_file.suffix_array_file.get(rank);
            assert!(res.is_some());
            assert_eq!(res.unwrap(), suffix);
        }

        let res = sufr_file.suffix_array_file.get_range(1..100);
        assert!(res.is_err());
        assert_eq!(
            res.as_ref().unwrap_err().to_string(),
            "Invalid range: 1..100"
        );

        let res = sufr_file.suffix_array_file.get_range(8..9);
        assert!(res.is_ok());
        assert_eq!(res.as_ref().unwrap(), &[13]);

        let res = sufr_file.suffix_array_file.get_range(8..13);
        assert!(res.is_ok());
        assert_eq!(res.as_ref().unwrap(), &[13, 11, 9, 2, 4]);

        let res = sufr_file.suffix_array_file.get_range(1..8);
        assert!(res.is_ok());
        assert_eq!(res.as_ref().unwrap(), &[0, 12, 10, 1, 3, 5, 7]);

        let all: Vec<_> = sufr_file.suffix_array_file.iter().collect();
        assert_eq!(all, &[14, 0, 12, 10, 1, 3, 5, 7, 13, 11, 9, 2, 4, 6, 8]);

        for (i, suffix) in sufr_file.suffix_array_file.iter().enumerate() {
            assert_eq!(suf_by_rank[i], suffix);
        }

        Ok(())
    }

    // --------------------------------------------------
    // The "compare" function is now deeply nested inside the SuffixSearch
    // which is created inside the "suffix_search" function and I'm lost
    // how to untangle and test this.
    //#[test]
    //fn test_compare() -> Result<()> {
    //    // 0  1  2  3  4  5  6  7  8  9 10 11 12 13 14
    //    // A  A  B  A  B  A  B  A  B  B  A  B  A  B  #
    //    let sufr_file: SufrFile<u32> = SufrFile::read("../data/inputs/abba.sufr")?;
    //
    //    // Compare to B to B with no skip
    //    let query = "B".as_bytes();
    //    let res = sufr_file.compare(query, 13, 0);
    //    assert_eq!(res.cmp, Ordering::Equal);
    //    assert_eq!(res.lcp, 1);
    //
    //    // Compare to B to B with skip = 1
    //    let query = "B".as_bytes();
    //    let res = sufr_file.compare(query, 13, 1);
    //    assert_eq!(res.cmp, Ordering::Equal);
    //    assert_eq!(res.lcp, 1);
    //
    //    // Compare to B to AB
    //    let query = "B".as_bytes();
    //    let res = sufr_file.compare(query, 12, 0);
    //    assert_eq!(res.cmp, Ordering::Greater);
    //    assert_eq!(res.lcp, 0);
    //
    //    // Compare to ABABA to ABBABAB#
    //    let query = "ABABA".as_bytes();
    //    let res = sufr_file.compare(query, 7, 2);
    //    assert_eq!(res.cmp, Ordering::Less);
    //    assert_eq!(res.lcp, 2);
    //
    //    // Compare to ABAB to ABABBABAB#
    //    let query = "ABABA".as_bytes();
    //    let res = sufr_file.compare(query, 5, 2);
    //    assert_eq!(res.cmp, Ordering::Less);
    //    assert_eq!(res.lcp, 4);
    //
    //    Ok(())
    //}
}
