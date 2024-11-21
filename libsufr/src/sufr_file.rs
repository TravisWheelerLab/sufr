use crate::{
    types::{
        Comparison, FromUsize, Int, LocateResult, LocateResultPosition, SearchOptions,
        SearchResult, SearchResultLocations,
    },
    util::{slice_u8_to_vec, usize_to_bytes},
};
use anyhow::{anyhow, bail, Result};
use home::home_dir;
use log::info;
use rayon::prelude::*;
use std::{
    cmp::{min, Ordering},
    fs::{self, File},
    io::{Read, Seek, SeekFrom, Write},
    mem,
    ops::Range,
    path::{Path, PathBuf},
    slice,
    time::Instant,
};

// --------------------------------------------------
#[derive(Debug)]
pub struct FileAccess<T>
where
    T: Int + FromUsize<T> + Sized + Send + Sync + serde::ser::Serialize,
{
    file: File,
    buffer: Vec<T>,
    buffer_size: usize,
    buffer_pos: usize,
    pub size: usize,
    start_position: u64,
    current_position: u64,
    end_position: u64,
    exhausted: bool,
}

impl<T> FileAccess<T>
where
    T: Int + FromUsize<T> + Sized + Send + Sync + serde::ser::Serialize,
{
    pub fn new(filename: &str, start: u64, num_elements: usize) -> Result<Self> {
        let file = File::open(filename)?;
        let size = num_elements * mem::size_of::<T>();
        Ok(FileAccess {
            file,
            buffer: vec![],
            buffer_size: 2usize.pow(30),
            buffer_pos: 0,
            size,
            start_position: start,
            current_position: start,
            end_position: start + size as u64,
            exhausted: false,
        })
    }

    pub fn reset(&mut self) {
        self.buffer = vec![];
        self.buffer_pos = 0;
        self.current_position = self.start_position;
        self.exhausted = false;
    }

    pub fn iter(&mut self) -> FileAccessIter<T> {
        FileAccessIter { file_access: self }
    }

    // --------------------------------------------------
    // TODO: Ignoring lots of Results to return Option
    pub fn get(&mut self, pos: usize) -> Option<T> {
        // Don't bother looking for something beyond the end
        let seek = self.start_position + (pos * mem::size_of::<T>()) as u64;
        if seek < self.end_position {
            let _ = self.file.seek(SeekFrom::Start(seek));
            let mut buffer: Vec<u8> = vec![0; mem::size_of::<T>()];
            let bytes_read = self.file.read(&mut buffer).unwrap();
            (bytes_read == mem::size_of::<T>()).then(|| {
                let res = unsafe {
                    std::slice::from_raw_parts(buffer.as_ptr() as *const _, 1)
                };
                res[0]
            })
        } else {
            None
        }
    }

    // --------------------------------------------------
    pub fn get_range(&mut self, range: Range<usize>) -> Result<Vec<T>> {
        let start = self.start_position as usize + (range.start * mem::size_of::<T>());
        let end = self.start_position as usize + (range.end * mem::size_of::<T>());
        let valid = self.start_position as usize..self.end_position as usize + 1;
        if valid.contains(&start) && valid.contains(&end) {
            self.file.seek(SeekFrom::Start(start as u64))?;
            let mut buffer: Vec<u8> = vec![0; end - start];
            let bytes_read = self.file.read(&mut buffer)?;
            let num_vals = bytes_read / mem::size_of::<T>();
            Ok(slice_u8_to_vec(&buffer, num_vals))
        } else {
            bail!("Invalid range: {range:?}")
        }
    }
}

// --------------------------------------------------
#[derive(Debug)]
pub struct FileAccessIter<'a, T>
where
    T: Int + FromUsize<T> + Sized + Send + Sync + serde::ser::Serialize,
{
    file_access: &'a mut FileAccess<T>,
}

impl<T> Iterator for FileAccessIter<'_, T>
where
    T: Int + FromUsize<T> + Sized + Send + Sync + serde::ser::Serialize,
{
    type Item = T;

    fn next(&mut self) -> Option<Self::Item> {
        if self.file_access.exhausted {
            None
        } else {
            // Fill the buffer
            if self.file_access.buffer.is_empty()
                || self.file_access.buffer_pos == self.file_access.buffer.len()
            {
                if self.file_access.current_position >= self.file_access.end_position {
                    self.file_access.exhausted = true;
                    return None;
                }

                self.file_access
                    .file
                    .seek(SeekFrom::Start(self.file_access.current_position))
                    .unwrap();

                let bytes_wanted = min(
                    self.file_access.buffer_size * mem::size_of::<T>(),
                    (self.file_access.end_position - self.file_access.current_position)
                        as usize,
                );

                let mut buffer: Vec<u8> = vec![0; bytes_wanted];
                self.file_access.file.read_exact(&mut buffer).unwrap();
                self.file_access.current_position =
                    self.file_access.file.stream_position().unwrap();

                let num_vals = bytes_wanted / mem::size_of::<T>();
                self.file_access.buffer = slice_u8_to_vec(&buffer, num_vals);
                self.file_access.buffer_pos = 0;
            }

            let val = self
                .file_access
                .buffer
                .get(self.file_access.buffer_pos)
                .copied();

            self.file_access.buffer_pos += 1;
            val
        }
    }
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
            num_suffixes: if sa.is_empty() {
                num_suffixes
            } else {
                sa.len()
            },
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
}

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
    pub query_low_memory: bool,
    pub text_pos: usize,
    pub suffix_array_pos: usize,
    pub lcp_pos: usize,
    pub max_query_len: T,
    pub text_len: T,
    pub num_suffixes: T,
    pub num_sequences: T,
    pub sequence_starts: Vec<T>,
    pub headers: Vec<String>,
    pub text: Vec<u8>,
    pub suffix_array_mem: Vec<T>,
    pub suffix_array_mem_mql: Option<usize>,
    pub suffix_array_rank_mem: Vec<usize>,
    pub suffix_array_file: FileAccess<T>,
    pub lcp_file: FileAccess<T>,
}

// --------------------------------------------------
impl<T> SufrFile<T>
where
    T: Int + FromUsize<T> + Sized + Send + Sync,
{
    // Read serialized ".sufr" file
    pub fn read(filename: &str) -> Result<SufrFile<T>> {
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
        let num_suffixes = usize::from_ne_bytes(buffer);

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

        // Text
        let mut text = vec![0; text_len];
        file.read_exact(&mut text)?;

        // Suffix Array
        let suffix_array_file: FileAccess<T> =
            FileAccess::new(filename, suffix_array_pos as u64, num_suffixes)?;
        file.seek_relative(suffix_array_file.size as i64)?;

        // LCP
        let lcp_file: FileAccess<T> =
            FileAccess::new(filename, lcp_pos as u64, num_suffixes)?;
        file.seek_relative(lcp_file.size as i64)?;

        // Headers are variable in length so they are at the end
        let mut buffer = vec![];
        file.read_to_end(&mut buffer)?;
        let headers: Vec<String> = bincode::deserialize(&buffer)?;

        Ok(SufrFile {
            filename: filename.to_string(),
            version,
            is_dna,
            allow_ambiguity,
            ignore_softmask,
            query_low_memory: false,
            text_pos,
            suffix_array_pos,
            lcp_pos,
            text_len: T::from_usize(text_len),
            num_suffixes: T::from_usize(num_suffixes),
            max_query_len,
            num_sequences,
            sequence_starts,
            headers,
            text,
            suffix_array_file,
            lcp_file,
            suffix_array_mem: vec![],
            suffix_array_mem_mql: None,
            suffix_array_rank_mem: vec![],
        })
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
        let mut previous: Option<usize> = None;
        let mut errors: Vec<String> = vec![];
        let text_len = self.text_len.to_usize();
        let num_suffixes = self.num_suffixes.to_usize();

        for i in 0..num_suffixes {
            if i > 0 && i % 1_000_000 == 0 {
                info!("Checked {i}");
            }
            let cur_sa = self.suffix_array_file.get(i).expect("sa").to_usize();
            let cur_lcp = self.lcp_file.get(i).expect("lcp").to_usize();

            if let Some(prev_sa) = previous {
                let check_lcp = self.find_lcp(cur_sa, prev_sa, text_len);
                if check_lcp != cur_lcp {
                    errors.push(format!(
                        "{cur_sa} (r. {i}): LCP {cur_lcp} should be {check_lcp}"
                    ));
                }

                let is_less = match (
                    self.text.get(prev_sa + cur_lcp),
                    self.text.get(cur_sa + cur_lcp),
                ) {
                    (Some(a), Some(b)) => a < b,
                    (None, Some(_)) => true,
                    _ => false,
                };

                if !is_less {
                    errors.push(format!("{cur_sa} (r. {i}): greater than previous"));
                }

                if !errors.is_empty() {
                    dbg!(errors);
                    panic!("blah");
                }
            }
            previous = Some(cur_sa);
        }
        Ok(errors)
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
    pub fn subsample_suffix_array(
        &mut self,
        max_query_len: usize,
    ) -> (Vec<T>, Vec<usize>) {
        let max_query_len = T::from_usize(max_query_len);

        // Ensure we start from the beginning of the SA/LCP files
        self.lcp_file.reset();
        self.suffix_array_file.reset();

        // 37s to process 15/hs1
        let max_len = self.num_suffixes.to_usize();
        let mut suffix_array: Vec<T> = Vec::with_capacity(max_len);
        let mut rank: Vec<usize> = Vec::with_capacity(max_len);

        //for (i, suffix) in self
        //    .lcp_file
        //    .iter()
        //    .zip(self.suffix_array_file.iter())
        //    .enumerate()
        //    .filter_map(|(i, (lcp, suffix))| {
        //        (lcp < max_query_len).then_some((i, suffix))
        //    })
        //{
        //    suffix_array.push(suffix);
        //    rank.push(i);
        //}

        for (i, (lcp, suffix)) in self
            .lcp_file
            .iter()
            .zip(self.suffix_array_file.iter())
            .enumerate()
        {
            if lcp < max_query_len {
                suffix_array.push(suffix);
                //rank.push(T::from_usize(i));
                rank.push(i);
            }
        }

        // 78s to process 15/hs1
        //let ranked_suffixes: Vec<(usize, T)> = self
        //    .lcp_file
        //    .iter()
        //    .zip(self.suffix_array_file.iter())
        //    .enumerate()
        //    .filter_map(|(rank, (lcp, suffix))| {
        //        (lcp < max_query_len).then_some((rank, suffix))
        //    })
        //    .collect();
        //let mut suffix_array: Vec<T> = Vec::with_capacity(ranked_suffixes.len());
        //let mut rank: Vec<usize> = Vec::with_capacity(ranked_suffixes.len());
        //for (suffix_rank, suffix) in ranked_suffixes {
        //    suffix_array.push(suffix);
        //    rank.push(suffix_rank);
        //}

        (suffix_array, rank)
    }

    // --------------------------------------------------
    pub fn set_suffix_array_mem(&mut self, mut max_query_len: usize) -> Result<()> {
        // If ".sufr" file was built with a nonzero (T::default) max_query_len
        // Then this is the value we must use
        if self.max_query_len > T::default() {
            if max_query_len > 0 {
                max_query_len = min(max_query_len, self.max_query_len.to_usize());
            } else {
                max_query_len = self.max_query_len.to_usize();
            }
        }

        if max_query_len == self.max_query_len.to_usize() {
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
                let num_elements = usize::from_ne_bytes(buffer);

                let mut buffer = vec![0; num_elements * mem::size_of::<T>()];
                file.read_exact(&mut buffer)?;
                self.suffix_array_mem = slice_u8_to_vec(&buffer, num_elements);

                let mut buffer = vec![];
                file.read_to_end(&mut buffer)?;
                self.suffix_array_rank_mem = unsafe {
                    std::slice::from_raw_parts(
                        buffer.as_ptr() as *const _,
                        num_elements,
                    )
                    .to_vec()
                };

                info!(
                    "Read compressed SA ({}/{}) from cache file {} in {:?}",
                    self.suffix_array_mem.len(),
                    self.num_suffixes,
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
                    self.num_suffixes,
                    now.elapsed()
                );

                // Write cache file
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

                let bytes = unsafe {
                    slice::from_raw_parts(
                        self.suffix_array_rank_mem.as_ptr() as *const u8,
                        self.suffix_array_rank_mem.len() * std::mem::size_of::<usize>(),
                    )
                };
                file.write_all(bytes)?;

                info!(
                    "Wrote to cache {} in {:?}",
                    cache_path.display(),
                    now.elapsed()
                );
            }
        }

        Ok(())
    }

    // --------------------------------------------------
    pub fn suffix_search(
        &mut self,
        args: SearchOptions,
    ) -> Result<Vec<SearchResult<T>>> {
        self.query_low_memory = args.low_memory;
        // TODO: I don't use this now because I'm always going to disk for the range.
        let _n = if self.query_low_memory {
            self.num_suffixes.to_usize()
        } else {
            let max_query_len =
                args.max_query_len.unwrap_or(self.max_query_len.to_usize());
            self.set_suffix_array_mem(max_query_len)?;
            self.suffix_array_mem.len()
        };

        let res: Vec<_> = args
            .queries
            .clone()
            .into_par_iter()
            .enumerate()
            .map(|(query_num, query)| -> Result<SearchResult<T>> {
                let search_file: FileAccess<T> = FileAccess::new(
                    &self.filename,
                    self.suffix_array_pos as u64,
                    self.num_suffixes.to_usize(),
                )?;
                let mut search = SuffixSearch::new(
                    &self.text,
                    search_file,
                    &self.suffix_array_mem,
                    &self.suffix_array_rank_mem,
                    args.low_memory,
                    self.num_suffixes.to_usize(),
                );
                search.search(query_num, &query, args.find_suffixes)
            })
            .flatten() // TODO: Do I throw away these errors?
            .collect();
        //dbg!(&res);
        //Ok(res.into_iter().flatten().collect())
        Ok(res)
    }

    // --------------------------------------------------
    pub fn locate(&mut self, args: SearchOptions) -> Result<Vec<LocateResult<T>>> {
        let search_result = &self.suffix_search(args)?;
        let seq_starts = self.sequence_starts.clone();
        let seq_names = self.headers.clone();
        let now = Instant::now();
        let mut locate_result: Vec<LocateResult<T>> = vec![];

        // Augment the search with relative sequence positions
        for res in search_result {
            let mut positions = vec![];
            if let Some(locs) = &res.locations {
                for (rank, suffix) in locs.ranks.clone().zip(locs.suffixes.clone()) {
                    let i = seq_starts.partition_point(|&val| val <= suffix) - 1;
                    positions.push(LocateResultPosition {
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

        info!("Search finished in {:?}", now.elapsed());

        Ok(locate_result)
    }
}

// --------------------------------------------------
#[cfg(test)]
mod test {
    use crate::{
        sufr_file::SufrFile,
        types::{LocateResult, LocateResultPosition, SearchOptions},
    };
    use anyhow::Result;

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

        let mut sufr_file: SufrFile<u32> = SufrFile::read("tests/inputs/abba.sufr")?;

        for val in &[true, false] {
            let args = SearchOptions {
                queries: vec!["A".to_string()],
                max_query_len: None,
                low_memory: *val,
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
                        LocateResultPosition {
                            rank: 1,
                            suffix: 0,
                            sequence_name: "1".to_string(),
                            sequence_position: 0,
                        },
                        LocateResultPosition {
                            rank: 2,
                            suffix: 12,
                            sequence_name: "1".to_string(),
                            sequence_position: 12,
                        },
                        LocateResultPosition {
                            rank: 3,
                            suffix: 10,
                            sequence_name: "1".to_string(),
                            sequence_position: 10,
                        },
                        LocateResultPosition {
                            rank: 4,
                            suffix: 1,
                            sequence_name: "1".to_string(),
                            sequence_position: 1,
                        },
                        LocateResultPosition {
                            rank: 5,
                            suffix: 3,
                            sequence_name: "1".to_string(),
                            sequence_position: 3,
                        },
                        LocateResultPosition {
                            rank: 6,
                            suffix: 5,
                            sequence_name: "1".to_string(),
                            sequence_position: 5,
                        },
                        LocateResultPosition {
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
                low_memory: *val,
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
                        LocateResultPosition {
                            rank: 8,
                            suffix: 13,
                            sequence_name: "1".to_string(),
                            sequence_position: 13,
                        },
                        LocateResultPosition {
                            rank: 9,
                            suffix: 11,
                            sequence_name: "1".to_string(),
                            sequence_position: 11,
                        },
                        LocateResultPosition {
                            rank: 10,
                            suffix: 9,
                            sequence_name: "1".to_string(),
                            sequence_position: 9,
                        },
                        LocateResultPosition {
                            rank: 11,
                            suffix: 2,
                            sequence_name: "1".to_string(),
                            sequence_position: 2,
                        },
                        LocateResultPosition {
                            rank: 12,
                            suffix: 4,
                            sequence_name: "1".to_string(),
                            sequence_position: 4,
                        },
                        LocateResultPosition {
                            rank: 13,
                            suffix: 6,
                            sequence_name: "1".to_string(),
                            sequence_position: 6,
                        },
                        LocateResultPosition {
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
                low_memory: *val,
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
                        LocateResultPosition {
                            rank: 3,
                            suffix: 10,
                            sequence_name: "1".to_string(),
                            sequence_position: 10,
                        },
                        LocateResultPosition {
                            rank: 4,
                            suffix: 1,
                            sequence_name: "1".to_string(),
                            sequence_position: 1,
                        },
                        LocateResultPosition {
                            rank: 5,
                            suffix: 3,
                            sequence_name: "1".to_string(),
                            sequence_position: 3,
                        },
                        LocateResultPosition {
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
                low_memory: *val,
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
                    positions: vec![LocateResultPosition {
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
                low_memory: *val,
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
        let input_file = "tests/inputs/abba.sufr";
        let mut sufr_file: SufrFile<u32> = SufrFile::read(input_file)?;
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
    //    let sufr_file: SufrFile<u32> = SufrFile::read("tests/inputs/abba.sufr")?;
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
