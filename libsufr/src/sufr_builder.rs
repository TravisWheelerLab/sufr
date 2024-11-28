use crate::{
    types::{FromUsize, Int, OUTFILE_VERSION, SENTINEL_CHARACTER},
    util::{slice_u8_to_vec, usize_to_bytes, vec_to_slice_u8},
};
use anyhow::{anyhow, bail, Result};
use log::info;
use rand::Rng;
use rayon::prelude::*;
use regex::Regex;
use std::{
    cmp::{max, min, Ordering},
    collections::HashSet,
    fs::{self, File, OpenOptions},
    io::{BufWriter, Seek, SeekFrom, Write},
    mem,
    path::PathBuf,
    sync::{Arc, Mutex},
    time::Instant,
};
use tempfile::NamedTempFile;

// --------------------------------------------------
#[derive(Debug)]
pub struct Partition<T>
where
    T: Int + FromUsize<T> + Sized + Send + Sync + serde::ser::Serialize,
{
    order: usize,
    len: usize,
    first_suffix: T,
    last_suffix: T,
    sa_path: PathBuf,
    lcp_path: PathBuf,
}

// --------------------------------------------------
struct PartitionBuildResult<T>
where
    T: Int + FromUsize<T> + Sized + Send + Sync + serde::ser::Serialize,
{
    builders: Vec<Arc<Mutex<PartitionBuilder<T>>>>,
    num_suffixes: usize,
}

// --------------------------------------------------
#[derive(Debug)]
struct PartitionBuilder<T>
where
    T: Int + FromUsize<T> + Sized + Send + Sync + serde::ser::Serialize,
{
    vals: Vec<T>,
    capacity: usize,
    len: usize,
    total_len: usize,
    path: PathBuf,
}

// --------------------------------------------------
impl<T> PartitionBuilder<T>
where
    T: Int + FromUsize<T> + Sized + Send + Sync + serde::ser::Serialize,
{
    fn new(capacity: usize) -> Result<Self> {
        let tmp = NamedTempFile::new()?;
        let (_, path) = tmp.keep()?;

        Ok(PartitionBuilder {
            vals: vec![T::default(); capacity],
            len: 0,
            total_len: 0,
            capacity,
            path,
        })
    }

    pub fn add(&mut self, val: T) -> Result<()> {
        self.vals[self.len] = val;
        self.len += 1;
        if self.len == self.capacity {
            self.write()?;
            // Reset
            for i in 0..self.len {
                self.vals[i] = T::default();
            }
            self.len = 0;
        }

        Ok(())
    }

    pub fn write(&mut self) -> Result<()> {
        if self.len > 0 {
            let mut file = OpenOptions::new()
                .create(true)
                .append(true)
                .open(&self.path)?;
            file.write_all(vec_to_slice_u8(&self.vals[0..self.len]))?;
            self.total_len += self.len;
        }
        Ok(())
    }
}

// --------------------------------------------------
#[derive(Debug)]
pub struct SufrBuilderArgs {
    pub text: Vec<u8>,
    pub max_query_len: Option<usize>,
    pub is_dna: bool,
    pub allow_ambiguity: bool,
    pub ignore_softmask: bool,
    pub sequence_starts: Vec<usize>,
    pub headers: Vec<String>,
    pub num_partitions: usize,
    pub sequence_delimiter: u8,
    pub seed_mask: Option<String>,
}

// --------------------------------------------------
#[derive(Debug)]
pub struct SufrBuilder<T>
where
    T: Int + FromUsize<T> + Sized + Send + Sync + serde::ser::Serialize,
{
    pub version: u8,
    pub is_dna: bool,
    pub allow_ambiguity: bool,
    pub ignore_softmask: bool,
    pub max_query_len: T,
    pub text_len: T,
    pub num_suffixes: T,
    pub num_sequences: T,
    pub sequence_starts: Vec<T>,
    pub headers: Vec<String>,
    pub text: Vec<u8>,
    pub partitions: Vec<Partition<T>>,
    pub sequence_delimiter: u8,
    pub seed_mask_orig: Vec<u8>,
    pub seed_mask: Vec<usize>,
}

// --------------------------------------------------
impl<T> SufrBuilder<T>
where
    T: Int + FromUsize<T> + Sized + Send + Sync,
{
    pub fn new(args: SufrBuilderArgs) -> Result<SufrBuilder<T>> {
        let text: Vec<_> = args
            .text
            .iter()
            .map(|b| {
                if (97..=122).contains(b) {
                    if args.ignore_softmask {
                        b'N'
                    } else {
                        // only shift lowercase ASCII
                        b & 0b1011111
                    }
                } else {
                    *b
                }
            })
            .collect();
        let text_len = T::from_usize(text.len());

        if args.seed_mask.is_some() && args.max_query_len.is_some() {
            bail!("Cannot use max_query_len and seed_mask together");
        }

        // Validate seed mask before max_query_len
        let (seed_mask_orig, seed_mask): (Vec<u8>, Vec<usize>) = match args.seed_mask {
            Some(mask) => {
                let seed_re = Regex::new("^[01]+$").unwrap();
                if !seed_re.is_match(&mask) {
                    bail!("Invalid mask: {mask}");
                }
                let bytes: Vec<u8> = mask.bytes().collect();
                let mask: Vec<usize> = bytes
                    .iter()
                    .enumerate()
                    .filter_map(|(i, &b)| (b == b'1').then_some(i))
                    .collect();
                if mask.is_empty() {
                    bail!("Seed mask must contain at least one 1");
                }
                (bytes, mask)
            }
            None => (vec![], vec![]),
        };

        // Having a seed mask implies max_query_len
        let max_query_len = if seed_mask.is_empty() {
            args.max_query_len.map_or(T::default(), T::from_usize)
        } else {
            T::from_usize(seed_mask_orig.len())
        };

        let mut sa = SufrBuilder {
            version: OUTFILE_VERSION,
            is_dna: args.is_dna,
            allow_ambiguity: args.allow_ambiguity,
            ignore_softmask: args.ignore_softmask,
            max_query_len,
            text_len,
            num_suffixes: T::default(),
            text,
            num_sequences: T::from_usize(args.sequence_starts.len()),
            sequence_starts: args
                .sequence_starts
                .into_iter()
                .map(T::from_usize)
                .collect::<Vec<_>>(),
            headers: args.headers,
            partitions: vec![],
            sequence_delimiter: args.sequence_delimiter,
            seed_mask_orig,
            seed_mask,
        };
        sa.sort(args.num_partitions)?;
        Ok(sa)
    }

    // --------------------------------------------------
    // Assumes pos is always found -- danger
    pub fn string_at(&self, pos: usize) -> String {
        self.text
            .get(pos..)
            .map(|v| String::from_utf8(v.to_vec()).unwrap())
            .unwrap()
    }

    // --------------------------------------------------
    #[inline(always)]
    pub fn find_lcp(&self, start1: T, start2: T, len: T) -> T {
        //println!("start1 {start1} start2 {start2} len {len}");
        let start1 = start1.to_usize();
        let start2 = start2.to_usize();
        let len = len.to_usize();
        let text_len = self.text_len.to_usize();
        let end1 = min(start1 + len, text_len);
        let end2 = min(start2 + len, text_len);

        // TODO: This is probably bad to collect potentially huge ranges
        let mut range1: Vec<_> = (start1..end1).collect();
        let mut range2: Vec<_> = (start2..end2).collect();

        //println!("BEFORE: range1 {range1:?} range2 {range2:?}");
        //let vals1: Vec<u8> = range1.iter().map(|v| self.text[*v]).collect();
        //let vals2: Vec<u8> = range2.iter().map(|v| self.text[*v]).collect();
        //println!(
        //    "BEFORE: vals1 {} vals2 {}",
        //    String::from_utf8(vals1).unwrap(),
        //    String::from_utf8(vals2).unwrap(),
        //);

        if !self.seed_mask.is_empty() {
            //println!("seed_mask {:?}", self.seed_mask);
            let mut tmp1 = vec![];
            for pos in &self.seed_mask {
                if let Some(val) = range1.get(*pos) {
                    tmp1.push(*val);
                } else {
                    break;
                }
            }
            range1 = tmp1.clone();

            let mut tmp2 = vec![];
            for pos in &self.seed_mask {
                if let Some(val) = range2.get(*pos) {
                    tmp2.push(*val);
                } else {
                    break;
                }
            }
            range2 = tmp2.clone();
            //println!("AFTER: range1 {tmp1:?} range2 {tmp2:?}");
        }
        //let vals1: Vec<u8> = range1.iter().map(|v| self.text[*v]).collect();
        //let vals2: Vec<u8> = range2.iter().map(|v| self.text[*v]).collect();

        // Original
        //unsafe {
        //    T::from_usize(
        //        (start1..end1)
        //            .zip(start2..end2)
        //            .take_while(|(a, b)| {
        //                self.text.get_unchecked(*a) == self.text.get_unchecked(*b)
        //            })
        //            .count(),
        //    )
        //}

        // Return count using subset
        unsafe {
            T::from_usize(
                range1
                    .into_iter()
                    .zip(range2.into_iter())
                    .take_while(|(a, b)| {
                        self.text.get_unchecked(*a) == self.text.get_unchecked(*b)
                    })
                    .count(),
            )
        }

        //let count = unsafe {
        //    T::from_usize(
        //        range1
        //            .into_iter()
        //            .zip(range2.into_iter())
        //            .take_while(|(a, b)| {
        //                self.text.get_unchecked(*a) == self.text.get_unchecked(*b)
        //            })
        //            .count(),
        //    )
        //};
        //
        //if self.seed_mask.is_empty() {
        //    count
        //} else {
        //    if count == self.seed_mask.len() {
        //        // Return last element
        //        self.seed_mask[count - 1]
        //    } else {
        //        // Return next mask position minus 1
        //        self.seed_mask[count] - 1
        //    }
        //}
    }

    // --------------------------------------------------
    #[inline(always)]
    pub fn is_less(&self, s1: T, s2: T) -> bool {
        if s1 == s2 {
            false
        } else {
            let max_query_len = if self.max_query_len > T::default() {
                self.max_query_len
            } else {
                self.text_len
            };

            // Make a "find_lcp_full_offset" that calls "find_lcp" and 
            // then uses mask to find full length of match
            let len_lcp = self.find_lcp(s1, s2, max_query_len).to_usize();

            // No need to look at next character
            //if len_lcp == max_query_len {
            //    what to return?
            //    false
            //}

            match (
                self.text.get(s1.to_usize() + len_lcp),
                self.text.get(s2.to_usize() + len_lcp),
            ) {
                (Some(a), Some(b)) => a < b,
                (None, Some(_)) => true,
                _ => false,
            }
        }
    }

    // --------------------------------------------------
    pub fn upper_bound(&self, target: T, pivots: &[T]) -> usize {
        // See if target is less than the first element
        if pivots.is_empty() || self.is_less(target, pivots[0]) {
            0
        } else {
            // Find where all the values are less than target
            let i = pivots.partition_point(|&p| self.is_less(p, target));

            // If the value at the partition is the same as the target
            if pivots.get(i).map_or(false, |&v| v == target) {
                // Then return the next value, which might be out of range
                i + 1
            } else {
                // Else return the partition point
                i
            }
        }
    }

    // --------------------------------------------------
    fn partition(
        &mut self,
        num_partitions: usize,
        //) -> Result<(Vec<Arc<Mutex<PartitionBuilder<T>>>>, usize)> {
    ) -> Result<PartitionBuildResult<T>> {
        // Create more partitions than requested because
        // we can't know how big they will end up being
        let max_partitions = self.text_len.to_usize() / 4;
        let num_partitions = if num_partitions * 10 < max_partitions {
            num_partitions * 10
        } else if num_partitions * 5 < max_partitions {
            num_partitions * 5
        } else if num_partitions * 2 < max_partitions {
            num_partitions * 2
        } else if num_partitions < max_partitions {
            num_partitions
        } else {
            max_partitions
        };

        // Randomly select some pivots
        let now = Instant::now();
        let pivot_sa = self.select_pivots(self.text.len(), num_partitions);
        let num_pivots = pivot_sa.len();
        info!(
            "Selected {num_pivots} pivot{} in {:?}",
            if num_pivots == 1 { "" } else { "s" },
            now.elapsed()
        );

        let capacity = 4096;
        let mut builders: Vec<_> = vec![];
        for _ in 0..num_partitions {
            let builder: PartitionBuilder<T> = PartitionBuilder::new(capacity)?;
            builders.push(Arc::new(Mutex::new(builder)));
        }

        let now = Instant::now();
        self.text
            .par_iter()
            .enumerate()
            .try_for_each(|(i, &val)| -> Result<()> {
                if val == SENTINEL_CHARACTER
                    || !self.is_dna // Allow anything else if not DNA
                    || (b"ACGT".contains(&val) || self.allow_ambiguity)
                {
                    let partition_num = self.upper_bound(T::from_usize(i), &pivot_sa);
                    match builders[partition_num].lock() {
                        Ok(mut guard) => {
                            guard.add(T::from_usize(i))?;
                        }
                        Err(e) => panic!("Failed to lock: {e}"),
                    }
                }
                Ok(())
            })?;

        // Flush out any remaining buffers
        let mut num_suffixes = 0;
        for builder in &builders {
            match builder.lock() {
                Ok(mut val) => {
                    val.write()?;
                    num_suffixes += val.total_len;
                }
                Err(e) => panic!("Failed to lock: {e}"),
            }
        }

        info!(
            "Wrote {num_suffixes} unsorted suffixes to partition{} in {:?}",
            if num_pivots == 1 { "" } else { "s" },
            now.elapsed()
        );

        //Ok((builders, num_suffixes))
        Ok(PartitionBuildResult {
            builders,
            num_suffixes,
        })
    }

    // --------------------------------------------------
    pub fn sort(&mut self, num_partitions: usize) -> Result<()> {
        let mut partition_build = self.partition(num_partitions)?;
        // Be sure to round up to get all the suffixes
        let num_per_partition = (partition_build.num_suffixes as f64
            / num_partitions as f64)
            .ceil() as usize;
        let total_sort_time = Instant::now();
        let mut num_taken = 0;
        let mut partition_inputs = vec![vec![]; num_partitions];

        // We (probably) have many more partitions than we need,
        // so here we accumulate the small partitions from the left
        // stopping when we reach a boundary like 1M/partition.
        // This evens out the workload to sort the partitions.
        #[allow(clippy::needless_range_loop)]
        for partition_num in 0..num_partitions {
            let boundary = num_per_partition * (partition_num + 1);
            while !partition_build.builders.is_empty() {
                let part = partition_build.builders.remove(0);
                match part.lock() {
                    Ok(builder) => {
                        if builder.total_len > 0 {
                            partition_inputs[partition_num]
                                .push((builder.path.clone(), builder.total_len));
                            num_taken += builder.total_len;
                        }
                    }
                    Err(e) => panic!("Can't get partition: {e}"),
                }

                // Let the last partition soak up the rest
                if partition_num < num_partitions - 1 && num_taken > boundary {
                    break;
                }
            }
        }

        // Ensure we got all the suffixes
        if num_taken != partition_build.num_suffixes {
            bail!(
                "Took {num_taken} but needed to take {}",
                partition_build.num_suffixes
            );
        }

        let mut partitions: Vec<Option<Partition<T>>> =
            (0..num_partitions).map(|_| None).collect();

        partitions.par_iter_mut().enumerate().try_for_each(
            |(partition_num, partition)| -> Result<()> {
                // Find the suffixes in this partition
                let mut part_sa = vec![];
                for (path, len) in &partition_inputs[partition_num] {
                    let buffer = fs::read(path)?;
                    let mut part: Vec<T> = slice_u8_to_vec(&buffer, *len);
                    part_sa.append(&mut part);
                }

                let len = part_sa.len();
                if len > 0 {
                    let mut sa_w = part_sa.clone();
                    let mut lcp = vec![T::default(); len];
                    let mut lcp_w = vec![T::default(); len];
                    self.merge_sort(&mut sa_w, &mut part_sa, len, &mut lcp, &mut lcp_w);

                    // Write to disk
                    let mut sa_file = NamedTempFile::new()?;
                    let _ = sa_file.write(vec_to_slice_u8(&part_sa))?;
                    let mut lcp_file = NamedTempFile::new()?;
                    let _ = lcp_file.write(vec_to_slice_u8(&lcp))?;
                    let (_, sa_path) = sa_file.keep()?;
                    let (_, lcp_path) = lcp_file.keep()?;

                    *partition = Some(Partition {
                        order: partition_num,
                        len,
                        first_suffix: *part_sa.first().unwrap(),
                        last_suffix: *part_sa.last().unwrap(),
                        sa_path,
                        lcp_path,
                    });
                }
                Ok(())
            },
        )?;

        // Get rid of None/unwrap Some, put in order
        let mut partitions: Vec<_> = partitions.into_iter().flatten().collect();
        partitions.sort_by_key(|p| p.order);

        let sizes: Vec<_> = partitions.iter().map(|p| p.len).collect();
        let total_size = sizes.iter().sum::<usize>();
        info!(
            "Sorted {total_size} suffixes in {num_partitions} partitions (avg {}) in {:?}",
            total_size / num_partitions,
            total_sort_time.elapsed()
        );
        self.num_suffixes = T::from_usize(sizes.iter().sum());
        self.partitions = partitions;

        Ok(())
    }

    // --------------------------------------------------
    pub fn merge_sort(
        &self,
        x: &mut [T],
        y: &mut [T],
        n: usize,
        lcp: &mut [T],
        lcp_w: &mut [T],
    ) {
        if n == 1 {
            lcp[0] = T::default();
        } else {
            let mid = n / 2;
            self.merge_sort(
                &mut y[..mid],
                &mut x[..mid],
                mid,
                &mut lcp_w[..mid],
                &mut lcp[..mid],
            );

            self.merge_sort(
                &mut y[mid..],
                &mut x[mid..],
                n - mid,
                &mut lcp_w[mid..],
                &mut lcp[mid..],
            );

            self.merge(x, mid, lcp_w, y, lcp);
        }
    }

    // --------------------------------------------------
    fn merge(
        &self,
        suffix_array: &mut [T],
        mid: usize,
        lcp_w: &mut [T],
        target_sa: &mut [T],
        target_lcp: &mut [T],
    ) {
        let (mut x, mut y) = suffix_array.split_at_mut(mid);
        let (mut lcp_x, mut lcp_y) = lcp_w.split_at_mut(mid);
        let mut len_x = x.len();
        let mut len_y = y.len();
        let mut m = T::default(); // Last LCP from left side (x)
        let mut idx_x = 0; // Index into x (left side)
        let mut idx_y = 0; // Index into y (right side)
        let mut idx_target = 0; // Index into target

        while idx_x < len_x && idx_y < len_y {
            let l_x = lcp_x[idx_x];

            match l_x.cmp(&m) {
                Ordering::Greater => {
                    target_sa[idx_target] = x[idx_x];
                    target_lcp[idx_target] = l_x;
                }
                Ordering::Less => {
                    target_sa[idx_target] = y[idx_y];
                    target_lcp[idx_target] = m;
                    m = l_x;
                }
                Ordering::Equal => {
                    // Length of shorter suffix
                    let max_n = self.text_len - max(x[idx_x], y[idx_y]);

                    // Prefix-context length for the suffixes
                    let context = if self.max_query_len > T::default() {
                        min(self.max_query_len, max_n)
                    } else {
                        max_n
                    };

                    // LCP(X_i, Y_j)
                    let len_lcp =
                        m + self.find_lcp(x[idx_x] + m, y[idx_y] + m, context - m);

                    //if !self.seed_mask.is_empty() {
                    //    if len_lcp == self.seed_mask.len() {
                    //        len_lcp = self.seed_mask[len_lcp.to_usize()]
                    //    } else {
                    //        len_lcp = self.seed_mask[len_lcp + 1] - 1
                    //    }
                    //}

                    // If the len of the LCP is the entire shorter
                    // sequence, take that.
                    if len_lcp == max_n {
                        target_sa[idx_target] = max(x[idx_x], y[idx_y])
                    }
                    // Else, look at the next char after the LCP
                    // to determine order.
                    else if self.text[(x[idx_x] + len_lcp).to_usize()]
                        < self.text[(y[idx_y] + len_lcp).to_usize()]
                    {
                        target_sa[idx_target] = x[idx_x]
                    }
                    // Else take from the right
                    else {
                        target_sa[idx_target] = y[idx_y]
                    }

                    // If we took from the right...
                    if target_sa[idx_target] == x[idx_x] {
                        target_lcp[idx_target] = l_x;
                    } else {
                        target_lcp[idx_target] = m
                    }
                    m = len_lcp;
                }
            }

            if target_sa[idx_target] == x[idx_x] {
                idx_x += 1;
            } else {
                idx_y += 1;
                mem::swap(&mut x, &mut y);
                mem::swap(&mut len_x, &mut len_y);
                mem::swap(&mut lcp_x, &mut lcp_y);
                mem::swap(&mut idx_x, &mut idx_y);
            }

            idx_target += 1;
        }

        // Copy rest of the data from X to Z.
        while idx_x < len_x {
            target_sa[idx_target] = x[idx_x];
            target_lcp[idx_target] = lcp_x[idx_x];
            idx_x += 1;
            idx_target += 1;
        }

        // Copy rest of the data from Y to Z.
        if idx_y < len_y {
            target_sa[idx_target] = y[idx_y];
            target_lcp[idx_target] = m;
            idx_y += 1;
            idx_target += 1;

            while idx_y < len_y {
                target_sa[idx_target] = y[idx_y];
                target_lcp[idx_target] = lcp_y[idx_y];
                idx_y += 1;
                idx_target += 1;
            }
        }
    }

    // --------------------------------------------------
    #[inline(always)]
    fn select_pivots(&self, text_len: usize, num_partitions: usize) -> Vec<T> {
        if num_partitions > 1 {
            // Use a HashMap because selecting pivots one-at-a-time
            // can result in duplicates.
            let num_pivots = num_partitions - 1;
            let rng = &mut rand::thread_rng();
            let mut pivot_sa = HashSet::<T>::new();
            loop {
                let pos = rng.gen_range(0..text_len);
                if self.is_dna && !b"ACGT$".contains(&self.text[pos]) {
                    continue;
                }
                let _ = pivot_sa.insert(T::from_usize(pos));
                if pivot_sa.len() == num_pivots {
                    break;
                }
            }
            // Sort the selected pivots
            let mut pivot_sa: Vec<T> = pivot_sa.iter().cloned().collect();
            let mut sa_w = pivot_sa.clone();
            let len = pivot_sa.len();
            let mut lcp = vec![T::default(); len];
            let mut lcp_w = vec![T::default(); len];
            self.merge_sort(&mut sa_w, &mut pivot_sa, len, &mut lcp, &mut lcp_w);
            pivot_sa
        } else {
            vec![]
        }
    }

    // --------------------------------------------------
    // Serialize data to a ".sufr" file
    pub fn write(&self, filename: &str) -> Result<usize> {
        let mut file = BufWriter::new(
            File::create(filename).map_err(|e| anyhow!("{filename}: {e}"))?,
        );

        let mut bytes_out: usize = 0;

        // Various metadata
        let is_dna: u8 = if self.is_dna { 1 } else { 0 };
        let allow_ambiguity: u8 = if self.allow_ambiguity { 1 } else { 0 };
        let ignore_softmask: u8 = if self.ignore_softmask { 1 } else { 0 };
        bytes_out +=
            file.write(&[OUTFILE_VERSION, is_dna, allow_ambiguity, ignore_softmask])?;

        // Text length
        bytes_out += file.write(&usize_to_bytes(self.text_len.to_usize()))?;

        // Locations of text, suffix array, and LCP
        // Will be corrected at the end
        let locs_pos = file.stream_position()?;
        bytes_out += file.write(&usize_to_bytes(0usize))?;
        bytes_out += file.write(&usize_to_bytes(0usize))?;
        bytes_out += file.write(&usize_to_bytes(0usize))?;

        // Number of suffixes
        bytes_out += file.write(&usize_to_bytes(self.num_suffixes.to_usize()))?;

        // Max query length
        bytes_out += file.write(&usize_to_bytes(self.max_query_len.to_usize()))?;

        // Number of sequences
        bytes_out += file.write(&usize_to_bytes(self.sequence_starts.len()))?;

        // Sequence starts
        bytes_out += file.write(vec_to_slice_u8(&self.sequence_starts))?;

        // Seed mask
        let seed_mask_len = self.seed_mask_orig.len();
        bytes_out += file.write(&usize_to_bytes(seed_mask_len))?;
        if seed_mask_len > 0 {
            file.write_all(&self.seed_mask_orig)?;
            bytes_out += seed_mask_len;
        }

        // Text
        let text_pos = bytes_out;
        file.write_all(&self.text)?;
        bytes_out += self.text.len();

        // Stitch partitioned suffix files together
        let sa_pos = bytes_out;
        for partition in &self.partitions {
            let buffer = fs::read(&partition.sa_path)?;
            bytes_out += &buffer.len();
            file.write_all(&buffer)?;
        }

        // Stitch partitioned LCP files together
        let lcp_pos = bytes_out;
        for (i, partition) in self.partitions.iter().enumerate() {
            let buffer = fs::read(&partition.lcp_path)?;
            bytes_out += &buffer.len();

            if i == 0 {
                file.write_all(&buffer)?;
            } else {
                // Fix LCP boundary
                let mut lcp: Vec<T> = slice_u8_to_vec(&buffer, partition.len);
                if let Some(val) = lcp.first_mut() {
                    *val = self.find_lcp(
                        self.partitions[i - 1].last_suffix,
                        partition.first_suffix,
                        self.text_len,
                    );
                }
                file.write_all(vec_to_slice_u8(&lcp))?;
            }
        }

        // Headers are variable in length so they are at the end
        bytes_out += file.write(&bincode::serialize(&self.headers)?)?;

        // Go back to header and record the locations
        file.seek(SeekFrom::Start(locs_pos))?;
        let _ = file.write(&usize_to_bytes(text_pos))?;
        let _ = file.write(&usize_to_bytes(sa_pos))?;
        let _ = file.write(&usize_to_bytes(lcp_pos))?;

        Ok(bytes_out)
    }
}

// --------------------------------------------------
#[cfg(test)]
mod test {
    use super::{SufrBuilder, SufrBuilderArgs};
    use anyhow::Result;

    #[test]
    fn test_seed_mask() -> Result<()> {
        let text = b"TTTAGC".to_vec();
        let args = SufrBuilderArgs {
            text: text.clone(),
            max_query_len: None,
            is_dna: false,
            allow_ambiguity: false,
            ignore_softmask: false,
            sequence_starts: vec![0],
            headers: vec!["1".to_string()],
            num_partitions: 2,
            sequence_delimiter: b'N',
            seed_mask: Some("1234".to_string()),
        };
        let res: Result<SufrBuilder<u32>> = SufrBuilder::new(args);
        assert!(res.is_err());
        assert_eq!(res.unwrap_err().to_string(), "Invalid mask: 1234");

        let args = SufrBuilderArgs {
            text: text.clone(),
            max_query_len: None,
            is_dna: false,
            allow_ambiguity: false,
            ignore_softmask: false,
            sequence_starts: vec![0],
            headers: vec!["1".to_string()],
            num_partitions: 2,
            sequence_delimiter: b'N',
            seed_mask: Some("0000".to_string()),
        };
        let res: Result<SufrBuilder<u32>> = SufrBuilder::new(args);
        assert!(res.is_err());
        assert_eq!(
            res.unwrap_err().to_string(),
            "Seed mask must contain at least one 1"
        );

        let args = SufrBuilderArgs {
            text: text.clone(),
            max_query_len: None,
            is_dna: false,
            allow_ambiguity: false,
            ignore_softmask: false,
            sequence_starts: vec![0],
            headers: vec!["1".to_string()],
            num_partitions: 2,
            sequence_delimiter: b'N',
            seed_mask: Some("1010".to_string()),
        };
        let res: Result<SufrBuilder<u32>> = SufrBuilder::new(args);
        assert!(res.is_ok());
        let res = res.unwrap();
        assert_eq!(res.seed_mask_orig, b"1010");
        assert_eq!(res.max_query_len, 4);

        Ok(())
    }

    #[test]
    fn test_upper_bound() -> Result<()> {
        //          012345
        let text = b"TTTAGC".to_vec();
        let args = SufrBuilderArgs {
            text,
            max_query_len: None,
            is_dna: false,
            allow_ambiguity: false,
            ignore_softmask: false,
            sequence_starts: vec![0],
            headers: vec!["1".to_string()],
            num_partitions: 2,
            sequence_delimiter: b'N',
            seed_mask: None,
        };
        let sufr: SufrBuilder<u32> = SufrBuilder::new(args)?;

        // The suffix "AGC$" is found before "GC$" and "C$
        assert_eq!(sufr.upper_bound(3, &[5, 4]), 0);

        // The suffix "TAGC$" is beyond all the values
        assert_eq!(sufr.upper_bound(2, &[3, 5, 4]), 3);

        // The "C$" is the last value
        assert_eq!(sufr.upper_bound(5, &[3, 5, 4]), 2);

        //           0123456789
        let text = b"ACGTNNACGT".to_vec();
        //let text_len = text.len();
        let args = SufrBuilderArgs {
            text,
            max_query_len: None,
            is_dna: false,
            allow_ambiguity: false,
            ignore_softmask: false,
            sequence_starts: vec![0],
            headers: vec!["1".to_string()],
            num_partitions: 2,
            sequence_delimiter: b'N',
            seed_mask: None,
        };
        let sufr: SufrBuilder<u32> = SufrBuilder::new(args)?;

        // ACGTNNACGT$ == ACGTNNACGT$
        assert_eq!(sufr.upper_bound(0, &[0]), 1);

        // ACGTNNACGT$ > ACGT$
        assert_eq!(sufr.upper_bound(0, &[6]), 1);

        // ACGT$ < ACGTNNACGT$
        assert_eq!(sufr.upper_bound(6, &[0]), 0);

        // ACGT$ == ACGT$
        assert_eq!(sufr.upper_bound(6, &[6]), 1);

        // Pivots = [CGT$, GT$]
        // ACGTNNACGT$ < CGT$ => p0
        assert_eq!(sufr.upper_bound(0, &[7, 8]), 0);

        // CGTNNACGT$ > CGT$  => p1
        assert_eq!(sufr.upper_bound(1, &[7, 8]), 1);

        // GT$ == GT$  => p1
        assert_eq!(sufr.upper_bound(1, &[7, 8]), 1);

        // T$ > GT$  => p2
        assert_eq!(sufr.upper_bound(9, &[7, 8]), 2);

        // T$ < TNNACGT$ => p0
        assert_eq!(sufr.upper_bound(9, &[3]), 0);

        Ok(())
    }
}
