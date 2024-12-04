use crate::{
    types::{FromUsize, Int, OUTFILE_VERSION, SENTINEL_CHARACTER},
    util::{
        seed_mask_difference, seed_mask_positions, slice_u8_to_vec, usize_to_bytes,
        valid_seed_mask, vec_to_slice_u8,
    },
};
use anyhow::{anyhow, bail, Result};
use log::info;
use rand::Rng;
use rayon::prelude::*;
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
#[derive(Debug)]
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
                if !valid_seed_mask(&mask) {
                    bail!("Invalid mask: {mask}");
                }
                let nums: Vec<u8> = mask
                    .as_bytes()
                    .iter()
                    .flat_map(|b| match b {
                        b'1' => Some(1),
                        b'0' => Some(0),
                        _ => None,
                    })
                    .collect();
                let mask_pos = seed_mask_positions(&nums);
                if mask_pos.is_empty() {
                    bail!("Seed mask must contain at least one 1");
                }
                (nums, mask_pos)
            }
            None => (vec![], vec![]),
        };

        // Having a seed mask implies max_query_len
        let max_query_len = if seed_mask.is_empty() {
            args.max_query_len.map_or(T::default(), T::from_usize)
        } else {
            // Number of 1s/care positions
            T::from_usize(seed_mask.iter().sum())
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
        let start1 = start1.to_usize();
        let start2 = start2.to_usize();
        let seed_diff: Vec<usize> = seed_mask_difference(&self.seed_mask);
        let text_len = self.text_len.to_usize();
        let (end1, end2) = if seed_diff.is_empty() {
            let len = len.to_usize();
            (min(start1 + len, text_len), min(start2 + len, text_len))
        } else {
            // Get the actual span of the spaced seeds
            let a_len = self
                .seed_mask
                .iter()
                .map(|offset| start1 + offset)
                .filter(|&v| v < text_len)
                .count();
            let b_len = self
                .seed_mask
                .iter()
                .map(|offset| start2 + offset)
                .filter(|&v| v < text_len)
                .count();
            (start1 + a_len, start2 + b_len)
        };

        unsafe {
            return T::from_usize(
                (start1..end1)
                    .zip(start2..end2)
                    .enumerate()
                    .take_while(|(i, (a, b))| {
                        let add = seed_diff.get(*i).unwrap_or(&0);
                        self.text.get_unchecked(*a + add)
                            == self.text.get_unchecked(*b + add)
                    })
                    .count(),
            );
        }
    }

    // --------------------------------------------------
    #[inline(always)]
    pub fn find_lcp_full_offset(&self, lcp: T) -> T {
        let lcp = lcp.to_usize();
        let mut full_offset = lcp;

        // If there were some found in common
        if full_offset > 0 {
            // Then add the offset from the seed difference
            // AND the *next* offset, if possible
            let seed_diff = seed_mask_difference(&self.seed_mask);
            for val in &[lcp - 1, lcp] {
                if let Some(add) = seed_diff.get(*val) {
                    full_offset += add;
                }
            }
        }

        T::from_usize(full_offset)
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

            let len_lcp =
                self.find_lcp_full_offset(self.find_lcp(s1, s2, max_query_len));
            let len_lcp = len_lcp.to_usize();

            if len_lcp == max_query_len.to_usize() {
                // The strings are equal
                false
            } else {
                // Look at the next character
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
    }

    // --------------------------------------------------
    pub fn upper_bound(&self, target: T, pivots: &[T]) -> usize {
        // If pivots is empty (no partitions) or
        // If target is less than the first element
        if pivots.first().map_or(false, |&v| self.is_less(target, v)) {
            0
        } else if pivots.last().map_or(false, |&v| !self.is_less(target, v)) {
            // If target is greater than the last element
            pivots.len()
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
                    let full_len_lcp = self.find_lcp_full_offset(len_lcp);

                    // If the len of the LCP is the entire shorter
                    // sequence, take that.
                    if full_len_lcp >= max_n {
                        target_sa[idx_target] = max(x[idx_x], y[idx_y])
                    }
                    // Else, look at the next char after the LCP
                    // to determine order.
                    else if self.text[(x[idx_x] + full_len_lcp).to_usize()]
                        < self.text[(y[idx_y] + full_len_lcp).to_usize()]
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
    use pretty_assertions::assert_eq;

    #[test]
    fn test_is_less() -> Result<()> {
        //           012345
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

        // 1: TTAGC
        // 0: TTTAGC
        assert!(sufr.is_less(1, 0));

        // 0: TTTAGC
        // 1: TTAGC
        assert!(!sufr.is_less(0, 1));

        // 2: TAGC
        // 3: AGC
        assert!(!sufr.is_less(2, 3));

        // 3: AGC
        // 0: TTTAGC
        assert!(sufr.is_less(3, 0));

        Ok(())
    }

    #[test]
    fn test_is_less_max_query_len() -> Result<()> {
        //           012345
        let text = b"TTTAGC".to_vec();
        let args = SufrBuilderArgs {
            text,
            max_query_len: Some(2),
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

        // 1: TTAGC
        // 0: TTTAGC
        // This is true w/o MQL 2 but here they are equal
        // ("TT" == "TT")
        assert!(!sufr.is_less(1, 0));

        // 0: TTTAGC
        // 1: TTAGC
        // ("TT" == "TT")
        assert!(!sufr.is_less(0, 1));

        // 2: TAGC
        // 3: AGC
        assert!(!sufr.is_less(2, 3));

        // 3: AGC
        // 0: TTTAGC
        assert!(sufr.is_less(3, 0));

        Ok(())
    }

    #[test]
    fn test_is_less_seed_mask() -> Result<()> {
        //           012345
        let text = b"TTTTAT".to_vec();
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
            seed_mask: Some("101".to_string()),
        };
        let sufr: SufrBuilder<u32> = SufrBuilder::new(args)?;

        // 0: TTTTAT
        // 1: TTTAT
        // "T-T" vs "T-T"
        assert!(!sufr.is_less(0, 1));

        // 1: TTTAT
        // 0: TTTTAT
        // "T-T" vs "T-T"
        //assert!(!sufr.is_less(1, 0));

        // 0: TTTTAT
        // 3: TAT
        // "T-T" vs "T-T"
        //assert!(!sufr.is_less(0, 3));

        // 3: TAT
        // 0: TTTTAT
        // "T-T" vs "T-T"
        //assert!(!sufr.is_less(3, 0));

        Ok(())
    }

    #[test]
    fn test_find_lcp_no_seed_mask() -> Result<()> {
        //           012345
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

        // 0: TTTAGC
        // 1:  TTAGC
        // 6: len of text
        assert_eq!(sufr.find_lcp(0, 1, 6), 2);

        // 0: TTTAGC
        // 2:   TAGC
        // 6: len of text
        assert_eq!(sufr.find_lcp(0, 2, 6), 1);

        // 0: TTTAGC
        // 1:  TTAGC
        // 1: max query len = 1
        assert_eq!(sufr.find_lcp(0, 1, 1), 1);

        // 0: TTTAGC
        // 3:    AGC
        // 6: len of text
        assert_eq!(sufr.find_lcp(0, 3, 6), 0);

        Ok(())
    }

    #[test]
    fn test_find_lcp_with_seed_mask() -> Result<()> {
        //           012345
        let text = b"TTTTTA".to_vec();
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
            seed_mask: Some("1101".to_string()),
        };
        let sufr: SufrBuilder<u32> = SufrBuilder::new(args)?;

        // 0: TTTTTA
        // 1:  TTTTA
        assert_eq!(sufr.find_lcp(0, 1, 3), 3);

        // 0: TTTTTA
        // 2:   TTTA
        assert_eq!(sufr.find_lcp(0, 2, 3), 2);

        // 0: TTTTTA
        // 5:      A
        assert_eq!(sufr.find_lcp(0, 5, 3), 0);

        Ok(())
    }

    #[test]
    fn test_find_lcp_full_offset() -> Result<()> {
        //           012345
        let text = b"TTTTTA".to_vec();
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
            seed_mask: Some("1101".to_string()),
        };
        let sufr: SufrBuilder<u32> = SufrBuilder::new(args)?;

        // 0: TTTTTA
        // 1:  TTTTA
        // There are only 3 "care" positions but they span 4 bp
        let lcp = sufr.find_lcp(0, 1, 3);
        assert_eq!(lcp, 3);
        assert_eq!(sufr.find_lcp_full_offset(lcp), 4);

        // 0: TTTTTA
        // 3:    TTA
        // There are only two Ts in common, but it should span 3 to the end
        let lcp = sufr.find_lcp(0, 3, 3);
        assert_eq!(lcp, 2);
        assert_eq!(sufr.find_lcp_full_offset(lcp), 3);

        // 0: TTTTTA
        // 4:     TA
        // There is only one T and the As are in a "care" position
        let lcp = sufr.find_lcp(0, 4, 3);
        assert_eq!(lcp, 1);
        assert_eq!(sufr.find_lcp_full_offset(lcp), 1);

        // 0: TTTTTA
        // 5:      A
        // The A is in a "care" position
        let lcp = sufr.find_lcp(0, 5, 3);
        assert_eq!(lcp, 0);
        assert_eq!(sufr.find_lcp_full_offset(lcp), 0);

        Ok(())
    }

    #[test]
    fn test_upper_bound_1() -> Result<()> {
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

        Ok(())
    }

    #[test]
    fn test_upper_bound_2() -> Result<()> {
        //           0123456789
        let text = b"ACGTNNACGT".to_vec();
        let args = SufrBuilderArgs {
            // 10 $
            text,                           //  6 ACGT$
            max_query_len: None,            //  0 ACGTNNACGT$
            is_dna: false,                  //  7 CGT$
            allow_ambiguity: false,         //  1 CGTNNACGT$
            ignore_softmask: false,         //  8 GT$
            sequence_starts: vec![0],       //  2 GTNNACGT$
            headers: vec!["1".to_string()], //  5 NACGT$
            num_partitions: 2,              //  4 NNACGT$
            sequence_delimiter: b'N',       //  9 T$
            seed_mask: None,                //  3 TNNACGT$
        };

        let sufr: SufrBuilder<u64> = SufrBuilder::new(args)?;

        // ACGTNNACGT$ == ACGTNNACGT$
        assert_eq!(sufr.upper_bound(0, &[0]), 1);

        // ACGTNNACGT$ (0) > ACGT$ (6)
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

    #[test]
    fn test_upper_bound_seed_mask() -> Result<()> {
        //           0123456789
        let text = b"ACGTNNACGT".to_vec();
        let args = SufrBuilderArgs {
            text,                               //  6 ACGT$
            max_query_len: None,                //  0 ACGTNNACGT$
            is_dna: false,                      //  7 CGT$
            allow_ambiguity: false,             //  1 CGTNNACGT$
            ignore_softmask: false,             //  8 GT$
            sequence_starts: vec![0],           //  2 GTNNACGT$
            headers: vec!["1".to_string()],     //  4 NNACGT$
            num_partitions: 2,                  //  5 NACGT$
            sequence_delimiter: b'N',           //  9 T$
            seed_mask: Some("101".to_string()), //  3 TNNACGT$
        };
        let sufr: SufrBuilder<u32> = SufrBuilder::new(args)?;

        // two bins
        // ACGTNNACGT$ == ACGTNNACGT$ (A-G)
        // so goes in bin 1
        assert_eq!(sufr.upper_bound(0, &[0]), 1);

        // two bins
        // ACGTNNACGT$ == ACGT$ (A-G)
        // so goes in bin 1
        assert_eq!(sufr.upper_bound(0, &[6]), 1);

        // two bins
        // ACGT$ == ACGTNNACGT$ (A-G)
        // so goes in bin 1
        assert_eq!(sufr.upper_bound(6, &[0]), 1);

        // ACGT$ == ACGT$ (A-G)
        assert_eq!(sufr.upper_bound(6, &[6]), 1);

        // Pivots = [CGT$, GT$]
        // ACGTNNACGT$ < CGT$ => p0
        assert_eq!(sufr.upper_bound(0, &[7, 8]), 0);

        // Pivots = [CGT$, GT$]
        // CGTNNACGT$ == CGT$ (C-T) => p1
        assert_eq!(sufr.upper_bound(1, &[7, 8]), 1);

        // Pivots = [CGT$, GT$]
        // GT$ == GT$ => p2
        assert_eq!(sufr.upper_bound(8, &[7, 8]), 2);

        // Pivots = [CGT$, GT$]
        // T$ > GT$  => p2
        assert_eq!(sufr.upper_bound(9, &[7, 8]), 2);

        // T$ == TNNACGT$ (only compare T) => p1
        assert_eq!(sufr.upper_bound(9, &[3]), 1);

        Ok(())
    }
}
