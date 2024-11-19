use anyhow::{anyhow, bail, Result};
use home::home_dir;
use log::info;
use needletail::parse_fastx_file;
use rand::Rng;
use rayon::prelude::*;
use std::{
    cmp::{max, min, Ordering},
    collections::HashSet,
    fmt::{Debug, Display},
    fs::{self, File, OpenOptions},
    hash::Hash,
    io::{BufWriter, Read, Seek, SeekFrom, Write},
    mem,
    ops::{Add, Div, Range, Sub},
    path::{Path, PathBuf},
    slice,
    sync::{Arc, Mutex},
    time::Instant,
};
use tempfile::NamedTempFile;

const OUTFILE_VERSION: u8 = 4;
const SENTINEL: u8 = b'$';

// --------------------------------------------------
#[derive(Debug)]
pub struct Comparison {
    cmp: Ordering,
    lcp: usize,
}

// --------------------------------------------------
#[derive(Debug)]
pub struct SequenceFileData {
    pub seq: Vec<u8>,
    pub start_positions: Vec<usize>,
    pub headers: Vec<String>,
}

// --------------------------------------------------
pub trait Int:
    Debug
    + Add<Output = Self>
    + Sub<Output = Self>
    + Div<Output = Self>
    + Copy
    + Default
    + Display
    + Ord
    + Hash
    + serde::ser::Serialize
{
    fn to_usize(&self) -> usize;
}

impl Int for u32 {
    fn to_usize(&self) -> usize {
        *self as usize
    }
}

impl Int for u64 {
    fn to_usize(&self) -> usize {
        *self as usize
    }
}

pub trait FromUsize<T> {
    fn from_usize(val: usize) -> T;
}

impl FromUsize<u32> for u32 {
    fn from_usize(val: usize) -> u32 {
        val as u32
    }
}

impl FromUsize<u64> for u64 {
    fn from_usize(val: usize) -> u64 {
        val as u64
    }
}

// --------------------------------------------------
#[derive(Debug)]
pub struct Locate {
    pub queries: Vec<String>,
    pub max_query_len: Option<usize>,
    pub low_memory: bool,
}

// --------------------------------------------------
#[derive(Debug, PartialEq)]
pub struct LocateResult<T>
where
    T: Int + FromUsize<T> + Sized + Send + Sync + serde::ser::Serialize,
{
    pub query: String,
    pub positions: Vec<LocateResultPosition<T>>,
    pub ranks: Range<usize>,
}

// --------------------------------------------------
#[derive(Debug, PartialEq)]
pub struct LocateResultPosition<T>
where
    T: Int + FromUsize<T> + Sized + Send + Sync + serde::ser::Serialize,
{
    pub suffix: T,
    pub sequence_name: String,
    pub sequence_position: T,
}

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

#[derive(Debug)]
pub struct PartitionBuilder<T>
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
            file.write_all(SufrBuilder::vec_to_slice_u8(&self.vals[0..self.len]))?;
            self.total_len += self.len;
        }
        Ok(())
    }
}

struct PartitionBuildResult<T>
where
    T: Int + FromUsize<T> + Sized + Send + Sync + serde::ser::Serialize,
{
    builders: Vec<Arc<Mutex<PartitionBuilder<T>>>>,
    num_suffixes: usize,
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
}

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

        let mut sa = SufrBuilder {
            version: OUTFILE_VERSION,
            is_dna: args.is_dna,
            allow_ambiguity: args.allow_ambiguity,
            ignore_softmask: args.ignore_softmask,
            max_query_len: args.max_query_len.map_or(T::default(), T::from_usize),
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
        let len = len.to_usize();
        let end1 = min(start1 + len, self.text_len.to_usize());
        let end2 = min(start2 + len, self.text_len.to_usize());
        unsafe {
            return T::from_usize(
                (start1..end1)
                    .zip(start2..end2)
                    .take_while(|(a, b)| {
                        self.text.get_unchecked(*a) == self.text.get_unchecked(*b)
                    })
                    .count(),
            );
        }
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
            let len_lcp = self.find_lcp(s1, s2, max_query_len).to_usize();

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
                if val == SENTINEL
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
                    let mut part: Vec<T> = SufrBuilder::slice_u8_to_vec(&buffer, *len);
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
                    let _ = sa_file.write(Self::vec_to_slice_u8(&part_sa))?;
                    let mut lcp_file = NamedTempFile::new()?;
                    let _ = lcp_file.write(Self::vec_to_slice_u8(&lcp))?;
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
        bytes_out += file.write(Self::vec_to_slice_u8(&self.sequence_starts))?;

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
                let mut lcp: Vec<T> = Self::slice_u8_to_vec(&buffer, partition.len);
                if let Some(val) = lcp.first_mut() {
                    *val = self.find_lcp(
                        self.partitions[i - 1].last_suffix,
                        partition.first_suffix,
                        self.text_len,
                    );
                }
                file.write_all(Self::vec_to_slice_u8(&lcp))?;
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

    // --------------------------------------------------
    pub fn vec_to_slice_u8(vec: &[T]) -> &[u8] {
        unsafe {
            slice::from_raw_parts(
                vec.as_ptr() as *const _,
                std::mem::size_of_val(vec), //vec.len() * mem::size_of::<T>(),
            )
        }
    }

    // --------------------------------------------------
    pub fn slice_u8_to_vec(buffer: &[u8], len: usize) -> Vec<T> {
        unsafe { std::slice::from_raw_parts(buffer.as_ptr() as *const _, len).to_vec() }
    }
}

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
    size: usize,
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
        let valid = self.start_position as usize..self.end_position as usize;
        if valid.contains(&start) && valid.contains(&end) {
            self.file.seek(SeekFrom::Start(start as u64))?;
            let mut buffer: Vec<u8> = vec![0; end - start];
            let bytes_read = self.file.read(&mut buffer)?;
            let num_vals = bytes_read / mem::size_of::<T>();
            Ok(SufrBuilder::slice_u8_to_vec(&buffer, num_vals))
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
                self.file_access.buffer =
                    SufrBuilder::slice_u8_to_vec(&buffer, num_vals);
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
            SufrBuilder::slice_u8_to_vec(&buffer, num_sequences.to_usize());

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
                self.suffix_array_mem =
                    SufrBuilder::slice_u8_to_vec(&buffer, num_elements);

                let mut buffer = vec![];
                file.read_to_end(&mut buffer)?;
                //self.suffix_array_rank_mem =
                //    SufrBuilder::slice_u8_to_vec(&buffer, num_elements);
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
    fn get_suffix(&mut self, pos: usize) -> Option<T> {
        if self.query_low_memory {
            self.suffix_array_file.get(pos)
        } else {
            self.suffix_array_mem.get(pos).copied()
        }
    }

    // --------------------------------------------------
    pub fn locate(&mut self, args: Locate) -> Result<Vec<Result<LocateResult<T>>>> {
        self.query_low_memory = args.low_memory;
        let n = if self.query_low_memory {
            self.num_suffixes.to_usize()
        } else {
            let max_query_len =
                args.max_query_len.unwrap_or(self.max_query_len.to_usize());
            self.set_suffix_array_mem(max_query_len)?;
            self.suffix_array_mem.len()
        };
        let seq_starts = self.sequence_starts.clone();
        let seq_names = self.headers.clone();

        let now = Instant::now();
        let res: Vec<_> = args
            .queries
            .into_iter()
            .map(|query| -> Result<LocateResult<T>> {
                let qry = query.as_bytes();

                if let Some(start) = self.suffix_search_first(qry, 0, n - 1, 0, 0) {
                    let end = self
                        .suffix_search_last(qry, start, n - 1, n, 0, 0)
                        .unwrap_or(start);

                    // Rank is empty when we have the full SA in memory
                    // AND when doing low-memory searches
                    let (suffixes, ranks) = if self.suffix_array_rank_mem.is_empty() {
                        let (start_rank, end_rank) = (start, end + 1);
                        // For low-memory, go to disk
                        let suffixes = if self.suffix_array_mem.is_empty() {
                            self.suffix_array_file.get_range(start_rank..end_rank)?
                        } else {
                            // Otherwise, get from memory
                            self.suffix_array_mem[start_rank..end_rank].to_vec()
                        };
                        (suffixes, start_rank..end_rank)
                    } else {
                        // This is the case for the compressed/in-memory SA
                        let start_rank = self.suffix_array_rank_mem[start];
                        let end_rank = if start == end {
                            if start == self.suffix_array_rank_mem.len() - 1 {
                                // We're on the last rank, so go to end
                                self.num_suffixes.to_usize()
                            } else {
                                // Use the next LCP rank
                                self.suffix_array_rank_mem[start + 1]
                            }
                        } else {
                            self.suffix_array_rank_mem[end] + 1
                        };

                        // I have to go to disk to get the actual suffixes
                        let suffixes =
                            self.suffix_array_file.get_range(start_rank..end_rank)?;
                        (suffixes, start_rank..end_rank)
                    };

                    let positions: Vec<_> = suffixes
                        .iter()
                        .map(|&suffix| {
                            let i =
                                seq_starts.partition_point(|&val| val <= suffix) - 1;
                            LocateResultPosition {
                                suffix,
                                sequence_name: seq_names[i].clone(),
                                sequence_position: suffix - seq_starts[i],
                            }
                        })
                        .collect();

                    Ok(LocateResult {
                        query: query.to_string(),
                        positions,
                        ranks,
                    })
                } else {
                    Err(anyhow!("{query}"))
                }
            })
            .collect();

        info!("Search finished in {:?}", now.elapsed());
        Ok(res)
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
}

// --------------------------------------------------
// Utility function to find length of the input text
// determine whether to build u32 or u64
pub fn read_text_length(filename: &str) -> Result<usize> {
    let mut file = File::open(filename).map_err(|e| anyhow!("{filename}: {e}"))?;

    // Meta (version, is_dna)
    let mut buffer = [0; 4];
    file.read_exact(&mut buffer)?;

    let outfile_version = buffer[0];
    if outfile_version == OUTFILE_VERSION {
        // Length of text is the next usize
        let mut buffer = [0; 8];
        file.read_exact(&mut buffer)?;
        Ok(usize::from_ne_bytes(buffer))
    } else {
        bail!("Unknown sufr version {outfile_version}");
    }
}

// --------------------------------------------------
// Utility function to read FASTA/Q file for sequence
// data needed by SufrBuilder
pub fn read_sequence_file(
    filename: &str,
    sequence_delimiter: u8,
) -> Result<SequenceFileData> {
    let mut reader = parse_fastx_file(filename)?;
    let mut seq: Vec<u8> = Vec::with_capacity(u32::MAX as usize);
    let mut headers: Vec<String> = vec![];
    let mut start_positions: Vec<usize> = vec![];
    let mut i = 0;
    while let Some(rec) = reader.next() {
        let rec = rec?;
        if i > 0 {
            seq.push(sequence_delimiter);
        }

        // Record current length as start position
        start_positions.push(seq.len());
        let mut tmp: Vec<u8> = rec.seq().iter().copied().collect();
        seq.append(&mut tmp);
        i += 1;

        headers.push(String::from_utf8(rec.id().to_vec())?);
    }

    // File delimiter
    seq.push(SENTINEL);

    Ok(SequenceFileData {
        seq,
        start_positions,
        headers,
    })
}

// --------------------------------------------------
// Utility function used by SufrBuilder
fn usize_to_bytes(value: usize) -> Vec<u8> {
    // Determine the size of usize in bytes
    let size = std::mem::size_of::<usize>();

    // Create a vector to hold the bytes
    let mut bytes = Vec::with_capacity(size);

    // Convert usize to bytes
    for i in 0..size {
        bytes.push((value >> (i * 8)) as u8);
    }

    bytes
}

// --------------------------------------------------
#[cfg(test)]
mod tests {
    use super::{
        read_sequence_file, read_text_length, usize_to_bytes, Locate, LocateResult,
        LocateResultPosition, SufrBuilder, SufrBuilderArgs, SufrFile,
    };
    use crate::OUTFILE_VERSION;
    use anyhow::Result;
    use std::cmp::Ordering;
    use tempfile::NamedTempFile;

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

    #[test]
    fn test_slice_u8_to_vec() -> Result<()> {
        let res: Vec<u32> = SufrBuilder::slice_u8_to_vec(&[0, 0, 0, 0], 1);
        assert_eq!(res, &[0u32]);

        let res: Vec<u64> = SufrBuilder::slice_u8_to_vec(&[0, 0, 0, 0, 0, 0, 0, 0], 1);
        assert_eq!(res, &[0u64]);

        let res: Vec<u32> = SufrBuilder::slice_u8_to_vec(&[1, 0, 0, 0], 1);
        assert_eq!(res, &[1u32]);

        let res: Vec<u64> = SufrBuilder::slice_u8_to_vec(&[1, 0, 0, 0, 0, 0, 0, 0], 1);
        assert_eq!(res, &[1u64]);

        let res: Vec<u32> = SufrBuilder::slice_u8_to_vec(&[255, 255, 255, 255], 1);
        assert_eq!(res, &[u32::MAX]);

        let res: Vec<u64> =
            SufrBuilder::slice_u8_to_vec(&[255, 255, 255, 255, 255, 255, 255, 255], 1);
        assert_eq!(res, &[u64::MAX]);

        let res: Vec<u32> =
            SufrBuilder::slice_u8_to_vec(&[0, 0, 0, 0, 255, 255, 255, 255], 2);
        assert_eq!(res, &[0u32, u32::MAX]);

        let res: Vec<u64> = SufrBuilder::slice_u8_to_vec(
            &[
                0, 0, 0, 0, 0, 0, 0, 0, 255, 255, 255, 255, 255, 255, 255, 255,
            ],
            2,
        );
        assert_eq!(res, &[0u64, u64::MAX]);

        Ok(())
    }

    #[test]
    fn test_vec_to_slice_u8() -> Result<()> {
        let res = SufrBuilder::vec_to_slice_u8(&[0u32]);
        assert_eq!(res, &[0, 0, 0, 0]);

        let res = SufrBuilder::vec_to_slice_u8(&[0u64]);
        assert_eq!(res, &[0, 0, 0, 0, 0, 0, 0, 0]);

        let res = SufrBuilder::vec_to_slice_u8(&[1u32]);
        assert_eq!(res, &[1, 0, 0, 0]);

        let res = SufrBuilder::vec_to_slice_u8(&[1u64]);
        assert_eq!(res, &[1, 0, 0, 0, 0, 0, 0, 0]);

        let res = SufrBuilder::vec_to_slice_u8(&[u32::MAX]);
        assert_eq!(res, &[255, 255, 255, 255]);

        let res = SufrBuilder::vec_to_slice_u8(&[u64::MAX]);
        assert_eq!(res, &[255, 255, 255, 255, 255, 255, 255, 255]);

        let res = SufrBuilder::vec_to_slice_u8(&[0u32, u32::MAX]);
        assert_eq!(res, &[0, 0, 0, 0, 255, 255, 255, 255]);

        let res = SufrBuilder::vec_to_slice_u8(&[0u64, u64::MAX]);
        assert_eq!(
            res,
            &[0, 0, 0, 0, 0, 0, 0, 0, 255, 255, 255, 255, 255, 255, 255, 255]
        );

        Ok(())
    }

    #[test]
    fn test_compare() -> Result<()> {
        // 0  1  2  3  4  5  6  7  8  9 10 11 12 13 14
        // A  A  B  A  B  A  B  A  B  B  A  B  A  B  #
        let sufr_file: SufrFile<u32> = SufrFile::read("tests/inputs/abba.sufr")?;

        // Compare to B to B with no skip
        let query = "B".as_bytes();
        let res = sufr_file.compare(query, 13, 0);
        assert_eq!(res.cmp, Ordering::Equal);
        assert_eq!(res.lcp, 1);

        // Compare to B to B with skip = 1
        let query = "B".as_bytes();
        let res = sufr_file.compare(query, 13, 1);
        assert_eq!(res.cmp, Ordering::Equal);
        assert_eq!(res.lcp, 1);

        // Compare to B to AB
        let query = "B".as_bytes();
        let res = sufr_file.compare(query, 12, 0);
        assert_eq!(res.cmp, Ordering::Greater);
        assert_eq!(res.lcp, 0);

        // Compare to ABABA to ABBABAB#
        let query = "ABABA".as_bytes();
        let res = sufr_file.compare(query, 7, 2);
        assert_eq!(res.cmp, Ordering::Less);
        assert_eq!(res.lcp, 2);

        // Compare to ABAB to ABABBABAB#
        let query = "ABABA".as_bytes();
        let res = sufr_file.compare(query, 5, 2);
        assert_eq!(res.cmp, Ordering::Less);
        assert_eq!(res.lcp, 4);

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

        let mut sufr_file: SufrFile<u32> = SufrFile::read("tests/inputs/abba.sufr")?;

        for val in &[true, false] {
            let args = Locate {
                queries: vec!["A".to_string()],
                max_query_len: None,
                low_memory: *val,
            };
            let res = sufr_file.locate(args);
            assert!(res.is_ok());
            let res = res.unwrap();
            assert_eq!(res.len(), 1);
            let val = &res[0];
            assert!(val.is_ok());
            let res = val.as_ref().unwrap();
            assert_eq!(
                res,
                &LocateResult {
                    query: "A".to_string(),
                    ranks: 1..8,
                    positions: vec![
                        LocateResultPosition {
                            suffix: 0,
                            sequence_name: "1".to_string(),
                            sequence_position: 0,
                        },
                        LocateResultPosition {
                            suffix: 12,
                            sequence_name: "1".to_string(),
                            sequence_position: 12,
                        },
                        LocateResultPosition {
                            suffix: 10,
                            sequence_name: "1".to_string(),
                            sequence_position: 10,
                        },
                        LocateResultPosition {
                            suffix: 1,
                            sequence_name: "1".to_string(),
                            sequence_position: 1,
                        },
                        LocateResultPosition {
                            suffix: 3,
                            sequence_name: "1".to_string(),
                            sequence_position: 3,
                        },
                        LocateResultPosition {
                            suffix: 5,
                            sequence_name: "1".to_string(),
                            sequence_position: 5,
                        },
                        LocateResultPosition {
                            suffix: 7,
                            sequence_name: "1".to_string(),
                            sequence_position: 7,
                        },
                    ]
                }
            );
        }

        for val in &[true, false] {
            let args = Locate {
                queries: vec!["B".to_string()],
                max_query_len: None,
                low_memory: *val,
            };
            let res = sufr_file.locate(args);
            assert!(res.is_ok());
            let res = res.unwrap();
            assert_eq!(res.len(), 1);
            let val = &res[0];
            assert!(val.is_ok());
            let res = val.as_ref().unwrap();
            assert_eq!(
                res,
                &LocateResult {
                    query: "B".to_string(),
                    ranks: 8..15,
                    positions: vec![
                        LocateResultPosition {
                            suffix: 13,
                            sequence_name: "1".to_string(),
                            sequence_position: 13,
                        },
                        LocateResultPosition {
                            suffix: 11,
                            sequence_name: "1".to_string(),
                            sequence_position: 11,
                        },
                        LocateResultPosition {
                            suffix: 9,
                            sequence_name: "1".to_string(),
                            sequence_position: 9,
                        },
                        LocateResultPosition {
                            suffix: 2,
                            sequence_name: "1".to_string(),
                            sequence_position: 2,
                        },
                        LocateResultPosition {
                            suffix: 4,
                            sequence_name: "1".to_string(),
                            sequence_position: 4,
                        },
                        LocateResultPosition {
                            suffix: 6,
                            sequence_name: "1".to_string(),
                            sequence_position: 6,
                        },
                        LocateResultPosition {
                            suffix: 8,
                            sequence_name: "1".to_string(),
                            sequence_position: 8,
                        },
                    ]
                }
            );
        }

        for val in &[true, false] {
            let args = Locate {
                queries: vec!["ABAB".to_string()],
                max_query_len: None,
                low_memory: *val,
            };
            let res = sufr_file.locate(args);
            assert!(res.is_ok());
            let res = res.unwrap();
            assert_eq!(res.len(), 1);
            let val = &res[0];
            assert!(val.is_ok());
            let res = val.as_ref().unwrap();
            assert_eq!(
                res,
                &LocateResult {
                    query: "ABAB".to_string(),
                    ranks: 3..7,
                    positions: vec![
                        LocateResultPosition {
                            suffix: 10,
                            sequence_name: "1".to_string(),
                            sequence_position: 10,
                        },
                        LocateResultPosition {
                            suffix: 1,
                            sequence_name: "1".to_string(),
                            sequence_position: 1,
                        },
                        LocateResultPosition {
                            suffix: 3,
                            sequence_name: "1".to_string(),
                            sequence_position: 3,
                        },
                        LocateResultPosition {
                            suffix: 5,
                            sequence_name: "1".to_string(),
                            sequence_position: 5,
                        },
                    ]
                }
            );
        }

        for val in &[true, false] {
            let args = Locate {
                queries: vec!["ABABB".to_string()],
                max_query_len: None,
                low_memory: *val,
            };
            let res = sufr_file.locate(args);
            assert!(res.is_ok());
            let res = res.unwrap();
            assert_eq!(res.len(), 1);
            let val = &res[0];
            assert!(val.is_ok());
            let res = val.as_ref().unwrap();
            assert_eq!(
                res,
                &LocateResult {
                    query: "ABABB".to_string(),
                    ranks: 6..7,
                    positions: vec![
                        LocateResultPosition {
                            suffix: 5,
                            sequence_name: "1".to_string(),
                            sequence_position: 5,
                        },
                    ]
                }
            );
        }

        for val in &[true, false] {
            let args = Locate {
                queries: vec!["BBBB".to_string()],
                max_query_len: None,
                low_memory: *val,
            };
            let res = sufr_file.locate(args);
            assert!(res.is_ok());
            assert!(res.is_ok());
            let res = res.unwrap();
            assert_eq!(res.len(), 1);
            let val = &res[0];
            assert!(val.is_err());
            assert_eq!(val.as_ref().unwrap_err().to_string(), "BBBB".to_string());
        }

        Ok(())
    }

    #[test]
    fn test_usize_to_bytes() -> Result<()> {
        assert_eq!(usize_to_bytes(usize::MIN), [0, 0, 0, 0, 0, 0, 0, 0]);
        assert_eq!(usize_to_bytes(1), [1, 0, 0, 0, 0, 0, 0, 0]);
        assert_eq!(usize_to_bytes(10), [10, 0, 0, 0, 0, 0, 0, 0]);
        assert_eq!(usize_to_bytes(100), [100, 0, 0, 0, 0, 0, 0, 0]);
        assert_eq!(usize_to_bytes(1000), [232, 3, 0, 0, 0, 0, 0, 0]);
        assert_eq!(usize_to_bytes(10000), [16, 39, 0, 0, 0, 0, 0, 0]);
        assert_eq!(usize_to_bytes(100000), [160, 134, 1, 0, 0, 0, 0, 0]);
        assert_eq!(
            usize_to_bytes(usize::MAX),
            [255, 255, 255, 255, 255, 255, 255, 255]
        );
        Ok(())
    }

    #[test]
    fn test_read_sequence_file() -> Result<()> {
        let file = "tests/inputs/2.fa";
        let sequence_delimiter = b'N';
        let res = read_sequence_file(file, sequence_delimiter);
        assert!(res.is_ok());
        let data = res.unwrap();
        assert_eq!(data.seq, b"ACGTacgtNacgtACGT$");
        assert_eq!(data.start_positions, [0, 9]);
        assert_eq!(data.headers, ["ABC", "DEF"]);
        Ok(())
    }

    #[test]
    fn test_read_text_length() -> Result<()> {
        let sufr_file = "tests/inputs/2.sufr";
        let res = read_text_length(sufr_file);
        assert!(res.is_ok());
        let len = res.unwrap();
        assert_eq!(len, 18);
        Ok(())
    }

    #[test]
    fn test_write_read_suffix_file_32() -> Result<()> {
        let seq_file = "tests/inputs/2.fa";
        let sequence_delimiter = b'N';
        let seq_data = read_sequence_file(seq_file, sequence_delimiter)?;
        let args = SufrBuilderArgs {
            text: seq_data.seq,
            max_query_len: None,
            is_dna: true,
            allow_ambiguity: false,
            ignore_softmask: false,
            sequence_starts: seq_data.start_positions.into_iter().collect(),
            headers: seq_data.headers,
            num_partitions: 2,
            sequence_delimiter,
        };
        let sufr_builder: SufrBuilder<u32> = SufrBuilder::new(args)?;
        // 17 13 9 0 4 14 10 1 5 15 11 2 6 16 12 3 7
        let sorted_sa = [17, 13, 9, 0, 4, 14, 10, 1, 5, 15, 11, 2, 6, 16, 12, 3, 7];
        let lcp = [0, 0, 4, 8, 4, 0, 3, 7, 3, 0, 2, 6, 2, 0, 1, 5, 1];
        let outfile = NamedTempFile::new()?;
        let outpath = &outfile.path().to_str().unwrap();
        let res = sufr_builder.write(outpath);
        assert!(res.is_ok());
        assert!(outfile.path().exists());

        let res: Result<SufrFile<u32>> = SufrFile::read(outpath);
        assert!(res.is_ok());

        let mut sufr_file = res.unwrap();
        assert_eq!(sufr_file.version, OUTFILE_VERSION);
        assert!(sufr_file.is_dna);
        assert_eq!(sufr_file.text_len, 18);
        assert_eq!(sufr_file.num_sequences, 2);
        assert_eq!(sufr_file.sequence_starts, [0, 9]);
        assert_eq!(sufr_file.headers, ["ABC", "DEF"]);
        assert_eq!(sufr_file.text, b"ACGTACGTNACGTACGT$");

        let file_sa: Vec<_> = sufr_file.suffix_array_file.iter().collect();
        assert_eq!(file_sa, sorted_sa);
        let file_lcp: Vec<_> = sufr_file.lcp_file.iter().collect();
        assert_eq!(file_lcp, lcp);
        Ok(())
    }

    #[test]
    fn test_write_read_suffix_file_64() -> Result<()> {
        let seq_file = "tests/inputs/1.fa";
        let sequence_delimiter = b'N';
        let seq_data = read_sequence_file(seq_file, sequence_delimiter)?;
        let args = SufrBuilderArgs {
            text: seq_data.seq,
            max_query_len: None,
            is_dna: true,
            allow_ambiguity: true,
            ignore_softmask: false,
            sequence_starts: seq_data.start_positions.into_iter().collect(),
            headers: seq_data.headers,
            num_partitions: 2,
            sequence_delimiter,
        };

        let suffix_array: SufrBuilder<u64> = SufrBuilder::new(args)?;
        let sorted_sa = [10, 6, 0, 7, 1, 8, 2, 5, 4, 9, 3];
        let lcp = [0, 0, 4, 0, 3, 0, 2, 0, 1, 0, 1];
        let outfile = NamedTempFile::new()?;
        let outpath = &outfile.path().to_str().unwrap();
        let res = suffix_array.write(outpath);
        assert!(res.is_ok());
        assert!(outfile.path().exists());

        let res: Result<SufrFile<u64>> = SufrFile::read(outpath);
        assert!(res.is_ok());

        let mut sufr_file = res.unwrap();
        assert_eq!(sufr_file.version, OUTFILE_VERSION);
        assert!(sufr_file.is_dna);
        assert_eq!(sufr_file.text_len, 11);
        assert_eq!(sufr_file.num_sequences, 1);
        assert_eq!(sufr_file.sequence_starts, [0]);
        assert_eq!(sufr_file.headers, ["1"]);
        assert_eq!(sufr_file.text, b"ACGTNNACGT$");

        let file_sa: Vec<_> = sufr_file.suffix_array_file.iter().collect();
        assert_eq!(file_sa, sorted_sa);
        let file_lcp: Vec<_> = sufr_file.lcp_file.iter().collect();
        assert_eq!(file_lcp, lcp);
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

    #[test]
    fn test_subsample_suffix_array() -> Result<()> {
        let seq_file = "tests/inputs/smol.fa";
        let sequence_delimiter = b'N';
        let seq_data = read_sequence_file(seq_file, sequence_delimiter)?;
        let builder_args = SufrBuilderArgs {
            text: seq_data.seq,
            max_query_len: None,
            is_dna: true,
            allow_ambiguity: false,
            ignore_softmask: false,
            sequence_starts: seq_data.start_positions,
            headers: seq_data.headers,
            num_partitions: 2,
            sequence_delimiter,
        };
        let sufr: SufrBuilder<u32> = SufrBuilder::new(builder_args)?;
        let outfile = NamedTempFile::new()?;
        let outpath = &outfile.path().to_str().unwrap();
        let res = sufr.write(outpath);
        assert!(res.is_ok());
        assert_eq!(sufr.num_suffixes, 364);

        let mut sufr_file: SufrFile<u32> = SufrFile::read(outpath)?;
        let full_sa: Vec<_> = sufr_file.suffix_array_file.iter().collect();
        let full_lcp: Vec<_> = sufr_file.lcp_file.iter().collect();
        assert_eq!(full_sa.len(), 364);
        assert_eq!(full_lcp.len(), 364);

        let max_query_len = 1;
        let (sub_sa, sub_rank) = sufr_file.subsample_suffix_array(max_query_len);
        assert_eq!(sub_sa.len(), 5);
        assert_eq!(sub_rank.len(), 5);
        // $, A, C, G, T
        assert_eq!(sub_sa, vec![365, 364, 92, 224, 363]);
        assert_eq!(sub_rank, vec![0, 1, 94, 191, 284]);

        let max_query_len = 2;
        let (sub_sa, sub_rank) = sufr_file.subsample_suffix_array(max_query_len);
        assert_eq!(sub_sa.len(), 20);
        assert_eq!(sub_rank.len(), 20);
        // $, A$, AA, AC, AG, AN, AT, CA, CC, CG, CT, GA, GC, GG, GN, GT, TA, TC, TG, TT
        assert_eq!(
            sub_sa,
            vec![
                365, 364, 358, 91, 341, 255, 362, 92, 339, 233, 296, 224, 88, 129, 110,
                96, 363, 217, 223, 356
            ]
        );
        assert_eq!(
            sub_rank,
            vec![
                0, 1, 2, 38, 49, 70, 71, 94, 112, 143, 170, 191, 216, 252, 269, 270,
                284, 298, 315, 343
            ]
        );

        let max_query_len = 3;
        let (sub_sa, sub_rank) = sufr_file.subsample_suffix_array(max_query_len);
        assert_eq!(sub_sa.len(), 71);
        assert_eq!(sub_rank.len(), 71);

        let max_query_len = 5;
        let (sub_sa, sub_rank) = sufr_file.subsample_suffix_array(max_query_len);
        assert_eq!(sub_sa.len(), 293);
        assert_eq!(sub_rank.len(), 293);

        let max_query_len = 3;
        let (sub_sa, sub_rank) = sufr_file.subsample_suffix_array(max_query_len);
        assert_eq!(sub_sa.len(), 71);
        assert_eq!(sub_rank.len(), 71);

        Ok(())
    }
}
