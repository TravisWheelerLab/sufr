use anyhow::{anyhow, bail, Result};
use log::info;
use needletail::parse_fastx_file;
//use range_minimum_query::Rmq;
//use rand::seq::SliceRandom;
use rand::Rng;
use rayon::prelude::*;
//use seq_io::fasta::{Reader, Record};
use std::{
    cmp::{max, min, Ordering},
    collections::HashSet,
    fmt::{Debug, Display},
    fs::{self, File, OpenOptions},
    hash::Hash,
    io::{Read, Seek, SeekFrom, Write},
    mem,
    //ops::Range,
    ops::{Add, Div, Range, Sub},
    path::PathBuf,
    slice,
    sync::{Arc, Mutex},
    time::Instant,
};
use tempfile::NamedTempFile;

const OUTFILE_VERSION: u8 = 2;

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
            let mut out = OpenOptions::new()
                .create(true)
                .append(true)
                .open(&self.path)?;
            out.write_all(SufrBuilder::vec_to_slice_u8(&self.vals[0..self.len]))?;
            self.total_len += self.len;
        }
        Ok(())
    }
}

// --------------------------------------------------
#[derive(Debug)]
pub struct SufrBuilder<T>
where
    T: Int + FromUsize<T> + Sized + Send + Sync + serde::ser::Serialize,
{
    pub version: u8,
    pub is_dna: bool,
    pub max_context: T,
    pub len: T,
    pub suffix_array_len: T,
    pub num_sequences: T,
    pub sequence_starts: Vec<T>,
    pub headers: Vec<String>,
    pub text: Vec<u8>,
    pub partitions: Vec<Partition<T>>,
}

impl<T> SufrBuilder<T>
where
    T: Int + FromUsize<T> + Sized + Send + Sync,
{
    pub fn new(
        text: Vec<u8>,
        len: T,
        max_context: Option<T>,
        is_dna: bool,
        sequence_starts: Vec<T>,
        headers: Vec<String>,
        num_partitions: usize,
    ) -> Result<SufrBuilder<T>> {
        let mut sa = SufrBuilder {
            version: OUTFILE_VERSION,
            is_dna,
            max_context: max_context.unwrap_or(len),
            len,
            suffix_array_len: T::default(),
            text,
            num_sequences: T::from_usize(sequence_starts.len()),
            sequence_starts,
            headers,
            partitions: vec![],
        };
        sa.sort(num_partitions)?;
        Ok(sa)
    }

    // --------------------------------------------------
    //pub fn check_lcp(&self) -> Vec<T> {
    //    let mut errors = vec![];
    //    for i in 1..self.suffix_array.len() {
    //        let len_lcp = self.find_lcp(
    //            self.suffix_array[i],
    //            self.suffix_array[i - 1],
    //            self.max_context,
    //        );
    //        if len_lcp != self.lcp[i] {
    //            errors.push(T::from_usize(i));
    //        }
    //    }
    //    errors
    //}

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
        let end1 = min(start1 + len, self.len.to_usize());
        let end2 = min(start2 + len, self.len.to_usize());
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
            let len_lcp = self.find_lcp(s1, s2, self.max_context).to_usize();

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
    pub fn upper_bound(&self, target: T, sa: &[T]) -> usize {
        // See if target is less than the first element
        if sa.is_empty() || self.is_less(target, sa[0]) {
            //T::default()
            0
        } else {
            // Find where all the values are less than target
            let i = sa.partition_point(|&p| self.is_less(p, target));

            // If the value at the partition is the same as the target
            if sa.get(i).map_or(false, |&v| v == target) {
                // Then return the next value, which might be out of range
                //T::from_usize(i + 1)
                i + 1
            } else {
                // Else return the partition point
                //T::from_usize(i)
                i
            }
        }
    }

    // --------------------------------------------------
    fn partition(
        &mut self,
        num_partitions: usize,
    ) -> Result<(Vec<Arc<Mutex<PartitionBuilder<T>>>>, usize)> {
        // Create more partitions than requested because
        // we can't know how big they will end up being
        let max_partitions = self.len.to_usize() / 4;
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
            .try_for_each(|(i, val)| -> Result<()> {
                if !self.is_dna || b"ACGT#".contains(val) {
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
            "Wrote unsorted partition{} in {:?}",
            if num_pivots == 1 { "" } else { "s" },
            now.elapsed()
        );

        Ok((builders, num_suffixes))
    }

    // --------------------------------------------------
    pub fn sort(&mut self, num_partitions: usize) -> Result<()> {
        // We will get more partition files than num_partitions
        let (mut part_files, num_suffixes) = self.partition(num_partitions)?;
        let num_per_partition = num_suffixes / num_partitions;
        let total_sort_time = Instant::now();
        let mut partition_inputs = vec![vec![]; num_partitions];

        for part_num in 0..num_partitions {
            let mut num_taken = 0;
            while !part_files.is_empty() {
                let part = part_files.remove(0);
                match part.lock() {
                    Ok(builder) => {
                        if builder.total_len > 0 {
                            partition_inputs[part_num]
                                .push((builder.path.clone(), builder.total_len));
                            num_taken += builder.total_len;
                        }
                    }
                    Err(e) => panic!("Can't get partition: {e}"),
                }

                if num_taken > num_per_partition {
                    break;
                }
            }
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
                    //let now = Instant::now();
                    let mut sa_w = part_sa.clone();
                    let mut lcp = vec![T::default(); len];
                    let mut lcp_w = vec![T::default(); len];
                    self.merge_sort(&mut sa_w, &mut part_sa, len, &mut lcp, &mut lcp_w);

                    // Write to disk
                    let mut sa_file = NamedTempFile::new()?;
                    let _ = sa_file.write(Self::vec_to_slice_u8(&part_sa))?;
                    let mut lcp_file = NamedTempFile::new()?;
                    let _ = lcp_file.write(Self::vec_to_slice_u8(&lcp))?;
                    //info!("LCP partition {partition_num}: {lcp:?}");
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
                    //info!(
                    //    "Partition {partition_num:4}: Sorted/wrote {len:8} in {:?}",
                    //    now.elapsed()
                    //);
                }
                Ok(())
            },
        )?;

        // Get rid of None/unwrap Some, put in order
        let mut partitions: Vec<_> = partitions.into_iter().flatten().collect();
        partitions.sort_by_key(|p| p.order);

        let sizes: Vec<_> = partitions.iter().map(|p| p.len).collect();
        info!(
            "Sorted {num_partitions} partitions (avg {}) in {:?}",
            sizes.iter().sum::<usize>() / num_partitions,
            total_sort_time.elapsed()
        );
        self.suffix_array_len = T::from_usize(sizes.iter().sum());
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
                    let max_n = self.len - max(x[idx_x], y[idx_y]);

                    // Prefix-context length for the suffixes
                    let context = min(self.max_context, max_n);

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
                if self.is_dna && !b"ACGT".contains(&self.text[pos]) {
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
    pub fn write(&self, mut out: impl Write) -> Result<usize> {
        let mut bytes_out = 0;

        // Header: version/is_dna
        let is_dna: u8 = if self.is_dna { 1 } else { 0 };
        bytes_out += out.write(&[OUTFILE_VERSION, is_dna])?;

        // Text length
        bytes_out += out.write(&usize_to_bytes(self.len.to_usize()))?;

        // SA length
        bytes_out += out.write(&usize_to_bytes(self.suffix_array_len.to_usize()))?;

        // Max context
        bytes_out += out.write(&usize_to_bytes(self.max_context.to_usize()))?;

        // Number of sequences
        bytes_out += out.write(&usize_to_bytes(self.sequence_starts.len()))?;

        // Sequence starts
        bytes_out += out.write(Self::vec_to_slice_u8(&self.sequence_starts))?;

        // Text
        bytes_out += out.write(&self.text)?;

        // Stitch partitioned suffix files together
        for partition in &self.partitions {
            let buffer = fs::read(&partition.sa_path)?;
            bytes_out += out.write(&buffer)?;
        }

        // Stitch partitioned LCP files together
        for (i, partition) in self.partitions.iter().enumerate() {
            let buffer = fs::read(&partition.lcp_path)?;

            if i == 0 {
                bytes_out += out.write(&buffer)?;
            } else {
                // Fix LCP boundary
                let mut lcp = Self::slice_u8_to_vec(&buffer, partition.len);
                if let Some(val) = lcp.first_mut() {
                    *val = self.find_lcp(
                        self.partitions[i - 1].last_suffix,
                        partition.first_suffix,
                        self.max_context,
                    );
                }
                bytes_out += out.write(Self::vec_to_slice_u8(&lcp))?;
            }
        }

        // Headers are variable in length so they are at the end
        bytes_out += out.write(&bincode::serialize(&self.headers)?)?;

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
#[derive(Debug, Clone)]
pub struct FileAccess<T>
where
    T: Int + FromUsize<T> + Sized + Send + Sync + serde::ser::Serialize,
{
    filename: String,
    buffer: Vec<T>,
    buffer_size: usize,
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
    pub fn new(filename: &str, start: u64, num_elements: usize) -> Self {
        let size = num_elements * mem::size_of::<T>();
        FileAccess {
            filename: filename.to_string(),
            buffer: vec![],
            buffer_size: 1048576,
            size,
            start_position: start,
            current_position: start,
            end_position: start + size as u64,
            exhausted: false,
        }
    }

    // --------------------------------------------------
    // TODO: Ignoring lots of Results to return Option
    pub fn get(&self, pos: usize) -> Option<T> {
        let mut file = File::open(&self.filename).unwrap();

        // Don't bother looking for something beyond the end
        let seek = self.start_position + (pos * mem::size_of::<T>()) as u64;
        if seek < self.end_position {
            let _ = file.seek(SeekFrom::Start(seek));
            let mut buffer: Vec<u8> = vec![0; mem::size_of::<T>()];
            let bytes_read = file.read(&mut buffer).unwrap();
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
    pub fn get_range(&self, range: Range<usize>) -> Result<Vec<T>> {
        let mut file = File::open(&self.filename).unwrap();
        let start = self.start_position as usize + (range.start * mem::size_of::<T>());
        let end = self.start_position as usize + (range.end * mem::size_of::<T>());
        let valid = self.start_position as usize..self.end_position as usize;
        if valid.contains(&start) && valid.contains(&end) {
            file.seek(SeekFrom::Start(start as u64))?;
            let mut buffer: Vec<u8> = vec![0; end - start];
            let bytes_read = file.read(&mut buffer)?;
            let num_vals = bytes_read / mem::size_of::<T>();
            Ok(SufrBuilder::slice_u8_to_vec(&buffer, num_vals))
        } else {
            bail!("Invalid range: {range:?}")
        }
    }
}

impl<T> Iterator for FileAccess<T>
where
    T: Int + FromUsize<T> + Sized + Send + Sync + serde::ser::Serialize,
{
    type Item = T;

    fn next(&mut self) -> Option<Self::Item> {
        if self.exhausted {
            None
        } else {
            // Fill the buffer
            if self.buffer.is_empty() {
                if self.current_position >= self.end_position {
                    self.exhausted = true;
                    return None;
                }

                let mut file = File::open(&self.filename).unwrap();
                file.seek(SeekFrom::Start(self.current_position)).unwrap();
                let mut bytes_wanted = self.buffer_size * mem::size_of::<T>();
                if self.current_position + bytes_wanted as u64 > self.end_position {
                    bytes_wanted = (self.end_position - self.current_position) as usize;
                }
                let mut buffer: Vec<u8> = vec![0; bytes_wanted];
                let bytes_read = file.read(&mut buffer).unwrap();
                self.current_position = file.stream_position().unwrap();
                let num_vals = bytes_read / mem::size_of::<T>();
                self.buffer = SufrBuilder::slice_u8_to_vec(&buffer, num_vals);
                self.buffer.reverse();
            }

            self.buffer.pop()
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
    pub max_context: T,
    pub len: T,
    pub num_suffixes: T,
    pub num_sequences: T,
    pub sequence_starts: Vec<T>,
    pub headers: Vec<String>,
    pub text: Vec<u8>,
    pub suffix_array: FileAccess<T>,
    pub lcp: FileAccess<T>,
}

impl<T> SufrFile<T>
where
    T: Int + FromUsize<T> + Sized + Send + Sync,
{
    // Read serialized ".sufr" file
    pub fn read(filename: &str) -> Result<SufrFile<T>> {
        let mut file = File::open(filename).map_err(|e| anyhow!("{filename}: {e}"))?;

        // Meta (version, is_dna)
        let mut buffer = [0; 2];
        file.read_exact(&mut buffer)?;
        let version = buffer[0];
        let is_dna = buffer[1] == 1;

        // Length of text
        let mut buffer = [0; 8];
        file.read_exact(&mut buffer)?;
        let text_len = usize::from_ne_bytes(buffer);

        // Length of SA
        let mut buffer = [0; 8];
        file.read_exact(&mut buffer)?;
        let num_suffixes = usize::from_ne_bytes(buffer);

        // Max context
        let mut buffer = [0; 8];
        file.read_exact(&mut buffer)?;
        let max_context = T::from_usize(usize::from_ne_bytes(buffer));

        // Number of sequences
        let mut buffer = [0; 8];
        file.read_exact(&mut buffer)?;
        let num_sequences = T::from_usize(usize::from_ne_bytes(buffer));

        // Sequence starts
        let mut buffer = vec![0u8; num_sequences.to_usize() * mem::size_of::<T>()];
        file.read_exact(&mut buffer)?;
        let sequence_starts: Vec<T> =
            SufrBuilder::slice_u8_to_vec(&buffer, num_sequences.to_usize());

        // Text
        let mut text = vec![0u8; text_len];
        file.read_exact(&mut text)?;

        // Suffix Array
        let suffix_array: FileAccess<T> = FileAccess::new(
            filename,
            file.stream_position()?,
            num_suffixes,
        );
        file.seek_relative(suffix_array.size as i64)?;

        // LCP
        let lcp: FileAccess<T> = FileAccess::new(
            filename,
            file.stream_position()?,
            num_suffixes,
        );
        file.seek_relative(lcp.size as i64)?;

        // Headers are variable in length so they are at the end
        let mut buffer = vec![];
        file.read_to_end(&mut buffer)?;
        let headers: Vec<String> = bincode::deserialize(&buffer)?;

        Ok(SufrFile {
            filename: filename.to_string(),
            version,
            is_dna,
            len: T::from_usize(text_len),
            num_suffixes: T::from_usize(num_suffixes),
            max_context,
            num_sequences,
            sequence_starts,
            headers,
            text,
            suffix_array,
            lcp,
        })
    }

    // --------------------------------------------------
    pub fn check_suffix_array(&self) -> Result<Vec<usize>> {
        let sa = self.suffix_array.clone();
        let lcp = self.lcp.clone();
        let mut previous: Option<usize> = None;
        let mut errors = vec![];

        for (cur_sa, len_lcp) in sa.into_iter().zip(lcp.into_iter()) {
            let cur_sa = cur_sa.to_usize();
            let len_lcp = len_lcp.to_usize();
            if let Some(prev_sa) = previous {
                let is_less = match (
                    self.text.get(prev_sa + len_lcp),
                    self.text.get(cur_sa + len_lcp),
                ) {
                    (Some(a), Some(b)) => a < b,
                    (None, Some(_)) => true,
                    _ => false,
                };
                if !is_less {
                    errors.push(cur_sa);
                }
            }
            previous = Some(cur_sa);
        }
        Ok(errors)
    }

    // --------------------------------------------------
    pub fn string_at(&self, pos: usize, len: Option<usize>) -> String {
        let end = len.map_or(self.len.to_usize(), |n| pos + n);
        self.text
            .get(pos..end)
            .map(|v| String::from_utf8(v.to_vec()).unwrap())
            .unwrap()
    }

    // --------------------------------------------------
    pub fn search(&self, query: &str) -> Option<(usize, usize)> {
        let qry = query.as_bytes();
        let n = self.num_suffixes.to_usize();
        self.suffix_search_first(qry, 0, n - 1, 0, 0).map(|first| {
            let last = self.suffix_search_last(qry, first, n - 1, n, 0, 0).unwrap();
            // TODO: Return zero-offset or 1-based?
            (first, last)
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

            let mid_cmp = self.compare(
                qry,
                self.suffix_array.get(mid)?.to_usize(),
                min(left_lcp, right_lcp),
            );

            if mid_cmp.cmp == Ordering::Equal
                && (mid == 0
                    || self
                        .compare(qry, self.suffix_array.get(mid - 1)?.to_usize(), 0)
                        .cmp
                        == Ordering::Greater)
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
            let mid_cmp = self.compare(
                qry,
                self.suffix_array.get(mid)?.to_usize(),
                min(left_lcp, right_lcp),
            );

            if mid_cmp.cmp == Ordering::Equal
                && (mid == n - 1
                    || self
                        .compare(qry, self.suffix_array.get(mid + 1)?.to_usize(), 0)
                        .cmp
                        == Ordering::Less)
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
            .skip(skip)
            .zip(self.text.get(suffix_pos + skip..).unwrap())
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
// Utility function to find suffix array length to
// determine whether to build u32 or u64
pub fn read_suffix_length(filename: &str) -> Result<usize> {
    let mut file = File::open(filename).map_err(|e| anyhow!("{filename}: {e}"))?;

    // Meta (version, is_dna)
    let mut buffer = [0; 2];
    file.read_exact(&mut buffer)?;

    let outfile_version = buffer[0];
    if outfile_version == OUTFILE_VERSION {
        // Length of SA is the next usize
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
pub fn read_sequence_file(filename: &str) -> Result<SequenceFileData> {
    let mut reader = parse_fastx_file(filename)?;
    //let mut seq : Vec<u8> = vec![];
    let mut seq = Vec::with_capacity(u32::MAX as usize);
    let mut headers: Vec<String> = vec![];
    let mut start_positions: Vec<usize> = vec![];
    let mut i = 0;
    while let Some(rec) = reader.next() {
        let rec = rec?;
        if i > 0 {
            // Sequence delimiter
            seq.push(b'$');
        }

        // Record current length as start position
        start_positions.push(seq.len());

        // Uppercase (mask w/32)
        let mut current: Vec<u8> = rec.seq().iter().map(|b| b & 0b1011111).collect();
        seq.append(&mut current);
        i += 1;

        headers.push(String::from_utf8(rec.id().to_vec())?);
    }

    // File delimiter
    seq.push(b'#');

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
        read_sequence_file, read_suffix_length, usize_to_bytes, SufrBuilder, SufrFile,
    };
    use anyhow::Result;
    use std::cmp::Ordering;

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
    fn test_search() -> Result<()> {
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
        //
        let sufr_file : SufrFile<u32> = SufrFile::read("tests/inputs/abba.sufr")?;

        let res = sufr_file.search("B");
        assert_eq!(res, Some((8, 14)));

        let res = sufr_file.search("AB");
        assert_eq!(res, Some((2, 7)));

        let res = sufr_file.search("BABAB");
        assert_eq!(res, Some((10, 12)));

        let res = sufr_file.search("ABAA");
        assert_eq!(res, None);

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
        let res = read_sequence_file(file);
        assert!(res.is_ok());
        let data = res.unwrap();
        assert_eq!(data.seq, b"ACGTACGT$ACGTACGT#");
        assert_eq!(data.start_positions, [0, 9]);
        assert_eq!(data.headers, ["ABC", "DEF"]);
        Ok(())
    }

    #[test]
    fn test_read_suffix_length() -> Result<()> {
        let sufr_file = "tests/inputs/2.sufr";
        let res = read_suffix_length(sufr_file);
        assert!(res.is_ok());
        let len = res.unwrap();
        assert_eq!(len, 18);
        Ok(())
    }

    //#[test]
    //fn test_write_read_suffix_file_32() -> Result<()> {
    //    let seq_file = "tests/inputs/2.fa";
    //    let seq_data = read_sequence_file(seq_file)?;
    //    let max_context: Option<u32> = None;
    //    let is_dna = true;
    //    let start_positions: Vec<_> =
    //        seq_data.start_positions.iter().map(|&v| v as u32).collect();
    //    let len = seq_data.seq.len() as u32;
    //    let num_partitions = 2;
    //    let suffix_array: SuffixArray<u32> = SuffixArray::new(
    //        seq_data.seq,
    //        len,
    //        max_context,
    //        is_dna,
    //        start_positions,
    //        seq_data.headers,
    //        num_partitions,
    //    )?;
    //
    //    let sorted_sa =
    //        [17, 13, 4, 9, 0, 14, 5, 10, 1, 15, 6, 11, 2, 16, 7, 12, 3];
    //    let lcp = [0, 0, 4, 4, 8, 0, 3, 3, 7, 0, 2, 2, 6, 0, 1, 1, 5];
    //    let outfile = NamedTempFile::new()?;
    //    let outpath = &outfile.path().to_str().unwrap();
    //    let out = BufWriter::new(
    //        File::create(outpath).map_err(|e| anyhow!("{outpath}: {e}"))?,
    //    );
    //    let res = suffix_array.write(out);
    //    assert!(res.is_ok());
    //    assert!(outfile.path().exists());
    //
    //    let res: Result<SuffixArray<u32>> = SuffixArray::read(outpath);
    //    assert!(res.is_ok());
    //
    //    let sa = res.unwrap();
    //    assert_eq!(sa.version, 1);
    //    assert!(sa.is_dna);
    //    assert_eq!(sa.len, 18);
    //    assert_eq!(sa.num_sequences, 2);
    //    assert_eq!(sa.sequence_starts, [0, 9]);
    //    assert_eq!(sa.headers, ["ABC", "DEF"]);
    //    assert_eq!(sa.suffix_array, sorted_sa);
    //    assert_eq!(sa.lcp, lcp);
    //    assert_eq!(sa.text, b"ACGTACGT$ACGTACGT#");
    //    Ok(())
    //}
    //
    //#[test]
    //fn test_write_read_suffix_file_64() -> Result<()> {
    //    let seq_file = "tests/inputs/2.fa";
    //    let seq_data = read_sequence_file(seq_file)?;
    //    let max_context: Option<u64> = None;
    //    let is_dna = true;
    //    let start_positions: Vec<_> =
    //        seq_data.start_positions.iter().map(|&v| v as u64).collect();
    //    let len = seq_data.seq.len() as u64;
    //    let num_partitions = 2;
    //    let suffix_array: SuffixArray<u64> = SuffixArray::new(
    //        seq_data.seq,
    //        len,
    //        max_context,
    //        is_dna,
    //        start_positions,
    //        seq_data.headers,
    //        num_partitions,
    //    )?;
    //
    //    let sorted_sa: &[u64] =
    //        &[17, 13, 4, 9, 0, 14, 5, 10, 1, 15, 6, 11, 2, 16, 7, 12, 3];
    //    let lcp: &[u64] =
    //        &[0, 0, 4, 4, 8, 0, 3, 3, 7, 0, 2, 2, 6, 0, 1, 1, 5];
    //    let outfile = NamedTempFile::new()?;
    //    let outpath = &outfile.path().to_str().unwrap();
    //    let out = BufWriter::new(
    //        File::create(outpath).map_err(|e| anyhow!("{outpath}: {e}"))?,
    //    );
    //    let res = suffix_array.write(out);
    //    assert!(res.is_ok());
    //    assert!(outfile.path().exists());
    //
    //    let res: Result<SuffixArray<u64>> = SuffixArray::read(outpath);
    //    assert!(res.is_ok());
    //
    //    let sa = res.unwrap();
    //    assert_eq!(sa.version, 1);
    //    assert!(sa.is_dna);
    //    assert_eq!(sa.len, 18);
    //    assert_eq!(sa.num_sequences, 2);
    //    assert_eq!(sa.sequence_starts, [0, 9]);
    //    assert_eq!(sa.headers, ["ABC", "DEF"]);
    //    assert_eq!(sa.suffix_array, sorted_sa);
    //    assert_eq!(sa.lcp, lcp);
    //    assert_eq!(sa.text, b"ACGTACGT$ACGTACGT#");
    //    Ok(())
    //}
    //
    //#[test]
    //fn test_upper_bound() -> Result<()> {
    //    //          012345
    //    let text = "TTTAGC".as_bytes().to_vec();
    //    let len = text.len();
    //    let max_context: Option<u32> = None;
    //    let is_dna = false;
    //    let sequence_starts = vec![0];
    //    let headers = vec!["1".to_string()];
    //    let num_partitions = 2;
    //    let sa: SuffixArray<u32> = SuffixArray::new(
    //        text,
    //        len as u32,
    //        max_context,
    //        is_dna,
    //        sequence_starts,
    //        headers,
    //        num_partitions,
    //    )?;
    //
    //    // The suffix "AGC$" is found before "GC$" and "C$
    //    assert_eq!(sa.upper_bound(3, &[5, 4]), None);
    //
    //    // The suffix "TAGC$" is beyond all the values
    //    assert_eq!(sa.upper_bound(2, &[3, 5, 4]), Some(3));
    //
    //    // The "C$" is the last value
    //    assert_eq!(sa.upper_bound(5, &[3, 5, 4]), Some(2));
    //
    //    Ok(())
    //}
}
