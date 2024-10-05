use anyhow::{anyhow, bail, Result};
use log::info;
use needletail::parse_fastx_file;
//use rand::seq::SliceRandom;
use rand::Rng;
use rayon::prelude::*;
//use seq_io::fasta::{Reader, Record};
use std::{
    cmp::{max, min, Ordering},
    collections::HashSet,
    fmt::{Debug, Display},
    fs::{self, File},
    hash::Hash,
    io::{Read, Write},
    mem,
    //ops::Range,
    ops::{Add, Div, Sub},
    path::PathBuf,
    slice,
    time::Instant,
};
use tempfile::NamedTempFile;

const OUTFILE_VERSION: u8 = 1;

#[derive(Debug)]
pub struct Comparison {
    cmp: Ordering,
    lcp: usize,
}

#[derive(Debug)]
pub struct SequenceFileData {
    pub seq: Vec<u8>,
    pub start_positions: Vec<usize>,
    pub headers: Vec<String>,
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
pub struct SuffixArrayBuilder<T>
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

impl<T> SuffixArrayBuilder<T>
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
    ) -> Result<SuffixArrayBuilder<T>> {
        let mut sa = SuffixArrayBuilder {
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
                        self.text.get_unchecked(*a)
                            == self.text.get_unchecked(*b)
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
    pub fn upper_bound(&self, target: T, sa: &[T]) -> T {
        // See if target is less than the first element
        if sa.is_empty() || self.is_less(target, sa[0]) {
            T::default()
        } else {
            // Find where all the values are less than target
            let i = sa.partition_point(|&p| self.is_less(p, target));

            // If the value at the partition is the same as the target
            if sa.get(i).map_or(false, |&v| v == target) {
                // Then return the next value, which might be out of range
                T::from_usize(i + 1)
            } else {
                // Else return the partition point
                T::from_usize(i)
            }
        }
    }

    //fn partition(
    //    &mut self,
    //    text: &[u8],
    //    num_partitions: usize,
    //) -> Result<()> {
    //    Ok(())
    //}

    // --------------------------------------------------
    pub fn sort(&mut self, num_partitions: usize) -> Result<()> {
        // Randomly select some pivots
        let now = Instant::now();
        let pivot_sa = self.select_pivots(self.text.len(), num_partitions);
        let num_pivots = pivot_sa.len();
        info!(
            "Selected {num_pivots} pivot{} in {:?}",
            if num_pivots == 1 { "" } else { "s" },
            now.elapsed()
        );

        // Closure to find suffixes belonging to a partition
        let partitioner = |wanted: T| -> Vec<T> {
            self.text
                .par_iter()
                .enumerate()
                .flat_map(|(i, val)| {
                    if !self.is_dna || b"ACGT#".contains(val) {
                        let part =
                            self.upper_bound(T::from_usize(i), &pivot_sa);
                        (part == wanted).then_some(T::from_usize(i))
                    } else {
                        None
                    }
                })
                .collect()
        };

        let now = Instant::now();
        let mut partitions: Vec<Option<Partition<T>>> =
            (0..num_partitions).map(|_| None).collect();
        //partitions.iter_mut().enumerate().try_for_each(
        partitions.par_iter_mut().enumerate().try_for_each(
            |(partition_num, partition)| -> Result<()> {
                //let now = Instant::now();
                // Find the suffixes in this partition
                let mut part_sa = partitioner(T::from_usize(partition_num));
                //info!(
                //    "Found {} suffixes for partition {partition_num} in {:?}",
                //    part_sa.len(),
                //    now.elapsed()
                //);
                let len = part_sa.len();
                if len > 0 {
                    //let now = Instant::now();
                    //info!("UNSORTED partition {partition_num}: {part_sa:?}");
                    let mut sa_w = part_sa.clone();
                    let mut lcp = vec![T::default(); len];
                    let mut lcp_w = vec![T::default(); len];
                    self.merge_sort(
                        &mut sa_w,
                        &mut part_sa,
                        len,
                        &mut lcp,
                        &mut lcp_w,
                    );
                    //info!("SORTED partition {partition_num}  : {part_sa:?}");

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
                }
                //info!(
                //    "Sorted {} in partition {partition_num} in {:?}",
                //    part_sa.len(),
                //    now.elapsed()
                //);
                Ok(())
            },
        )?;

        // Get rid of None/unwrap Some, put in order
        let mut partitions: Vec<_> =
            partitions.into_iter().flatten().collect();
        partitions.sort_by_key(|p| p.order);

        let sizes: Vec<_> = partitions.iter().map(|p| p.len).collect();
        info!(
            "Split/sorted {num_partitions} partitions (avg {}) in {:?}",
            sizes.iter().sum::<usize>() / num_partitions,
            now.elapsed()
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
                    let len_lcp = m + self.find_lcp(
                        x[idx_x] + m,
                        y[idx_y] + m,
                        context - m,
                    );

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
    // Search SuffixArray
    //pub fn search(&self, query: &str) -> Option<Range<usize>> {
    //    //let suffixes: Vec<_> = self
    //    //    .suffix_array
    //    //    .iter()
    //    //    .zip(self.suffix_array.iter().map(|&p| self.string_at(p.to_usize())))
    //    //    .map(|(p,s)| format!("{p:3}: {s}"))
    //    //    .collect();
    //    //println!("{}", suffixes.join("\n"));
    //    let qry = query.as_bytes();
    //    let n = self.suffix_array.len();
    //    self.suffix_search_first(qry, 0, n - 1, 0, 0).map(|first| {
    //        let last =
    //            self.suffix_search_last(qry, first, n - 1, n, 0, 0).unwrap();
    //        first..last
    //    })
    //}

    // --------------------------------------------------
    //fn suffix_search_first(
    //    &self,
    //    qry: &[u8],
    //    low: usize,
    //    high: usize,
    //    left_lcp: usize,
    //    right_lcp: usize,
    //) -> Option<usize> {
    //    if high >= low {
    //        let mid = low + ((high - low) / 2);
    //
    //        let mid_cmp = self.compare(
    //            qry,
    //            self.suffix_array[mid].to_usize(),
    //            min(left_lcp, right_lcp),
    //        );
    //
    //        if mid_cmp.cmp == Ordering::Equal
    //            && (mid == 0
    //                || self
    //                    .compare(
    //                        qry,
    //                        self.suffix_array[mid - 1].to_usize(),
    //                        0,
    //                    )
    //                    .cmp
    //                    == Ordering::Greater)
    //        {
    //            Some(mid)
    //        } else if mid_cmp.cmp == Ordering::Greater {
    //            self.suffix_search_first(
    //                qry,
    //                mid + 1,
    //                high,
    //                mid_cmp.lcp,
    //                right_lcp,
    //            )
    //        } else {
    //            // Ordering::Less
    //            self.suffix_search_first(
    //                qry,
    //                low,
    //                mid - 1,
    //                left_lcp,
    //                mid_cmp.lcp,
    //            )
    //        }
    //    } else {
    //        None
    //    }
    //}

    // --------------------------------------------------
    //fn suffix_search_last(
    //    &self,
    //    qry: &[u8],
    //    low: usize,
    //    high: usize,
    //    n: usize,
    //    left_lcp: usize,
    //    right_lcp: usize,
    //) -> Option<usize> {
    //    if high >= low {
    //        let mid = low + ((high - low) / 2);
    //        let mid_cmp = self.compare(
    //            qry,
    //            self.suffix_array[mid].to_usize(),
    //            min(left_lcp, right_lcp),
    //        );
    //
    //        if mid_cmp.cmp == Ordering::Equal
    //            && (mid == n - 1
    //                || self
    //                    .compare(
    //                        qry,
    //                        self.suffix_array[mid + 1].to_usize(),
    //                        0,
    //                    )
    //                    .cmp
    //                    == Ordering::Less)
    //        {
    //            Some(mid)
    //        } else if mid_cmp.cmp == Ordering::Less {
    //            self.suffix_search_last(
    //                qry,
    //                low,
    //                mid - 1,
    //                n,
    //                left_lcp,
    //                mid_cmp.lcp,
    //            )
    //        } else {
    //            self.suffix_search_last(
    //                qry,
    //                mid + 1,
    //                high,
    //                n,
    //                mid_cmp.lcp,
    //                right_lcp,
    //            )
    //        }
    //    } else {
    //        None
    //    }
    //}

    // --------------------------------------------------
    pub fn compare(
        &self,
        query: &[u8],
        suffix_pos: usize,
        skip: usize,
    ) -> Comparison {
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

        //let q = String::from_utf8(query.to_vec()).unwrap();
        //let s =
        //    String::from_utf8(self.text.get(suffix_pos..).unwrap().to_vec())
        //        .unwrap();
        //println!(
        //    "q '{q}' s at {suffix_pos} '{s}' skip {skip} = {:?}",
        //    Comparison { lcp, cmp }
        //);
        Comparison { lcp, cmp }
    }

    // --------------------------------------------------
    #[inline(always)]
    fn select_pivots(
        &self,
        text_len: usize,
        num_partitions: usize,
    ) -> Vec<T> {
        if num_partitions > 1 {
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
            let mut pivot_sa: Vec<T> = pivot_sa.iter().cloned().collect();
            let mut sa_w = pivot_sa.clone();
            let len = pivot_sa.len();
            let mut lcp = vec![T::default(); len];
            let mut lcp_w = vec![T::default(); len];
            self.merge_sort(
                &mut sa_w,
                &mut pivot_sa,
                len,
                &mut lcp,
                &mut lcp_w,
            );
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
        bytes_out +=
            out.write(&usize_to_bytes(self.suffix_array_len.to_usize()))?;

        // Max context
        bytes_out +=
            out.write(&usize_to_bytes(self.max_context.to_usize()))?;

        // Number of sequences
        bytes_out +=
            out.write(&usize_to_bytes(self.sequence_starts.len()))?;

        // Sequence starts
        bytes_out +=
            out.write(Self::vec_to_slice_u8(&self.sequence_starts))?;

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

        // Sequence
        bytes_out += out.write(&self.text)?;
        //let bytes_text = out.write(&self.text)?;
        //info!("Wrote {bytes_text} '{:?}'", self.text);
        //bytes_out += bytes_text;

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
        unsafe {
            std::slice::from_raw_parts(buffer.as_ptr() as *const _, len)
                .to_vec()
        }
    }
}

// --------------------------------------------------
#[derive(Debug)]
pub struct SuffixArray<T>
where
    T: Int + FromUsize<T> + Sized + Send + Sync + serde::ser::Serialize,
{
    pub version: u8,
    pub is_dna: bool,
    pub max_context: T,
    pub len: T,
    pub num_sequences: T,
    pub sequence_starts: Vec<T>,
    pub headers: Vec<String>,
    pub text: Vec<u8>,
    pub suffix_array: Vec<T>,
    pub lcp: Vec<T>,
}

impl<T> SuffixArray<T>
where
    T: Int + FromUsize<T> + Sized + Send + Sync,
{
    pub fn string_at(&self, pos: usize) -> String {
        self.text
            .get(pos..)
            .map(|v| String::from_utf8(v.to_vec()).unwrap())
            .unwrap()
    }

    //pub fn check_order(&self) -> Vec<T> {
    //    let mut errors = vec![];
    //    for window in self.suffix_array.windows(2) {
    //        if let [prev, cur] = window {
    //            if !self.is_less(*prev, *cur) {
    //                errors.push(*prev);
    //            }
    //        }
    //    }
    //    errors
    //}

    // Read serialized ".sufr" file
    pub fn read(filename: &str) -> Result<SuffixArray<T>> {
        println!("Reading '{filename}'");
        let mut file =
            File::open(filename).map_err(|e| anyhow!("{filename}: {e}"))?;

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
        let sa_len = usize::from_ne_bytes(buffer);

        // Max context
        let mut buffer = [0; 8];
        file.read_exact(&mut buffer)?;
        let max_context = T::from_usize(usize::from_ne_bytes(buffer));

        // Number of sequences
        let mut buffer = [0; 8];
        file.read_exact(&mut buffer)?;
        let num_sequences = T::from_usize(usize::from_ne_bytes(buffer));

        // Sequence starts
        let mut buffer =
            vec![0u8; num_sequences.to_usize() * mem::size_of::<T>()];
        file.read_exact(&mut buffer)?;
        let sequence_starts: Vec<T> = SuffixArrayBuilder::slice_u8_to_vec(
            &buffer,
            num_sequences.to_usize(),
        );

        // Suffix Array
        let mut buffer = vec![0u8; sa_len * mem::size_of::<T>()];
        file.read_exact(&mut buffer)?;
        let suffix_array: Vec<T> =
            SuffixArrayBuilder::slice_u8_to_vec(&buffer, sa_len);

        // LCP Array
        let mut buffer = vec![0u8; sa_len * mem::size_of::<T>()];
        file.read_exact(&mut buffer)?;
        let lcp: Vec<T> =
            SuffixArrayBuilder::slice_u8_to_vec(&buffer, sa_len);

        // Sequence/text
        let mut text = vec![0u8; text_len];
        file.read_exact(&mut text)?;

        // Headers are variable in length so they are at the end
        let mut buffer = vec![];
        file.read_to_end(&mut buffer)?;
        let headers: Vec<String> = bincode::deserialize(&buffer)?;

        Ok(SuffixArray {
            version,
            is_dna,
            len: T::from_usize(text_len),
            max_context,
            num_sequences,
            sequence_starts,
            headers,
            suffix_array,
            lcp,
            text,
        })
    }
}

// --------------------------------------------------
// Utility function to find suffix array length to
// determine whether to build u32 or u64
pub fn read_suffix_length(filename: &str) -> Result<usize> {
    let mut file =
        File::open(filename).map_err(|e| anyhow!("{filename}: {e}"))?;

    // Meta (version, is_dna)
    let mut buffer = [0; 2];
    file.read_exact(&mut buffer)?;

    let outfile_version = buffer[0];
    if outfile_version == 1 {
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
// data needed by SuffixArrayBuilder
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
        let mut current: Vec<u8> =
            rec.seq().iter().map(|b| b & 0b1011111).collect();
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
// Utility function used by SuffixArrayBuilder
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
        read_sequence_file, read_suffix_length, usize_to_bytes,
        SuffixArrayBuilder,
    };
    use anyhow::{anyhow, Result};
    use std::{cmp::Ordering, fs::File, io::BufWriter};
    use tempfile::NamedTempFile;

    #[test]
    fn test_slice_u8_to_vec() -> Result<()> {
        let res: Vec<u32> =
            SuffixArrayBuilder::slice_u8_to_vec(&[0, 0, 0, 0], 1);
        assert_eq!(res, &[0u32]);

        let res: Vec<u64> =
            SuffixArrayBuilder::slice_u8_to_vec(&[0, 0, 0, 0, 0, 0, 0, 0], 1);
        assert_eq!(res, &[0u64]);

        let res: Vec<u32> =
            SuffixArrayBuilder::slice_u8_to_vec(&[1, 0, 0, 0], 1);
        assert_eq!(res, &[1u32]);

        let res: Vec<u64> =
            SuffixArrayBuilder::slice_u8_to_vec(&[1, 0, 0, 0, 0, 0, 0, 0], 1);
        assert_eq!(res, &[1u64]);

        let res: Vec<u32> =
            SuffixArrayBuilder::slice_u8_to_vec(&[255, 255, 255, 255], 1);
        assert_eq!(res, &[u32::MAX]);

        let res: Vec<u64> = SuffixArrayBuilder::slice_u8_to_vec(
            &[255, 255, 255, 255, 255, 255, 255, 255],
            1,
        );
        assert_eq!(res, &[u64::MAX]);

        let res: Vec<u32> = SuffixArrayBuilder::slice_u8_to_vec(
            &[0, 0, 0, 0, 255, 255, 255, 255],
            2,
        );
        assert_eq!(res, &[0u32, u32::MAX]);

        let res: Vec<u64> = SuffixArrayBuilder::slice_u8_to_vec(
            &[
                0, 0, 0, 0, 0, 0, 0, 0, 255, 255, 255, 255, 255, 255, 255,
                255,
            ],
            2,
        );
        assert_eq!(res, &[0u64, u64::MAX]);

        Ok(())
    }

    #[test]
    fn test_vec_to_slice_u8() -> Result<()> {
        let res = SuffixArrayBuilder::vec_to_slice_u8(&[0u32]);
        assert_eq!(res, &[0, 0, 0, 0]);

        let res = SuffixArrayBuilder::vec_to_slice_u8(&[0u64]);
        assert_eq!(res, &[0, 0, 0, 0, 0, 0, 0, 0]);

        let res = SuffixArrayBuilder::vec_to_slice_u8(&[1u32]);
        assert_eq!(res, &[1, 0, 0, 0]);

        let res = SuffixArrayBuilder::vec_to_slice_u8(&[1u64]);
        assert_eq!(res, &[1, 0, 0, 0, 0, 0, 0, 0]);

        let res = SuffixArrayBuilder::vec_to_slice_u8(&[u32::MAX]);
        assert_eq!(res, &[255, 255, 255, 255]);

        let res = SuffixArrayBuilder::vec_to_slice_u8(&[u64::MAX]);
        assert_eq!(res, &[255, 255, 255, 255, 255, 255, 255, 255]);

        let res = SuffixArrayBuilder::vec_to_slice_u8(&[0u32, u32::MAX]);
        assert_eq!(res, &[0, 0, 0, 0, 255, 255, 255, 255]);

        let res = SuffixArrayBuilder::vec_to_slice_u8(&[0u64, u64::MAX]);
        assert_eq!(
            res,
            &[
                0, 0, 0, 0, 0, 0, 0, 0, 255, 255, 255, 255, 255, 255, 255,
                255
            ]
        );

        Ok(())
    }

    #[test]
    fn test_compare() -> Result<()> {
        // 0  1  2  3  4  5  6  7  8  9 10 11 12 13 14
        // A  A  B  A  B  A  B  A  B  B  A  B  A  B  #
        let seq_file = "tests/inputs/abba.fa";
        let seq_data = read_sequence_file(seq_file)?;
        let max_context: Option<u32> = None;
        let is_dna = false;
        let start_positions: Vec<_> =
            seq_data.start_positions.iter().map(|&v| v as u32).collect();
        let len = seq_data.seq.len() as u32;
        let num_partitions = 1;
        let sa: SuffixArrayBuilder<u32> = SuffixArrayBuilder::new(
            seq_data.seq,
            len,
            max_context,
            is_dna,
            start_positions,
            seq_data.headers,
            num_partitions,
        )?;

        // Compare to B to B with no skip
        let query = "B".as_bytes();
        let res = sa.compare(query, 13, 0);
        assert_eq!(res.cmp, Ordering::Equal);
        assert_eq!(res.lcp, 1);

        // Compare to B to B with skip = 1
        let query = "B".as_bytes();
        let res = sa.compare(query, 13, 1);
        assert_eq!(res.cmp, Ordering::Equal);
        assert_eq!(res.lcp, 1);

        // Compare to B to AB
        let query = "B".as_bytes();
        let res = sa.compare(query, 12, 0);
        assert_eq!(res.cmp, Ordering::Greater);
        assert_eq!(res.lcp, 0);

        // Compare to ABABA to ABBABAB#
        let query = "ABABA".as_bytes();
        let res = sa.compare(query, 7, 2);
        assert_eq!(res.cmp, Ordering::Less);
        assert_eq!(res.lcp, 2);

        // Compare to ABAB to ABABBABAB#
        let query = "ABABA".as_bytes();
        let res = sa.compare(query, 5, 2);
        assert_eq!(res.cmp, Ordering::Less);
        assert_eq!(res.lcp, 4);

        Ok(())
    }

    #[test]
    fn test_search() -> Result<()> {
        // 14: #
        //  0: AABABABABBABAB#
        // 12: AB#
        // 10: ABAB#
        //  1: ABABABABBABAB#
        //  3: ABABABBABAB#
        //  5: ABABBABAB#
        //  7: ABBABAB#
        // 13: B#
        // 11: BAB#
        //  9: BABAB#
        //  2: BABABABBABAB#
        //  4: BABABBABAB#
        //  6: BABBABAB#
        //  8: BBABAB#
        let seq_file = "tests/inputs/abba.fa";
        let seq_data = read_sequence_file(seq_file)?;
        let max_context: Option<u32> = None;
        let is_dna = false;
        let start_positions: Vec<_> =
            seq_data.start_positions.iter().map(|&v| v as u32).collect();
        let len = seq_data.seq.len() as u32;
        let num_partitions = 1;
        let sa: SuffixArrayBuilder<u32> = SuffixArrayBuilder::new(
            seq_data.seq,
            len,
            max_context,
            is_dna,
            start_positions,
            seq_data.headers,
            num_partitions,
        )?;

        let res = sa.search("B");
        assert_eq!(res, Some(8..14));

        let res = sa.search("AB");
        assert_eq!(res, Some(2..7));

        let res = sa.search("BABAB");
        assert_eq!(res, Some(10..12));

        let res = sa.search("ABAA");
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
        assert_eq!(len, 16);
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
