use anyhow::{anyhow, Result};
use log::info;
use rand::seq::SliceRandom;
use rayon::prelude::*;
use std::{
    cmp::{max, min, Ordering},
    fmt::Debug,
    fmt::Display,
    fs::File,
    io::{BufWriter, Write},
    mem,
    ops::{Add, Sub},
    slice,
    time::Instant,
};

const OUTFILE_VERSION: u8 = 1;

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
    + Copy
    + Default
    + Display
    + Ord
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

// --------------------------------------------------
#[derive(Clone, Debug)]
pub struct SuffixArrayBuilder<T>
where
    T: Int + FromUsize<T> + Sized + Send + Sync,
{
    // The original text stored as bytes
    pub text: Vec<u8>,

    // The length of the original text
    pub len: T,

    // The maximum length when comparing suffixes
    pub max_context: T,

    // Whether or not to skip suffixes that start with N
    pub is_dna: bool,
}

impl<T> SuffixArrayBuilder<T>
where
    T: Int + FromUsize<T> + Sized + Send + Sync,
{
    // --------------------------------------------------
    pub fn new(
        text: Vec<u8>,
        len: T,
        max_context: Option<T>,
        is_dna: bool,
    ) -> SuffixArrayBuilder<T> {
        SuffixArrayBuilder {
            text,
            len,
            max_context: max_context.unwrap_or(len),
            is_dna,
        }
    }

    // --------------------------------------------------
    pub fn check_order(&self, sa: &[T]) -> Vec<T> {
        let mut errors = vec![];
        for window in sa.windows(2) {
            if let [prev, cur] = window {
                if !self.is_less(*prev, *cur) {
                    errors.push(*prev);
                }
            }
        }
        errors
    }

    // --------------------------------------------------
    pub fn check_lcp(&self, sa: &[T], lcp: &[T]) -> Vec<T> {
        let mut errors = vec![];
        for i in 1..sa.len() {
            let len_lcp = self.find_lcp(sa[i], sa[i - 1], self.max_context);
            if len_lcp != lcp[i] {
                errors.push(T::from_usize(i));
            }
        }
        errors
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
    //#[inline(always)]
    //pub fn find_lcp_unchecked(s1: &[u8], s2: &[u8], len: T) -> T {
    //    for i in 0..len.to_usize() {
    //        unsafe {
    //            if s1.get_unchecked(i) != s2.get_unchecked(i) {
    //                return T::from_usize(i);
    //            }
    //        }
    //    }

    //    len
    //}

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
    pub fn upper_bound(&self, target: T, sa: &[T]) -> Option<T> {
        // See if target is less than the first element
        if self.is_less(target, sa[0]) {
            None
        } else {
            // Find where all the values are less than target
            let i = sa.partition_point(|&p| self.is_less(p, target));

            // If the value at the partition is the same as the target
            if sa.get(i).map_or(false, |&v| v == target) {
                // Then return the next value, which might be out of range
                Some(T::from_usize(i + 1))
            } else {
                // Else return the partition point
                Some(T::from_usize(i))
            }
        }
    }

    // --------------------------------------------------
    pub fn sort(&self, num_partitions: usize) -> (Vec<T>, Vec<T>) {
        let now = Instant::now();
        let suffixes: Vec<T> = if self.is_dna {
            // Text will already be "uppercase" (but bytes)
            self.text
                .iter()
                .enumerate()
                .filter_map(|(i, c)| {
                    b"ACGT".contains(c).then_some(T::from_usize(i))
                })
                .collect()
        } else {
            (0..self.text.len()).map(T::from_usize).collect()
        };
        info!(
            "Constructed unsorted suffix array of len {} in {:?}",
            suffixes.len(),
            now.elapsed()
        );

        // Select random pivots
        let now = Instant::now();
        let pivot_sa = self.select_pivots(&suffixes, num_partitions - 1);
        info!(
            "Selected/sorted {} pivots in {:?}",
            pivot_sa.len(),
            now.elapsed()
        );

        // Partition the suffixes using the pivots
        let now = Instant::now();
        let partitions = self.partition(suffixes, pivot_sa);
        let sizes: Vec<_> = partitions.iter().map(|p| p.len()).collect();
        info!(
            "Split into {num_partitions} partitions (avg {}) in {:?}",
            sizes.iter().sum::<usize>() / sizes.len(),
            now.elapsed()
        );

        // Sort partitions
        let now = Instant::now();
        let sorted: Vec<_> = partitions
            .into_par_iter()
            .filter(|v| !v.is_empty())
            .map(|mut part_sa| {
                let len = part_sa.len();
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
                (part_sa, lcp)
            })
            .collect();
        info!("Sorted partitions in {:?}", now.elapsed());

        // Concatenate the partitioned SA/LCP vectors
        let now = Instant::now();
        let sa: Vec<_> = sorted.iter().flat_map(|t| t.0.clone()).collect();
        let mut lcp: Vec<_> =
            sorted.into_iter().flat_map(|t| t.1.clone()).collect();
        info!("Concatenated partitions in {:?}", now.elapsed());

        // Fix the LCP boundaries
        let now = Instant::now();
        let mut pos = 0;
        for size in sizes.into_iter().filter(|&v| v > 0) {
            pos += size;

            if pos >= sa.len() {
                break;
            }
            if lcp[pos] == T::default() {
                lcp[pos] =
                    self.find_lcp(sa[pos - 1], sa[pos], self.max_context);
            }
        }
        info!("Fixed LCP boundaries in {:?}", now.elapsed());

        (sa, lcp)
    }

    // --------------------------------------------------
    #[inline(always)]
    fn select_pivots(&self, suffixes: &[T], num_pivots: usize) -> Vec<T> {
        let mut rng = &mut rand::thread_rng();
        let mut pivot_sa: Vec<_> = suffixes
            .choose_multiple(&mut rng, num_pivots)
            .cloned()
            .collect();
        let mut sa_w = pivot_sa.clone();
        let len = pivot_sa.len();
        let mut lcp = vec![T::default(); len];
        let mut lcp_w = vec![T::default(); len];
        self.merge_sort(&mut sa_w, &mut pivot_sa, len, &mut lcp, &mut lcp_w);
        pivot_sa
    }

    // --------------------------------------------------
    #[inline(always)]
    fn partition(&self, suffixes: Vec<T>, pivot_sa: Vec<T>) -> Vec<Vec<T>> {
        // Find the highest partition for each suffix
        let parts: Vec<_> = suffixes
            .par_iter()
            .map(|pos| self.upper_bound(*pos, &pivot_sa).unwrap_or_default())
            .collect();

        // Copy the suffixes into the correct partition
        let mut partitions: Vec<Vec<T>> = vec![vec![]; pivot_sa.len() + 1];
        for (suffix, part) in suffixes.iter().zip(&parts) {
            partitions[part.to_usize()].push(*suffix);
        }

        partitions
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
    pub fn write(&self, sa: &[T], lcp: &[T], outfile: &str) -> Result<()> {
        let mut out = BufWriter::new(
            File::create(outfile).map_err(|e| anyhow!("{outfile}: {e}"))?,
        );

        // Header: version/is_dna
        let is_dna: u8 = if self.is_dna { 1 } else { 0 };
        let _ = out.write(&[OUTFILE_VERSION, is_dna]);

        // Suffix array length
        let _ = out.write(&usize_to_bytes(sa.len()))?;

        // Write out suffix array/LCP as raw bytes
        let slice_sa: &[u8] = unsafe {
            slice::from_raw_parts(
                sa.as_ptr() as *const _,
                sa.len() * mem::size_of::<T>(),
            )
        };
        out.write_all(slice_sa)?;

        let slice_lcp: &[u8] = unsafe {
            slice::from_raw_parts(
                lcp.as_ptr() as *const _,
                lcp.len() * mem::size_of::<T>(),
            )
        };
        out.write_all(slice_lcp)?;

        Ok(())
    }
}

// --------------------------------------------------
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
    use super::SuffixArrayBuilder;
    use anyhow::Result;

    #[test]
    fn test_upper_bound() -> Result<()> {
        //          012345
        let text = "TTTAGC".as_bytes().to_vec();
        let len = text.len();
        let max_context: Option<u32> = None;
        let is_dna = false;
        let sa: SuffixArrayBuilder<u32> =
            SuffixArrayBuilder::new(text, len as u32, max_context, is_dna);

        // The suffix "AGC$" is found before "GC$" and "C$
        assert_eq!(sa.upper_bound(3, &[5, 4]), None);

        // The suffix "TAGC$" is beyond all the values
        assert_eq!(sa.upper_bound(2, &[3, 5, 4]), Some(3));

        // The "C$" is the last value
        assert_eq!(sa.upper_bound(5, &[3, 5, 4]), Some(2));

        Ok(())
    }

    //#[test]
    //fn test_sort_subarrays() -> Result<()> {
    //    let text = Cursor::new("AACTGCGGAT$");
    //    let max_context = None;
    //    let ignore_start_n = false;
    //    let suf_arr = SuffixArray::new(text, max_context, ignore_start_n);

    //    // Ensure we get two subarrays
    //    let subs = suf_arr.sort_subarrays(2);
    //    assert_eq!(subs.len(), 2);

    //    for (sa, _lcp, _) in subs {
    //        let suffixes: Vec<String> =
    //            sa.iter().map(|&p| suf_arr._string_at(p)).collect();

    //        // Check that suffixes are correctly ordered
    //        for pair in suffixes.windows(2) {
    //            if let [a, b] = pair {
    //                assert!(a < b);
    //            }
    //        }
    //    }

    //    Ok(())
    //}

    #[test]
    fn test_find_lcp() -> Result<()> {
        //          012345
        let text = "TTTAGC".as_bytes().to_vec();
        let len = text.len();
        let max_context: Option<u32> = None;
        let is_dna = true;
        let sa: SuffixArrayBuilder<u32> =
            SuffixArrayBuilder::new(text, len as u32, max_context, is_dna);

        assert_eq!(sa.find_lcp(0 as u32, 1 as u32, 1 as u32), 1);
        assert_eq!(sa.find_lcp(0 as u32, 1 as u32, 10 as u32), 2);

        //    assert_eq!(
        //        find_lcp(
        //            "A".to_string().as_bytes(),
        //            "A".to_string().as_bytes(),
        //            1
        //        ),
        //        1
        //    );

        //    assert_eq!(
        //        find_lcp(
        //            "A".to_string().as_bytes(),
        //            "AA".to_string().as_bytes(),
        //            1
        //        ),
        //        1
        //    );

        //    assert_eq!(
        //        find_lcp(
        //            "AA".to_string().as_bytes(),
        //            "AAC".to_string().as_bytes(),
        //            3
        //        ),
        //        2
        //    );

        //    assert_eq!(
        //        find_lcp(
        //            "AC".to_string().as_bytes(),
        //            "ACA".to_string().as_bytes(),
        //            2
        //        ),
        //        2
        //    );

        Ok(())
    }
}
