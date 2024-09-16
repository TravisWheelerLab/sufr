use format_num::NumberFormat;
use log::info;
use rayon::prelude::*;
use std::{
    cmp::{max, min, Ordering},
    fmt::Display,
    mem,
    ops::{Add, Range, Sub},
};

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
    Add<Output = Self> + Sub<Output = Self> + Copy + Default + Display + Ord
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
pub struct SuffixArray<T>
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
    pub ignore_start_n: bool,
}

impl<T> SuffixArray<T>
where
    T: Int + FromUsize<T> + Sized + Send + Sync,
{
    // --------------------------------------------------
    pub fn new(
        text: Vec<u8>,
        len: T,
        max_context: Option<T>,
        ignore_start_n: bool,
    ) -> SuffixArray<T> {
        SuffixArray {
            text,
            len,
            max_context: max_context.unwrap_or(len),
            ignore_start_n,
        }
    }

    // --------------------------------------------------
    #[allow(dead_code)]
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
    #[inline(always)]
    fn find_lcp(s1: &[u8], s2: &[u8], len: T) -> T {
        for i in 0..len.to_usize() {
            unsafe {
                if s1.get_unchecked(i) != s2.get_unchecked(i) {
                    return T::from_usize(i);
                }
            }
        }

        len
    }
    // --------------------------------------------------
    fn convert_slices_to_vecs(vec_of_slices: Vec<&[T]>) -> Vec<Vec<T>> {
        vec_of_slices
            .into_iter() // Convert Vec<&[usize]> into an iterator
            .map(|slice| slice.to_vec()) // Convert each slice into a Vec<usize>
            .collect() // Collect into a Vec<Vec<usize>>
    }

    // --------------------------------------------------
    fn transpose(matrix: Vec<Vec<Option<&[T]>>>) -> Vec<Vec<&[T]>> {
        // Determine the number of columns (max length of rows)
        let num_cols = matrix.iter().map(|row| row.len()).max().unwrap_or(0);

        // Collect transposed rows, unpack the Options to remove None values
        (0..num_cols)
            .map(|col| {
                matrix
                    .iter()
                    .filter_map(|row| row.get(col))
                    .cloned()
                    .flatten()
                    .filter(|v| !v.is_empty())
                    .collect()
            })
            .collect()
    }

    // --------------------------------------------------
    // Assumes pos is always found -- danger
    // TODO: Remove?
    pub fn _string_at(&self, pos: usize) -> String {
        self.text
            .get(pos..)
            .map(|v| String::from_utf8(v.to_vec()).unwrap())
            .unwrap()
    }

    // --------------------------------------------------
    pub fn is_less(&self, s1: T, s2: T) -> bool {
        let lcp = (s1.to_usize()..self.max_context.to_usize())
            .zip(s2.to_usize()..self.max_context.to_usize())
            .take_while(|(a, b)| self.text[*a] == self.text[*b])
            .count();

        match (
            self.text.get(s1.to_usize() + lcp),
            self.text.get(s2.to_usize() + lcp),
        ) {
            (Some(a), Some(b)) => a < b,
            (None, Some(_)) => true,
            _ => false,
        }
    }

    // --------------------------------------------------
    fn upper_bound(&self, target: T, sa: &[T]) -> Option<T> {
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
    pub fn locate_pivots(
        &self,
        sub_sas: &Vec<Vec<T>>,
        pivots: Vec<T>,
    ) -> Vec<Vec<Option<Range<usize>>>> {
        // The pivots divide the SA into
        // sections, so we may have a section
        // before the first pivot and one after
        // the last such that two pivots may
        // divide the SA into 3 parts like so:
        // 0 ----- P1 -- P2 ------ END
        // We could also have no values that fall
        // into a pivot range so that we get just
        // something at the end:
        // 0 <None> P1 <None> P2 ------ END
        // Or something at the start:
        // 0 ----- P1 <None> P2 <None> END
        // Or nothing at all:
        // 0 <None> P1 <None> P2 <None> END
        sub_sas
            //.iter()
            .into_par_iter()
            .map(|sub_sa| {
                let mut sub_locs: Vec<Option<Range<usize>>> = vec![];
                let mut prev_end: Option<T> = None;
                let mut exhausted = false;

                for &pivot in &pivots {
                    if exhausted {
                        // Add None once a pivot has consumed the suffixes
                        sub_locs.push(None);
                    } else {
                        let found = self.upper_bound(pivot, sub_sa);
                        match found {
                            Some(i) => {
                                sub_locs.push(Some(
                                    prev_end.map_or(0, |p| p.to_usize())
                                        ..i.to_usize(),
                                ));

                                // Check if sub SA is exhausted
                                exhausted = i == T::from_usize(sub_sa.len());
                            }
                            _ => sub_locs.push(None),
                        }
                        prev_end = found;
                    }
                }

                let last_index = prev_end.unwrap_or(T::from_usize(0));
                if last_index < T::from_usize(sub_sa.len()) {
                    sub_locs.push(Some(last_index.to_usize()..sub_sa.len()));
                } else {
                    sub_locs.push(None);
                }

                sub_locs
            })
            .collect()
    }

    // --------------------------------------------------
    // TODO: Take SA/LCPs as mutable
    #[allow(unused_assignments)]
    pub fn merge_part_subs(
        &self,
        part_sas: &[Vec<&[T]>],
        part_lcps: &[Vec<&[T]>],
    ) -> Vec<T> {
        //let pairs: Vec<_> = part_sas.iter().zip(part_lcps).collect();
        //let merged_subs: Vec<_> = pairs
        let merged_subs: Vec<_> = (0..part_sas.len())
            .into_par_iter()
            .map(|i| {
                //.map(|(target_sa, target_lcp)| {
                // TODO: Avoid allocation here?
                //let lens: Vec<_> =
                //    part_sas[i].iter().map(|v| v.len()).collect();
                //dbg!(lens);
                //let tmp: Vec<usize> = *part_sas[i].iter().concat();
                //dbg!(tmp);

                let mut target_sa =
                    Self::convert_slices_to_vecs(part_sas[i].to_vec());
                let mut target_lcp =
                    Self::convert_slices_to_vecs(part_lcps[i].to_vec());

                // Iteratively merge in pairs
                while target_sa.len() > 1 {
                    let mut i = 0;
                    let mut tmp_sa = vec![];
                    let mut tmp_lcp = vec![];

                    while (2 * i) < target_sa.len() {
                        let first = 2 * i;
                        let second = first + 1;
                        if second >= target_sa.len() {
                            tmp_sa.push(target_sa[first].to_vec());
                            tmp_lcp.push(target_lcp[first].to_vec());
                        } else {
                            let mut source_sa = target_sa[first].to_vec();
                            let mid = source_sa.len();
                            let mut sa2 = target_sa[second].to_vec();
                            source_sa.append(&mut sa2);

                            let mut source_lcp = target_lcp[first].to_vec();
                            let mut lcp2 = target_lcp[second].to_vec();
                            source_lcp.append(&mut lcp2);

                            // Zero out first LCP
                            if let Some(v) = source_lcp.get_mut(0) {
                                *v = T::from_usize(0);
                            }

                            // Create working copies for merge
                            let target_sa = source_sa.clone();
                            let target_lcp = source_lcp.clone();
                            let to = source_sa.len() - 1;
                            let mut sa = vec![source_sa, target_sa];
                            let mut lcp = vec![source_lcp, target_lcp];
                            let source = 0;
                            let target = 1;
                            let from = 0;
                            self.iter_merge(
                                &mut sa, &mut lcp, source, target, from, mid,
                                to,
                            );
                            tmp_sa.push(sa[target].to_vec());
                            tmp_lcp.push(lcp[target].to_vec());
                        }
                        i += 1;
                    }
                    //(target_sa, target_lcp, tmp_sa, tmp_lcp) =
                    //    (tmp_sa, tmp_lcp, target_sa, target_lcp);
                    mem::swap(&mut target_sa, &mut tmp_sa);
                    mem::swap(&mut target_lcp, &mut tmp_lcp);
                    //target_sa = mem::take(&mut tmp_sa);
                }
                target_sa.first().map_or(vec![], |v| v.to_vec())
            })
            .filter(|v| !v.is_empty())
            .collect();

        // Concatenate all the merged subs into a single vector
        merged_subs.into_iter().flatten().collect()
    }

    // --------------------------------------------------
    pub fn partition_subarrays<'a>(
        &self,
        sa: &'a [Vec<T>],
        lcp: &'a [Vec<T>],
        pivot_locs: Vec<Vec<Option<Range<usize>>>>,
    ) -> (Vec<Vec<&'a [T]>>, Vec<Vec<&'a [T]>>) {
        // Given the SA/LCP arrays for each of the subarrays
        // and the pivot locations, subdivide the SA/LCP
        // subarrays further.
        let lcp_pairs: Vec<_> = lcp.iter().zip(&pivot_locs).collect();
        let ret_lcp: Vec<_> = lcp_pairs
            .into_par_iter()
            .map(|(sub_lcp, sub_pivots)| {
                let mut tmp = vec![];
                for range in sub_pivots {
                    let res = range.clone().map(|r| sub_lcp.get(r).unwrap());
                    tmp.push(res);
                }
                tmp
            })
            .collect();

        let sa_pairs: Vec<_> = sa.iter().zip(&pivot_locs).collect();
        let ret_sa: Vec<_> = sa_pairs
            .into_par_iter()
            .map(|(sub_sa, sub_pivots)| {
                let mut tmp = vec![];
                for range in sub_pivots {
                    let res = range.clone().map(|r| sub_sa.get(r).unwrap());
                    tmp.push(res);
                }
                tmp
            })
            .collect();

        (Self::transpose(ret_sa), Self::transpose(ret_lcp))
    }

    // --------------------------------------------------
    fn gen_lcp(&self, sa: &[T]) -> Vec<T> {
        (0..sa.len())
            .map(|i| {
                if i == 0 {
                    T::from_usize(0)
                } else {
                    let s1 = &self.text[sa[i - 1].to_usize()..];
                    let s2 = &self.text[sa[i].to_usize()..];

                    T::from_usize(
                        s1.iter()
                            .take(self.max_context.to_usize())
                            .zip(s2)
                            .take_while(|(a, b)| a == b)
                            .count(),
                    )
                }
            })
            .collect()
    }

    // --------------------------------------------------
    // TODO: Is it possible to take the sub_pivots mutably?
    // As written, I make temporary copies to use mem::swap
    pub fn select_pivots(&self, mut sub_pivots: Vec<Vec<T>>) -> Vec<T> {
        let num_partitions = sub_pivots.len();

        // The pivots were plucked from their sorted
        // suffixes and have no LCPs, so construct.
        // TODO: make sure we actually need to do this
        let mut sub_lcps: Vec<Vec<T>> = sub_pivots
            .iter()
            .map(|pivots| self.gen_lcp(pivots))
            .collect();

        // This uses a sliding function to merge
        // Given 0, 1, 2, 3, 4, 5, it will
        // Merge 0/1, 2/3, 5 => 0, 1, 2
        // Merge 0/1, 2      => 0, 1
        // Merge 0/1         => 0
        while sub_pivots.len() > 1 {
            let mut i = 0;
            let mut tmp_sa = vec![];
            let mut tmp_lcp = vec![];

            while (2 * i) < sub_pivots.len() {
                let first = 2 * i;
                let second = first + 1;
                if second >= sub_pivots.len() {
                    // TODO: This is the base case, down to just one
                    // partition. Do I have to make copies of the final
                    // pivots or could I do something clever here?
                    tmp_sa.push(sub_pivots[first].to_vec());
                    tmp_lcp.push(sub_lcps[first].to_vec());
                } else {
                    // I have to make copies of the data to pass
                    // mutably into the merge function
                    let mut source_sa = sub_pivots[first].to_vec();
                    let mid = source_sa.len();
                    let mut sa2 = sub_pivots[second].to_vec();
                    source_sa.append(&mut sa2);

                    let mut source_lcp = sub_lcps[first].to_vec();
                    let mut lcp2 = sub_lcps[second].to_vec();
                    source_lcp.append(&mut lcp2);

                    // Create working copies for merge
                    let target_sa = source_sa.clone();
                    let target_lcp = source_lcp.clone();
                    let from = 0;
                    let to = source_sa.len() - 1;
                    let mut sa = vec![source_sa, target_sa];
                    let mut lcp = vec![source_lcp, target_lcp];
                    let source = 0;
                    let target = 1;
                    self.iter_merge(
                        &mut sa, &mut lcp, source, target, from, mid, to,
                    );

                    tmp_sa.push(sa[target].to_vec());
                    tmp_lcp.push(lcp[target].to_vec());
                }
                i += 1;
            }
            mem::swap(&mut sub_pivots, &mut tmp_sa);
            mem::swap(&mut sub_lcps, &mut tmp_lcp);
        }

        sub_pivots
            .first()
            .map(|pivots| self.sample_pivots(pivots, num_partitions - 1))
            .unwrap_or_default()
    }

    // --------------------------------------------------
    pub fn sample_pivots(&self, sa: &[T], num_pivots: usize) -> Vec<T> {
        let step = sa.len() / if num_pivots > 0 { num_pivots } else { 1 };

        (0..num_pivots).map(|i| sa[(i + 1) * step - 1]).collect()
    }

    // --------------------------------------------------
    // TODO: return type struct
    pub fn sort_subarrays(
        &self,
        suggested_num_partitions: usize,
    ) -> Vec<(Vec<T>, Vec<T>, Vec<T>)> {
        let mut suffixes: Vec<T> = if self.ignore_start_n {
            // Text will already be "uppercase" (but bytes), N = 78
            self.text
                .iter()
                .enumerate()
                .filter_map(|(i, &c)| (c != 78).then_some(T::from_usize(i)))
                .collect()
        } else {
            (0..self.text.len()).map(T::from_usize).collect()
        };
        let num_suffixes = suffixes.len();
        let num_parts = if suggested_num_partitions < num_suffixes {
            suggested_num_partitions
        } else if num_suffixes > 500 {
            16
        } else {
            1
        };
        let subset_size = num_suffixes / num_parts;
        let len = num_suffixes as f64;
        let mut pivots_per_part = (32.0 * len.log10()).ceil() as usize;
        if (pivots_per_part > subset_size) || (pivots_per_part < 1) {
            pivots_per_part =
                if subset_size == 1 { 1 } else { subset_size / 2 };
        }

        let num_fmt = NumberFormat::new();
        info!(
            "{num_parts} partition{} with {} element{} and {pivots_per_part} pivot{} each",
            if subset_size == 1 { "" } else { "s" },
            num_fmt.format(",.0", subset_size as f64),
            if num_parts == 1 { "" } else { "s" },
            if pivots_per_part == 1 { "" } else { "s" },
        );

        //let now = Instant::now();
        let partitions: Vec<_> = (0..num_parts)
            .map(|i| {
                let len = subset_size
                    + if i == num_parts - 1 {
                        num_suffixes % num_parts
                    } else {
                        0
                    };
                suffixes.drain(0..len).collect::<Vec<_>>()
            })
            .collect();
        //println!("Finished making parts in {:?}", now.elapsed());

        partitions
            .into_par_iter()
            .map(|source_sa| {
                //let now = Instant::now();
                let len = source_sa.len();
                let target_sa = source_sa.clone();
                let source_lcp = vec![T::default(); len];
                let target_lcp = vec![T::default(); len];
                let mut sa = vec![source_sa, target_sa];
                let mut lcp = vec![source_lcp, target_lcp];
                //println!("Finished allocs in {:?}", now.elapsed());

                let source = 0;
                let target = 1;
                //let now = Instant::now();
                let final_target = self.iter_merge_sort(
                    &mut sa,
                    &mut lcp,
                    source,
                    target,
                    len - 1,
                );

                //println!("Finished merge in {:?}", now.elapsed());
                let sub_sa = sa[final_target].clone();
                let sub_lcp = lcp[final_target].clone();
                let pivots = self.sample_pivots(&sub_sa, pivots_per_part);

                //if let Ok(mut num) = counter.lock() {
                //    *num += 1;
                //    if *num % 1000 == 0 {
                //        info!("  Sorted {num} subarrays...");
                //    }
                //}

                (sub_sa, sub_lcp, pivots)
            })
            .collect()
    }

    // --------------------------------------------------
    fn iter_merge_sort(
        &self,
        sa: &mut [Vec<T>],
        lcp: &mut [Vec<T>],
        mut source: usize,
        mut target: usize,
        high: usize,
    ) -> usize {
        let mut m = 1;

        // divide the array into blocks of size `m`
        while m <= high {
            // Don't swap the first time
            if m > 1 {
                source = (source + 1) % 2;
                target = (target + 1) % 2;
            }

            for i in (0..=high).step_by(2 * m) {
                let from = i;
                let mut to = i + 2 * m - 1;
                if to > high {
                    to = high;
                }
                let mut mid = i + m;
                if mid > high {
                    mid = high;
                }
                self.iter_merge(sa, lcp, source, target, from, mid, to);
            }
            // m = [1, 2, 4, 8, 16...]
            m *= 2;
        }

        // Return the index of the final target
        target
    }

    // --------------------------------------------------
    fn iter_merge(
        &self,
        sa: &mut [Vec<T>],
        lcp: &mut [Vec<T>],
        source_idx: usize,
        target_idx: usize,
        from: usize,
        mid: usize,
        to: usize,
    ) {
        //println!("\n>>> MERGE FROM {from} MID {mid} TO {to} source");
        let mut m = T::from_usize(0); // Last LCP from left side (x)
        let mut idx_x = from; // Index into x
        let mut idx_y = mid; // Index into y
        let mut k = from; // Index into target
        let mut take_x = mid - from;
        let mut take_y = to - mid + 1;
        //println!(
        //    "idx_x {idx_x} idx_y {idx_y} take_x {take_x} take_y {take_y}"
        //);
        //println!(
        //    "x {}-{} {:?} y {}-{} {:?}",
        //    idx_x,
        //    idx_x + take_x,
        //    &sa[source][idx_x..idx_x + take_x],
        //    idx_y,
        //    idx_y + take_y,
        //    &sa[source][idx_y..idx_y + take_y],
        //);

        while take_x > 0 && take_y > 0 {
            let l_x = lcp[source_idx][idx_x];

            //println!("k {k} end_x {end_x} end_y {end_y} l_x {l_x} m {m}");
            //println!(
            //    "idx_x {} = {:?}",
            //    source_sa[idx_x],
            //    String::from_utf8((&self.text[source_sa[idx_x]..]).to_vec())
            //);
            //println!(
            //    "idx_y {} = {:?}",
            //    source_sa[idx_y],
            //    String::from_utf8((&self.text[source_sa[idx_y]..]).to_vec())
            //);
            match l_x.cmp(&m) {
                Ordering::Greater => {
                    sa[target_idx][k] = sa[source_idx][idx_x];
                    lcp[target_idx][k] = l_x;
                }
                Ordering::Less => {
                    sa[target_idx][k] = sa[source_idx][idx_y];
                    lcp[target_idx][k] = m;
                    m = l_x;
                }
                Ordering::Equal => {
                    // Find the length of shorter suffix
                    let max_n = self.len
                        - max(sa[source_idx][idx_x], sa[source_idx][idx_y]);

                    // Prefix-context length for the suffixes
                    let context = min(self.max_context, max_n);

                    // LCP(X_i, Y_j)
                    let len_lcp = m + Self::find_lcp(
                        &self.text[(sa[source_idx][idx_x] + m).to_usize()..],
                        &self.text[(sa[source_idx][idx_y] + m).to_usize()..],
                        context - m,
                    );

                    // If the len of the LCP is the entire shorter
                    // sequence, take that (the one with the higher SA value)
                    if len_lcp == max_n {
                        sa[target_idx][k] =
                            max(sa[source_idx][idx_x], sa[source_idx][idx_y]);
                    }
                    // Else, look at the next char after the LCP
                    // to determine order.
                    else if self.text
                        [(sa[source_idx][idx_x] + len_lcp).to_usize()]
                        < self.text
                            [(sa[source_idx][idx_y] + len_lcp).to_usize()]
                    {
                        sa[target_idx][k] = sa[source_idx][idx_x];
                    }
                    // Else take from the right
                    else {
                        sa[target_idx][k] = sa[source_idx][idx_y];
                    }

                    // If we took from the right...
                    if sa[target_idx][k] == sa[source_idx][idx_x] {
                        lcp[target_idx][k] = l_x;
                    } else {
                        lcp[target_idx][k] = m
                    }
                    m = len_lcp;
                }
            }

            if sa[target_idx][k] == sa[source_idx][idx_x] {
                idx_x += 1;
                take_x -= 1;
            } else {
                idx_y += 1;
                take_y -= 1;
                // Swap X/Y
                (idx_x, idx_y) = (idx_y, idx_x);
                (take_x, take_y) = (take_y, take_x);
            }

            //println!(
            //    "END LOOP: idx_x {idx_x} idx_y {idx_y} take_x {take_x} take_y {take_y}"
            //);

            k += 1;
        }

        // Copy rest of the data from X to Z.
        while take_x > 0 {
            sa[target_idx][k] = sa[source_idx][idx_x];
            lcp[target_idx][k] = lcp[source_idx][idx_x];
            idx_x += 1;
            take_x -= 1;
            k += 1;
        }

        // Copy rest of the data from Y to Z.
        if take_y > 0 {
            sa[target_idx][k] = sa[source_idx][idx_y];
            lcp[target_idx][k] = m;
            idx_y += 1;
            take_y -= 1;
            k += 1;

            while take_y > 0 {
                sa[target_idx][k] = sa[source_idx][idx_y];
                lcp[target_idx][k] = lcp[source_idx][idx_y];
                idx_y += 1;
                take_y -= 1;
                k += 1;
            }
        }
    }
}

// --------------------------------------------------
//#[inline]
//fn byte_to_base(byte: &u8) -> Option<u8> {
//    match byte {
//        65 | 97 => Some(0),  // A
//        67 | 99 => Some(1),  // C
//        71 | 103 => Some(3), // G
//        84 | 116 => Some(2), // T
//        _ => None,
//    }
//}

// --------------------------------------------------
/// Reads input file into an array of 2-bit encoded values
/// E.g., A = 0, C = 1, G = 2, T = 3
//fn _read_input(mut input: impl BufRead) -> Result<Vec<u8>> {
//    let mut ret = vec![];
//    let start = Instant::now();
//    let buffer_size = 1024 * 128;

//    // Mask with 32 to uppercase ASCII values
//    let ret: Vec<_> = input
//        .bytes()
//        .map_while(Result::ok)
//        .map(|b| b & 0b1011111)
//        .collect();

// 1.05537875s
//let mut buffer = vec![0; buffer_size];
//loop {
//    let bytes_read = input.read(&mut buffer)?;
//    if bytes_read == 0 {
//        break;
//    }
//    ret.extend(
//        buffer[..bytes_read].iter().filter_map(|b| byte_to_base(b)),
//    );
//}

// 1.699109792s in debug, 1.031828417s in release
//let mut buffer = vec![0; buffer_size];
//loop {
//    let bytes_read = input.read(&mut buffer)?;
//    if bytes_read == 0 {
//        break;
//    }
//    for byte in buffer[..bytes_read].iter() {
//        if let Some(checked) = match byte {
//            65 | 97 => Some(0),
//            67 | 99 => Some(1),
//            71 | 103 => Some(3),
//            84 | 116 => Some(2),
//            _ => None,
//        } {
//            ret.push(checked);
//        }
//    }
//}

// 53.436009s in debug, 1.731694375s in release
//let lookup = HashMap::<u8, u8>::from([
//    (65, 0),  // A
//    (97, 0),  // a
//    (67, 1),  // C
//    (99, 1),  // c
//    (71, 3),  // G
//    (103, 3), // g
//    (84, 2),  // T
//    (116, 2), // t
//]);
//let ret: Vec<u8> = input
//    .bytes()
//    .filter_map(|b| b.ok())
//    .filter_map(|b| lookup.get(&b))
//    .cloned()
//    .collect();

// 35.509148416s in debug, 1.198498916s in release
//let lookup = HashMap::from([
//    (65, 0),  // A
//    (97, 0),  // a
//    (67, 1),  // C
//    (99, 1),  // c
//    (71, 3),  // G
//    (103, 3), // g
//    (84, 2),  // T
//    (116, 2), // t
//]);

//let mut buffer = vec![0; 65536];
//loop {
//    let bytes_read = input.read(&mut buffer)?;
//    if bytes_read == 0 {
//        break;
//    }
//    for byte in buffer[..bytes_read].iter() {
//        if let Some(val) = lookup.get(byte) {
//            ret.push(*val);
//        }
//    }
//}

// 10.935326709s in debug, 1.262571875s in release
//for byte in input.bytes() {
//    if let Some(val) = match byte.unwrap() {
//        65 | 97 => Some(0),
//        67 | 99 => Some(1),
//        71 | 103 => Some(3),
//        84 | 116 => Some(2),
//        _ => None,
//    } {
//        ret.push(val);
//    }
//}

// 10.899836542s in debug, 675.702666ms in release
//for byte in input.bytes() {
//    let base = BASE_LOOKUP[byte.unwrap() as usize];
//    if base != 255 {
//        ret.push(base);
//    }
//}

// 18.955611916s in debug, 857.590666ms in release
//let ret: Vec<u8> = input
//    .bytes()
//    .filter_map(|b| b.ok())
//    .map(|b| BASE_LOOKUP[b as usize])
//    .filter(|&b| b < 255)
//    .collect();

//debug!("Read input in {:?}", start.elapsed());

//Ok(ret)
//}

// --------------------------------------------------
#[cfg(test)]
mod tests {
    use super::{find_lcp, SuffixArray};
    use anyhow::Result;
    use std::io::Cursor;

    #[test]
    fn test_upper_bound() {
        //                      012345
        let text = Cursor::new("TTTAGC");
        let max_context = None;
        let ignore_start_n = false;
        let sa = SuffixArray::new(text, max_context, ignore_start_n);

        // The suffix "AGC$" is found before "GC$" and "C$
        assert_eq!(sa.upper_bound(3, &[5, 4]), None);

        // The suffix "TAGC$" is beyond all the values
        assert_eq!(sa.upper_bound(2, &[3, 5, 4]), Some(3));

        // The "C$" is the last value
        assert_eq!(sa.upper_bound(5, &[3, 5, 4]), Some(2));
    }

    #[test]
    fn test_sort_subarrays() -> Result<()> {
        let text = Cursor::new("AACTGCGGAT$");
        let max_context = None;
        let ignore_start_n = false;
        let suf_arr = SuffixArray::new(text, max_context, ignore_start_n);

        // Ensure we get two subarrays
        let subs = suf_arr.sort_subarrays(2);
        assert_eq!(subs.len(), 2);

        for (sa, _lcp, _) in subs {
            let suffixes: Vec<String> =
                sa.iter().map(|&p| suf_arr._string_at(p)).collect();

            // Check that suffixes are correctly ordered
            for pair in suffixes.windows(2) {
                if let [a, b] = pair {
                    assert!(a < b);
                }
            }
        }

        Ok(())
    }

    //#[test]
    //fn test_generate() -> Result<()> {
    //    let max_context = None;

    //    let cursor = Cursor::new("AACTGCGGAT");
    //    let s = SuffixArray::new(cursor, max_context);
    //    let (sa, _lcp) = s.generate();
    //    assert_eq!(sa, &[10, 0, 1, 8, 5, 2, 7, 4, 6, 9, 3,]);

    //    let cursor = Cursor::new("CTCACC");
    //    let s = SuffixArray::new(cursor, max_context);
    //    let (sa, _lcp) = s.generate();
    //    assert_eq!(sa, &[6, 3, 5, 2, 4, 0, 1]);
    //    Ok(())
    //}

    #[test]
    fn test_find_lcp() -> Result<()> {
        assert_eq!(
            find_lcp(
                "A".to_string().as_bytes(),
                "C".to_string().as_bytes(),
                1
            ),
            0
        );

        assert_eq!(
            find_lcp(
                "A".to_string().as_bytes(),
                "A".to_string().as_bytes(),
                1
            ),
            1
        );

        assert_eq!(
            find_lcp(
                "A".to_string().as_bytes(),
                "AA".to_string().as_bytes(),
                1
            ),
            1
        );

        assert_eq!(
            find_lcp(
                "AA".to_string().as_bytes(),
                "AAC".to_string().as_bytes(),
                3
            ),
            2
        );

        assert_eq!(
            find_lcp(
                "AC".to_string().as_bytes(),
                "ACA".to_string().as_bytes(),
                2
            ),
            2
        );

        Ok(())
    }
}
