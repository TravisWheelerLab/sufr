use log::info;
use rand::seq::SliceRandom;
use rayon::prelude::*;
use std::{
    cmp::{max, min, Ordering},
    fmt::Debug,
    fmt::Display,
    mem,
    ops::{Add, Sub},
    time::Instant,
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
    //pub fn check_order(&self, sa: &[T]) -> Vec<T> {
    //    let mut errors = vec![];
    //    for window in sa.windows(2) {
    //        if let [prev, cur] = window {
    //            if !self.is_less(*prev, *cur) {
    //                errors.push(*prev);
    //            }
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
    pub fn find_lcp_unchecked(s1: &[u8], s2: &[u8], len: T) -> T {
        for i in 0..len.to_usize() {
            unsafe {
                if s1.get_unchecked(i) != s2.get_unchecked(i) {
                    return T::from_usize(i);
                }
            }
        }

        len
    }

    #[inline(always)]
    pub fn find_lcp(&self, start1: T, start2: T, len: T) -> T {
        let start1 = start1.to_usize();
        let start2 = start2.to_usize();
        let len = len.to_usize();
        let end1 = min(start1 + len, self.len.to_usize());
        let end2 = min(start2 + len, self.len.to_usize());
        //T::from_usize(
        //    (start1..end1)
        //        .zip(start2..end2)
        //        .take_while(|(a, b)| self.text[*a] == self.text[*b])
        //        .count(),
        //)
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
    pub fn partition_by_pivot(&self, num_partitions: usize) -> Vec<T> {
        // Read the input
        let now = Instant::now();
        let suffixes: Vec<T> = if self.ignore_start_n {
            // Text will already be "uppercase" (but bytes), N = 78
            self.text
                .iter()
                .enumerate()
                .filter_map(|(i, &c)| (c != 78).then_some(T::from_usize(i)))
                .collect()
        } else {
            (0..self.text.len()).map(T::from_usize).collect()
        };
        info!("Read input in {:?}", now.elapsed());

        // Select random pivots
        let now = Instant::now();
        let pivot_sa =
            self.select_random_pivots(&suffixes, num_partitions - 1);
        info!(
            "Selected/sorted {} pivots in {:?}",
            pivot_sa.len(),
            now.elapsed()
        );

        // Partition the suffixes using the pivots
        let now = Instant::now();
        let partitions = self.partition_by_pivots(suffixes, pivot_sa);
        let sizes: Vec<_> = partitions.iter().map(|p| p.len()).collect();
        info!(
            "Split into {num_partitions} partitions (avg {}) in {:?}",
            sizes.iter().sum::<usize>() / sizes.len(),
            now.elapsed()
        );

        // Sort the partitions
        let now = Instant::now();
        let res: Vec<_> = partitions
            .into_par_iter()
            .filter(|v| !v.is_empty())
            //.map(|part_sa| {
            //    // Iterative
            //    let len = part_sa.len();
            //    let target_sa = part_sa.clone();
            //    let source_lcp = vec![T::default(); len];
            //    let target_lcp = vec![T::default(); len];
            //    let mut sa = vec![part_sa, target_sa];
            //    let mut lcp = vec![source_lcp, target_lcp];
            //    let source = 0;
            //    let target = 1;
            //    let final_target = self.iter_merge_sort(
            //        &mut sa,
            //        &mut lcp,
            //        source,
            //        target,
            //        len - 1,
            //    );
            //    sa[final_target].clone()
            //})
            .map(|mut part_sa| {
                // Recursive
                let len = part_sa.len();
                let mut sa_w = part_sa.clone();
                let mut lcp = vec![T::default(); len];
                let mut lcp_w = vec![T::default(); len];
                self.recursive_merge_sort(
                    &mut sa_w,
                    &mut part_sa,
                    len,
                    &mut lcp,
                    &mut lcp_w,
                );
                part_sa
            })
            .flatten()
            .collect();

        info!("Sorted/merged partitions in {:?}", now.elapsed());

        res
    }

    // --------------------------------------------------
    #[inline(always)]
    fn select_random_pivots(
        &self,
        suffixes: &[T],
        num_pivots: usize,
    ) -> Vec<T> {
        let mut rng = &mut rand::thread_rng();
        let mut pivot_sa: Vec<_> = suffixes
            .choose_multiple(&mut rng, num_pivots)
            .cloned()
            .collect();

        // Iterative
        //let len = pivot_sa.len();
        //let target_sa = pivot_sa.clone();
        //let source_lcp = vec![T::default(); len];
        //let target_lcp = vec![T::default(); len];
        //let mut sa = vec![pivot_sa, target_sa];
        //let mut lcp = vec![source_lcp, target_lcp];
        //let source = 0;
        //let target = 1;
        //let final_target =
        //self.iter_merge_sort(&mut sa, &mut lcp, source, target, len - 1);
        //sa[final_target].clone()

        // Recursive
        let len = pivot_sa.len();
        let mut sa_w = pivot_sa.clone();
        let mut lcp = vec![T::default(); len];
        let mut lcp_w = vec![T::default(); len];
        self.recursive_merge_sort(
            &mut sa_w,
            &mut pivot_sa,
            len,
            &mut lcp,
            &mut lcp_w,
        );
        pivot_sa
    }

    // --------------------------------------------------
    #[inline(always)]
    fn partition_by_pivots(
        &self,
        suffixes: Vec<T>,
        pivot_sa: Vec<T>,
    ) -> Vec<Vec<T>> {
        let parts: Vec<_> = suffixes
            .par_iter()
            .map(|pos| {
                self.upper_bound(*pos, &pivot_sa).unwrap_or(T::default())
            })
            .collect();

        let mut partitions: Vec<Vec<T>> = vec![vec![]; pivot_sa.len() + 1];
        for (pos, part) in suffixes.iter().zip(&parts) {
            partitions[part.to_usize()].push(*pos);
        }
        partitions
    }

    // --------------------------------------------------
    pub fn recursive_merge_sort(
        self: &Self,
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
            self.recursive_merge_sort(
                &mut y[..mid],
                &mut x[..mid],
                mid,
                &mut lcp_w[..mid],
                &mut lcp[..mid],
            );

            self.recursive_merge_sort(
                &mut y[mid..],
                &mut x[mid..],
                n - mid,
                &mut lcp_w[mid..],
                &mut lcp[mid..],
            );

            self.recursive_merge(x, mid, lcp_w, y, lcp);
        }
    }

    // --------------------------------------------------
    fn recursive_merge<'a>(
        self: &Self,
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
        let mut idx_x = 0; // Index into x
        let mut idx_y = 0; // Index into y
        let mut k = 0; // Index into z

        while idx_x < len_x && idx_y < len_y {
            let l_x = lcp_x[idx_x];

            if l_x > m {
                target_sa[k] = x[idx_x];
                target_lcp[k] = l_x;
            } else if l_x < m {
                target_sa[k] = y[idx_y];
                target_lcp[k] = m;
                m = l_x;
            } else {
                // Length of shorter suffix
                let max_n = self.len - max(x[idx_x], y[idx_y]);

                // Prefix-context length for the suffixes
                let context = min(self.max_context, max_n);

                // LCP(X_i, Y_j)
                //let len_lcp = m + Self::find_lcp(
                //    &self.text[(x[idx_x] + m).to_usize()..],
                //    &self.text[(y[idx_y] + m).to_usize()..],
                //    context - m,
                //);
                let len_lcp = m + self.find_lcp(
                    x[idx_x] + m,
                    y[idx_y] + m,
                    context - m,
                );

                // If the len of the LCP is the entire shorter
                // sequence, take that.
                if len_lcp == max_n {
                    target_sa[k] = max(x[idx_x], y[idx_y])
                }
                // Else, look at the next char after the LCP
                // to determine order.
                else if self.text[(x[idx_x] + len_lcp).to_usize()]
                    < self.text[(y[idx_y] + len_lcp).to_usize()]
                {
                    target_sa[k] = x[idx_x]
                }
                // Else take from the right
                else {
                    target_sa[k] = y[idx_y]
                }

                // If we took from the right...
                if target_sa[k] == x[idx_x] {
                    target_lcp[k] = l_x;
                } else {
                    target_lcp[k] = m
                }
                m = len_lcp;
            }

            if target_sa[k] == x[idx_x] {
                idx_x += 1;
            } else {
                idx_y += 1;
                mem::swap(&mut x, &mut y);
                mem::swap(&mut len_x, &mut len_y);
                mem::swap(&mut lcp_x, &mut lcp_y);
                mem::swap(&mut idx_x, &mut idx_y);
            }

            k += 1;
        }

        // Copy rest of the data from X to Z.
        while idx_x < len_x {
            target_sa[k] = x[idx_x];
            target_lcp[k] = lcp_x[idx_x];
            idx_x += 1;
            k += 1;
        }

        // Copy rest of the data from Y to Z.
        if idx_y < len_y {
            target_sa[k] = y[idx_y];
            target_lcp[k] = m;
            idx_y += 1;
            k += 1;

            while idx_y < len_y {
                target_sa[k] = y[idx_y];
                target_lcp[k] = lcp_y[idx_y];
                idx_y += 1;
                k += 1;
            }
        }
    }

    // --------------------------------------------------
    fn _iter_merge_sort(
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
                self._iter_merge(sa, lcp, source, target, from, mid, to);
            }
            // m = [1, 2, 4, 8, 16...]
            m *= 2;
        }

        // Return the index of the final target
        target
    }

    // --------------------------------------------------
    fn _iter_merge(
        &self,
        sa: &mut [Vec<T>],
        lcp: &mut [Vec<T>],
        source_idx: usize,
        target_idx: usize,
        from: usize,
        mid: usize,
        to: usize,
    ) {
        let mut m = T::from_usize(0); // Last LCP from left side (x)
        let mut idx_x = from; // Index into x
        let mut idx_y = mid; // Index into y
        let mut k = from; // Index into target
        let mut take_x = mid - from;
        let mut take_y = to - mid + 1;

        while take_x > 0 && take_y > 0 {
            let l_x = lcp[source_idx][idx_x];
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
                    //let len_lcp = m + Self::find_lcp(
                    //    &self.text[(sa[source_idx][idx_x] + m).to_usize()..],
                    //    &self.text[(sa[source_idx][idx_y] + m).to_usize()..],
                    //    context - m,
                    //);
                    let len_lcp = m + self.find_lcp(
                        sa[source_idx][idx_x] + m,
                        sa[source_idx][idx_y] + m,
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
    use super::SuffixArray;
    use anyhow::Result;

    #[test]
    fn test_upper_bound() -> Result<()> {
        //          012345
        let text = "TTTAGC".as_bytes().to_vec();
        let len = text.len();
        let max_context: Option<u32> = None;
        let ignore_start_n = false;
        let sa: SuffixArray<u32> =
            SuffixArray::new(text, len as u32, max_context, ignore_start_n);

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
        let ignore_start_n = false;
        let sa: SuffixArray<u32> =
            SuffixArray::new(text, len as u32, max_context, ignore_start_n);

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
