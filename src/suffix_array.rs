use anyhow::{bail, Result};
use log::{debug, info};
use rayon::prelude::*;
use std::{
    cmp::{max, min},
    io::BufRead,
    mem,
    ops::Range,
};

// --------------------------------------------------
#[derive(Clone, Debug)]
pub struct SuffixArray {
    pub text: Vec<u8>,
    pub len: usize,
    pub max_context: usize,
}

impl SuffixArray {
    pub fn new(
        input: impl BufRead,
        max_context: Option<usize>,
    ) -> SuffixArray {
        // TODO: Is this good and proper? Any use of the results
        // will need to similarly filter the inputs.
        let mut text: Vec<_> = input
            .bytes()
            .map_while(Result::ok)
            .map(|b| b & 0b1011111) // Uppercase
            .filter(|&b| b >= 0b1000001 && b <= 0b1011010) // A-Z
            .collect();
        text.push(b'$');
        let len = text.len();

        SuffixArray {
            text,
            len,
            max_context: max_context.unwrap_or(len),
        }
    }

    pub fn string_at(self: &Self, pos: usize) -> String {
        String::from_utf8(self.text[pos..].to_vec()).unwrap()
    }

    fn upper_bound(self: &Self, target: &str, sa: &[usize]) -> Option<usize> {
        if sa
            .first()
            .map_or(false, |&p| self.string_at(p).as_str() > target)
        {
            None
        } else {
            let i =
                sa.partition_point(|&p| self.string_at(p).as_str() <= target);
            if sa
                .get(i)
                .map_or(false, |&p| self.string_at(p).as_str() == target)
            {
                Some(i + 1)
            } else {
                Some(i)
            }
        }
    }

    // --------------------------------------------------
    pub fn locate_pivots(
        self: &Self,
        sub_sas: &Vec<Vec<usize>>,
        pivot_suffixes: Vec<String>,
    ) -> Vec<Vec<Option<Range<usize>>>> {
        sub_sas
            .into_par_iter()
            .map(|sub_sa| {
                let mut sub_locs = vec![];
                let mut prev_end: Option<usize> = None;

                for suffix in &pivot_suffixes {
                    let found = self.upper_bound(suffix, sub_sa);
                    //let sub: Vec<_> =
                    //    sub_sa.iter().map(|&p| self.string_at(p)).collect();
                    //println!("suffix '{suffix}' in {sub:?} = {found:?}");
                    match found {
                        Some(i) => {
                            sub_locs.push(Some(prev_end.unwrap_or(0)..i))
                        }
                        _ => sub_locs.push(None),
                    }
                    prev_end = found;
                }

                // The last partition
                let last_index = prev_end.unwrap_or(0);
                if last_index < sub_sa.len() {
                    sub_locs.push(Some(last_index..sub_sa.len()));
                } else {
                    sub_locs.push(None);
                }

                sub_locs
            })
            .collect()
    }

    // --------------------------------------------------
    pub fn merge_part_subs(
        self: &Self,
        part_sas: Vec<Vec<&[usize]>>,
        part_lcps: Vec<Vec<&[usize]>>,
    ) -> Vec<usize> {
        let pairs: Vec<_> = part_sas.iter().zip(part_lcps).collect();
        let merged_subs: Vec<_> = pairs
            //.into_par_iter()
            .iter()
            .map(|(part_sa, part_lcp)| {
                //let vals: Vec<_> = part_sa
                //    .iter()
                //    .map(|v| {
                //        v.iter()
                //            .map(|&p| self.string_at(p))
                //            .collect::<Vec<_>>()
                //    })
                //    .collect();
                //dbg!(vals);

                //println!("\n>>> PART_SA {part_sa:?}");
                let mut target_sa = convert_slices_to_vecs(part_sa.to_vec());
                let mut target_lcp =
                    convert_slices_to_vecs(part_lcp.to_vec());

                while target_sa.len() > 1 {
                    let mut i = 0;
                    //println!("i {i} target_sa len {}", target_sa.len());
                    let mut tmp_sa = vec![];
                    let mut tmp_lcp = vec![];
                    while (2 * i) < target_sa.len() {
                        let first = 2 * i;
                        let second = first + 1;
                        if second >= target_sa.len() {
                            //println!(
                            //    "pushing {:?} to target",
                            //    target_sa[first]
                            //);
                            tmp_sa.push(target_sa[first].to_vec());
                            tmp_lcp.push(target_lcp[first].to_vec());
                            //println!("tmp_sa {tmp_sa:?}");
                        } else {
                            //println!("Merging {first}/{second}");
                            let mut sa = target_sa[first].to_vec();
                            let mid = sa.len();
                            let mut sa2 = target_sa[second].to_vec();
                            sa.append(&mut sa2);

                            let mut lcp = target_lcp[first].to_vec();
                            let mut lcp2 = target_lcp[second].to_vec();
                            lcp.append(&mut lcp2);

                            // Zero out first LCP
                            if let Some(v) = lcp.get_mut(0) {
                                *v = 0;
                            }

                            // Create working copies for merge
                            let mut sa_w = sa.clone();
                            let mut lcp_w = lcp.clone();

                            //println!(
                            //    "sa = {sa:?} mid = {mid}, lcp = {lcp:?}"
                            //);

                            self.merge(
                                &mut sa, mid, &mut lcp_w, &mut sa_w, &mut lcp,
                            );

                            //println!("pushing {sa_w:?} to tmp");
                            tmp_sa.push(sa_w);
                            tmp_lcp.push(lcp);
                            //println!("tmp_sa {tmp_sa:?}");
                        }
                        i += 1;
                    }
                    mem::swap(&mut target_sa, &mut tmp_sa);
                    mem::swap(&mut target_lcp, &mut tmp_lcp);
                }
                //println!("target_sa = {target_sa:?}");
                //println!("target_lcp = {target_lcp:?}");
                target_sa.get(0).map_or(vec![], |v| v.to_vec())

                // Make a mutable copy of the first SA/LCP
                // Handle possibly empty
                //let mut sa: Vec<_> =
                //    part_sa.get(0).map_or(vec![], |v| v.to_vec());
                //let mut lcp = part_lcp.get(0).map_or(vec![], |v| v.to_vec());

                //// Zero out the first position as it's no longer useful
                //if let Some(first) = lcp.get_mut(0) {
                //    *first = 0;
                //}

                //debug!("First SA {sa:#?}");
                ////let first_sa: Vec<_> =
                ////    sa.iter().map(|&p| self.string_at(p)).collect();
                ////dbg!(first_sa);

                //// Process the rest, appending to the already done
                //for i in 1..part_sa.len() {
                //    // The "midpoint" is the length of the SA before appending
                //    let mid = sa.len();
                //    debug!("Using mid {mid}");

                //    // Merge the first two SA/LCPs
                //    let mut next_sa = part_sa[i].to_vec();
                //    let mut next_lcp = part_lcp[i].to_vec();
                //    if let Some(first) = next_lcp.get_mut(0) {
                //        *first = 0;
                //    }
                //    debug!("Adding next_sa {next_sa:#?}");
                //    //let next_sa_vals: Vec<_> =
                //    //    next_sa.iter().map(|&p| self.string_at(p)).collect();
                //    //dbg!(next_sa_vals);

                //    sa.append(&mut next_sa);
                //    lcp.append(&mut next_lcp);

                //    // Create working copies for merge
                //    let mut sa_w = sa.clone();
                //    let mut lcp_w = lcp.clone();

                //    //println!("sa = {sa:?} mid = {mid}, lcp = {lcp:?}");
                //    self.merge(&mut sa, mid, &mut lcp_w, &mut sa_w, &mut lcp);
                //    sa = sa_w;
                //    debug!("After merge sa = {sa:#?}");
                //    //let apres_merge: Vec<_> =
                //    //    sa.iter().map(|&p| self.string_at(p)).collect();
                //    //dbg!(&apres_merge);
                //}
                //sa
            })
            .filter(|v| !v.is_empty())
            .collect();

        debug!("merged_subs = {merged_subs:#?}");

        // Concatenate all the merged subs into a single vector
        merged_subs.into_iter().flatten().collect()
    }

    pub fn check(vals: &Vec<String>) -> Result<()> {
        let mut prev: Option<String> = None;
        for (i, cur) in vals.iter().enumerate() {
            if let Some(p) = prev {
                if p > *cur {
                    bail!(">>> {i} not sorted <<<")
                }
            }
            prev = Some(cur.clone());
        }
        Ok(())
    }

    pub fn partition_subarrays<'a>(
        self: &Self,
        sa: &'a Vec<Vec<usize>>,
        lcp: &'a Vec<Vec<usize>>,
        pivot_locs: Vec<Vec<Option<Range<usize>>>>,
    ) -> (Vec<Vec<&'a [usize]>>, Vec<Vec<&'a [usize]>>) {
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

        (transpose(ret_sa), transpose(ret_lcp))
    }

    pub fn select_pivots(self: &Self, sa: &Vec<Vec<usize>>) -> Vec<String> {
        let n = sa.len();
        let len = self.len as f64;
        let pivots_per_part =
            max(1, min((32.0 * len.log10()).ceil() as usize, self.len / n));
        //let pivots_per_part = 2;
        debug!("pivots_per_part {pivots_per_part}");

        let mut pivots: Vec<_> = sa
            .into_par_iter()
            .flat_map(|s| self.sample_pivots(s, pivots_per_part))
            .collect();
        //dbg!(&pivots);

        let num_pivots = pivots.len();
        let mut pivots_w = pivots.clone();
        let mut temp1 = vec![0; num_pivots];
        let mut temp2 = vec![0; num_pivots];
        self.merge_sort(
            &mut pivots,
            &mut pivots_w,
            num_pivots,
            &mut temp1,
            &mut temp2,
        );

        //println!("pivots after merge_sort = {:?}", &pivots);
        //println!("pivots_w after merge_sort = {:?}", &pivots_w);

        let num_partitions = sa.len();
        let final_pivots = self.sample_pivots(&pivots_w, num_partitions - 1);
        let suffixes: Vec<String> =
            final_pivots.iter().map(|&p| self.string_at(p)).collect();
        suffixes
    }

    pub fn sample_pivots(
        self: &Self,
        sa: &[usize],
        pivots_per_part: usize,
    ) -> Vec<usize> {
        let gap = sa.len()
            / if pivots_per_part > 0 {
                pivots_per_part
            } else {
                1
            };
        (0..pivots_per_part)
            .map(|i| sa[(i + 1) * gap - 1])
            .collect()
    }

    // TODO: Make this part of initialization? I don't like having
    // to worry about this when reusing SuffixArray for sort_subarrays
    fn calc_num_partitions(self: &Self, suggested: usize) -> usize {
        if suggested < self.len {
            suggested
        } else if self.len > 500 {
            100
        } else {
            1
        }
    }

    pub fn sort_subarrays(
        self: &Self,
        num_partitions: usize,
    ) -> Vec<(Vec<usize>, Vec<usize>)> {
        // Handle very small inputs
        // TODO: Better to have a max partition size?
        let text_len = self.len;
        let n = self.calc_num_partitions(num_partitions);
        let subset_size = text_len / n;
        info!(
            "Using {n} partition{} of {subset_size} elements",
            if n == 1 { "" } else { "s" }
        );

        // This seems weird to reuse SuffixArray in this loop
        // but I need access to the "text" inside the merge_sort
        let par_iter = (0..n).into_par_iter().map(|i| {
            let start = i * subset_size;
            let len = subset_size + if i == n - 1 { text_len % n } else { 0 };
            let tmp = self.clone();
            let (sub_sa, sub_lcp) = tmp.generate(start..start + len);
            (sub_sa, sub_lcp)
        });

        let subparts: Vec<_> = par_iter.collect();
        // It's not important for the parts to be in any order (?)
        //subparts.sort_by(|a, b| a.0.cmp(&b.0));
        subparts
    }

    // This is the standalone generator for a suffix array
    pub fn generate(
        self: &Self,
        range: Range<usize>,
    ) -> (Vec<usize>, Vec<usize>) {
        //let mut sa: Vec<usize> = (0..self.len).collect();
        let mut sa: Vec<usize> = range.collect();
        let len = sa.len();
        let mut sa_w = sa.clone();
        let mut lcp = vec![0; len];
        let mut lcp_w = vec![0; len];
        //self.merge_sort(&mut sa_w, &mut sa, self.len, &mut lcp, &mut lcp_w);
        self.merge_sort(&mut sa_w, &mut sa, len, &mut lcp, &mut lcp_w);
        (sa, lcp)
    }

    pub fn merge_sort(
        self: &Self,
        x: &mut [usize],
        y: &mut [usize],
        n: usize,
        lcp: &mut [usize],
        lcp_w: &mut [usize],
    ) {
        if n == 1 {
            lcp[0] = 0;
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

    fn merge(
        self: &Self,
        suffix_array: &mut [usize],
        mid: usize,
        lcp_w: &mut [usize],
        z: &mut [usize],
        lcp_z: &mut [usize],
    ) {
        let mut x = suffix_array[..mid].to_vec();
        let mut y = suffix_array[mid..].to_vec();
        let mut lcp_x = lcp_w[..mid].to_vec();
        let mut lcp_y = lcp_w[mid..].to_vec();
        //println!("merge x {x:?} y {y:?} lcp_x {lcp_x:?} lcp_y {lcp_y:?}");
        let mut len_x = x.len();
        let mut len_y = y.len();

        // Bookkeeping
        let mut m = 0;
        let mut i = 0; // Index into x
        let mut j = 0; // Index into y
        let mut k = 0; // Index into z

        while i < len_x && j < len_y {
            //println!("i {i} {} j {j} {}", x[i], y[j]);
            let l_x = lcp_x[i];
            //println!("l_x {l_x} m {m}");

            if l_x > m {
                z[k] = x[i];
                lcp_z[k] = l_x;
            } else if l_x < m {
                z[k] = y[j];
                lcp_z[k] = m;
                m = l_x;
            } else {
                // Length of shorter suffix
                let max_n = self.len - max(x[i], y[j]);

                // Prefix-context length for the suffixes
                let context = min(self.max_context, max_n);

                //println!(
                //    "left  = {:?}",
                //    String::from_utf8(self.text[(x[i] + m)..].to_vec())
                //);
                //println!(
                //    "right = {:?}",
                //    String::from_utf8(self.text[(y[j] + m)..].to_vec())
                //);
                //println!("context = {}", context - m);

                // LCP(X_i, Y_j)
                let n = m + lcp(
                    &self.text[(x[i] + m)..],
                    &self.text[(y[j] + m)..],
                    context - m,
                );
                //println!("m = {m}, n = {n}, max_n = {max_n}");
                //println!(
                //    "next char left '{}' vs right '{}'",
                //    self.text[x[i] + n],
                //    self.text[y[j] + n]
                //);

                // If the len of the LCP is the entire shorter
                // sequence, take that.
                if n == max_n {
                    z[k] = max(x[i], y[j])
                }
                // Else, look at the next char after the LCP
                // to determine order.
                else if self.text[x[i] + n] < self.text[y[j] + n] {
                    z[k] = x[i]
                }
                // Else take from the right
                else {
                    z[k] = y[j]
                }

                // If we took from the right...
                if z[k] == x[i] {
                    lcp_z[k] = l_x;
                } else {
                    lcp_z[k] = m
                }
                m = n;
            }

            if z[k] == x[i] {
                i += 1;
            } else {
                j += 1;
                mem::swap(&mut x, &mut y);
                mem::swap(&mut len_x, &mut len_y);
                mem::swap(&mut lcp_x, &mut lcp_y);
                mem::swap(&mut i, &mut j);
            }

            k += 1;
        }

        // Copy rest of the data from X to Z.
        while i < len_x {
            z[k] = x[i];
            lcp_z[k] = lcp_x[i];
            i += 1;
            k += 1;
        }

        // Copy rest of the data from Y to Z.
        if j < len_y {
            z[k] = y[j];
            lcp_z[k] = m;
            j += 1;
            k += 1;

            while j < len_y {
                z[k] = y[j];
                lcp_z[k] = lcp_y[j];
                j += 1;
                k += 1;
            }
        }
    }
}

// --------------------------------------------------
fn convert_slices_to_vecs(vec_of_slices: Vec<&[usize]>) -> Vec<Vec<usize>> {
    vec_of_slices
        .into_iter() // Convert Vec<&[usize]> into an iterator
        .map(|slice| slice.to_vec()) // Convert each slice into a Vec<usize>
        .collect() // Collect into a Vec<Vec<usize>>
}

// --------------------------------------------------
fn transpose(matrix: Vec<Vec<Option<&[usize]>>>) -> Vec<Vec<&[usize]>> {
    // Determine the number of columns (max length of rows)
    let num_cols = matrix.iter().map(|row| row.len()).max().unwrap_or(0);

    // Collect transposed rows, unpack the Options to remove None values
    (0..num_cols)
        .map(|col| {
            matrix
                .iter()
                .filter_map(|row| row.get(col))
                .cloned()
                .flat_map(|v| v)
                .filter(|v| !v.is_empty())
                .collect()
        })
        .collect()
}

// --------------------------------------------------
fn lcp(s1: &[u8], s2: &[u8], len: usize) -> usize {
    for i in 0..len {
        if s1[i] != s2[i] {
            return i;
        }
    }
    len
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
    //use std::io::Cursor;

    #[test]
    fn test_upper_bound() {
        //          012345
        let text = "TTTAGC".as_bytes().to_vec();
        let suf_arr = SuffixArray {
            text: text.clone(),
            len: text.len(),
            max_context: text.len(),
        };

        assert_eq!(suf_arr.upper_bound("A$", &[3, 5, 4]), None);
        assert_eq!(suf_arr.upper_bound("C$", &[3, 5, 4]), Some(2));
        assert_eq!(suf_arr.upper_bound("G$", &[3, 5, 4]), Some(1));
        //assert_eq!(suf_arr.upper_bound("T$", &[0, 1, 2]), 3);
        //assert_eq!(suf_arr.upper_bound("ACG$", &[0, 1, 2]), Some(0));
        //assert_eq!(suf_arr.upper_bound("G$", &[0, 1, 2]), Some(3));
        //assert_eq!(suf_arr.upper_bound("A$", &[1, 2]), Some(0));
        //assert_eq!(suf_arr.upper_bound("C$", &[0, 1, 2]), Some(1));
        //assert_eq!(suf_arr.upper_bound("G$", &[0, 1, 2]), Some(2));
        //assert_eq!(suf_arr.upper_bound("T$", &[0, 1, 2]), Some(3));
    }

    #[test]
    fn test_sort_subarrays() -> Result<()> {
        let text = "AACTGCGGAT$".as_bytes().to_vec();
        let suf_arr = SuffixArray {
            text: text.clone(),
            len: text.len(),
            max_context: text.len(),
        };

        // Ensure we get two subarrays
        let subs = suf_arr.sort_subarrays(2);
        assert_eq!(subs.len(), 2);

        for (sa, _lcp) in subs {
            let suffixes: Vec<String> = sa
                .iter()
                .flat_map(|&p| String::from_utf8(text[p..].to_vec()))
                .collect();

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

    //#[test]
    //fn test_lcp() -> Result<()> {
    //    assert_eq!(lcp(&[b'A'], &[b'C'], 1), 0);

    //    assert_eq!(lcp(&[b'A'], &[b'A'], 1), 1);

    //    assert_eq!(lcp(&[b'A', b'A'], &[b'A', b'A', b'C'], 2), 2);

    //    assert_eq!(lcp(&[b'A', b'C'], &[b'A', b'A', b'C'], 2), 1);

    //    assert_eq!(lcp(&[b'A', b'A'], &[b'A', b'A', b'C'], 1), 1);
    //    Ok(())
    //}
}
