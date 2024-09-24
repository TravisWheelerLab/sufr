use anyhow::{anyhow, bail, Result};
use log::info;
use needletail::parse_fastx_file;
use rand::seq::SliceRandom;
use rayon::prelude::*;
use std::{
    cmp::{max, min, Ordering},
    fmt::Debug,
    fmt::Display,
    fs::File,
    io::{BufWriter, Read, Write},
    mem,
    ops::{Add, Sub},
    slice,
    time::Instant,
};

const OUTFILE_VERSION: u8 = 1;

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
#[derive(Debug)]
pub struct SuffixArray<T>
where
    T: Int + FromUsize<T> + Sized + Send + Sync,
{
    pub version: u8,
    pub is_dna: bool,
    pub max_context: T,
    pub len: T,
    pub num_sequences: T,
    pub sequence_starts: Vec<T>,
    pub headers: Vec<String>,
    pub suffix_array: Vec<T>,
    pub lcp: Vec<T>,
    pub text: Vec<u8>,
}

impl<T> SuffixArray<T>
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
    ) -> SuffixArray<T> {
        let total_time = Instant::now();
        let now = Instant::now();
        let suffix_array: Vec<T> = if is_dna {
            // Text will already be "uppercase" (but bytes)
            text.iter()
                .enumerate()
                .filter_map(|(i, c)| {
                    b"ACGT".contains(c).then_some(T::from_usize(i))
                })
                .collect()
        } else {
            (0..text.len()).map(T::from_usize).collect()
        };
        let lcp = vec![T::default(); suffix_array.len()];
        info!(
            "Created unsorted suffix array of len {} in {:?}",
            suffix_array.len(),
            now.elapsed()
        );

        let mut sa = SuffixArray {
            version: OUTFILE_VERSION,
            is_dna,
            max_context: max_context.unwrap_or(len),
            len,
            text,
            num_sequences: T::from_usize(sequence_starts.len()),
            sequence_starts,
            headers,
            suffix_array,
            lcp,
        };
        sa.sort(num_partitions);
        info!(
            "Total time to create suffix array: {:?}",
            total_time.elapsed()
        );

        sa
    }

    // --------------------------------------------------
    pub fn check_order(&self) -> Vec<T> {
        let mut errors = vec![];
        for window in self.suffix_array.windows(2) {
            if let [prev, cur] = window {
                if !self.is_less(*prev, *cur) {
                    errors.push(*prev);
                }
            }
        }
        errors
    }

    // --------------------------------------------------
    pub fn check_lcp(&self) -> Vec<T> {
        let mut errors = vec![];
        for i in 1..self.suffix_array.len() {
            let len_lcp = self.find_lcp(
                self.suffix_array[i],
                self.suffix_array[i - 1],
                self.max_context,
            );
            if len_lcp != self.lcp[i] {
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
    pub fn sort(&mut self, num_partitions: usize) {
        // Select random pivots
        let now = Instant::now();
        let pivot_sa =
            self.select_pivots(&self.suffix_array, num_partitions - 1);
        info!(
            "Selected/sorted {} pivots in {:?}",
            pivot_sa.len(),
            now.elapsed()
        );

        // Partition the suffixes using the pivots
        let now = Instant::now();
        let partitions = self.partition(&self.suffix_array, pivot_sa);
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

        self.suffix_array = sa;
        self.lcp = lcp;
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
    fn partition(&self, suffixes: &Vec<T>, pivot_sa: Vec<T>) -> Vec<Vec<T>> {
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
    // Read serialized ".sufr" file
    pub fn read(filename: &str) -> Result<SuffixArray<T>> {
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
        let sequence_starts: Vec<T> = unsafe {
            std::slice::from_raw_parts(
                buffer.as_ptr() as *const _,
                num_sequences.to_usize(),
            )
            .to_vec()
        };

        // Suffix Array
        let mut buffer = vec![0u8; sa_len * mem::size_of::<T>()];
        file.read_exact(&mut buffer)?;
        let suffix_array: Vec<T> = unsafe {
            std::slice::from_raw_parts(buffer.as_ptr() as *const _, sa_len)
                .to_vec()
        };

        // LCP Array
        let mut buffer = vec![0u8; sa_len * mem::size_of::<T>()];
        file.read_exact(&mut buffer)?;
        let lcp: Vec<T> = unsafe {
            std::slice::from_raw_parts(buffer.as_ptr() as *const _, sa_len)
                .to_vec()
        };

        // Sequence -- add 2 for final "#" character delimiting text
        let mut buffer = vec![0u8; text_len * mem::size_of::<u8>()];
        file.read_exact(&mut buffer)?;
        let text: Vec<u8> = unsafe {
            std::slice::from_raw_parts(buffer.as_ptr(), text_len).to_vec()
        };

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

    // --------------------------------------------------
    // Serialize SuffixArray to a ".sufr" file
    pub fn write(&self, outfile: &str) -> Result<usize> {
        let mut out = BufWriter::new(
            File::create(outfile).map_err(|e| anyhow!("{outfile}: {e}"))?,
        );
        let mut bytes_out = 0;

        // Header: version/is_dna
        let is_dna: u8 = if self.is_dna { 1 } else { 0 };
        bytes_out += out.write(&[OUTFILE_VERSION, is_dna])?;

        // Text length
        bytes_out += out.write(&usize_to_bytes(self.len.to_usize()))?;

        // SA length
        bytes_out += out.write(&usize_to_bytes(self.suffix_array.len()))?;

        // Max context
        bytes_out +=
            out.write(&usize_to_bytes(self.max_context.to_usize()))?;

        // Number of sequences
        bytes_out +=
            out.write(&usize_to_bytes(self.sequence_starts.len()))?;

        // Sequence starts
        let slice_starts: &[u8] = unsafe {
            slice::from_raw_parts(
                self.sequence_starts.as_ptr() as *const _,
                self.sequence_starts.len() * mem::size_of::<T>(),
            )
        };
        bytes_out += out.write(slice_starts)?;

        // Suffix array
        let slice_sa: &[u8] = unsafe {
            slice::from_raw_parts(
                self.suffix_array.as_ptr() as *const _,
                self.suffix_array.len() * mem::size_of::<T>(),
            )
        };
        bytes_out += out.write(slice_sa)?;

        // LCP array
        let slice_lcp: &[u8] = unsafe {
            slice::from_raw_parts(
                self.lcp.as_ptr() as *const _,
                self.lcp.len() * mem::size_of::<T>(),
            )
        };
        bytes_out += out.write(slice_lcp)?;

        // Sequence
        bytes_out += out.write(&self.text)?;

        // Headers are variable in length so they are at the end
        bytes_out += out.write(&bincode::serialize(&self.headers)?)?;

        Ok(bytes_out)
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
    let mut seq = vec![];
    let mut headers = vec![];
    let mut start_positions = vec![];
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
        read_sequence_file, read_suffix_length, usize_to_bytes, SuffixArray,
    };
    use anyhow::Result;
    use tempfile::NamedTempFile;

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

    #[test]
    fn test_write_read_suffix_file_32() -> Result<()> {
        let seq_file = "tests/inputs/2.fa";
        let seq_data = read_sequence_file(seq_file)?;
        let max_context: Option<u32> = None;
        let is_dna = true;
        let start_positions: Vec<_> =
            seq_data.start_positions.iter().map(|&v| v as u32).collect();
        let len = seq_data.seq.len() as u32;
        let num_partitions = 2;
        let suffix_array: SuffixArray<u32> = SuffixArray::new(
            seq_data.seq,
            len,
            max_context,
            is_dna,
            start_positions,
            seq_data.headers,
            num_partitions,
        );

        let sorted_sa =
            [13, 4, 9, 0, 14, 5, 10, 1, 15, 6, 11, 2, 16, 7, 12, 3];
        let lcp = [0, 4, 4, 8, 0, 3, 3, 7, 0, 2, 2, 6, 0, 1, 1, 5];
        let outfile = NamedTempFile::new()?;
        let outpath = &outfile.path().to_str().unwrap();
        let res = suffix_array.write(outpath);
        assert!(res.is_ok());
        assert!(outfile.path().exists());

        let res: Result<SuffixArray<u32>> = SuffixArray::read(outpath);
        assert!(res.is_ok());

        let sa = res.unwrap();
        assert_eq!(sa.version, 1);
        assert_eq!(sa.is_dna, true);
        assert_eq!(sa.len, 18);
        assert_eq!(sa.num_sequences, 2);
        assert_eq!(sa.sequence_starts, [0, 9]);
        assert_eq!(sa.headers, ["ABC", "DEF"]);
        assert_eq!(sa.suffix_array, sorted_sa);
        assert_eq!(sa.lcp, lcp);
        assert_eq!(sa.text, b"ACGTACGT$ACGTACGT#");
        Ok(())
    }

    #[test]
    fn test_write_read_suffix_file_64() -> Result<()> {
        let seq_file = "tests/inputs/2.fa";
        let seq_data = read_sequence_file(seq_file)?;
        let max_context: Option<u64> = None;
        let is_dna = true;
        let start_positions: Vec<_> =
            seq_data.start_positions.iter().map(|&v| v as u64).collect();
        let len = seq_data.seq.len() as u64;
        let num_partitions = 2;
        let suffix_array: SuffixArray<u64> = SuffixArray::new(
            seq_data.seq,
            len,
            max_context,
            is_dna,
            start_positions,
            seq_data.headers,
            num_partitions,
        );

        let sorted_sa: Vec<u64> =
            vec![13, 4, 9, 0, 14, 5, 10, 1, 15, 6, 11, 2, 16, 7, 12, 3]
                .iter()
                .map(|&v| v as u64)
                .collect();
        let lcp: Vec<u64> =
            vec![0, 4, 4, 8, 0, 3, 3, 7, 0, 2, 2, 6, 0, 1, 1, 5]
                .iter()
                .map(|&v| v as u64)
                .collect();
        let outfile = NamedTempFile::new()?;
        let outpath = &outfile.path().to_str().unwrap();
        let res = suffix_array.write(outpath);
        assert!(res.is_ok());
        assert!(outfile.path().exists());

        let res: Result<SuffixArray<u64>> = SuffixArray::read(outpath);
        assert!(res.is_ok());

        let sa = res.unwrap();
        assert_eq!(sa.version, 1);
        assert_eq!(sa.is_dna, true);
        assert_eq!(sa.len, 18);
        assert_eq!(sa.num_sequences, 2);
        assert_eq!(sa.sequence_starts, [0, 9]);
        assert_eq!(sa.headers, ["ABC", "DEF"]);
        assert_eq!(sa.suffix_array, sorted_sa);
        assert_eq!(sa.lcp, lcp);
        assert_eq!(sa.text, b"ACGTACGT$ACGTACGT#");
        Ok(())
    }

    #[test]
    fn test_upper_bound() -> Result<()> {
        //          012345
        let text = "TTTAGC".as_bytes().to_vec();
        let len = text.len();
        let max_context: Option<u32> = None;
        let is_dna = false;
        let sequence_starts = vec![0];
        let headers = vec!["1".to_string()];
        let num_partitions = 2;
        let sa: SuffixArray<u32> = SuffixArray::new(
            text,
            len as u32,
            max_context,
            is_dna,
            sequence_starts,
            headers,
            num_partitions,
        );

        // The suffix "AGC$" is found before "GC$" and "C$
        assert_eq!(sa.upper_bound(3, &[5, 4]), None);

        // The suffix "TAGC$" is beyond all the values
        assert_eq!(sa.upper_bound(2, &[3, 5, 4]), Some(3));

        // The "C$" is the last value
        assert_eq!(sa.upper_bound(5, &[3, 5, 4]), Some(2));

        Ok(())
    }
}
