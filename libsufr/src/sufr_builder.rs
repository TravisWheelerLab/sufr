//! Create on-disk suffix/LCP arrays
//!
//! Sufr builds the suffix and LCP (longest common prefix) arrays on disk.
//! The suffixes are partitioned into temporary files, and these are sorted
//! in parallel.
//! Call the `write` method to serialized the data structures to
//! a file (preferably with the _.sufr_ extension) that can be read
//! by `sufr_file`.
//!

use crate::{
    types::{
        Int, SeedMask, SuffixSortType, SufrBuilderArgs, OUTFILE_VERSION,
        SENTINEL_CHARACTER,
    },
    util::{find_lcp_full_offset, slice_u8_to_slice_int, vec_to_slice_u8},
};
use anyhow::{anyhow, bail, Result};
use log::info;
use rayon::prelude::*;
use std::{
    borrow::Cow,
    cmp::{max, min, Ordering},
    fs::{self, File, OpenOptions},
    io::{BufWriter, Seek, SeekFrom, Write},
    mem,
    ops::Range,
    path::PathBuf,
    sync::Mutex,
    time::Instant,
};
use tempfile::NamedTempFile;

// --------------------------------------------------
/// A struct for partitioning, sorting, and writing suffixes to disk
#[derive(Debug)]
pub struct SufrBuilder<T: Int, B: ScratchBuffer = DiskScratchBuffer<T>> {
    /// The serialization version.
    pub version: u8,

    /// Whether or not the sequence is nucleotide.
    pub is_dna: bool,

    /// Whether or not the nucleotide sequence allows characters
    /// other than A, C, G, or T.
    pub allow_ambiguity: bool,

    /// Whether or not the nucleotide sequence ignores
    /// softmasked/lowercase bases.
    pub ignore_softmask: bool,

    /// The length of the given text.
    pub text_len: T,

    /// The number of suffixes that were indexed from the text, which
    /// could be less than `text_len` when ambiguity/softmasked values
    /// are ignored.
    pub num_suffixes: T,

    /// The number of sequences represented in the text.
    pub num_sequences: T,

    /// The positions in the text where each sequence starts.
    pub sequence_starts: Vec<T>,

    /// The names of the sequences in the text. Should be the same length
    /// as `sequence_starts`.
    pub sequence_names: Vec<String>,

    /// The text that was indexed.
    pub text: Vec<u8>,

    /// Whether the text is sorted fully, using a maximum query length,
    /// or a seed mask.
    pub sort_type: SuffixSortType,

    /// The number of partitions to use when building.
    partitions: Vec<Partition<B>>,

    /// The locations of long runs of Ns in nucleotide text.
    pub n_ranges: Vec<Range<usize>>,

    /// The name of the output file
    pub path: String,
}

// --------------------------------------------------
impl<T: Int, B: ScratchBuffer<Item = T> + Send + Sync> SufrBuilder<T, B> {
    /// Create a new suffix/LCP array.
    /// The results will live in temporary files on-disk.
    /// The integer values representing the positions of each suffix will
    /// be `u32` when the length of the text is less than 2^32 and `u64`,
    /// otherwise.
    ///
    /// ```
    /// use anyhow::Result;
    /// use std::{fs, path::Path};
    /// use libsufr::{
    ///     sufr_builder::SufrBuilder,
    ///     types::SufrBuilderArgs,
    ///     util::read_sequence_file,
    /// };
    ///
    /// fn main() -> Result<()> {
    ///     let path = Path::new("../data/inputs/1.fa");
    ///     let sequence_delimiter = b'%';
    ///     let seq_data = read_sequence_file(path, sequence_delimiter)?;
    ///     let text_len = seq_data.seq.len() as u64;
    ///     let outfile = "1.sufr";
    ///     let builder_args = SufrBuilderArgs {
    ///         text: seq_data.seq,
    ///         low_memory: false,
    ///         path: Some(outfile.to_string()),
    ///         max_query_len: None,
    ///         is_dna: true,
    ///         allow_ambiguity: false,
    ///         ignore_softmask: true,
    ///         sequence_starts: seq_data.start_positions.into_iter().collect(),
    ///         sequence_names: seq_data.sequence_names,
    ///         num_partitions: 1024,
    ///         seed_mask: None,
    ///     };
    ///
    ///     if text_len < u32::MAX as u64 {
    ///         let sufr_builder: SufrBuilder<u32> = SufrBuilder::new(builder_args)?;
    ///     } else {
    ///         let sufr_builder: SufrBuilder<u64> = SufrBuilder::new(builder_args)?;
    ///     }
    ///
    ///     fs::remove_file(&outfile)?;
    ///
    ///     Ok(())
    /// }
    /// ```
    pub fn new(args: SufrBuilderArgs) -> Result<SufrBuilder<T, B>> {
        let mut text = args.text;

        // Count the byte occurrences so that the best alphabet can be selected for partitioning...
        let mut occupancy = [0u8; 256];
        text.iter_mut().for_each(|b| {
            // Check for lowercase
            if (97..=122).contains(b) {
                if args.ignore_softmask {
                    *b = b'N'
                } else {
                    // only shift lowercase ASCII
                    *b &= 0b1011111
                }
            }
            occupancy[*b as usize] |= 1;
        });
        let text_len = T::from_usize(text.len());

        if args.seed_mask.is_some() && args.max_query_len.is_some() {
            bail!("Cannot use max_query_len and seed_mask together");
        }

        let sort_type = if let Some(mask) = args.seed_mask {
            let seed_mask = SeedMask::new(&mask)?;
            SuffixSortType::Mask(seed_mask)
        } else {
            SuffixSortType::MaxQueryLen(args.max_query_len.unwrap_or(0))
        };

        // Check for long runs of Ns when ambiguous bases are allowed.
        let mut n_ranges: Vec<Range<usize>> = vec![];
        if args.allow_ambiguity {
            let mut n_start: Option<usize> = None;
            let min_n = 1000;
            let now = Instant::now();
            for (i, &byte) in text.iter().enumerate() {
                if byte == b'N' {
                    if n_start.is_none() {
                        n_start = Some(i);
                    }
                } else {
                    if let Some(prev) = n_start {
                        if i - prev >= min_n {
                            n_ranges.push(prev..i);
                        }
                    }
                    n_start = None;
                }
            }
            info!("Scanned for runs of Ns in {:?}", now.elapsed());
        }

        let mut sa = SufrBuilder {
            version: OUTFILE_VERSION,
            is_dna: args.is_dna,
            allow_ambiguity: args.allow_ambiguity,
            ignore_softmask: args.ignore_softmask,
            sort_type,
            text_len,
            num_suffixes: T::default(),
            text,
            num_sequences: T::from_usize(args.sequence_starts.len()),
            sequence_starts: args
                .sequence_starts
                .into_iter()
                .map(T::from_usize)
                .collect::<Vec<_>>(),
            sequence_names: args.sequence_names,
            partitions: vec![],
            n_ranges,
            path: args.path.unwrap_or("out.sufr".to_string()),
        };
        sa.sort(args.num_partitions, occupancy)?;
        sa.write()?;
        Ok(sa)
    }

    // --------------------------------------------------
    // TODO: Remove? Only useful during debugging
    // Return the string at a given suffix position
    // Warning: Assumes pos is always found.
    //
    // Args:
    // * `pos`: the suffix position
    //pub(crate) fn string_at(&self, pos: usize) -> String {
    //    self.text
    //        .get(pos..)
    //        .map(|v| String::from_utf8(v.to_vec()).unwrap())
    //        .unwrap()
    //}

    // --------------------------------------------------
    /// If a suffix is in a long run of Ns, return the position of the final N
    ///
    /// Args:
    /// * `suffix`: suffix position
    fn find_n_run(&self, suffix: usize) -> Option<usize> {
        self.n_ranges
            .binary_search_by(|range| {
                if range.contains(&suffix) {
                    Ordering::Equal
                } else if range.start < suffix {
                    Ordering::Less
                } else {
                    Ordering::Greater
                }
            })
            .ok()
            .map(|i| self.n_ranges[i].end)
    }

    fn lcp(a: &[u8], b: &[u8]) -> usize {
        #[cfg(target_arch = "x86_64")]
        {
            if is_x86_feature_detected!("avx2") {
                // SAFETY: avx2 is available
                return unsafe { Self::lcp_avx2(a, b) };
            }
        }
        Self::lcp_scalar(a, b)
    }

    #[cfg(target_arch = "x86_64")]
    #[target_feature(enable = "avx2")]
    /// AVX2 implementation of LCP with an unrolled loop
    ///
    /// SAFETY: This function must only be called after validating that avx2 is actually available
    unsafe fn lcp_avx2(a: &[u8], b: &[u8]) -> usize {
        use std::arch::x86_64::*;
        let n = a.len().min(b.len());
        let mut i = 0;
        const UNROLL: usize = 4;
        while i + UNROLL * 32 <= n {
            for j in 0..UNROLL {
                let off = i + j * 32;
                // SAFETY: the bounds of a and b (256bits = 32 bytes; 4 iterations) are checked at the top of the while loop
                let va = unsafe { _mm256_loadu_si256(a.as_ptr().add(off) as *const _) };
                let vb = unsafe { _mm256_loadu_si256(b.as_ptr().add(off) as *const _) };
                let eq = _mm256_cmpeq_epi8(va, vb);
                let mask = _mm256_movemask_epi8(eq) as u32;
                if mask != 0xFFFF_FFFF {
                    return off + (!mask).trailing_zeros() as usize;
                }
            }
            i += UNROLL * 32;
        }
        i + Self::lcp_scalar(&a[i..n], &b[i..n])
    }

    fn lcp_scalar(a: &[u8], b: &[u8]) -> usize {
        std::iter::zip(a, b).take_while(|(a, b)| a == b).count()
    }

    // --------------------------------------------------
    /// Find the longest common prefix between two suffixes.
    ///
    /// Args:
    /// * `start1`: position of first suffix
    /// * `start2`: position of second suffix
    /// * `len`: the maximum length to check, e.g., the maximum query
    ///   length or the weight of the seed mask (number of 1/"care" positions)
    /// * `skip`: skip over the this many characters at the beginning.
    ///   Because of the incremental way the LCPs are calculated, we may know
    ///   that two suffixes share the `skip` number in common already.
    #[inline(always)]
    fn find_lcp(&self, start1: usize, start2: usize, len: T, skip: usize) -> T {
        // TODO: Could we use traits for SortType, parameterize the builder
        // on initialization and avoid conditionals here?
        match &self.sort_type {
            SuffixSortType::Mask(mask) => {
                // Use the seed diff vector to select only the
                // "care" positions up to the length of the text
                let a_vals = mask
                    .positions
                    .iter()
                    .skip(skip)
                    .map(|&offset| start1 + offset)
                    .filter(|&v| v < self.text_len.to_usize());

                let b_vals = mask
                    .positions
                    .iter()
                    .skip(skip)
                    .map(|&offset| start2 + offset)
                    .filter(|&v| v < self.text_len.to_usize());

                unsafe {
                    T::from_usize(
                        skip + a_vals
                            .zip(b_vals)
                            .take_while(|(a, b)| {
                                self.text.get_unchecked(*a)
                                    == self.text.get_unchecked(*b)
                            })
                            .count(),
                    )
                }
            }
            SuffixSortType::MaxQueryLen(max_query_len) => {
                match (&self.find_n_run(start1), &self.find_n_run(start2)) {
                    // If the two suffixes start in long stretches of Ns
                    // Then use the min of the end positions
                    (Some(end1), Some(end2)) => {
                        T::from_usize(min(end1 - start1, end2 - start2))
                    }
                    _ => {
                        let text_len = self.text_len.to_usize();
                        let len = if max_query_len > &0 {
                            *max_query_len
                        } else {
                            len.to_usize()
                        };
                        let start1 = start1 + skip;
                        let start2 = start2 + skip;
                        let end1 = min(start1 + len, text_len);
                        let end2 = min(start2 + len, text_len);
                        T::from_usize(
                            skip + Self::lcp(
                                &self.text[start1..end1],
                                &self.text[start2..end2],
                            ),
                        )
                    }
                }
            }
        }
    }

    // --------------------------------------------------
    /// Determine whether or not the first suffix position is lexicographically
    /// less than the second.
    /// This function is used to place suffixes into the highest partition
    /// for sorting.
    ///
    /// Args:
    /// * `start1`: the position of the first suffix
    /// * `start2`: the position of the second suffix
    #[cfg(test)]
    #[inline(always)]
    fn is_less(&self, start1: T, start2: T) -> bool {
        if start1 == start2 {
            false
        } else {
            let max_query_len = match &self.sort_type {
                SuffixSortType::Mask(seed_mask) => T::from_usize(seed_mask.weight),
                SuffixSortType::MaxQueryLen(max_query_len) => {
                    if max_query_len > &0 {
                        T::from_usize(*max_query_len)
                    } else {
                        self.text_len
                    }
                }
            };

            let len_lcp = find_lcp_full_offset(
                self.find_lcp(start1.to_usize(), start2.to_usize(), max_query_len, 0)
                    .to_usize(),
                &self.sort_type,
            );

            if len_lcp >= max_query_len.to_usize() {
                // The strings are equal(ish)
                false
            } else {
                // Look at the next character
                match (
                    self.text.get(start1.to_usize() + len_lcp),
                    self.text.get(start2.to_usize() + len_lcp),
                ) {
                    (Some(a), Some(b)) => a < b,
                    (None, Some(_)) => true,
                    _ => false,
                }
            }
        }
    }

    // --------------------------------------------------
    /// Find the highest partition to place a suffix for sorting.
    ///
    /// Args:
    /// * `suffix`: a suffix position
    /// * `pivots`: randomly selected suffix positions sorted lexicographically
    #[cfg(test)]
    #[inline(always)]
    fn upper_bound(&self, suffix: T, pivots: &[T]) -> usize {
        // Returns 0 when pivots is empty
        pivots.partition_point(|&p| self.is_less(p, suffix))
    }

    // --------------------------------------------------
    /// Write the suffixes into temporary files for sorting
    fn partition<A: PartitioningAlphabet + Send + Sync>(
        &mut self,
    ) -> Result<PartitionBuildResult<B>> {
        let alphabet = A::init();
        let mut buffers: Vec<_> = vec![];
        for _ in 0..A::NUM_PARTITIONS {
            let buffer = B::default();
            buffers.push(Mutex::new(buffer));
        }

        let now = Instant::now();
        const CHUNK_SIZE: usize = 1 << 19;
        (0..self.text.len().div_ceil(CHUNK_SIZE))
            .into_par_iter()
            .try_for_each(|chunk| -> Result<()> {
                let start = chunk * CHUNK_SIZE;
                let end = (start + CHUNK_SIZE).min(self.text.len());

                let mut suffixes = Vec::with_capacity(end - start);
                let mut keys = Vec::with_capacity(end - start);
                let mut counts = vec![0; A::NUM_PARTITIONS];

                for i in start..end {
                    let val = self.text.get(i).copied().unwrap_or(0);
                    if !(val == SENTINEL_CHARACTER
                        || !self.is_dna // Allow anything else if not DNA
                        || (b"ACGT".contains(&val) || self.allow_ambiguity))
                    {
                        continue;
                    }

                    let mut key = 0;
                    match &self.sort_type {
                        SuffixSortType::Mask(mask) => {
                            for (j, &off) in
                                mask.positions.iter().enumerate().take(A::COUNT)
                            {
                                let c = self.text.get(i + off).copied().unwrap_or(0);
                                key |= alphabet.lookup(c)
                                    << ((A::COUNT - j - 1) * A::BITS);
                            }
                        }
                        SuffixSortType::MaxQueryLen(_) => {
                            for (j, &c) in
                                self.text[i..].iter().enumerate().take(A::COUNT)
                            {
                                key |= alphabet.lookup(c)
                                    << ((A::COUNT - j - 1) * A::BITS);
                            }
                        }
                    }

                    keys.push(key);
                    counts[key as usize] += 1;
                    suffixes.push(T::from_usize(i));
                }

                let mut offsets = Vec::with_capacity(A::NUM_PARTITIONS);
                let mut total = 0;
                for &n in counts.iter() {
                    offsets.push(total);
                    total += n;
                }

                let mut ordered = vec![T::default(); keys.len()];
                let mut cursors = offsets.clone();
                for (&key, &suf) in keys.iter().zip(suffixes.iter()) {
                    let offset = &mut cursors[key as usize];
                    ordered[*offset] = suf;
                    *offset += 1;
                }

                for (partition_num, (&offset, &count)) in
                    offsets.iter().zip(counts.iter()).enumerate()
                {
                    if count > 0 {
                        let slice = &ordered[offset..offset + count];
                        match buffers[partition_num].lock() {
                            Ok(mut buf) => {
                                if buf.extend_from_slice(slice).is_err() {
                                    bail!("Unable to write data to disk")
                                }
                            }
                            Err(e) => bail!("{e}"),
                        }
                    }
                }
                Ok(())
            })?;

        // Flush out any remaining buffers
        let mut num_suffixes = 0;
        let buffers = buffers
            .into_iter()
            .map(|buffer| match buffer.into_inner() {
                Ok(mut buf) => {
                    buf.flush()?;
                    num_suffixes += buf.count();
                    Ok(buf)
                }
                Err(e) => panic!("Failed to lock: {e}"),
            })
            .collect::<Result<Vec<_>>>()?;

        info!(
            "Wrote {num_suffixes} unsorted suffixes to partition{} in {:?}",
            if A::NUM_PARTITIONS == 1 { "" } else { "s" },
            now.elapsed()
        );

        //Ok((builders, num_suffixes))
        Ok(PartitionBuildResult {
            buffers,
            num_suffixes,
        })
    }

    // --------------------------------------------------
    /// Sort the suffixes.
    ///
    /// Args:
    /// * `num_partitions`: the number of partitions to used
    /// * `occupancy`: occurrences of bytes, used to select partitioning strategy
    fn sort(&mut self, num_partitions: usize, mut occupancy: [u8; 256]) -> Result<()> {
        // ... and if there are none outside of the nucleic alphabet then use it.
        for &c in NucleicAcidAlphabet::ALPHABET.iter() {
            occupancy[c as usize] = 0;
        }
        let partition_build = if occupancy.iter().map(|&b| b as usize).sum::<usize>()
            == 0
        {
            info!("Detected Nucleic Acid alphabet: using 5-character prefix for partitioning");
            self.partition::<NucleicAcidAlphabet>()?
        } else {
            // Otherwise check the same for the amino acid alphabet
            for &c in AminoAcidAlphabet::ALPHABET.iter() {
                occupancy[c as usize] = 0;
            }
            if occupancy.iter().map(|&b| b as usize).sum::<usize>() == 0 {
                info!("Detected Amino Acid alphabet: using 3-character prefix for partitioning");
                self.partition::<AminoAcidAlphabet>()?
            } else {
                info!("Detected non-Nucleic/non-Amino Acid: using 2-byte prefix for partitioning");
                self.partition::<BytesAlphabet>()?
            }
        };

        // Be sure to round up to get all the suffixes
        let num_per_partition = (partition_build.num_suffixes as f64
            / num_partitions as f64)
            .ceil() as usize;
        let total_sort_time = Instant::now();
        let mut num_taken = 0;
        let mut partition_inputs = std::iter::repeat_with(|| Vec::new())
            .take(num_partitions)
            .collect::<Vec<_>>();

        // We (probably) have many more partitions than we need,
        // so here we accumulate the small partitions from the left
        // stopping when we reach a boundary like 1M/partition.
        // This evens out the workload to sort the partitions.
        let mut part_buffers = partition_build.buffers.into_iter();
        for (partition_num, partition_input) in partition_inputs.iter_mut().enumerate()
        {
            let boundary = num_per_partition * (partition_num + 1);
            for buffer in part_buffers.by_ref() {
                if buffer.count() > 0 {
                    num_taken += buffer.count();
                    partition_input.push(buffer);
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

        let mut partitions: Vec<Option<Partition<B>>> =
            (0..num_partitions).map(|_| None).collect();

        partitions
            .par_iter_mut()
            .zip(partition_inputs)
            .enumerate()
            .try_for_each(
                |(partition_num, (partition, partition_input))| -> Result<()> {
                    // Find the suffixes in this partition
                    let mut part_sa = vec![];
                    for buffer in partition_input {
                        part_sa.extend_from_slice(&buffer.consume()?);
                    }

                    let len = part_sa.len();
                    if len > 0 {
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

                        // Write to disk
                        let mut sa_buffer = B::default();
                        sa_buffer.extend_from_slice(&part_sa)?;
                        let mut lcp_buffer = B::default();
                        lcp_buffer.extend_from_slice(&lcp)?;

                        *partition = Some(Partition {
                            order: partition_num,
                            len,
                            first_suffix: part_sa.first().unwrap().to_usize(),
                            last_suffix: part_sa.last().unwrap().to_usize(),
                            sa_buffer,
                            lcp_buffer,
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
    fn merge_sort(
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
            const NESTED_PAR_GRAIN_SIZE: usize = 1 << 13;
            let mid = n / 2;
            let (xl, xr) = x.split_at_mut(mid);
            let (yl, yr) = y.split_at_mut(mid);
            let (lcpw_l, lcpw_r) = lcp_w.split_at_mut(mid);
            let (lcp_l, lcp_r) = lcp.split_at_mut(mid);

            if mid < NESTED_PAR_GRAIN_SIZE {
                self.merge_sort(yl, xl, mid, lcpw_l, lcp_l);
                self.merge_sort(yr, xr, n - mid, lcpw_r, lcp_r);
            } else {
                rayon::join(
                    || self.merge_sort(yl, xl, mid, lcpw_l, lcp_l),
                    || self.merge_sort(yr, xr, n - mid, lcpw_r, lcp_r),
                );
            }

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
                    let shorter_suffix = max(x[idx_x], y[idx_y]);
                    let max_n = self.text_len - shorter_suffix;

                    let context = match &self.sort_type {
                        SuffixSortType::Mask(seed_mask) => T::from_usize(
                            seed_mask
                                .positions
                                .iter()
                                .filter(|&i| *i < max_n.to_usize())
                                .count(),
                        ),
                        SuffixSortType::MaxQueryLen(max_query_len) => {
                            if max_query_len > &0 {
                                min(T::from_usize(*max_query_len), max_n)
                            } else {
                                max_n
                            }
                        }
                    };

                    // LCP(X_i, Y_j)
                    let (len_lcp, full_len_lcp) = if m < context {
                        let lcp = self.find_lcp(
                            x[idx_x].to_usize(),
                            y[idx_y].to_usize(),
                            context - m,
                            m.to_usize(), // skip
                        );
                        let full_lcp =
                            find_lcp_full_offset(lcp.to_usize(), &self.sort_type);
                        (lcp, T::from_usize(full_lcp))
                    } else {
                        (context, context)
                    };

                    // If full LCP equals context/MQL, take shorter suffix
                    if len_lcp >= context {
                        target_sa[idx_target] = shorter_suffix;
                    }
                    // Else, look at the next char after the LCP to determine order.
                    else {
                        let cmp = self.text[(x[idx_x] + full_len_lcp).to_usize()]
                            .cmp(&self.text[(y[idx_y] + full_len_lcp).to_usize()]);

                        match cmp {
                            Ordering::Equal => {
                                target_sa[idx_target] = shorter_suffix;
                            }
                            Ordering::Less => {
                                target_sa[idx_target] = x[idx_x];
                            }
                            Ordering::Greater => {
                                target_sa[idx_target] = y[idx_y];
                            }
                        }
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
    /// Serialize contents of the sorted partitions to a _.sufr_ file.
    /// Returns the number of bytes written to disk.
    ///
    /// Args:
    /// * `filename`: the name of the output file.
    fn write(&mut self) -> Result<usize> {
        let filename = &self.path;
        let mut file = BufWriter::new(
            File::create(filename).map_err(|e| anyhow!("{filename}: {e}"))?,
        );

        let mut bytes_out: usize = 0;

        // TODO (throughout this method): write() can write less than the whole buffer

        // Various metadata
        let is_dna: u8 = if self.is_dna { 1 } else { 0 };
        let allow_ambiguity: u8 = if self.allow_ambiguity { 1 } else { 0 };
        let ignore_softmask: u8 = if self.ignore_softmask { 1 } else { 0 };
        bytes_out +=
            file.write(&[OUTFILE_VERSION, is_dna, allow_ambiguity, ignore_softmask])?;

        // Text length
        bytes_out += file.write(&self.text_len.to_usize().to_le_bytes())?;

        // Locations of text, suffix array, and LCP
        // Will be corrected at the end
        let locs_pos = file.stream_position()?;
        bytes_out += file.write(&0usize.to_le_bytes())?;
        bytes_out += file.write(&0usize.to_le_bytes())?;
        bytes_out += file.write(&0usize.to_le_bytes())?;

        // Number of suffixes
        bytes_out += file.write(&self.num_suffixes.to_usize().to_le_bytes())?;

        // Max query length
        let max_query_len = if let SuffixSortType::MaxQueryLen(val) = &self.sort_type {
            *val
        } else {
            0
        };
        bytes_out += file.write(&max_query_len.to_le_bytes())?;

        // Number of sequences
        bytes_out += file.write(&self.sequence_starts.len().to_le_bytes())?;

        // Sequence starts
        bytes_out += file.write(vec_to_slice_u8(&self.sequence_starts))?;

        // Seed mask
        match &self.sort_type {
            SuffixSortType::Mask(seed_mask) => {
                bytes_out += file.write(&seed_mask.bytes.len().to_le_bytes())?;
                file.write_all(&seed_mask.bytes)?;
                bytes_out += seed_mask.bytes.len();
            }
            _ => bytes_out += file.write(&0usize.to_le_bytes())?,
        }

        // Text
        let text_pos = bytes_out;
        file.write_all(&self.text)?;
        bytes_out += self.text.len();

        // Stitch partitioned suffix files together
        // NB: Using mem::take() and consume() to free memory/disk space along the way
        let sa_pos = bytes_out;
        for partition in &mut self.partitions {
            let sa_buffer = mem::take(&mut partition.sa_buffer).consume()?;
            let sa_bytes = vec_to_slice_u8(&sa_buffer);
            bytes_out += sa_bytes.len();
            file.write_all(sa_bytes)?;
        }

        let lcp_pos = bytes_out;

        // Stitch partitioned LCP files together
        for i in 0..self.partitions.len() {
            let mut lcp = mem::take(&mut self.partitions[i].lcp_buffer).consume()?;

            if i != 0 {
                // Fix LCP boundary
                if let Some(val) = lcp.first_mut() {
                    *val = self.find_lcp(
                        self.partitions[i - 1].last_suffix,
                        self.partitions[i].first_suffix,
                        self.text_len,
                        0, // start at beginning
                    );
                }
            }

            let lcp_bytes = vec_to_slice_u8(&lcp);
            bytes_out += lcp_bytes.len();
            file.write_all(lcp_bytes)?;
        }

        // Sequence names are variable in length so they are at the end
        bytes_out += file.write(&bincode::serialize(&self.sequence_names)?)?;

        // Go back to header and record the locations
        file.seek(SeekFrom::Start(locs_pos))?;
        let _ = file.write(&text_pos.to_le_bytes())?;
        let _ = file.write(&sa_pos.to_le_bytes())?;
        let _ = file.write(&lcp_pos.to_le_bytes())?;

        Ok(bytes_out)
    }
}

/// Alphabet for partitioning by a prefix key,
/// i.e. placing suffixes into partitions based
/// on the first `COUNT` characters. The "key" is a packed integer;
/// `BITS`, `COUNT`, and `NUM_PARTITIONS` constants describe the size.
///
/// Implementations need only provide the value for `ALPHABET`; the rest
/// are computed.
///
/// Consumers should use only the `init()` and `lookup()` functions,
/// in particular to support `BytesAlphabet` which does not use a lookup
/// table.
trait PartitioningAlphabet {
    /// The alphabet used; must be in ascending sorted order
    const ALPHABET: &'static [u8];

    /// Number of bits required to represent a character in ALPHABET
    const BITS: usize = {
        let mut bits = 0;
        while (1usize << bits) < Self::ALPHABET.len() {
            bits += 1;
        }
        bits
    };

    /// Number of characters (packed into BITS) that fit into the accumulator (u16)
    const COUNT: usize = u16::BITS as usize / Self::BITS;

    /// Number of resulting partitions: `2 ^ BITS ^ COUNT`
    const NUM_PARTITIONS: usize = 2usize.pow(Self::BITS as u32).pow(Self::COUNT as u32);

    /// Returns the lookup table from a UTF-8/ASCII byte input to its rank.
    fn make_lookup_table() -> [u16; 256] {
        let mut lookup = [0u16; 256];
        for (i, &c) in Self::ALPHABET.iter().enumerate() {
            lookup[c as usize] = i as u16;
        }
        lookup
    }

    /// Initialize the alphabet, usually by creating a lookup table.
    fn init() -> Self;

    /// Get the value of `c` in the lower `BITS` bits of a `u16`.
    fn lookup(&self, c: u8) -> u16;
}

struct AminoAcidAlphabet([u16; 256]);
impl PartitioningAlphabet for AminoAcidAlphabet {
    const ALPHABET: &'static [u8] = b"$%*-ABCDEFGHIJKLMNOPQRSTUVWXYZ";

    fn init() -> Self {
        Self(Self::make_lookup_table())
    }

    fn lookup(&self, c: u8) -> u16 {
        self.0[c as usize]
    }
}

struct NucleicAcidAlphabet([u16; 256]);
impl PartitioningAlphabet for NucleicAcidAlphabet {
    const ALPHABET: &'static [u8] = b"$%ACGNT";

    fn init() -> Self {
        Self(Self::make_lookup_table())
    }

    fn lookup(&self, c: u8) -> u16 {
        self.0[c as usize]
    }
}

struct BytesAlphabet;
impl PartitioningAlphabet for BytesAlphabet {
    const ALPHABET: &'static [u8] = b"";
    const BITS: usize = 8;

    fn init() -> Self {
        Self
    }

    fn lookup(&self, c: u8) -> u16 {
        c as u16
    }
}

// --------------------------------------------------
/// Represents the partition values written to disk
#[derive(Debug)]
struct Partition<B: ScratchBuffer> {
    /// The sorted position of this parition.
    order: usize,

    /// The number of suffixes/LCP values contained in this partition.
    len: usize,

    /// The value of the first suffix. Used in stitching together the LCPs.
    first_suffix: usize,

    /// The value of the last suffix. Used in stitching together the LCPs.
    last_suffix: usize,

    /// The buffer containing the suffix array.
    sa_buffer: B,

    /// The buffer containing the LCP array.
    lcp_buffer: B,
}

// --------------------------------------------------
/// This struct provides access to the on-disk partitions.
#[derive(Debug)]
struct PartitionBuildResult<B: ScratchBuffer<Item: Int>> {
    /// A thread-safe vector of `PartitionBuilder` values
    buffers: Vec<B>,

    /// The total number of suffixes that were written to disk.
    num_suffixes: usize,
}

/// Trait representing scratch buffers that may be backed by disk or memory.
pub trait ScratchBuffer: Default {
    /// The type of items in the buffer
    type Item: Int;

    /// Append `val` to `self`.
    fn push(&mut self, val: Self::Item) -> Result<()>;

    /// Append a slice of `vals` to `self`.
    fn extend_from_slice(&mut self, vals: &[Self::Item]) -> Result<()>;

    /// Flush any unwritten data from `self` to the backing store.
    fn flush(&mut self) -> Result<()>;

    /// Return the number of items written to `self`.
    fn count(&self) -> usize;

    /// Return the full contents of `self`, either already in memory or loaded from disk.
    fn read(&self) -> Result<Cow<'_, [Self::Item]>>;

    /// Load the full contents of `self` into memory for the last time.
    fn consume(self) -> Result<Vec<Self::Item>>;
}

/// Implements a scratch buffer that writes its items to disk, then reads them back from disk later
#[derive(Debug)]
pub struct DiskScratchBuffer<T: Int> {
    buf: Vec<T>,
    count: usize,
    path: Option<PathBuf>,
}

impl<T: Int> DiskScratchBuffer<T> {
    const FLUSH_AT: usize = 4096;

    fn open(&mut self) -> Result<File> {
        if self.path.is_none() {
            let tmp = NamedTempFile::new()?;
            let (_, path) = tmp.keep()?;
            self.path = Some(path)
        }
        Ok(OpenOptions::new()
            .create(true)
            .append(true)
            .open(self.path.as_ref().unwrap())?)
    }
}

impl<T: Int> Default for DiskScratchBuffer<T> {
    fn default() -> Self {
        Self {
            buf: vec![],
            count: 0,
            path: None,
        }
    }
}

impl<T: Int> ScratchBuffer for DiskScratchBuffer<T> {
    type Item = T;
    fn push(&mut self, val: T) -> Result<()> {
        self.buf.push(val);
        if self.buf.len() >= Self::FLUSH_AT {
            self.flush()?;
        }

        Ok(())
    }

    fn extend_from_slice(&mut self, vals: &[T]) -> Result<()> {
        self.buf.extend_from_slice(vals);
        if self.buf.len() >= Self::FLUSH_AT {
            self.flush()?;
        }
        Ok(())
    }

    fn flush(&mut self) -> Result<()> {
        if !self.buf.is_empty() {
            let mut file = self.open()?;
            file.write_all(vec_to_slice_u8(&self.buf))?;
            self.count += self.buf.len();
            self.buf.clear();
            self.buf.shrink_to_fit();
        }
        Ok(())
    }

    fn count(&self) -> usize {
        self.count + self.buf.len()
    }

    fn read(&self) -> Result<Cow<'_, [T]>> {
        match &self.path {
            Some(path) => {
                let mut buffer = fs::read(path)?;
                let count = buffer.len() / std::mem::size_of::<T>();

                let mut data = Vec::with_capacity(self.buf.len() + count);
                data.extend_from_slice(slice_u8_to_slice_int(&mut buffer, count));
                data.extend_from_slice(&self.buf);
                Ok(Cow::Owned(data))
            }
            None => Ok(Cow::Borrowed(&self.buf)),
        }
    }

    fn consume(mut self) -> Result<Vec<Self::Item>> {
        let buf = match self.read()? {
            Cow::Owned(b) => b,
            Cow::Borrowed(b) => b.to_vec(),
        };
        if let Some(path) = self.path.take() {
            fs::remove_file(path)?;
        }
        Ok(buf)
    }
}

impl<T: Int> Drop for DiskScratchBuffer<T> {
    fn drop(&mut self) {
        if let Some(path) = self.path.take() {
            // best-effort since this is in Drop, and the file may have been consume()d already
            let _ = fs::remove_file(path);
        }
    }
}

/// Implements a scratch buffer that always keeps all items in memory
#[derive(Debug)]
pub struct MemoryScratchBuffer<T: Int> {
    data: Vec<T>,
}

impl<T: Int> Default for MemoryScratchBuffer<T> {
    fn default() -> Self {
        Self { data: vec![] }
    }
}

impl<T: Int> ScratchBuffer for MemoryScratchBuffer<T> {
    type Item = T;

    fn push(&mut self, val: T) -> Result<()> {
        self.data.push(val);
        Ok(())
    }

    fn extend_from_slice(&mut self, vals: &[T]) -> Result<()> {
        self.data.extend_from_slice(vals);
        Ok(())
    }

    fn flush(&mut self) -> Result<()> {
        Ok(())
    }

    fn count(&self) -> usize {
        self.data.len()
    }

    fn read(&self) -> Result<Cow<'_, [T]>> {
        Ok(Cow::Borrowed(self.data.as_ref()))
    }

    fn consume(self) -> Result<Vec<Self::Item>> {
        Ok(self.data)
    }
}

// --------------------------------------------------
#[cfg(test)]
mod test {
    use super::{SufrBuilder, SufrBuilderArgs};
    use anyhow::Result;
    use pretty_assertions::assert_eq;
    use std::fs;
    use tempfile::NamedTempFile;

    #[test]
    fn test_is_less() -> Result<()> {
        //           012345
        let text = b"TTTAGC".to_vec();
        let outfile = NamedTempFile::new()?;
        let args = SufrBuilderArgs {
            text,
            low_memory: true,
            path: Some(outfile.path().to_string_lossy().to_string()),
            max_query_len: None,
            is_dna: false,
            allow_ambiguity: false,
            ignore_softmask: false,
            sequence_starts: vec![0],
            sequence_names: vec!["1".to_string()],
            num_partitions: 2,
            seed_mask: None,
        };
        let sufr = SufrBuilder::<u32>::new(args)?;

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

        fs::remove_file(outfile)?;

        Ok(())
    }

    #[test]
    fn test_is_less_max_query_len() -> Result<()> {
        //           012345
        let text = b"TTTAGC".to_vec();
        let outfile = NamedTempFile::new()?;
        let args = SufrBuilderArgs {
            text,
            low_memory: true,
            path: Some(outfile.path().to_string_lossy().to_string()),
            max_query_len: Some(2),
            is_dna: false,
            allow_ambiguity: false,
            ignore_softmask: false,
            sequence_starts: vec![0],
            sequence_names: vec!["1".to_string()],
            num_partitions: 2,
            seed_mask: None,
        };
        let sufr = SufrBuilder::<u32>::new(args)?;

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

        fs::remove_file(outfile)?;

        Ok(())
    }

    #[test]
    fn test_is_less_seed_mask() -> Result<()> {
        //           012345
        let text = b"TTTTAT".to_vec();
        let outfile = NamedTempFile::new()?;
        let args = SufrBuilderArgs {
            text,
            low_memory: true,
            path: Some(outfile.path().to_string_lossy().to_string()),
            max_query_len: None,
            is_dna: false,
            allow_ambiguity: false,
            ignore_softmask: false,
            sequence_starts: vec![0],
            sequence_names: vec!["1".to_string()],
            num_partitions: 2,
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
        assert!(!sufr.is_less(1, 0));

        // 0: TTTTAT
        // 3: TAT
        // "T-T" vs "T-T"
        assert!(!sufr.is_less(0, 3));

        // 3: TAT
        // 0: TTTTAT
        // "T-T" vs "T-T"
        assert!(!sufr.is_less(3, 0));

        fs::remove_file(outfile)?;

        Ok(())
    }

    #[test]
    fn test_find_lcp_no_seed_mask() -> Result<()> {
        //           012345
        let text = b"TTTAGC".to_vec();
        let outfile = NamedTempFile::new()?;
        let args = SufrBuilderArgs {
            text,
            low_memory: true,
            path: Some(outfile.path().to_string_lossy().to_string()),
            max_query_len: None,
            is_dna: false,
            allow_ambiguity: false,
            ignore_softmask: false,
            sequence_starts: vec![0],
            sequence_names: vec!["1".to_string()],
            num_partitions: 2,
            seed_mask: None,
        };
        let sufr: SufrBuilder<u32> = SufrBuilder::new(args)?;

        // 0: TTTAGC
        // 1:  TTAGC
        // 6: len of text
        assert_eq!(sufr.find_lcp(0, 1, 6, 0), 2);

        // 0: TTTAGC
        // 2:   TAGC
        // 6: len of text
        assert_eq!(sufr.find_lcp(0, 2, 6, 0), 1);

        // 0: TTTAGC
        // 1:  TTAGC
        // 1: max query len = 1
        assert_eq!(sufr.find_lcp(0, 1, 1, 0), 1);

        // 0: TTTAGC
        // 3:    AGC
        // 6: len of text
        assert_eq!(sufr.find_lcp(0, 3, 6, 0), 0);

        // TODO: Add a test with skip

        fs::remove_file(outfile)?;

        Ok(())
    }

    #[test]
    fn test_find_lcp_with_seed_mask() -> Result<()> {
        //           012345
        let text = b"TTTTTA".to_vec();
        let outfile = NamedTempFile::new()?;
        let args = SufrBuilderArgs {
            text,
            low_memory: true,
            path: Some(outfile.path().to_string_lossy().to_string()),
            max_query_len: None,
            is_dna: false,
            allow_ambiguity: false,
            ignore_softmask: false,
            sequence_starts: vec![0],
            sequence_names: vec!["1".to_string()],
            num_partitions: 2,
            seed_mask: Some("1101".to_string()),
        };
        let sufr: SufrBuilder<u32> = SufrBuilder::new(args)?;

        // 0: TTTTTA
        // 1:  TTTTA
        assert_eq!(sufr.find_lcp(0, 1, 3, 0), 3);

        // 0: TTTTTA
        // 2:   TTTA
        assert_eq!(sufr.find_lcp(0, 2, 3, 0), 2);

        // 0: TTTTTA
        // 5:      A
        assert_eq!(sufr.find_lcp(0, 5, 3, 0), 0);

        fs::remove_file(outfile)?;

        Ok(())
    }

    #[test]
    fn test_upper_bound_1() -> Result<()> {
        //          012345
        let text = b"TTTAGC".to_vec();
        let outfile = NamedTempFile::new()?;
        let args = SufrBuilderArgs {
            text,
            low_memory: true,
            path: Some(outfile.path().to_string_lossy().to_string()),
            max_query_len: None,
            is_dna: false,
            allow_ambiguity: false,
            ignore_softmask: false,
            sequence_starts: vec![0],
            sequence_names: vec!["1".to_string()],
            num_partitions: 2,
            seed_mask: None,
        };
        let sufr: SufrBuilder<u32> = SufrBuilder::new(args)?;

        // The suffix "AGC$" is found before "GC$" and "C$
        assert_eq!(sufr.upper_bound(3, &[5, 4]), 0);

        // The suffix "TAGC$" is beyond all the values
        assert_eq!(sufr.upper_bound(2, &[3, 4, 5]), 3);

        // The "C$" is the last value
        assert_eq!(sufr.upper_bound(5, &[3, 4, 5]), 1);

        fs::remove_file(outfile)?;

        Ok(())
    }

    #[test]
    fn test_upper_bound_2() -> Result<()> {
        //           0123456789
        let text = b"ACGTNNACGT".to_vec();
        let outfile = NamedTempFile::new()?;
        let args = SufrBuilderArgs {
            text,
            low_memory: true,
            path: Some(outfile.path().to_string_lossy().to_string()),
            max_query_len: None,
            is_dna: false,
            allow_ambiguity: false,
            ignore_softmask: false,
            sequence_starts: vec![0],
            sequence_names: vec!["1".to_string()],
            num_partitions: 2,
            seed_mask: None,
        };

        let sufr: SufrBuilder<u64> = SufrBuilder::new(args)?;

        // ACGTNNACGT$ == ACGTNNACGT$
        assert_eq!(sufr.upper_bound(0, &[0]), 0);

        // ACGTNNACGT$ (0) > ACGT$ (6)
        assert_eq!(sufr.upper_bound(0, &[6]), 1);

        // ACGT$ < ACGTNNACGT$
        assert_eq!(sufr.upper_bound(6, &[0]), 0);

        // ACGT$ == ACGT$
        assert_eq!(sufr.upper_bound(6, &[6]), 0);

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

        fs::remove_file(outfile)?;

        Ok(())
    }

    #[test]
    fn test_upper_bound_seed_mask() -> Result<()> {
        //           0123456789
        let text = b"ACGTNNACGT".to_vec();
        let outfile = NamedTempFile::new()?;
        let args = SufrBuilderArgs {
            text,
            low_memory: true,
            path: Some(outfile.path().to_string_lossy().to_string()),
            max_query_len: None,
            is_dna: false,
            allow_ambiguity: false,
            ignore_softmask: false,
            sequence_starts: vec![0],
            sequence_names: vec!["1".to_string()],
            num_partitions: 2,
            seed_mask: Some("101".to_string()),
        };
        let sufr: SufrBuilder<u32> = SufrBuilder::new(args)?;

        // ACGTNNACGT$ == ACGTNNACGT$ (A-G)
        assert_eq!(sufr.upper_bound(0, &[0]), 0);

        // ACGTNNACGT$ == ACGT$ (A-G)
        assert_eq!(sufr.upper_bound(0, &[6]), 0);

        // ACGT$ == ACGTNNACGT$ (A-G)
        assert_eq!(sufr.upper_bound(6, &[0]), 0);

        // ACGT$ == ACGT$ (A-G)
        assert_eq!(sufr.upper_bound(6, &[6]), 0);

        // Pivots = [CGT$, GT$]
        // ACGTNNACGT$ < CGT$
        assert_eq!(sufr.upper_bound(0, &[7, 8]), 0);

        // Pivots = [CGT$, GT$]
        // CGTNNACGT$ == CGT$ (C-T)
        assert_eq!(sufr.upper_bound(1, &[7, 8]), 0);

        // Pivots = [CGT$, GT$]
        // GT$ == GT$
        assert_eq!(sufr.upper_bound(8, &[7, 8]), 1);

        // Pivots = [CGT$, GT$]
        // T$ > GT$  => p2
        assert_eq!(sufr.upper_bound(9, &[7, 8]), 2);

        // T$ == TNNACGT$ (only compare T)
        assert_eq!(sufr.upper_bound(9, &[3]), 0);

        fs::remove_file(outfile)?;

        Ok(())
    }

    #[test]
    fn test_alphabet_properties() {
        use super::{AminoAcidAlphabet, NucleicAcidAlphabet, PartitioningAlphabet};

        assert_eq!(NucleicAcidAlphabet::BITS, 3);
        assert_eq!(NucleicAcidAlphabet::COUNT, 5);

        assert_eq!(AminoAcidAlphabet::BITS, 5);
        assert_eq!(AminoAcidAlphabet::COUNT, 3);
    }
}
