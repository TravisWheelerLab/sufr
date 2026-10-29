use crate::{
    file_access::write_at,
    lcp::{lcp, LcpCache},
    suffix_array::SuffixArray,
    types::{Int, SeedMask, SuffixSortType, OUTFILE_VERSION, SENTINEL_CHARACTER},
    util::{self, find_lcp_full_offset, slice_int_to_slice_u8, vec_to_slice_u8},
};
use anyhow::{anyhow, bail, Result};
use log::info;
use rayon::prelude::*;
use std::{
    borrow::Cow,
    cell::RefCell,
    cmp::{max, min, Ordering},
    fs::{self, File, OpenOptions},
    io::{Read, Write},
    mem,
    ops::Range,
    path::{Path, PathBuf},
    sync::Mutex,
    time::Instant,
};
use tempfile::NamedTempFile;
use thread_local::ThreadLocal;

/// Configuration and construction of [`SuffixArray`]s.
///
/// Suffix Arrays can be constructed from a sequence file with [`SuffixArrayBuilder::from_sequence_file`],
/// or from in-memory sequence data with [`SuffixArrayBuilder::from_data`].
///
/// The input, processing, and output options can be configured with the various
/// `with_*` methods on this structure.
///
/// After configuration, call [`SuffixArrayBuilder::build_to_file`] to construct the suffix array
/// and write it out to a file.
///
/// # Examples
///
/// ```
/// use anyhow::Result;
/// use libsufr::SuffixArrayBuilder;
///
/// fn main() -> Result<()> {
///     let outfile = "3.sufr";
///     let suffix_array = SuffixArrayBuilder::from_sequence_file("../data/inputs/3.fa", b'%')?
///         .with_dna_options(true, false, true)
///         .build_to_file(outfile, false)?;
///
///     let meta = suffix_array.metadata()?;
///     assert_eq!(meta.text_len, 113);
///     assert_eq!(meta.len_suffixes, 101);
///
///     std::fs::remove_file(&outfile)?;
///     Ok(())
/// }
/// ```
///
pub struct SuffixArrayBuilder {
    // INPUT
    text: Vec<u8>,
    sequences: Vec<(usize, String)>,

    // INPUT OPTIONS
    is_dna: bool,
    allow_ambiguity: bool,
    ignore_softmask: bool,

    // PROCESSING OPTIONS
    sort_type: SuffixSortType,
    num_partitions: usize,
    memory_scratch: bool,

    // OUTPUT OPTIONS
    force_64bit: bool,
    write_lcp: bool,

    // BUILD STATE
    n_ranges: Vec<Range<usize>>,
}

impl SuffixArrayBuilder {
    /// Creates a new `SuffixArrayBuilder` from the contents of a sequence file.
    /// `sequence_delimiter` must be within the alphabet or sort below all characters -
    /// `%` is suggested for biological sequence.
    // TODO: "must be" - should either be nonconfigurable or should support other values.
    pub fn from_sequence_file(
        path: impl AsRef<Path>,
        sequence_delimiter: u8,
    ) -> Result<Self> {
        let now = Instant::now();
        let data = util::read_sequence_file(path.as_ref(), sequence_delimiter)?;
        info!(
            "Read input of len {} in {:?}",
            data.seq.len(),
            now.elapsed()
        );
        let sequences = data
            .start_positions
            .into_iter()
            .zip(data.sequence_names)
            .collect();
        Ok(Self::from_data(data.seq, sequences))
    }

    /// Creates a new `SuffixArrayBuilder` from in-memory `text`. `sequences`
    /// describes the start positions and names for subsequences of `text`,
    /// which should have a `%` inserted between them.
    pub fn from_data(text: Vec<u8>, sequences: Vec<(usize, String)>) -> Self {
        // TODO: validate text and sequences?
        // - sequences should be sorted, first starting at 0, length < text.len(), no duplicates. If empty, default?
        Self {
            text,
            sequences,
            is_dna: false,
            allow_ambiguity: false,
            ignore_softmask: false,
            sort_type: SuffixSortType::MaxQueryLen(0),
            num_partitions: 1024,
            memory_scratch: false,
            force_64bit: false,
            write_lcp: false,
            n_ranges: vec![],
        }
    }

    /// Sets whether the sequence text `is_dna` (default: `false`),
    /// and if so whether to index `N`s (`allow_ambiguity`, default: `false`).
    /// Additionally sets whether softmasked characters should be treated as `N`
    /// (default: `false`, treat soft-masked as unmasked by converting back to uppercase).
    pub fn with_dna_options(
        mut self,
        is_dna: bool,
        allow_ambiguity: bool,
        ignore_softmask: bool,
    ) -> Self {
        self.is_dna = is_dna;
        self.allow_ambiguity = allow_ambiguity;
        self.ignore_softmask = ignore_softmask;
        self
    }

    /// Sets the `sort_type` of this builder.
    pub fn with_sort_type(mut self, sort_type: SuffixSortType) -> Self {
        self.sort_type = sort_type;
        self
    }
    /// Sets the `max_query_len` of this builder.
    pub fn with_max_query_len(self, max_query_len: usize) -> Self {
        self.with_sort_type(SuffixSortType::MaxQueryLen(max_query_len))
    }
    /// Sets the `seed_mask` of this builder.
    pub fn with_seed_mask(self, mask: SeedMask) -> Self {
        self.with_sort_type(SuffixSortType::Mask(mask))
    }
    /// Sets the number of partitions to be used during the construction (default: `1024`).
    pub fn with_num_partitions(mut self, num_partitions: usize) -> Self {
        self.num_partitions = num_partitions;
        self
    }
    /// Sets whether to keep the partitions in memory during construction
    /// instead of in temporary files (default: `false`).
    pub fn with_memory_scratch(mut self, memory_scratch: bool) -> Self {
        self.memory_scratch = memory_scratch;
        self
    }

    // Sets whether or not to force 64-bit output even for small sequences (default: `false`).
    // Files must be read via SufrFile::<u64> directly so this is intended for testing only.
    #[cfg(test)]
    pub(crate) fn with_force_64bit(mut self, force_64bit: bool) -> Self {
        self.force_64bit = force_64bit;
        self
    }
    /// Sets whether or not to write the LCP array after construction (default: `false`).
    pub fn with_write_lcp(mut self, write_lcp: bool) -> Self {
        self.write_lcp = write_lcp;
        self
    }

    /// Builds the suffix array, writing it to `path` and loading it into memory.
    pub fn build_to_file(
        mut self,
        path: impl AsRef<Path>,
        low_memory: bool,
    ) -> Result<SuffixArray> {
        let path = path.as_ref();
        let occupancy = self.normalize_text();
        if self.allow_ambiguity {
            self.n_ranges = find_n_ranges(&self.text);
        }

        let fits_u32 = (self.text.len() as u64) < u32::MAX as u64 && !self.force_64bit;
        match (fits_u32, self.memory_scratch) {
            (true, false) => {
                self.build::<u32, DiskScratchBuffer<u32>>(path, occupancy)?
            }
            (true, true) => {
                self.build::<u32, MemoryScratchBuffer<u32>>(path, occupancy)?
            }
            (false, false) => {
                self.build::<u64, DiskScratchBuffer<u64>>(path, occupancy)?
            }
            (false, true) => {
                self.build::<u64, MemoryScratchBuffer<u64>>(path, occupancy)?
            }
        }
        SuffixArray::read(&path.to_string_lossy(), low_memory)
    }

    /// Upper-case soft-masked characters (or replace them with `N` when
    /// `ignore_softmask` is set) and return which bytes occur in the text.
    fn normalize_text(&mut self) -> [u8; 256] {
        let mut occupancy = [0u8; 256];
        for b in self.text.iter_mut() {
            if b.is_ascii_lowercase() {
                if self.ignore_softmask {
                    *b = b'N'
                } else {
                    *b = b.to_ascii_uppercase()
                }
            }
            occupancy[*b as usize] |= 1;
        }
        occupancy
    }

    fn build<T: Int, B: ScratchBuffer<Item = T> + Send + Sync>(
        &self,
        path: &Path,
        occupancy: [u8; 256],
    ) -> Result<()> {
        let partitions = self.sort::<T, B>(occupancy)?;
        self.write::<T, B>(path, partitions)?;
        Ok(())
    }

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
    /// * `lcp_cache`: optional LcpCache to store/retrieve from, for optimization
    #[inline(always)]
    fn find_lcp<T: Int>(
        &self,
        start1: usize,
        start2: usize,
        len: T,
        skip: usize,
        lcp_cache: Option<&mut LcpCache>,
    ) -> T {
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
                    .filter(|&v| v < self.text.len());

                let b_vals = mask
                    .positions
                    .iter()
                    .skip(skip)
                    .map(|&offset| start2 + offset)
                    .filter(|&v| v < self.text.len());

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
                        let text_len = self.text.len();
                        let len = if max_query_len > &0 {
                            *max_query_len
                        } else {
                            len.to_usize()
                        };
                        let start1 = start1 + skip;
                        let start2 = start2 + skip;
                        let end1 = min(start1 + len, text_len);
                        let end2 = min(start2 + len, text_len);
                        let lcp = lcp(
                            &self.text[start1..end1],
                            &self.text[start2..end2],
                            lcp_cache.as_deref().map(|c| (c, start1, start2)),
                        );
                        // IMPORTANT: Only cache if lcp < len (uncapped, "true" LCP)
                        if skip + lcp >= LcpCache::MIN_CACHEABLE && lcp < len {
                            if let Some(c) = lcp_cache {
                                c.set(start1 - skip, start2 - skip, skip + lcp);
                            }
                        }
                        T::from_usize(skip + lcp)
                    }
                }
            }
        }
    }

    // --------------------------------------------------
    /// Write the suffixes into temporary files for sorting
    fn partition<T, B, A>(&self) -> Result<PartitionBuildResult<B>>
    where
        T: Int,
        B: ScratchBuffer<Item = T> + Send + Sync,
        A: PartitioningAlphabet + Send + Sync,
    {
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
    fn sort<T: Int, B: ScratchBuffer<Item = T> + Send + Sync>(
        &self,
        mut occupancy: [u8; 256],
    ) -> Result<Vec<Partition<B>>> {
        let num_partitions = self.num_partitions;
        // ... and if there are none outside of the nucleic alphabet then use it.
        for &c in NucleicAcidAlphabet::ALPHABET.iter() {
            occupancy[c as usize] = 0;
        }
        let partition_build = if occupancy.iter().map(|&b| b as usize).sum::<usize>()
            == 0
        {
            info!("Detected Nucleic Acid alphabet: using 5-character prefix for partitioning");
            self.partition::<T, B, NucleicAcidAlphabet>()?
        } else {
            // Otherwise check the same for the amino acid alphabet
            for &c in AminoAcidAlphabet::ALPHABET.iter() {
                occupancy[c as usize] = 0;
            }
            if occupancy.iter().map(|&b| b as usize).sum::<usize>() == 0 {
                info!("Detected Amino Acid alphabet: using 3-character prefix for partitioning");
                self.partition::<T, B, AminoAcidAlphabet>()?
            } else {
                info!("Detected non-Nucleic/non-Amino Acid: using 2-byte prefix for partitioning");
                self.partition::<T, B, BytesAlphabet>()?
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

        let lcp_cache: ThreadLocal<RefCell<LcpCache>> = ThreadLocal::new();
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
                            &lcp_cache,
                        );

                        // Write to disk
                        let mut sa_buffer = B::default();
                        sa_buffer.extend_from_slice(&part_sa)?;
                        let lcp_buffer = if self.write_lcp {
                            let mut lcp_buffer = B::default();
                            lcp_buffer.extend_from_slice(&lcp)?;
                            Some(lcp_buffer)
                        } else {
                            None
                        };

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
        Ok(partitions)
    }

    // --------------------------------------------------
    fn merge_sort<T: Int>(
        &self,
        x: &mut [T],
        y: &mut [T],
        n: usize,
        lcp: &mut [T],
        lcp_w: &mut [T],
        lcp_cache: &ThreadLocal<RefCell<LcpCache>>,
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
                self.merge_sort(yl, xl, mid, lcpw_l, lcp_l, lcp_cache);
                self.merge_sort(yr, xr, n - mid, lcpw_r, lcp_r, lcp_cache);
            } else {
                rayon::join(
                    || self.merge_sort(yl, xl, mid, lcpw_l, lcp_l, lcp_cache),
                    || self.merge_sort(yr, xr, n - mid, lcpw_r, lcp_r, lcp_cache),
                );
            }

            self.merge(
                x,
                mid,
                lcp_w,
                y,
                lcp,
                &mut lcp_cache.get_or_default().borrow_mut(),
            );
        }
    }

    // --------------------------------------------------
    fn merge<T: Int>(
        &self,
        suffix_array: &mut [T],
        mid: usize,
        lcp_w: &mut [T],
        target_sa: &mut [T],
        target_lcp: &mut [T],
        lcp_cache: &mut LcpCache,
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
                    let max_n = T::from_usize(self.text.len()) - shorter_suffix;

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
                            Some(lcp_cache),
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
    /// The header and text are written first, then each partition's
    /// suffix array (and LCP array, when requested) is written in
    /// parallel at its precomputed offset, then the header is patched
    /// with the array positions.
    /// Returns the number of bytes written to disk.
    fn write<T: Int, B: ScratchBuffer<Item = T> + Send + Sync>(
        &self,
        path: &Path,
        partitions: Vec<Partition<B>>,
    ) -> Result<usize> {
        let now = Instant::now();
        let num_suffixes = partitions.iter().map(|p| p.len).sum();
        let out = self.open_output::<T>(path, num_suffixes)?;
        let int_size = mem::size_of::<T>();

        // Offset (in elements) of each partition in the arrays, and the
        // last suffix of the preceding partition for fixing the LCP boundary
        let mut offsets = Vec::with_capacity(partitions.len());
        let mut total = 0;
        let mut prev_last_suffix = None;
        for partition in &partitions {
            offsets.push((total, prev_last_suffix));
            total += partition.len;
            prev_last_suffix = Some(partition.last_suffix);
        }

        // NB: Consuming the partitions' buffers frees
        // memory/disk space along the way
        partitions.into_par_iter().zip(offsets).try_for_each(
            |(partition, (offset, prev_last_suffix))| -> Result<()> {
                let sa = partition.sa_buffer.consume()?;
                write_all_at(
                    &out.file,
                    vec_to_slice_u8(&sa),
                    (out.sa_pos + offset * int_size) as u64,
                )?;
                drop(sa);

                if let Some(lcp_buffer) = partition.lcp_buffer {
                    let mut lcp = lcp_buffer.consume()?;
                    if let Some(prev_last_suffix) = prev_last_suffix {
                        // Fix LCP boundary
                        if let Some(val) = lcp.first_mut() {
                            *val = self.find_lcp(
                                prev_last_suffix,
                                partition.first_suffix,
                                T::from_usize(self.text.len()),
                                0, // start at beginning
                                None,
                            );
                        }
                    }
                    write_all_at(
                        &out.file,
                        vec_to_slice_u8(&lcp),
                        (out.lcp_pos + offset * int_size) as u64,
                    )?;
                }
                Ok(())
            },
        )?;

        let bytes_out = self.finish_output(&out)?;
        info!(
            "Wrote suffix array to '{}' in {:?}",
            path.display(),
            now.elapsed()
        );
        Ok(bytes_out)
    }

    // --------------------------------------------------
    /// Create the output file and write the header and text. The suffix
    /// (and optional LCP) arrays are written afterwards at the positions
    /// recorded in the returned `Output`.
    fn open_output<T: Int>(&self, path: &Path, num_suffixes: usize) -> Result<Output> {
        let now = Instant::now();
        let file =
            File::create(path).map_err(|e| anyhow!("{}: {e}", path.display()))?;

        // Various metadata
        let is_dna: u8 = if self.is_dna { 1 } else { 0 };
        let allow_ambiguity: u8 = if self.allow_ambiguity { 1 } else { 0 };
        let ignore_softmask: u8 = if self.ignore_softmask { 1 } else { 0 };
        let mut header: Vec<u8> =
            vec![OUTFILE_VERSION, is_dna, allow_ambiguity, ignore_softmask];

        // Text length
        header.extend(self.text.len().to_le_bytes());

        // Locations of text, suffix array, and LCP; filled in by finish_output
        let locs_pos = header.len() as u64;
        header.extend(0usize.to_le_bytes());
        header.extend(0usize.to_le_bytes());
        header.extend(0usize.to_le_bytes());

        // Number of suffixes
        header.extend(num_suffixes.to_le_bytes());

        // Max query length
        let max_query_len = if let SuffixSortType::MaxQueryLen(val) = &self.sort_type {
            *val
        } else {
            0
        };
        header.extend(max_query_len.to_le_bytes());

        // Number of sequences
        header.extend(self.sequences.len().to_le_bytes());

        // Sequence starts
        let sequence_starts: Vec<T> = self
            .sequences
            .iter()
            .map(|&(start, _)| T::from_usize(start))
            .collect();
        header.extend_from_slice(vec_to_slice_u8(&sequence_starts));

        // Seed mask
        match &self.sort_type {
            SuffixSortType::Mask(seed_mask) => {
                header.extend(seed_mask.bytes.len().to_le_bytes());
                header.extend_from_slice(&seed_mask.bytes);
            }
            _ => header.extend(0usize.to_le_bytes()),
        }

        write_all_at(&file, &header, 0)?;

        // Text
        let text_pos = header.len();
        write_all_at(&file, &self.text, text_pos as u64)?;

        let int_size = mem::size_of::<T>();
        let sa_pos = text_pos + self.text.len();
        let sa_end = sa_pos + num_suffixes * int_size;
        let (lcp_pos, end_pos) = if self.write_lcp {
            (sa_end, sa_end + num_suffixes * int_size)
        } else {
            (0, sa_end)
        };
        info!("Wrote header and text in {:?}", now.elapsed());

        Ok(Output {
            file,
            locs_pos,
            text_pos,
            sa_pos,
            lcp_pos,
            end_pos,
        })
    }

    // --------------------------------------------------
    /// Write the sequence names after the arrays and record the array
    /// positions in the header. Returns the total number of bytes in
    /// the file.
    fn finish_output(&self, out: &Output) -> Result<usize> {
        // Sequence names are variable in length so they are at the end
        let names: Vec<&str> = self
            .sequences
            .iter()
            .map(|(_, name)| name.as_str())
            .collect();
        let names = bincode::serialize(&names)?;
        write_all_at(&out.file, &names, out.end_pos as u64)?;

        // Go back to header and record the locations
        let mut locs = out.text_pos.to_le_bytes().to_vec();
        locs.extend(out.sa_pos.to_le_bytes());
        locs.extend(out.lcp_pos.to_le_bytes());
        write_all_at(&out.file, &locs, out.locs_pos)?;

        Ok(out.end_pos + names.len())
    }
}

// --------------------------------------------------
/// The output file after the header and text have been written.
/// The suffix and LCP arrays are written into it at known offsets,
/// possibly from several threads.
struct Output {
    file: File,

    /// Where the header records the text/SA/LCP positions
    locs_pos: u64,

    /// Byte position of the text
    text_pos: usize,

    /// Byte position of the suffix array
    sa_pos: usize,

    /// Byte position of the LCP array, or 0 when none is written
    lcp_pos: usize,

    /// Byte position after the arrays, where the sequence names go
    end_pos: usize,
}

/// Write all of `buf` at `offset`, in parallel 64 MB chunks.
fn write_all_at(file: &File, buf: &[u8], offset: u64) -> Result<()> {
    let chunk_len = 1 << 26;
    buf.par_chunks(chunk_len)
        .enumerate()
        .try_for_each(|(i, mut chunk)| -> Result<()> {
            let mut pos = offset + (i * chunk_len) as u64;
            while !chunk.is_empty() {
                let n = write_at(file, chunk, pos)?;
                if n == 0 {
                    bail!("Failed to write to output file");
                }
                chunk = &chunk[n..];
                pos += n as u64;
            }
            Ok(())
        })
}

/// Find the runs of at least 1000 `N`s in `text`.
fn find_n_ranges(text: &[u8]) -> Vec<Range<usize>> {
    let mut n_ranges: Vec<Range<usize>> = vec![];
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
    n_ranges
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

    /// The buffer containing the LCP array, when one is written.
    lcp_buffer: Option<B>,
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
trait ScratchBuffer: Default {
    /// The type of items in the buffer
    type Item: Int;

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
struct DiskScratchBuffer<T: Int> {
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
                let mut file = File::open(path)?;
                let size = file.metadata()?.len() as usize;

                let file_count = size / std::mem::size_of::<T>();
                let mut data = Vec::with_capacity(file_count + self.buf.len());
                data.resize(file_count, T::from_usize(0));
                file.read_exact(slice_int_to_slice_u8(&mut data))?;

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
struct MemoryScratchBuffer<T: Int> {
    data: Vec<T>,
}

impl<T: Int> Default for MemoryScratchBuffer<T> {
    fn default() -> Self {
        Self { data: vec![] }
    }
}

impl<T: Int> ScratchBuffer for MemoryScratchBuffer<T> {
    type Item = T;

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
    use super::SuffixArrayBuilder;
    use crate::types::SeedMask;
    use anyhow::Result;
    use pretty_assertions::assert_eq;
    use std::fs;
    use tempfile::NamedTempFile;

    #[test]
    fn test_find_lcp_no_seed_mask() -> Result<()> {
        //           012345
        let text = b"TTTAGC".to_vec();
        let sufr = SuffixArrayBuilder::from_data(text, vec![(0, "1".to_string())]);

        // 0: TTTAGC
        // 1:  TTAGC
        // 6: len of text
        assert_eq!(sufr.find_lcp::<u32>(0, 1, 6, 0, None), 2);

        // 0: TTTAGC
        // 2:   TAGC
        // 6: len of text
        assert_eq!(sufr.find_lcp::<u32>(0, 2, 6, 0, None), 1);

        // 0: TTTAGC
        // 1:  TTAGC
        // 1: max query len = 1
        assert_eq!(sufr.find_lcp::<u32>(0, 1, 1, 0, None), 1);

        // 0: TTTAGC
        // 3:    AGC
        // 6: len of text
        assert_eq!(sufr.find_lcp::<u32>(0, 3, 6, 0, None), 0);

        // TODO: Add a test with skip

        Ok(())
    }

    #[test]
    fn test_find_lcp_with_seed_mask() -> Result<()> {
        //           012345
        let text = b"TTTTTA".to_vec();
        let sufr = SuffixArrayBuilder::from_data(text, vec![(0, "1".to_string())])
            .with_seed_mask(SeedMask::new("1101")?);

        // 0: TTTTTA
        // 1:  TTTTA
        assert_eq!(sufr.find_lcp::<u32>(0, 1, 3, 0, None), 3);

        // 0: TTTTTA
        // 2:   TTTA
        assert_eq!(sufr.find_lcp::<u32>(0, 2, 3, 0, None), 2);

        // 0: TTTTTA
        // 5:      A
        assert_eq!(sufr.find_lcp::<u32>(0, 5, 3, 0, None), 0);

        Ok(())
    }

    #[test]
    fn test_memory_scratch_matches_disk() -> Result<()> {
        let build = |memory_scratch| -> Result<Vec<u8>> {
            let outfile = NamedTempFile::new()?;
            SuffixArrayBuilder::from_sequence_file("../data/inputs/2.fa", b'%')?
                .with_memory_scratch(memory_scratch)
                .with_write_lcp(true)
                .build_to_file(outfile.path(), false)?;
            Ok(fs::read(outfile.path())?)
        };
        assert_eq!(build(true)?, build(false)?);
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
