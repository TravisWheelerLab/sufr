//! # Creating Suffix Arrays
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
        FromUsize, Int, SeedMask, SuffixSortType, SufrBuilderArgs, OUTFILE_VERSION,
        SENTINEL_CHARACTER,
    },
    util::{find_lcp_full_offset, slice_u8_to_vec, usize_to_bytes, vec_to_slice_u8},
};
use anyhow::{anyhow, bail, Result};
use log::info;
use rand::{rngs::StdRng, Rng, RngCore, SeedableRng};
use rayon::prelude::*;
use std::{
    cmp::{max, min, Ordering},
    collections::HashSet,
    fs::{self, File, OpenOptions},
    io::{BufWriter, Seek, SeekFrom, Write},
    mem,
    ops::Range,
    path::PathBuf,
    sync::{Arc, Mutex},
    time::Instant,
};
use tempfile::NamedTempFile;

// --------------------------------------------------
/// A struct for partitioning, sorting, and writing suffixes to disk
#[derive(Debug)]
pub struct SufrBuilder<T>
where
    T: Int + FromUsize<T> + Sized + Send + Sync + serde::ser::Serialize,
{
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
    partitions: Vec<Partition>,

    /// The locations of long runs of Ns in nucleotide text.
    pub n_ranges: Vec<Range<usize>>,

    /// The name of the output file
    pub path: String,
}

// --------------------------------------------------
impl<T> SufrBuilder<T>
where
    T: Int + FromUsize<T> + Sized + Send + Sync,
{
    /// Create a new suffix/LCP array.
    /// The results will live in temporary files on-disk.
    /// The integer values representing the positions of each suffix will
    /// be `u32` when the length of the text is less than 2^32 and `u64`,
    /// otherwise.
    ///
    /// ```
    /// use anyhow::Result;
    /// use std::path::Path;
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
    ///     let builder_args = SufrBuilderArgs {
    ///         text: seq_data.seq,
    ///         low_memory: false,
    ///         path: Some("1.sufr".to_string()),
    ///         max_query_len: None,
    ///         is_dna: true,
    ///         allow_ambiguity: false,
    ///         ignore_softmask: true,
    ///         sequence_starts: seq_data.start_positions.into_iter().collect(),
    ///         sequence_names: seq_data.sequence_names,
    ///         num_partitions: 1024,
    ///         seed_mask: None,
    ///         random_seed: 42,
    ///     };
    ///
    ///     if text_len < u32::MAX as u64 {
    ///         let sufr_builder: SufrBuilder<u32> = SufrBuilder::new(builder_args)?;
    ///     } else {
    ///         let sufr_builder: SufrBuilder<u64> = SufrBuilder::new(builder_args)?;
    ///     }
    ///
    ///     Ok(())
    /// }
    /// ```
    pub fn new(args: SufrBuilderArgs) -> Result<SufrBuilder<T>> {
        let text: Vec<_> = args
            .text
            .iter()
            .map(|b| {
                // Check for lowercase
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
        sa.sort(args.num_partitions, args.random_seed)?;
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
                        unsafe {
                            T::from_usize(
                                skip + (start1..end1)
                                    .zip(start2..end2)
                                    .take_while(|(a, b)| {
                                        self.text.get_unchecked(*a)
                                            == self.text.get_unchecked(*b)
                                    })
                                    .count(),
                            )
                        }
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
    #[inline(always)]
    fn upper_bound(&self, suffix: T, pivots: &[T]) -> usize {
        // Returns 0 when pivots is empty
        pivots.partition_point(|&p| self.is_less(p, suffix))
    }

    // --------------------------------------------------
    /// Write the suffixes into temporary files for sorting
    ///
    /// Args:
    /// * `num_partitions`: how many partitions to create
    /// * `random_seed`: a value for initializing the pseudo-random number
    ///   generator for reproducibility when selecting the suffixes used
    ///   for partitioning
    fn partition(
        &mut self,
        num_partitions: usize,
        random_seed: u64,
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
        let pivot_sa = self.select_pivots(self.text.len(), num_partitions, random_seed);
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
                if val == SENTINEL_CHARACTER
                    || !self.is_dna // Allow anything else if not DNA
                    || (b"ACGT".contains(&val) || self.allow_ambiguity)
                {
                    let suffix = T::from_usize(i);
                    let partition_num = self.upper_bound(suffix, &pivot_sa);
                    match builders[partition_num].lock() {
                        Ok(mut partition) => {
                            if partition.add(suffix).is_err() {
                                bail!("Unable to write data to disk")
                            }
                        }
                        Err(e) => bail!("{e}"),
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
    /// Sort the suffixes.
    ///
    /// Args:
    /// * `num_partitions`: the number of partitions to used
    /// * `random_seed`: a value for initializing the RNG
    fn sort(&mut self, num_partitions: usize, random_seed: u64) -> Result<()> {
        let mut partition_build = self.partition(num_partitions, random_seed)?;

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

        let mut partitions: Vec<Option<Partition>> =
            (0..num_partitions).map(|_| None).collect();

        partitions.par_iter_mut().enumerate().try_for_each(
            |(partition_num, partition)| -> Result<()> {
                // Find the suffixes in this partition
                let mut part_sa = vec![];
                for (path, len) in &partition_inputs[partition_num] {
                    let buffer = fs::read(path)?;
                    let mut part: Vec<T> = slice_u8_to_vec(&buffer, *len);
                    part_sa.append(&mut part);
                    fs::remove_file(path)?;
                }

                let len = part_sa.len();
                if len > 0 {
                    let mut sa_w = part_sa.clone();
                    let mut lcp = vec![T::default(); len];
                    let mut lcp_w = vec![T::default(); len];
                    self.merge_sort(&mut sa_w, &mut part_sa, len, &mut lcp, &mut lcp_w);

                    // Write to disk
                    let mut sa_file = NamedTempFile::new()?;
                    let _ = sa_file.write(vec_to_slice_u8(&part_sa))?;
                    let mut lcp_file = NamedTempFile::new()?;
                    let _ = lcp_file.write(vec_to_slice_u8(&lcp))?;
                    let (_, sa_path) = sa_file.keep()?;
                    let (_, lcp_path) = lcp_file.keep()?;

                    *partition = Some(Partition {
                        order: partition_num,
                        len,
                        first_suffix: part_sa.first().unwrap().to_usize(),
                        last_suffix: part_sa.last().unwrap().to_usize(),
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
    #[inline(always)]
    fn select_pivots(
        &self,
        text_len: usize,
        num_partitions: usize,
        random_seed: u64,
    ) -> Vec<T> {
        if num_partitions > 1 {
            // Use a HashMap because selecting pivots one-at-a-time
            // can result in duplicates.
            let num_pivots = num_partitions - 1;
            //let rng = &mut rand::thread_rng();
            let mut rng: Box<dyn RngCore> = if random_seed > 0 {
                Box::new(StdRng::seed_from_u64(random_seed))
            } else {
                Box::new(rand::rng())
            };
            let mut pivot_sa = HashSet::<T>::new();
            loop {
                let pos = rng.random_range(0..text_len);
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
    /// Serialize contents of the sorted partitions to a _.sufr_ file.
    /// Returns the number of bytes written to disk.
    ///
    /// Args:
    /// * `filename`: the name of the output file.
    fn write(&self) -> Result<()> {
        let filename = &self.path;
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
        let max_query_len = if let SuffixSortType::MaxQueryLen(val) = &self.sort_type {
            *val
        } else {
            0
        };
        bytes_out += file.write(&usize_to_bytes(max_query_len))?;

        // Number of sequences
        bytes_out += file.write(&usize_to_bytes(self.sequence_starts.len()))?;

        // Sequence starts
        bytes_out += file.write(vec_to_slice_u8(&self.sequence_starts))?;

        // Seed mask
        match &self.sort_type {
            SuffixSortType::Mask(seed_mask) => {
                bytes_out += file.write(&usize_to_bytes(seed_mask.bytes.len()))?;
                file.write_all(&seed_mask.bytes)?;
                bytes_out += seed_mask.bytes.len();
            }
            _ => bytes_out += file.write(&usize_to_bytes(0))?,
        }

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
            fs::remove_file(&partition.sa_path)?;
        }

        let lcp_pos = bytes_out;

        // Stitch partitioned LCP files together
        for (i, partition) in self.partitions.iter().enumerate() {
            let buffer = fs::read(&partition.lcp_path)?;
            bytes_out += &buffer.len();

            if i == 0 {
                file.write_all(&buffer)?;
            } else {
                // Fix LCP boundary
                let mut lcp: Vec<T> = slice_u8_to_vec(&buffer, partition.len);
                if let Some(val) = lcp.first_mut() {
                    *val = self.find_lcp(
                        self.partitions[i - 1].last_suffix,
                        partition.first_suffix,
                        self.text_len,
                        0, // start at beginning
                    );
                }
                file.write_all(vec_to_slice_u8(&lcp))?;
            }
            fs::remove_file(&partition.lcp_path)?;
        }

        // Sequence names are variable in length so they are at the end
        _ = file.write(&bincode::serialize(&self.sequence_names)?)?;

        // Go back to header and record the locations
        file.seek(SeekFrom::Start(locs_pos))?;
        let _ = file.write(&usize_to_bytes(text_pos))?;
        let _ = file.write(&usize_to_bytes(sa_pos))?;
        let _ = file.write(&usize_to_bytes(lcp_pos))?;

        Ok(())
    }
}

// --------------------------------------------------
/// Represents the partition values written to disk
#[derive(Debug)]
struct Partition {
    /// The sorted position of this parition.
    order: usize,

    /// The number of suffixes/LCP values contained in this partition.
    len: usize,

    /// The value of the first suffix. Used in stitching together the LCPs.
    first_suffix: usize,

    /// The value of the last suffix. Used in stitching together the LCPs.
    last_suffix: usize,

    /// The path to the file containing the suffix array.
    sa_path: PathBuf,

    /// The path to the file containing the LCP array.
    lcp_path: PathBuf,
}

// --------------------------------------------------
/// This struct provides access to the on-disk partitions.
#[derive(Debug)]
struct PartitionBuildResult<T>
where
    T: Int + FromUsize<T> + Sized + Send + Sync + serde::ser::Serialize,
{
    /// A thread-safe vector of `PartitionBuilder` values
    builders: Vec<Arc<Mutex<PartitionBuilder<T>>>>,

    /// The total number of suffixes that were written to disk.
    num_suffixes: usize,
}

// --------------------------------------------------
/// A struct for writing suffixes to disk.
#[derive(Debug)]
struct PartitionBuilder<T>
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
    /// Create a new `PartitionBuilder`. This struct is used to write it's
    /// suffix positions to a temporary file.
    ///
    /// Args:
    /// * `capacity`: the number of suffixes to hold in memory until the
    ///   writing to disk. This minimizes the number of times we access the disk
    ///   while also limiting the amount of memory used. Currently set to 4096
    ///   but it might be worth tuning this, perhaps use more memory to hit
    ///   disk less? Or if memory use is too high, lower and take a performance
    ///   hit for disk access?
    fn new(capacity: usize) -> Result<Self> {
        let tmp = NamedTempFile::new()?;
        let (_, path) = tmp.keep()?;

        Ok(PartitionBuilder {
            // Re-use a static vector to avoid repeated allocations
            vals: vec![T::default(); capacity],
            len: 0,
            total_len: 0,
            capacity,
            path,
        })
    }

    /// Add a suffix to the partition. When the internal array of values hits
    /// `capacity`, then write all the values to disk.
    ///
    /// Args:
    /// * `val`: the suffix position to add
    pub fn add(&mut self, val: T) -> Result<()> {
        self.vals[self.len] = val;
        self.len += 1;
        if self.len == self.capacity {
            self.write()?;
            self.len = 0;
        }

        Ok(())
    }

    /// Write the suffixes to disk. This must be called at the end to flush
    /// any remaining values after the last call(s) from `add`.
    pub fn write(&mut self) -> Result<()> {
        if self.len > 0 {
            let mut file = OpenOptions::new()
                .create(true)
                .append(true)
                .open(&self.path)?;
            file.write_all(vec_to_slice_u8(&self.vals[0..self.len]))?;
            self.total_len += self.len;
        }
        Ok(())
    }
}

// --------------------------------------------------
#[cfg(test)]
mod test {
    use super::{SufrBuilder, SufrBuilderArgs};
    use anyhow::Result;
    use pretty_assertions::assert_eq;

    #[test]
    fn test_is_less() -> Result<()> {
        //           012345
        let text = b"TTTAGC".to_vec();
        let args = SufrBuilderArgs {
            text,
            low_memory: true,
            path: None,
            max_query_len: None,
            is_dna: false,
            allow_ambiguity: false,
            ignore_softmask: false,
            sequence_starts: vec![0],
            sequence_names: vec!["1".to_string()],
            num_partitions: 2,
            seed_mask: None,
            random_seed: 0,
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

        Ok(())
    }

    #[test]
    fn test_is_less_max_query_len() -> Result<()> {
        //           012345
        let text = b"TTTAGC".to_vec();
        let args = SufrBuilderArgs {
            text,
            low_memory: true,
            path: None,
            max_query_len: Some(2),
            is_dna: false,
            allow_ambiguity: false,
            ignore_softmask: false,
            sequence_starts: vec![0],
            sequence_names: vec!["1".to_string()],
            num_partitions: 2,
            seed_mask: None,
            random_seed: 0,
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

        Ok(())
    }

    #[test]
    fn test_is_less_seed_mask() -> Result<()> {
        //           012345
        let text = b"TTTTAT".to_vec();
        let args = SufrBuilderArgs {
            text,
            low_memory: true,
            path: None,
            max_query_len: None,
            is_dna: false,
            allow_ambiguity: false,
            ignore_softmask: false,
            sequence_starts: vec![0],
            sequence_names: vec!["1".to_string()],
            num_partitions: 2,
            seed_mask: Some("101".to_string()),
            random_seed: 0,
        };
        let sufr: SufrBuilder<u32> = SufrBuilder::new(args)?;

        // 0: TTTTAT
        // 1: TTTAT
        // "T-T" vs "T-T"
        assert!(!sufr.is_less(0, 1));

        // 1: TTTAT
        // 0: TTTTAT
        // "T-T" vs "T-T"
        //assert!(!sufr.is_less(1, 0));

        // 0: TTTTAT
        // 3: TAT
        // "T-T" vs "T-T"
        //assert!(!sufr.is_less(0, 3));

        // 3: TAT
        // 0: TTTTAT
        // "T-T" vs "T-T"
        //assert!(!sufr.is_less(3, 0));

        Ok(())
    }

    #[test]
    fn test_find_lcp_no_seed_mask() -> Result<()> {
        //           012345
        let text = b"TTTAGC".to_vec();
        let args = SufrBuilderArgs {
            text,
            low_memory: true,
            path: None,
            max_query_len: None,
            is_dna: false,
            allow_ambiguity: false,
            ignore_softmask: false,
            sequence_starts: vec![0],
            sequence_names: vec!["1".to_string()],
            num_partitions: 2,
            seed_mask: None,
            random_seed: 0,
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

        Ok(())
    }

    #[test]
    fn test_find_lcp_with_seed_mask() -> Result<()> {
        //           012345
        let text = b"TTTTTA".to_vec();
        let args = SufrBuilderArgs {
            text,
            low_memory: true,
            path: None,
            max_query_len: None,
            is_dna: false,
            allow_ambiguity: false,
            ignore_softmask: false,
            sequence_starts: vec![0],
            sequence_names: vec!["1".to_string()],
            num_partitions: 2,
            seed_mask: Some("1101".to_string()),
            random_seed: 42,
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

        Ok(())
    }

    #[test]
    fn test_upper_bound_1() -> Result<()> {
        //          012345
        let text = b"TTTAGC".to_vec();
        let args = SufrBuilderArgs {
            text,
            low_memory: true,
            path: None,
            max_query_len: None,
            is_dna: false,
            allow_ambiguity: false,
            ignore_softmask: false,
            sequence_starts: vec![0],
            sequence_names: vec!["1".to_string()],
            num_partitions: 2,
            seed_mask: None,
            random_seed: 42,
        };
        let sufr: SufrBuilder<u32> = SufrBuilder::new(args)?;

        // The suffix "AGC$" is found before "GC$" and "C$
        assert_eq!(sufr.upper_bound(3, &[5, 4]), 0);

        // The suffix "TAGC$" is beyond all the values
        assert_eq!(sufr.upper_bound(2, &[3, 4, 5]), 3);

        // The "C$" is the last value
        assert_eq!(sufr.upper_bound(5, &[3, 4, 5]), 1);

        Ok(())
    }

    #[test]
    fn test_upper_bound_2() -> Result<()> {
        //           0123456789
        let text = b"ACGTNNACGT".to_vec();
        let args = SufrBuilderArgs {
            text,
            low_memory: true,
            path: None,
            max_query_len: None,
            is_dna: false,
            allow_ambiguity: false,
            ignore_softmask: false,
            sequence_starts: vec![0],
            sequence_names: vec!["1".to_string()],
            num_partitions: 2,
            seed_mask: None,
            random_seed: 42,
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

        Ok(())
    }

    #[test]
    fn test_upper_bound_seed_mask() -> Result<()> {
        //           0123456789
        let text = b"ACGTNNACGT".to_vec();
        let args = SufrBuilderArgs {
            text,
            low_memory: true,
            path: None,
            max_query_len: None,
            is_dna: false,
            allow_ambiguity: false,
            ignore_softmask: false,
            sequence_starts: vec![0],
            sequence_names: vec!["1".to_string()],
            num_partitions: 2,
            seed_mask: Some("101".to_string()),
            random_seed: 42,
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

        Ok(())
    }
}
