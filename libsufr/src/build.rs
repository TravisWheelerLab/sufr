use crate::{
    format::{Header, SeedMask, Sort, OUTFILE_VERSION},
    io::{as_bytes, as_bytes_mut, write_all_at, Int},
    lcp::{mask_offset, masked_lcp, suffix_lcp, LcpCache},
};
use anyhow::{anyhow, bail, Result};
use log::info;
use needletail::parse_fastx_file;
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

/// The character placed at the end of the text, which must not occur in
/// the text.
const SENTINEL_CHARACTER: u8 = b'$';

/// Configuration and construction of a _.sufr_ file.
///
/// Start from a sequence file with [`Builder::from_sequence_file`] or
/// from text in memory with [`Builder::from_text`], set options, then
/// write the file with [`Builder::write`].
pub struct Builder {
    text: Vec<u8>,
    sequences: Vec<(usize, String)>,
    is_dna: bool,
    allow_ambiguity: bool,
    ignore_softmask: bool,

    /// The seed mask the suffixes are sorted by, if any.
    seed_mask: Option<SeedMask>,

    /// The number of characters the suffixes are sorted by, or 0 for all.
    max_query_len: usize,

    num_partitions: usize,
    memory_scratch: bool,
    force_64bit: bool,
    write_lcp: bool,

    /// The runs of at least 1000 `N`s in the text.
    n_ranges: Vec<Range<usize>>,
}

impl Builder {
    /// A builder for the sequences in a FASTA or FASTQ file, joined by
    /// `sequence_delimiter`, which must be in the alphabet or sort below
    /// all of it: `%` suits biological sequence.
    pub fn from_sequence_file(
        path: impl AsRef<Path>,
        sequence_delimiter: u8,
    ) -> Result<Self> {
        // TODO: "must be" - should either be nonconfigurable or should support other values.
        let now = Instant::now();
        let data = read_sequence_file(path.as_ref(), sequence_delimiter)?;
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
        Ok(Self::from_text(data.seq, sequences))
    }

    /// A builder for `text`, whose `sequences` are given by start position
    /// and name, with a `%` between each two.
    pub fn from_text(text: Vec<u8>, sequences: Vec<(usize, String)>) -> Self {
        // TODO: validate text and sequences?
        // - sequences should be sorted, first starting at 0, length < text.len(), no duplicates. If empty, default?
        Self {
            text,
            sequences,
            is_dna: false,
            allow_ambiguity: false,
            ignore_softmask: false,
            seed_mask: None,
            max_query_len: 0,
            num_partitions: 1024,
            memory_scratch: false,
            force_64bit: false,
            write_lcp: false,
            n_ranges: vec![],
        }
    }

    /// Set whether the text `is_dna`, and if so whether to index `N`s
    /// (`allow_ambiguity`) and whether to replace softmasked (lowercase)
    /// characters with `N` rather than upper-case them (`ignore_softmask`).
    /// All three default to `false`.
    pub fn dna(
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

    /// Set how the suffixes are sorted (default: `Sort::MaxQueryLen(0)`,
    /// by all of their characters).
    pub fn sort(mut self, sort: Sort) -> Self {
        (self.seed_mask, self.max_query_len) = match sort {
            Sort::Mask(seed_mask) => (Some(seed_mask), 0),
            Sort::MaxQueryLen(max_query_len) => (None, max_query_len),
        };
        self
    }

    /// Set the number of partitions the suffixes are sorted in
    /// (default: 1024).
    pub fn partitions(mut self, num_partitions: usize) -> Self {
        self.num_partitions = num_partitions;
        self
    }

    /// Set whether to keep the partitions in memory during construction
    /// instead of in temporary files (default: `false`).
    pub fn memory_scratch(mut self, memory_scratch: bool) -> Self {
        self.memory_scratch = memory_scratch;
        self
    }

    // only Index::<u64> can read a small text written this
    // way, so this is for tests
    #[cfg(test)]
    pub(crate) fn force_64bit(mut self, force_64bit: bool) -> Self {
        self.force_64bit = force_64bit;
        self
    }

    /// Set whether to write the LCP array (default: `false`).
    pub fn write_lcp(mut self, write_lcp: bool) -> Self {
        self.write_lcp = write_lcp;
        self
    }

    /// Build the suffix array and write it to `path`.
    pub fn write(mut self, path: impl AsRef<Path>) -> Result<()> {
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
        Ok(())
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
        let partitions = self.sort_partitions::<T, B>(occupancy)?;
        self.write_partitions::<T, B>(path, partitions)?;
        Ok(())
    }

    /// The end of the run of `N`s holding `suffix`, if it is in one.
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

    /// The LCP of the suffixes at `start1` and `start2`, which share their
    /// first `skip` characters, comparing at most `len` characters.
    #[inline(always)]
    fn find_lcp<T: Int>(
        &self,
        start1: usize,
        start2: usize,
        len: T,
        skip: usize,
        lcp_cache: Option<&mut LcpCache>,
    ) -> T {
        // with a mask, `len` is not consulted: the comparison
        // stops at the end of the mask or of the text
        if let Some(mask) = &self.seed_mask {
            return T::from_usize(
                skip + masked_lcp(
                    &self.text[start1..],
                    &self.text[start2..],
                    &mask.positions[skip..],
                ),
            );
        }

        if let (Some(end1), Some(end2)) =
            (self.find_n_run(start1), self.find_n_run(start2))
        {
            return T::from_usize(min(end1 - start1, end2 - start2));
        }

        // TODO: a nonzero max_query_len replaces `len` rather than
        // capping it, so past a nonzero `skip` the comparison runs
        // beyond max_query_len and the LCP array can hold values
        // above it
        let len = if self.max_query_len > 0 {
            self.max_query_len
        } else {
            len.to_usize()
        };
        T::from_usize(suffix_lcp(&self.text, start1, start2, skip, len, lcp_cache))
    }

    /// Place each suffix in the scratch buffer for its `A::COUNT`
    /// character prefix.
    fn partition<T, B, A>(&self) -> Result<PartitionedSuffixes<B>>
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
                        || !self.is_dna
                        || (b"ACGT".contains(&val) || self.allow_ambiguity))
                    {
                        continue;
                    }

                    let mut key = 0;
                    match &self.seed_mask {
                        Some(mask) => {
                            for (j, &off) in
                                mask.positions.iter().enumerate().take(A::COUNT)
                            {
                                let c = self.text.get(i + off).copied().unwrap_or(0);
                                key |= alphabet.lookup(c)
                                    << ((A::COUNT - j - 1) * A::BITS);
                            }
                        }
                        None => {
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

        Ok(PartitionedSuffixes {
            buffers,
            num_suffixes,
        })
    }

    /// Sort the suffixes in `num_partitions` partitions, choosing the
    /// prefix they are first placed by from the bytes that `occupancy`
    /// marks as present.
    fn sort_partitions<T: Int, B: ScratchBuffer<Item = T> + Send + Sync>(
        &self,
        occupancy: [u8; 256],
    ) -> Result<Vec<Partition<B>>> {
        let num_partitions = self.num_partitions;
        let partitioned = self.partition_by_alphabet::<T, B>(occupancy)?;
        let total_sort_time = Instant::now();
        let partition_inputs = group_buffers(partitioned, num_partitions)?;

        let mut partitions: Vec<Option<Partition<B>>> =
            (0..num_partitions).map(|_| None).collect();

        let lcp_cache: ThreadLocal<RefCell<LcpCache>> = ThreadLocal::new();
        partitions
            .par_iter_mut()
            .zip(partition_inputs)
            .enumerate()
            .try_for_each(|(order, (partition, buffers))| -> Result<()> {
                *partition = self.sort_partition(order, buffers, &lcp_cache)?;
                Ok(())
            })?;

        // drop the empty partitions and put the rest in order
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

    /// Partition the suffixes by the nucleic acid alphabet when the text
    /// holds only its bytes, else by the amino acid alphabet when it holds
    /// only those, else by bytes.
    fn partition_by_alphabet<T: Int, B: ScratchBuffer<Item = T> + Send + Sync>(
        &self,
        mut occupancy: [u8; 256],
    ) -> Result<PartitionedSuffixes<B>> {
        for &c in NucleicAcidAlphabet::ALPHABET.iter() {
            occupancy[c as usize] = 0;
        }
        if occupancy.iter().map(|&b| b as usize).sum::<usize>() == 0 {
            info!("Detected Nucleic Acid alphabet: using 5-character prefix for partitioning");
            return self.partition::<T, B, NucleicAcidAlphabet>();
        }

        for &c in AminoAcidAlphabet::ALPHABET.iter() {
            occupancy[c as usize] = 0;
        }
        if occupancy.iter().map(|&b| b as usize).sum::<usize>() == 0 {
            info!("Detected Amino Acid alphabet: using 3-character prefix for partitioning");
            self.partition::<T, B, AminoAcidAlphabet>()
        } else {
            info!("Detected non-Nucleic/non-Amino Acid: using 2-byte prefix for partitioning");
            self.partition::<T, B, BytesAlphabet>()
        }
    }

    /// Sort the suffixes in `buffers` into the partition at `order`, or
    /// `None` when there are none.
    fn sort_partition<T: Int, B: ScratchBuffer<Item = T>>(
        &self,
        order: usize,
        buffers: Vec<B>,
        lcp_cache: &ThreadLocal<RefCell<LcpCache>>,
    ) -> Result<Option<Partition<B>>> {
        let mut part_sa = vec![];
        for buffer in buffers {
            part_sa.extend_from_slice(&buffer.consume()?);
        }

        let len = part_sa.len();
        if len == 0 {
            return Ok(None);
        }

        let mut sa_w = part_sa.clone();
        let mut lcp = vec![T::default(); len];
        let mut lcp_w = vec![T::default(); len];
        self.merge_sort(
            Arrays {
                sa: &mut sa_w,
                lcp: &mut lcp_w,
            },
            Arrays {
                sa: &mut part_sa,
                lcp: &mut lcp,
            },
            lcp_cache,
        );

        let mut sa_buffer = B::default();
        sa_buffer.extend_from_slice(&part_sa)?;
        let lcp_buffer = if self.write_lcp {
            let mut lcp_buffer = B::default();
            lcp_buffer.extend_from_slice(&lcp)?;
            Some(lcp_buffer)
        } else {
            None
        };

        Ok(Some(Partition {
            order,
            len,
            first_suffix: part_sa.first().unwrap().to_usize(),
            last_suffix: part_sa.last().unwrap().to_usize(),
            sa_buffer,
            lcp_buffer,
        }))
    }

    /// Sort the suffixes in `out` with their LCPs, using `work`, which
    /// starts with the same suffixes, as scratch space.
    fn merge_sort<T: Int>(
        &self,
        mut work: Arrays<'_, T>,
        mut out: Arrays<'_, T>,
        lcp_cache: &ThreadLocal<RefCell<LcpCache>>,
    ) {
        let n = work.sa.len();
        if n == 1 {
            out.lcp[0] = T::default();
        } else {
            const NESTED_PAR_GRAIN_SIZE: usize = 1 << 13;
            let mid = n / 2;
            let (work_l, work_r) = work.split_at_mut(mid);
            let (out_l, out_r) = out.split_at_mut(mid);

            // each half is sorted into `work`, swapping the
            // roles, so the merge below reads `work`
            if mid < NESTED_PAR_GRAIN_SIZE {
                self.merge_sort(out_l, work_l, lcp_cache);
                self.merge_sort(out_r, work_r, lcp_cache);
            } else {
                rayon::join(
                    || self.merge_sort(out_l, work_l, lcp_cache),
                    || self.merge_sort(out_r, work_r, lcp_cache),
                );
            }

            self.merge(work, mid, out, &mut lcp_cache.get_or_default().borrow_mut());
        }
    }

    /// Merge the sorted runs `from[..mid]` and `from[mid..]` into `into`.
    fn merge<T: Int>(
        &self,
        from: Arrays<'_, T>,
        mid: usize,
        into: Arrays<'_, T>,
        lcp_cache: &mut LcpCache,
    ) {
        let Arrays {
            sa: target_sa,
            lcp: target_lcp,
        } = into;
        let (mut x, mut y) = from.sa.split_at_mut(mid);
        let (mut lcp_x, mut lcp_y) = from.lcp.split_at_mut(mid);
        let mut len_x = x.len();
        let mut len_y = y.len();
        // the last LCP taken from x
        let mut m = T::default();
        let mut idx_x = 0;
        let mut idx_y = 0;
        let mut idx_target = 0;

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

                    let context = match &self.seed_mask {
                        Some(seed_mask) => T::from_usize(
                            seed_mask
                                .positions
                                .iter()
                                .filter(|&i| *i < max_n.to_usize())
                                .count(),
                        ),
                        None if self.max_query_len > 0 => {
                            min(T::from_usize(self.max_query_len), max_n)
                        }
                        None => max_n,
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
                        let full_lcp = match &self.seed_mask {
                            Some(seed_mask) => mask_offset(lcp.to_usize(), seed_mask),
                            None => lcp.to_usize(),
                        };
                        (lcp, T::from_usize(full_lcp))
                    } else {
                        (context, context)
                    };

                    // an LCP reaching the context sorts the shorter suffix
                    // first, and otherwise the byte after the LCP decides
                    if len_lcp >= context {
                        target_sa[idx_target] = shorter_suffix;
                    } else {
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

        while idx_x < len_x {
            target_sa[idx_target] = x[idx_x];
            target_lcp[idx_target] = lcp_x[idx_x];
            idx_x += 1;
            idx_target += 1;
        }

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

    /// Write the text, then each partition's suffix array (and LCP array,
    /// when requested) in parallel at its offset, then the header and
    /// sequence names. Returns the number of bytes written.
    fn write_partitions<T: Int, B: ScratchBuffer<Item = T> + Send + Sync>(
        &self,
        path: &Path,
        partitions: Vec<Partition<B>>,
    ) -> Result<usize> {
        let now = Instant::now();
        let num_suffixes = partitions.iter().map(|p| p.len).sum();
        let int_size = mem::size_of::<T>();
        let header = self.header(int_size, num_suffixes);
        let file =
            File::create(path).map_err(|e| anyhow!("{}: {e}", path.display()))?;
        write_all_at(&file, &self.text, header.text_pos as u64)?;

        let offsets = partition_offsets(&partitions);

        // consuming the partitions' buffers frees memory
        // and disk space along the way
        partitions.into_par_iter().zip(offsets).try_for_each(
            |(partition, (offset, prev_last_suffix))| -> Result<()> {
                let sa = partition.sa_buffer.consume()?;
                write_all_at(
                    &file,
                    as_bytes(&sa),
                    (header.suffix_array_pos + offset * int_size) as u64,
                )?;
                drop(sa);

                if let (Some(lcp_buffer), Some(lcp_pos)) =
                    (partition.lcp_buffer, header.lcp_pos)
                {
                    let mut lcp = lcp_buffer.consume()?;
                    // the first LCP of a partition is with the last
                    // suffix of the one before
                    if let Some(prev_last_suffix) = prev_last_suffix {
                        if let Some(val) = lcp.first_mut() {
                            *val = self.find_lcp(
                                prev_last_suffix,
                                partition.first_suffix,
                                T::from_usize(self.text.len()),
                                0,
                                None,
                            );
                        }
                    }
                    write_all_at(
                        &file,
                        as_bytes(&lcp),
                        (lcp_pos + offset * int_size) as u64,
                    )?;
                }
                Ok(())
            },
        )?;

        let bytes_out = header.write(&file)?;
        info!(
            "Wrote suffix array to '{}' in {:?}",
            path.display(),
            now.elapsed()
        );
        Ok(bytes_out)
    }

    /// The header of the output file, with `int_size` bytes per
    /// suffix array value.
    fn header(&self, int_size: usize, num_suffixes: usize) -> Header {
        let mut header = Header {
            version: OUTFILE_VERSION,
            is_dna: self.is_dna,
            allow_ambiguity: self.allow_ambiguity,
            ignore_softmask: self.ignore_softmask,
            text_len: self.text.len(),
            num_suffixes,
            sort: match &self.seed_mask {
                Some(seed_mask) => Sort::Mask(seed_mask.clone()),
                None => Sort::MaxQueryLen(self.max_query_len),
            },
            sequence_starts: self.sequences.iter().map(|(start, _)| *start).collect(),
            sequence_names: self
                .sequences
                .iter()
                .map(|(_, name)| name.clone())
                .collect(),
            int_size,
            text_pos: 0,
            suffix_array_pos: 0,
            lcp_pos: None,
        };
        header.set_positions(self.write_lcp);
        header
    }
}

/// The offset (in elements) of each partition in the arrays, with the
/// last suffix of the partition before it.
fn partition_offsets<B: ScratchBuffer>(
    partitions: &[Partition<B>],
) -> Vec<(usize, Option<usize>)> {
    let mut offsets = Vec::with_capacity(partitions.len());
    let mut total = 0;
    let mut prev_last_suffix = None;
    for partition in partitions {
        offsets.push((total, prev_last_suffix));
        total += partition.len;
        prev_last_suffix = Some(partition.last_suffix);
    }
    offsets
}

/// Gather the buffers, which are in prefix order, into `num_partitions`
/// groups of about the same number of suffixes.
fn group_buffers<B: ScratchBuffer>(
    partitioned: PartitionedSuffixes<B>,
    num_partitions: usize,
) -> Result<Vec<Vec<B>>> {
    // round up to take every suffix
    let num_per_partition =
        (partitioned.num_suffixes as f64 / num_partitions as f64).ceil() as usize;
    let mut num_taken = 0;
    let mut partition_inputs = std::iter::repeat_with(|| Vec::new())
        .take(num_partitions)
        .collect::<Vec<_>>();

    // there are (probably) many more buffers than partitions,
    // so each partition takes buffers from the left until it
    // passes its share, which evens out the sorting work
    let mut part_buffers = partitioned.buffers.into_iter();
    for (partition_num, partition_input) in partition_inputs.iter_mut().enumerate() {
        let boundary = num_per_partition * (partition_num + 1);
        for buffer in part_buffers.by_ref() {
            if buffer.count() > 0 {
                num_taken += buffer.count();
                partition_input.push(buffer);
            }

            // the last partition takes the rest
            if partition_num < num_partitions - 1 && num_taken > boundary {
                break;
            }
        }
    }

    if num_taken != partitioned.num_suffixes {
        bail!(
            "Took {num_taken} but needed to take {}",
            partitioned.num_suffixes
        );
    }
    Ok(partition_inputs)
}

/// The text and sequences read from a FASTA/FASTQ file.
struct SequenceFileData {
    /// The sequences joined by the delimiter, ending in the sentinel.
    seq: Vec<u8>,

    /// The start of each sequence in `seq`.
    start_positions: Vec<usize>,

    /// The name of each sequence.
    sequence_names: Vec<String>,
}

/// Read the sequences in the FASTA/FASTQ file at `path`, joined by
/// `sequence_delimiter`.
fn read_sequence_file(path: &Path, sequence_delimiter: u8) -> Result<SequenceFileData> {
    let mut reader = parse_fastx_file(path)?;
    let mut seq: Vec<u8> = Vec::with_capacity(u32::MAX as usize);
    let mut sequence_names: Vec<String> = vec![];
    let mut start_positions: Vec<usize> = vec![];
    let mut i = 0;
    while let Some(rec) = reader.next() {
        let rec = rec?;
        if i > 0 {
            seq.push(sequence_delimiter);
        }

        start_positions.push(seq.len());
        seq.extend_from_slice(&rec.seq());
        i += 1;

        // the name is the ID up to its first whitespace
        let id = std::str::from_utf8(rec.id())?
            .split_whitespace()
            .next()
            .map_or((i + 1).to_string(), |v| v.to_string());

        sequence_names.push(id);
    }

    seq.push(SENTINEL_CHARACTER);

    Ok(SequenceFileData {
        seq,
        start_positions,
        sequence_names,
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

/// The sorted suffixes and LCPs of one partition.
#[derive(Debug)]
struct Partition<B: ScratchBuffer> {
    /// The position of this partition in the suffix array.
    order: usize,

    /// The number of suffixes/LCP values contained in this partition.
    len: usize,

    /// The first suffix, for the LCP at the boundary with the partition
    /// before.
    first_suffix: usize,

    /// The last suffix, for the LCP at the boundary with the partition
    /// after.
    last_suffix: usize,

    /// The buffer containing the suffix array.
    sa_buffer: B,

    /// The buffer containing the LCP array, when one is written.
    lcp_buffer: Option<B>,
}

/// The unsorted suffixes in one scratch buffer per prefix, in prefix
/// order.
#[derive(Debug)]
struct PartitionedSuffixes<B: ScratchBuffer<Item: Int>> {
    buffers: Vec<B>,
    num_suffixes: usize,
}

/// A slice of the suffix array and the LCPs beside it.
struct Arrays<'a, T> {
    sa: &'a mut [T],
    lcp: &'a mut [T],
}

impl<T> Arrays<'_, T> {
    fn split_at_mut(&mut self, mid: usize) -> (Arrays<'_, T>, Arrays<'_, T>) {
        let (sa_l, sa_r) = self.sa.split_at_mut(mid);
        let (lcp_l, lcp_r) = self.lcp.split_at_mut(mid);
        (
            Arrays {
                sa: sa_l,
                lcp: lcp_l,
            },
            Arrays {
                sa: sa_r,
                lcp: lcp_r,
            },
        )
    }
}

/// A buffer of suffixes or LCPs during construction, in memory or in a
/// temporary file.
trait ScratchBuffer: Default {
    /// The type of the values.
    type Item: Int;

    /// Append `vals`.
    fn extend_from_slice(&mut self, vals: &[Self::Item]) -> Result<()>;

    /// Write any values still in memory to the backing store.
    fn flush(&mut self) -> Result<()>;

    /// The number of values appended.
    fn count(&self) -> usize;

    /// All the values, borrowed when they are in memory.
    fn read(&self) -> Result<Cow<'_, [Self::Item]>>;

    /// All the values, freeing the backing store.
    fn consume(self) -> Result<Vec<Self::Item>>;
}

/// A scratch buffer in a temporary file.
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
            file.write_all(as_bytes(&self.buf))?;
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
                file.read_exact(as_bytes_mut(&mut data))?;

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
            // best effort, since this is Drop and consume()
            // may have removed the file already
            let _ = fs::remove_file(path);
        }
    }
}

/// A scratch buffer in memory.
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

#[cfg(test)]
mod test {
    use super::{read_sequence_file, Builder};
    use crate::{
        format::{SeedMask, Sort},
        SuffixArray,
    };
    use anyhow::Result;
    use pretty_assertions::assert_eq;
    use std::{fs, path::Path};
    use tempfile::NamedTempFile;

    #[test]
    fn test_read_sequence_file() -> Result<()> {
        let file = Path::new("../data/inputs/2.fa");
        let sequence_delimiter = b'N';
        let res = read_sequence_file(file, sequence_delimiter);
        assert!(res.is_ok());
        let data = res.unwrap();
        assert_eq!(data.seq, b"ACGTacgtNacgtACGT$");
        assert_eq!(data.start_positions, [0, 9]);
        assert_eq!(data.sequence_names, ["ABC", "DEF"]);
        Ok(())
    }

    #[test]
    fn test_find_lcp_no_seed_mask() -> Result<()> {
        //           012345
        let text = b"TTTAGC".to_vec();
        let sufr = Builder::from_text(text, vec![(0, "1".to_string())]);

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
        let sufr = Builder::from_text(text, vec![(0, "1".to_string())])
            .sort(Sort::Mask(SeedMask::new("1101")?));

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
            Builder::from_sequence_file("../data/inputs/2.fa", b'%')?
                .memory_scratch(memory_scratch)
                .write_lcp(true)
                .write(outfile.path())?;
            Ok(fs::read(outfile.path())?)
        };
        assert_eq!(build(true)?, build(false)?);
        Ok(())
    }

    #[test]
    fn test_write() -> Result<()> {
        let outfile = NamedTempFile::new()?;
        Builder::from_sequence_file("../data/inputs/3.fa", b'%')?
            .dna(true, false, true)
            .write(outfile.path())?;
        let suffix_array = SuffixArray::open(&outfile.path().to_string_lossy())?;
        let meta = suffix_array.metadata()?;
        assert_eq!(meta.text_len, 113);
        assert_eq!(meta.num_suffixes, 101);
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
