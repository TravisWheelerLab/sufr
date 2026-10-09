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
    file_access::{read_exact_at, write_at},
    lcp::{lcp, LcpCache},
    radix::{sort_items, Item, KeyEncoder},
    types::{
        Int, SeedMask, SortStrategy, SuffixSortType, SufrBuilderArgs,
        OUTFILE_VERSION, SENTINEL_CHARACTER,
    },
    util::{find_lcp_full_offset, slice_int_to_slice_u8, vec_to_slice_u8},
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
    path::PathBuf,
    sync::Mutex,
    time::Instant,
};
use tempfile::NamedTempFile;
use thread_local::ThreadLocal;

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

    /// Whether to write the LCP array
    pub write_lcp: bool,
}

// --------------------------------------------------
/// The output file after the header and text have been written.
/// The suffix and LCP arrays are written into it at known offsets,
/// possibly from several threads.
struct Output {
    file: File,

    /// Where the header records the text/SA/LCP positions
    locs_pos: u64,

    /// Where the header records the number of suffixes
    num_suffixes_pos: u64,

    /// Byte position of the text
    text_pos: usize,

    /// Byte position of the suffix array
    sa_pos: usize,

    /// Byte position of the LCP array, or 0 when none is written.
    /// Valid after `set_num_suffixes`.
    lcp_pos: usize,

    /// Byte position after the arrays, where the sequence names go.
    /// Valid after `set_num_suffixes`.
    end_pos: usize,

    /// Bytes per suffix array element
    int_size: usize,

    /// Whether an LCP array follows the suffix array
    write_lcp: bool,
}

impl Output {
    /// Fix the array layout once the number of suffixes is known.
    fn set_num_suffixes(&mut self, num_suffixes: usize) {
        let sa_end = self.sa_pos + num_suffixes * self.int_size;
        if self.write_lcp {
            self.lcp_pos = sa_end;
            self.end_pos = sa_end + num_suffixes * self.int_size;
        } else {
            self.lcp_pos = 0;
            self.end_pos = sa_end;
        }
    }
}

/// The first and last sorted item of each radix partition, for the LCP
/// at each partition boundary
type PartitionEnds<T> = Mutex<Vec<Option<(Item<T>, Item<T>)>>>;

/// Fill `buf` from the file starting at `offset`, in parallel 64 MB chunks.
fn read_all_at(file: &File, buf: &mut [u8], offset: u64) -> Result<()> {
    let chunk_len = 1 << 26;
    buf.par_chunks_mut(chunk_len)
        .enumerate()
        .try_for_each(|(i, chunk)| -> Result<()> {
            let pos = offset + (i * chunk_len) as u64;
            read_exact_at(file, chunk, pos)?;
            Ok(())
        })
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
    ///     types::{SortStrategy, SufrBuilderArgs},
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
    ///         sort_strategy: SortStrategy::Merge,
    ///         write_lcp: true,
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

        // Normalize lowercase in place, in parallel, and count the byte
        // occurrences so that the best alphabet can be selected for
        // partitioning
        let ignore_softmask = args.ignore_softmask;
        let occupancy: [u8; 256] = text
            .par_chunks_mut(1 << 22)
            .map(|chunk| {
                let mut occupancy = [0u8; 256];
                for b in chunk.iter_mut() {
                    // Check for lowercase
                    if (97..=122).contains(b) {
                        if ignore_softmask {
                            *b = b'N'
                        } else {
                            // only shift lowercase ASCII
                            *b &= 0b1011111
                        }
                    }
                    occupancy[*b as usize] |= 1;
                }
                occupancy
            })
            .reduce(
                || [0u8; 256],
                |mut a, b| {
                    for (x, y) in a.iter_mut().zip(b) {
                        *x |= y;
                    }
                    a
                },
            );
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
            write_lcp: args.write_lcp,
        };
        match args.sort_strategy {
            SortStrategy::Merge => {
                sa.sort(args.num_partitions, occupancy)?;
                sa.write()?;
            }
            SortStrategy::RadixInMemory => sa.sort_radix_in_memory()?,
            SortStrategy::RadixPartitioned => {
                sa.sort_radix_partitioned(args.num_partitions)?
            }
        }
        Ok(sa)
    }

    // --------------------------------------------------
    /// Whether the suffix starting with byte `val` is indexed
    fn is_suffix_start(&self, val: u8) -> bool {
        val == SENTINEL_CHARACTER
            || !self.is_dna
            || (b"ACGT".contains(&val) || self.allow_ambiguity)
    }

    // --------------------------------------------------
    /// The offsets that make up a suffix's packed radix key: the care
    /// positions of the seed mask, or `0..max_query_len`.
    fn radix_key_encoder(&self) -> Result<KeyEncoder> {
        let positions: Vec<usize> = match &self.sort_type {
            SuffixSortType::Mask(seed_mask) => seed_mask.positions.clone(),
            SuffixSortType::MaxQueryLen(0) => {
                bail!("Radix sort requires a seed mask or a max query length")
            }
            SuffixSortType::MaxQueryLen(max_query_len) => {
                if !self.n_ranges.is_empty() {
                    bail!(
                        "Radix sort does not support long runs of Ns with \
                         ambiguity allowed in max-query-len mode"
                    );
                }
                (0..*max_query_len).collect()
            }
        };
        KeyEncoder::new(&self.text, &positions)
    }

    // --------------------------------------------------
    /// Which byte ranks start an indexed suffix (for a text already
    /// converted to ranks).
    fn start_ranks(&self, encoder: &KeyEncoder) -> [bool; 256] {
        encoder.rank_table(|b| self.is_suffix_start(b))
    }

    // --------------------------------------------------
    /// All indexed suffix positions in descending order, from the
    /// rank-converted text. Equal keys keep this order through the sort:
    /// the shorter suffix (larger position) comes first, as in `merge`.
    fn suffix_positions_descending(&self, start_rank: &[bool; 256]) -> Vec<T> {
        let chunk_len = 1 << 22;
        let mut chunks: Vec<Vec<T>> = self
            .text
            .par_chunks(chunk_len)
            .enumerate()
            .map(|(chunk_num, chunk)| {
                let base = chunk_num * chunk_len;
                chunk
                    .iter()
                    .enumerate()
                    .rev()
                    .filter(|&(_, &val)| start_rank[val as usize])
                    .map(|(i, _)| T::from_usize(base + i))
                    .collect()
            })
            .collect();
        chunks.reverse();
        chunks.concat()
    }

    // --------------------------------------------------
    /// Write sorted items into the output as the suffix-array elements
    /// starting at `offset`, plus their LCPs when requested. `prev` is the
    /// last item of the preceding partition, for the LCP at the boundary.
    fn emit_radix_partition(
        &self,
        out: &Output,
        encoder: &KeyEncoder,
        items: &[Item<T>],
        offset: usize,
        prev: Option<Item<T>>,
    ) -> Result<()> {
        if items.is_empty() {
            return Ok(());
        }
        let int_size = mem::size_of::<T>();
        let chunk_len = 1 << 22;

        // Each chunk builds its slice of the suffix (and LCP) array and
        // writes it at its own offset, so no partition-sized copy is made.
        let now = Instant::now();
        items
            .par_chunks(chunk_len)
            .enumerate()
            .try_for_each(|(i, chunk)| -> Result<()> {
                let chunk_offset = offset + i * chunk_len;
                let sa: Vec<T> = chunk.iter().map(|item| item.pos).collect();
                write_all_at(
                    &out.file,
                    vec_to_slice_u8(&sa),
                    (out.sa_pos + chunk_offset * int_size) as u64,
                )?;
                if self.write_lcp {
                    let lcp: Vec<T> = chunk
                        .iter()
                        .enumerate()
                        .map(|(j, b)| {
                            let a = if j > 0 {
                                Some(&chunk[j - 1])
                            } else if i > 0 {
                                Some(&items[i * chunk_len - 1])
                            } else {
                                prev.as_ref()
                            };
                            match a {
                                Some(a) => T::from_usize(encoder.lcp(
                                    a.key,
                                    a.pos.to_usize(),
                                    b.key,
                                    b.pos.to_usize(),
                                )),
                                None => T::default(),
                            }
                        })
                        .collect();
                    write_all_at(
                        &out.file,
                        vec_to_slice_u8(&lcp),
                        (out.lcp_pos + chunk_offset * int_size) as u64,
                    )?;
                }
                Ok(())
            })?;
        info!(
            "Wrote {} suffixes{} at offset {offset} in {:?}",
            items.len(),
            if self.write_lcp { " and LCPs" } else { "" },
            now.elapsed()
        );
        Ok(())
    }

    // --------------------------------------------------
    /// Sort every suffix in memory by packed key with an LSD radix sort
    /// and write the output.
    fn sort_radix_in_memory(&mut self) -> Result<()> {
        let total_time = Instant::now();
        let encoder = self.radix_key_encoder()?;
        info!(
            "Radix keys: {} symbols x {} bits = {} bits, pext {}",
            encoder.positions.len(),
            encoder.bits,
            encoder.total_bits,
            if encoder.uses_pext() { "on" } else { "off" }
        );
        let start_rank = self.start_ranks(&encoder);

        // The output gets the original text; after that the in-memory
        // copy is converted to ranks for key extraction.
        let mut out = self.open_output()?;
        let now = Instant::now();
        encoder.convert_to_ranks(&mut self.text);
        info!("Converted text to ranks in {:?}", now.elapsed());

        let now = Instant::now();
        let positions = self.suffix_positions_descending(&start_rank);
        let num_suffixes = positions.len();
        out.set_num_suffixes(num_suffixes);
        let mut items: Vec<Item<T>> = positions
            .into_par_iter()
            .map(|pos| Item {
                key: encoder.key(&self.text, pos.to_usize()),
                pos,
            })
            .collect();
        info!("Computed {num_suffixes} keys in {:?}", now.elapsed());

        let now = Instant::now();
        let mut scratch: Vec<Item<T>> = Vec::new();
        sort_items(&mut items, &mut scratch, encoder.total_bits);
        drop(scratch);
        info!("Sorted {num_suffixes} suffixes in {:?}", now.elapsed());

        self.emit_radix_partition(&out, &encoder, &items, 0, None)?;
        drop(items);
        self.finish_output(&out, num_suffixes)?;
        self.num_suffixes = T::from_usize(num_suffixes);
        info!(
            "Sorted and wrote {num_suffixes} suffixes (in-memory radix) in {:?}",
            total_time.elapsed()
        );
        Ok(())
    }

    // --------------------------------------------------
    /// Bucket suffix positions on disk by the high bits of their packed
    /// key, grouping buckets into roughly `num_partitions` partitions of
    /// similar size, then radix-sort the partitions in memory, a few at
    /// a time (`SUFR_RADIX_CONCURRENCY`, default 2).
    fn sort_radix_partitioned(&mut self, num_partitions: usize) -> Result<()> {
        let total_time = Instant::now();
        let encoder = self.radix_key_encoder()?;
        let num_partitions = num_partitions.max(1);
        let digit_bits = encoder.total_bits.min(16);
        let shift = encoder.total_bits - digit_bits;
        let num_buckets = 1usize << digit_bits;
        info!(
            "Radix keys: {} symbols x {} bits = {} bits, pext {}; bucketing on top {digit_bits} bits",
            encoder.positions.len(),
            encoder.bits,
            encoder.total_bits,
            if encoder.uses_pext() { "on" } else { "off" }
        );
        let start_rank = self.start_ranks(&encoder);

        // The output gets the original text; after that the in-memory
        // copy is converted to ranks for key extraction.
        let mut out = self.open_output()?;
        let now = Instant::now();
        encoder.convert_to_ranks(&mut self.text);
        info!("Converted text to ranks in {:?}", now.elapsed());

        // Pass 1: histogram of the top digit over all indexed suffixes
        let now = Instant::now();
        let chunk_len = 1 << 22;
        let histogram: Vec<usize> = self
            .text
            .par_chunks(chunk_len)
            .enumerate()
            .map(|(chunk_num, chunk)| {
                let base = chunk_num * chunk_len;
                let mut hist = vec![0usize; num_buckets];
                for (i, &val) in chunk.iter().enumerate() {
                    if start_rank[val as usize] {
                        let key = encoder.key(&self.text, base + i);
                        hist[(key >> shift) as usize] += 1;
                    }
                }
                hist
            })
            .reduce(
                || vec![0usize; num_buckets],
                |mut a, b| {
                    for (x, y) in a.iter_mut().zip(b) {
                        *x += y;
                    }
                    a
                },
            );
        let num_suffixes: usize = histogram.iter().sum();
        out.set_num_suffixes(num_suffixes);
        info!(
            "Bucket histogram of {num_suffixes} suffixes in {:?}",
            now.elapsed()
        );

        // Group consecutive buckets into partitions of about equal size
        let target = num_suffixes.div_ceil(num_partitions).max(1);
        let mut bucket_to_partition = vec![0usize; num_buckets];
        let mut partition_num = 0;
        let mut running = 0;
        for (bucket, &count) in histogram.iter().enumerate() {
            bucket_to_partition[bucket] = partition_num;
            running += count;
            if running >= target && partition_num + 1 < num_partitions {
                partition_num += 1;
                running = 0;
            }
        }
        let num_partitions = partition_num + 1;

        // Pass 2a: count, per 1M-byte sub-chunk of the text, how many
        // positions go to each partition, so that every sub-chunk can write
        // its positions at a known offset in each partition file. Offsets
        // are assigned from the end of the text so each file is in
        // descending position order (the tie rule for equal keys).
        let now = Instant::now();
        let text_len = self.text.len();
        let int_size = mem::size_of::<T>();
        let sub_chunk = 1 << 20;
        let num_sub = text_len.div_ceil(sub_chunk);
        let text = &self.text;
        let sub_counts: Vec<Vec<u32>> = (0..num_sub)
            .into_par_iter()
            .map(|s| {
                let from = s * sub_chunk;
                let to = min(from + sub_chunk, text_len);
                let mut counts = vec![0u32; num_partitions];
                for pos in from..to {
                    if start_rank[text[pos] as usize] {
                        let key = encoder.key(text, pos);
                        counts[bucket_to_partition[(key >> shift) as usize]] += 1;
                    }
                }
                counts
            })
            .collect();
        let mut sub_offsets = vec![vec![0usize; num_partitions]; num_sub];
        let mut counts = vec![0usize; num_partitions];
        for s in (0..num_sub).rev() {
            for (p, count) in counts.iter_mut().enumerate() {
                sub_offsets[s][p] = *count;
                *count += sub_counts[s][p] as usize;
            }
        }
        info!("Counted positions per partition in {:?}", now.elapsed());

        // Pass 2b: scatter positions to one file per partition with
        // positioned writes from every thread
        let now = Instant::now();
        let mut files: Vec<(PathBuf, File)> = Vec::with_capacity(num_partitions);
        for _ in 0..num_partitions {
            let (file, path) = NamedTempFile::new()?.keep()?;
            files.push((path, file));
        }
        (0..num_sub)
            .into_par_iter()
            .try_for_each(|s| -> Result<()> {
                let from = s * sub_chunk;
                let to = min(from + sub_chunk, text_len);
                let mut buckets: Vec<Vec<T>> = (0..num_partitions)
                    .map(|p| Vec::with_capacity(sub_counts[s][p] as usize))
                    .collect();
                for pos in (from..to).rev() {
                    if start_rank[text[pos] as usize] {
                        let key = encoder.key(text, pos);
                        let partition = bucket_to_partition[(key >> shift) as usize];
                        buckets[partition].push(T::from_usize(pos));
                    }
                }
                for (p, vals) in buckets.iter().enumerate() {
                    if !vals.is_empty() {
                        write_all_at(
                            &files[p].1,
                            vec_to_slice_u8(vals),
                            (sub_offsets[s][p] * int_size) as u64,
                        )?;
                    }
                }
                Ok(())
            })?;
        let paths: Vec<PathBuf> = files.into_iter().map(|(path, _)| path).collect();
        info!(
            "Wrote {} unsorted suffixes to {num_partitions} partition{} in {:?}",
            counts.iter().sum::<usize>(),
            if num_partitions == 1 { "" } else { "s" },
            now.elapsed()
        );

        // Sort partitions in memory, several at a time so that reading,
        // keying, sorting and writing overlap. Each partition is written
        // straight into the output file at its offset; the LCP at each
        // partition boundary is patched afterwards.
        let concurrency = std::env::var("SUFR_RADIX_CONCURRENCY")
            .ok()
            .and_then(|v| v.parse::<usize>().ok())
            .unwrap_or(2)
            .clamp(1, num_partitions);
        info!(
            "Sorting {concurrency} partition{} at a time",
            if concurrency == 1 { "" } else { "s" }
        );
        let mut offsets = Vec::with_capacity(num_partitions);
        let mut total = 0;
        for &count in &counts {
            offsets.push(total);
            total += count;
        }
        let next = std::sync::atomic::AtomicUsize::new(0);
        let ends: PartitionEnds<T> = Mutex::new(vec![None; num_partitions]);
        std::thread::scope(|scope| -> Result<()> {
            let workers: Vec<_> = (0..concurrency)
                .map(|_| {
                    scope.spawn(|| -> Result<()> {
                        // Buffers reused across this worker's partitions so
                        // the pages are faulted in once
                        let mut positions: Vec<T> = Vec::new();
                        let mut items: Vec<Item<T>> = Vec::new();
                        let mut scratch: Vec<Item<T>> = Vec::new();
                        loop {
                            let p = next
                                .fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                            if p >= num_partitions {
                                return Ok(());
                            }
                            let len = counts[p];
                            let path = &paths[p];
                            if len == 0 {
                                fs::remove_file(path)?;
                                continue;
                            }
                            let now = Instant::now();
                            positions.resize(len, T::default());
                            read_all_at(
                                &File::open(path)?,
                                slice_int_to_slice_u8(&mut positions),
                                0,
                            )?;
                            fs::remove_file(path)?;
                            positions
                                .par_iter()
                                .map(|&pos| Item {
                                    key: encoder.key(&self.text, pos.to_usize()),
                                    pos,
                                })
                                .collect_into_vec(&mut items);
                            info!(
                                "Partition {p}: read and keyed {len} suffixes in {:?}",
                                now.elapsed()
                            );

                            let now = Instant::now();
                            sort_items(&mut items, &mut scratch, encoder.total_bits);
                            info!(
                                "Partition {p}: sorted {len} suffixes in {:?}",
                                now.elapsed()
                            );

                            self.emit_radix_partition(
                                &out,
                                &encoder,
                                &items,
                                offsets[p],
                                None,
                            )?;
                            let first_last = (items[0], items[len - 1]);
                            match ends.lock() {
                                Ok(mut ends) => ends[p] = Some(first_last),
                                Err(e) => bail!("{e}"),
                            }
                        }
                    })
                })
                .collect();
            for worker in workers {
                worker
                    .join()
                    .map_err(|_| anyhow!("A partition sorting thread panicked"))??;
            }
            Ok(())
        })?;

        // Fix the LCP at each partition boundary
        if self.write_lcp {
            let ends = ends.into_inner().map_err(|e| anyhow!("{e}"))?;
            let mut prev: Option<Item<T>> = None;
            for (p, entry) in ends.iter().enumerate() {
                if let Some((first, last)) = entry {
                    if let Some(prev) = prev {
                        let lcp = T::from_usize(encoder.lcp(
                            prev.key,
                            prev.pos.to_usize(),
                            first.key,
                            first.pos.to_usize(),
                        ));
                        write_all_at(
                            &out.file,
                            vec_to_slice_u8(std::slice::from_ref(&lcp)),
                            (out.lcp_pos + offsets[p] * int_size) as u64,
                        )?;
                    }
                    prev = Some(*last);
                }
            }
        }
        self.finish_output(&out, num_suffixes)?;

        self.num_suffixes = T::from_usize(num_suffixes);
        info!(
            "Sorted and wrote {num_suffixes} suffixes in {num_partitions} partitions (partitioned radix) in {:?}",
            total_time.elapsed()
        );
        Ok(())
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
    /// * `lcp_cache`: optional LcpCache to store/retrieve from, for optimization
    #[inline(always)]
    fn find_lcp(
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
                self.find_lcp(
                    start1.to_usize(),
                    start2.to_usize(),
                    max_query_len,
                    0,
                    None,
                )
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
    fn merge(
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
    fn write(&mut self) -> Result<usize> {
        let now = Instant::now();
        let mut out = self.open_output()?;
        out.set_num_suffixes(self.num_suffixes.to_usize());
        let int_size = mem::size_of::<T>();

        // Offset (in elements) of each partition in the arrays, and the
        // last suffix of the preceding partition for fixing the LCP boundary
        let mut offsets = Vec::with_capacity(self.partitions.len());
        let mut total = 0;
        let mut prev_last_suffix = None;
        for partition in &self.partitions {
            offsets.push((total, prev_last_suffix));
            total += partition.len;
            prev_last_suffix = Some(partition.last_suffix);
        }

        // NB: Taking the partitions and consuming their buffers frees
        // memory/disk space along the way
        let partitions = mem::take(&mut self.partitions);
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
                                self.text_len,
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

        let bytes_out = self.finish_output(&out, self.num_suffixes.to_usize())?;
        info!(
            "Wrote suffix array to '{}' in {:?}",
            self.path,
            now.elapsed()
        );
        Ok(bytes_out)
    }

    // --------------------------------------------------
    /// Create the output file and write the header and text. The suffix
    /// (and optional LCP) arrays are written afterwards at the positions
    /// recorded in the returned `Output`, once `set_num_suffixes` has
    /// fixed the layout.
    fn open_output(&self) -> Result<Output> {
        let now = Instant::now();
        let filename = &self.path;
        let file = File::create(filename).map_err(|e| anyhow!("{filename}: {e}"))?;

        // Various metadata
        let is_dna: u8 = if self.is_dna { 1 } else { 0 };
        let allow_ambiguity: u8 = if self.allow_ambiguity { 1 } else { 0 };
        let ignore_softmask: u8 = if self.ignore_softmask { 1 } else { 0 };
        let mut header: Vec<u8> =
            vec![OUTFILE_VERSION, is_dna, allow_ambiguity, ignore_softmask];

        // Text length
        header.extend(self.text_len.to_usize().to_le_bytes());

        // Locations of text, suffix array, and LCP; filled in by finish_output
        let locs_pos = header.len() as u64;
        header.extend(0usize.to_le_bytes());
        header.extend(0usize.to_le_bytes());
        header.extend(0usize.to_le_bytes());

        // Number of suffixes; filled in by finish_output
        let num_suffixes_pos = header.len() as u64;
        header.extend(0usize.to_le_bytes());

        // Max query length
        let max_query_len = if let SuffixSortType::MaxQueryLen(val) = &self.sort_type {
            *val
        } else {
            0
        };
        header.extend(max_query_len.to_le_bytes());

        // Number of sequences
        header.extend(self.sequence_starts.len().to_le_bytes());

        // Sequence starts
        header.extend_from_slice(vec_to_slice_u8(&self.sequence_starts));

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

        let sa_pos = text_pos + self.text.len();
        info!("Wrote header and text in {:?}", now.elapsed());

        Ok(Output {
            file,
            locs_pos,
            num_suffixes_pos,
            text_pos,
            sa_pos,
            lcp_pos: 0,
            end_pos: sa_pos,
            int_size: mem::size_of::<T>(),
            write_lcp: self.write_lcp,
        })
    }

    // --------------------------------------------------
    /// Write the sequence names after the arrays and record the number of
    /// suffixes and the array positions in the header. Returns the total
    /// number of bytes in the file.
    fn finish_output(&self, out: &Output, num_suffixes: usize) -> Result<usize> {
        // Sequence names are variable in length so they are at the end
        let names = bincode::serialize(&self.sequence_names)?;
        write_all_at(&out.file, &names, out.end_pos as u64)?;

        // Go back to header and record the number of suffixes and the
        // locations
        write_all_at(&out.file, &num_suffixes.to_le_bytes(), out.num_suffixes_pos)?;
        let mut locs = out.text_pos.to_le_bytes().to_vec();
        locs.extend(out.sa_pos.to_le_bytes());
        locs.extend(out.lcp_pos.to_le_bytes());
        write_all_at(&out.file, &locs, out.locs_pos)?;

        Ok(out.end_pos + names.len())
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
    use super::{SortStrategy, SufrBuilder, SufrBuilderArgs};
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
            sort_strategy: SortStrategy::Merge,
            write_lcp: true,
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
            sort_strategy: SortStrategy::Merge,
            write_lcp: true,
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
            sort_strategy: SortStrategy::Merge,
            write_lcp: true,
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
            sort_strategy: SortStrategy::Merge,
            write_lcp: true,
        };
        let sufr: SufrBuilder<u32> = SufrBuilder::new(args)?;

        // 0: TTTAGC
        // 1:  TTAGC
        // 6: len of text
        assert_eq!(sufr.find_lcp(0, 1, 6, 0, None), 2);

        // 0: TTTAGC
        // 2:   TAGC
        // 6: len of text
        assert_eq!(sufr.find_lcp(0, 2, 6, 0, None), 1);

        // 0: TTTAGC
        // 1:  TTAGC
        // 1: max query len = 1
        assert_eq!(sufr.find_lcp(0, 1, 1, 0, None), 1);

        // 0: TTTAGC
        // 3:    AGC
        // 6: len of text
        assert_eq!(sufr.find_lcp(0, 3, 6, 0, None), 0);

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
            sort_strategy: SortStrategy::Merge,
            write_lcp: true,
        };
        let sufr: SufrBuilder<u32> = SufrBuilder::new(args)?;

        // 0: TTTTTA
        // 1:  TTTTA
        assert_eq!(sufr.find_lcp(0, 1, 3, 0, None), 3);

        // 0: TTTTTA
        // 2:   TTTA
        assert_eq!(sufr.find_lcp(0, 2, 3, 0, None), 2);

        // 0: TTTTTA
        // 5:      A
        assert_eq!(sufr.find_lcp(0, 5, 3, 0, None), 0);

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
            sort_strategy: SortStrategy::Merge,
            write_lcp: true,
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
            sort_strategy: SortStrategy::Merge,
            write_lcp: true,
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
            sort_strategy: SortStrategy::Merge,
            write_lcp: true,
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
