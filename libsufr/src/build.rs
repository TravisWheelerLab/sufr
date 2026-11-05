use crate::{
    format::{Header, SeedMask, SuffixSortType, OUTFILE_VERSION},
    io::{as_bytes, write_all_at, Int},
    lcp::{mask_offset, masked_lcp, suffix_lcp, LcpCache},
};
use anyhow::{anyhow, Result};
use log::info;
use needletail::parse_fastx_file;
use rayon::prelude::*;
use std::{
    cell::RefCell,
    cmp::{max, min, Ordering},
    fs::File,
    mem,
    ops::Range,
    path::Path,
    sync::Mutex,
    time::Instant,
};
use thread_local::ThreadLocal;

/// The character placed at the end of the text, which must not occur in
/// the text.
const SENTINEL_CHARACTER: u8 = b'$';

/// The kind of text a suffix array is built over.
pub enum Alphabet {
    /// Nucleotides.
    Dna {
        /// Whether to index suffixes starting with a byte other than `ACGT`.
        allow_ambiguity: bool,

        /// Whether to replace lowercase bytes with `N` rather than
        /// upper-case them.
        ignore_softmask: bool,
    },

    /// Amino acids.
    Amino {
        /// Whether to index suffixes starting with an unknown residue.
        index_unknown: bool,
    },

    /// Any bytes.
    Bytes,
}

/// Configuration and construction of a _.sufr_ file.
///
/// Start from a sequence file with [`SuffixArrayBuilder::from_sequence_file`] or
/// from text in memory with [`SuffixArrayBuilder::from_data`], set options, then
/// write the file with [`SuffixArrayBuilder::build_to_file`].
pub struct SuffixArrayBuilder {
    text: Vec<u8>,
    sequences: Vec<(usize, String)>,
    alphabet: Alphabet,

    /// The seed mask the suffixes are sorted by, if any.
    seed_mask: Option<SeedMask>,

    /// The number of characters the suffixes are sorted by, or 0 for all.
    max_query_len: usize,

    num_partitions: usize,

    /// The most bytes of suffixes to hold in memory while sorting, or
    /// `None` for the length of the text.
    memory_budget: Option<usize>,
    force_64bit: bool,
    write_lcp: bool,
}

impl SuffixArrayBuilder {
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
        Ok(Self::from_data(data.seq, sequences))
    }

    /// A builder for `text`, whose `sequences` are given by start position
    /// and name, with a `%` between each two.
    pub fn from_data(text: Vec<u8>, sequences: Vec<(usize, String)>) -> Self {
        // TODO: validate text and sequences?
        // - sequences should be sorted, first starting at 0, length < text.len(), no duplicates. If empty, default?
        Self {
            text,
            sequences,
            alphabet: Alphabet::Bytes,
            seed_mask: None,
            max_query_len: 0,
            num_partitions: 1024,
            memory_budget: None,
            force_64bit: false,
            write_lcp: false,
        }
    }

    /// Set the alphabet of the text (default: `Alphabet::Bytes`).
    pub fn with_alphabet(mut self, alphabet: Alphabet) -> Self {
        self.alphabet = alphabet;
        self
    }

    /// Set how the suffixes are sorted (default: `SuffixSortType::MaxQueryLen(0)`,
    /// by all of their characters).
    pub fn with_sort_type(mut self, sort_type: SuffixSortType) -> Self {
        (self.seed_mask, self.max_query_len) = match sort_type {
            SuffixSortType::Mask(seed_mask) => (Some(seed_mask), 0),
            SuffixSortType::MaxQueryLen(max_query_len) => (None, max_query_len),
        };
        self
    }

    /// Sort the suffixes by their first `max_query_len` characters.
    pub fn with_max_query_len(self, max_query_len: usize) -> Self {
        self.with_sort_type(SuffixSortType::MaxQueryLen(max_query_len))
    }

    /// Sort the suffixes by the care positions of `mask`.
    pub fn with_seed_mask(self, mask: SeedMask) -> Self {
        self.with_sort_type(SuffixSortType::Mask(mask))
    }

    /// Set the number of partitions the suffixes are sorted in
    /// (default: 1024).
    pub fn with_num_partitions(mut self, num_partitions: usize) -> Self {
        self.num_partitions = num_partitions;
        self
    }

    // only TypedSuffixArray::<u64> can read a small text written this
    // way, so this is for tests
    #[cfg(test)]
    pub(crate) fn with_force_64bit(mut self, force_64bit: bool) -> Self {
        self.force_64bit = force_64bit;
        self
    }

    /// Set the most bytes of suffixes to hold in memory while sorting
    /// (default: the length of the text); a text whose suffixes need
    /// more is scanned once per that many.
    pub fn with_memory_budget(mut self, bytes: usize) -> Self {
        self.memory_budget = Some(bytes);
        self
    }

    /// Set whether to write the LCP array (default: `false`).
    pub fn with_write_lcp(mut self, write_lcp: bool) -> Self {
        self.write_lcp = write_lcp;
        self
    }

    /// Build the suffix array and write it to `path`.
    pub fn build_to_file(self, path: impl AsRef<Path>) -> Result<()> {
        let path = path.as_ref();
        let fits_u32 = (self.text.len() as u64) < u32::MAX as u64 && !self.force_64bit;
        if fits_u32 {
            self.construct::<u32>(path)
        } else {
            self.construct::<u64>(path)
        }
    }

    /// Build the suffix array of the text with `W` per suffix and LCP and
    /// write it to `path`.
    fn construct<W: Int>(mut self, path: &Path) -> Result<()> {
        let mut text = mem::take(&mut self.text);
        let occupancy = normalize(&mut text, &self.alphabet);
        let text = text.as_slice();
        let keep = keep_table(&self.alphabet);
        let n_ranges = match self.alphabet {
            Alphabet::Dna {
                allow_ambiguity: true,
                ..
            } => find_n_ranges(text),
            _ => vec![],
        };
        let encoder = PrefixEncoder::from_occupancy(occupancy, self.max_query_len);
        let mask = self.seed_mask.as_ref();
        let max_query_len = self.max_query_len;

        // ---
        // count the suffixes of each bucket, those whose first byte
        // `keep` marks, and cut the buckets into partitions

        let now = Instant::now();
        let num_buckets = encoder.num_buckets();
        let counts = (0..text.len().div_ceil(CHUNK_SIZE))
            .into_par_iter()
            .fold(
                || vec![0; num_buckets],
                |mut counts, chunk| {
                    let start = chunk * CHUNK_SIZE;
                    let end = (start + CHUNK_SIZE).min(text.len());
                    match mask {
                        Some(mask) => {
                            for i in start..end {
                                if keep[text[i] as usize] {
                                    let key = encoder.encode(text, i, Some(mask));
                                    counts[encoder.bucket(encoder.capped(key))] += 1;
                                }
                            }
                        }
                        None => {
                            let mut key = encoder.encode(text, start, None);
                            for i in start..end {
                                if keep[text[i] as usize] {
                                    counts[encoder.bucket(encoder.capped(key))] += 1;
                                }
                                key = encoder.roll(text, i, key);
                            }
                        }
                    }
                    counts
                },
            )
            .reduce(
                || vec![0; num_buckets],
                |mut total, counts| {
                    for (total, count) in total.iter_mut().zip(counts) {
                        *total += count;
                    }
                    total
                },
            );
        let (ranges, buckets) = cut_partitions(&counts, self.num_partitions);
        let mut partition_of = vec![0; num_buckets];
        for (p, range) in buckets.iter().enumerate() {
            partition_of[range.clone()].fill(p);
        }
        let num_suffixes = ranges.last().map_or(0, |range| range.end);
        info!(
            "Counted {num_suffixes} suffixes into {num_buckets} buckets and {} partitions in {:?}",
            ranges.len(),
            now.elapsed()
        );

        // ---
        // the header, the file and the text

        let (is_dna, allow_ambiguity, ignore_softmask) = match self.alphabet {
            Alphabet::Dna {
                allow_ambiguity,
                ignore_softmask,
            } => (true, allow_ambiguity, ignore_softmask),
            Alphabet::Amino { .. } | Alphabet::Bytes => (false, false, false),
        };
        let mut header = Header {
            version: OUTFILE_VERSION,
            is_dna,
            allow_ambiguity,
            ignore_softmask,
            text_len: text.len(),
            num_suffixes,
            sort_type: match mask {
                Some(seed_mask) => SuffixSortType::Mask(seed_mask.clone()),
                None => SuffixSortType::MaxQueryLen(max_query_len),
            },
            sequence_starts: self.sequences.iter().map(|(start, _)| *start).collect(),
            sequence_names: self
                .sequences
                .iter()
                .map(|(_, name)| name.clone())
                .collect(),
            int_size: mem::size_of::<W>(),
            text_pos: 0,
            suffix_array_pos: 0,
            lcp_pos: None,
        };
        header.set_positions(self.write_lcp);

        let now = Instant::now();
        let file =
            File::create(path).map_err(|e| anyhow!("{}: {e}", path.display()))?;
        write_all_at(&file, text, header.text_pos as u64)?;

        // ---
        // each turn scans the text for the suffixes of as many
        // partitions as the budget holds, and sorts them from memory

        let per_suffix = mem::size_of::<W>() + mem::size_of::<u64>();
        let budget = self.memory_budget.unwrap_or(text.len());
        let mut first_suffixes = vec![];
        let mut last_suffixes = vec![];
        for turn in turns(&ranges, budget / per_suffix) {
            // gather the suffixes of this turn's partitions with
            // their keys, from one scan of the text
            let now = Instant::now();
            let first = turn.start;
            let count = turn.len();
            // one compare for "in this turn's buckets", so the branch is
            // rarely taken whichever turn this is; as lo <= b && b < hi
            // the first test is a coin flip in the middle turns and the
            // scan ran three times slower there
            let lo = buckets[first].start;
            let width = buckets[turn.end - 1].end - lo;
            let gathered: Vec<Mutex<(Vec<W>, Vec<u64>)>> = turn
                .map(|p| {
                    let len = ranges[p].len();
                    Mutex::new((Vec::with_capacity(len), Vec::with_capacity(len)))
                })
                .collect();

            // the buffers are sized once per thread and filled by index,
            // so the loop holds no call and keeps its state in registers
            (0..text.len().div_ceil(CHUNK_SIZE))
                .into_par_iter()
                .try_for_each_init(
                    || (vec![W::default(); CHUNK_SIZE], vec![0u64; CHUNK_SIZE]),
                    |(suffixes, keys), chunk| -> Result<()> {
                        let start = chunk * CHUNK_SIZE;
                        let end = (start + CHUNK_SIZE).min(text.len());
                        let mut n = 0;
                        match mask {
                            Some(mask) => {
                                for i in start..end {
                                    let key = encoder.capped(encoder.encode(
                                        text,
                                        i,
                                        Some(mask),
                                    ));
                                    let bucket = encoder.bucket(key);
                                    if keep[text[i] as usize]
                                        && bucket.wrapping_sub(lo) < width
                                    {
                                        suffixes[n] = W::from_usize(i);
                                        keys[n] = key;
                                        n += 1;
                                    }
                                }
                            }
                            None => {
                                let mut key = encoder.encode(text, start, None);
                                for i in start..end {
                                    let capped = encoder.capped(key);
                                    let bucket = encoder.bucket(capped);
                                    if keep[text[i] as usize]
                                        && bucket.wrapping_sub(lo) < width
                                    {
                                        suffixes[n] = W::from_usize(i);
                                        keys[n] = capped;
                                        n += 1;
                                    }
                                    key = encoder.roll(text, i, key);
                                }
                            }
                        }
                        let suffixes = &suffixes[..n];
                        let keys = &keys[..n];

                        let partition =
                            |i: usize| partition_of[encoder.bucket(keys[i])] - first;
                        let (suffixes, starts) =
                            counting_sort(suffixes, count, partition);
                        let (keys, _) = counting_sort(keys, count, partition);
                        for (p, run) in starts.windows(2).enumerate() {
                            let run = run[0]..run[1];
                            if !run.is_empty() {
                                let mut g =
                                    gathered[p].lock().map_err(|e| anyhow!("{e}"))?;
                                g.0.extend_from_slice(&suffixes[run.clone()]);
                                g.1.extend_from_slice(&keys[run]);
                            }
                        }
                        Ok(())
                    },
                )?;

            let gathered = gathered
                .into_iter()
                .map(|g| g.into_inner().map_err(|e| anyhow!("{e}")))
                .collect::<Result<Vec<_>>>()?;
            info!(
                "Gathered {} suffixes of partitions {first}..{} in {:?}",
                gathered.iter().map(|g| g.0.len()).sum::<usize>(),
                first + count,
                now.elapsed()
            );

            // sort each partition by key, then merge sort only the
            // ties, and write it with its LCPs
            let now = Instant::now();
            let lcp_cache: ThreadLocal<RefCell<LcpCache>> = ThreadLocal::new();
            let sorted = gathered
                .into_par_iter()
                .enumerate()
                .map(|(i, (suffixes, keys))| -> Result<(usize, usize)> {
                    let range = &ranges[first + i];
                    let len = suffixes.len();

                    // by key is by the first key_chars characters, so only a
                    // run of equal keys is left for the merge sort, and the
                    // LCP across unequal keys is read off the keys
                    let mut keyed: Vec<(u64, W)> =
                        keys.into_iter().zip(suffixes).collect();
                    keyed.sort_unstable_by_key(|&(key, _)| key);
                    let mut sa: Vec<W> = keyed.iter().map(|&(_, s)| s).collect();
                    let mut lcp = vec![W::default(); len];
                    let mut run = 0;
                    for i in 1..=len {
                        if i < len && keyed[i].0 == keyed[run].0 {
                            continue;
                        }
                        if i - run > 1 {
                            merge_sort(
                                text,
                                &mut sa[run..i].to_vec(),
                                &mut lcp[run..i].to_vec(),
                                &mut sa[run..i],
                                &mut lcp[run..i],
                                mask,
                                max_query_len,
                                &n_ranges,
                                &lcp_cache,
                            );
                        }
                        if run > 0 {
                            lcp[run] = W::from_usize(
                                encoder.shared(keyed[run - 1].0, keyed[run].0),
                            );
                        }
                        run = i;
                    }
                    drop(keyed);

                    let offset = range.start * mem::size_of::<W>();
                    write_all_at(
                        &file,
                        as_bytes(&sa),
                        (header.suffix_array_pos + offset) as u64,
                    )?;
                    if let Some(lcp_pos) = header.lcp_pos {
                        write_all_at(&file, as_bytes(&lcp), (lcp_pos + offset) as u64)?;
                    }
                    Ok((sa[0].to_usize(), sa[len - 1].to_usize()))
                })
                .collect::<Result<Vec<_>>>()?;

            info!(
                "Sorted partitions {first}..{} in {:?}",
                first + count,
                now.elapsed()
            );

            for (first_suffix, last_suffix) in sorted {
                first_suffixes.push(first_suffix);
                last_suffixes.push(last_suffix);
            }
        }

        // ---
        // the LCP at each partition's start is against the partition
        // before, which is known once both are sorted

        if let Some(lcp_pos) = header.lcp_pos {
            let first_lcps: Vec<W> = (1..first_suffixes.len())
                .into_par_iter()
                .map(|i| {
                    find_lcp(
                        text,
                        last_suffixes[i - 1],
                        first_suffixes[i],
                        W::from_usize(text.len()),
                        0,
                        mask,
                        max_query_len,
                        &n_ranges,
                        None,
                    )
                })
                .collect();
            for (range, first_lcp) in ranges.iter().skip(1).zip(first_lcps) {
                let at = lcp_pos + range.start * mem::size_of::<W>();
                write_all_at(&file, as_bytes(&[first_lcp]), at as u64)?;
            }
        }

        header.write(&file)?;
        info!(
            "Wrote suffix array to '{}' in {:?}",
            path.display(),
            now.elapsed()
        );
        Ok(())
    }
}

/// Upper-case soft-masked characters, or replace them with `N` when the
/// alphabet is `Dna` with `ignore_softmask`, and return which bytes occur.
fn normalize(text: &mut [u8], alphabet: &Alphabet) -> [bool; 256] {
    let ignore_softmask = matches!(
        alphabet,
        Alphabet::Dna {
            ignore_softmask: true,
            ..
        }
    );
    let mut occupancy = [false; 256];
    for b in text.iter_mut() {
        if b.is_ascii_lowercase() {
            if ignore_softmask {
                *b = b'N'
            } else {
                *b = b.to_ascii_uppercase()
            }
        }
        occupancy[*b as usize] = true;
    }
    occupancy
}

/// The first bytes of the suffixes that are indexed: `ACGT` and the
/// sentinel for `Dna` without `allow_ambiguity`, else every byte.
fn keep_table(alphabet: &Alphabet) -> [bool; 256] {
    match alphabet {
        Alphabet::Dna {
            allow_ambiguity: false,
            ..
        } => {
            let mut keep = [false; 256];
            for &c in b"ACGT" {
                keep[c as usize] = true;
            }
            keep[SENTINEL_CHARACTER as usize] = true;
            keep
        }
        // TODO: index_unknown is recorded but not acted on; dropping
        // `*`, `-` and `X` is out of scope for now
        Alphabet::Dna { .. } | Alphabet::Amino { .. } | Alphabet::Bytes => [true; 256],
    }
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

/// The suffixes one task of a pass over the text takes.
const CHUNK_SIZE: usize = 1 << 19;

/// The rank of each byte of the text in byte order, packed into a key:
/// the ranks of a suffix's first `key_chars` characters at `bits` each
/// from the top of a `u64`, so keys compare as the characters do and
/// two keys share `leading_zeros(a ^ b) / bits` characters.
pub(super) struct PrefixEncoder {
    table: [u8; 256],

    /// The bits per rank.
    bits: u32,

    /// The characters in a key.
    key_chars: usize,

    /// The characters in a bucket number, the key's top ones.
    bucket_chars: usize,

    /// The unused low bits of a key.
    shift: u32,

    /// Keeps the first `max_query_len` characters of a key, or all.
    cap: u64,
}

impl PrefixEncoder {
    /// The encoder for a text whose bytes `occupancy` marks, sorted by
    /// its first `max_query_len` characters, or all of them for 0.
    pub(super) fn from_occupancy(occupancy: [bool; 256], max_query_len: usize) -> Self {
        let mut table = [0u8; 256];
        let mut radix: usize = 0;
        for (c, &seen) in occupancy.iter().enumerate() {
            if seen {
                table[c] = radix as u8;
                radix += 1;
            }
        }
        let bits = radix.next_power_of_two().trailing_zeros().max(1);
        let key_chars = 64 / bits as usize;
        let bucket_chars = 16 / bits as usize;
        let shift = 64 - bits * key_chars as u32;
        let cap = if max_query_len > 0 && max_query_len < key_chars {
            u64::MAX << (64 - bits * max_query_len as u32)
        } else {
            u64::MAX
        };
        info!(
            "{radix} distinct bytes at {bits} bits: {key_chars}-character keys, \
             {bucket_chars}-character buckets"
        );
        Self {
            table,
            bits,
            key_chars,
            bucket_chars,
            shift,
            cap,
        }
    }

    /// The number of distinct bucket numbers.
    pub(super) fn num_buckets(&self) -> usize {
        1 << (self.bits * self.bucket_chars as u32)
    }

    /// `key` with the characters past `max_query_len` cleared.
    #[inline(always)]
    pub(super) fn capped(&self, key: u64) -> u64 {
        key & self.cap
    }

    /// The bucket number of a suffix with the capped key `key`.
    #[inline(always)]
    pub(super) fn bucket(&self, key: u64) -> usize {
        (key >> (64 - self.bits * self.bucket_chars as u32)) as usize
    }

    /// The characters two suffixes with the unequal keys `a` and `b`
    /// share.
    #[inline(always)]
    pub(super) fn shared(&self, a: u64, b: u64) -> usize {
        ((a ^ b).leading_zeros() / self.bits) as usize
    }

    /// The key of the suffix at `i` before the cap: its first
    /// `key_chars` characters, or those at the care positions of `mask`.
    pub(super) fn encode(&self, text: &[u8], i: usize, mask: Option<&SeedMask>) -> u64 {
        let mut key = 0u64;
        for j in 0..self.key_chars {
            // a character past the end of the text or the mask
            // ranks 0, so a short suffix keeps every digit
            let c = match mask {
                Some(mask) => mask.positions.get(j).and_then(|&off| text.get(i + off)),
                None => text.get(i + j),
            };
            let rank = c.map_or(0, |&c| self.table[c as usize] as u64);
            key = (key << self.bits) | rank;
        }
        key << self.shift
    }

    /// The key of the suffix after the one at `i` with the key `key`
    /// before the cap, without a mask: the first character shifts out
    /// of the top and the character `key_chars` ahead comes in at the
    /// bottom.
    #[inline(always)]
    pub(super) fn roll(&self, text: &[u8], i: usize, key: u64) -> u64 {
        let incoming = text
            .get(i + self.key_chars)
            .map_or(0, |&c| self.table[c as usize] as u64);
        (key << self.bits) | (incoming << self.shift)
    }
}

/// The suffixes of each partition as a range of the suffix array, and
/// its buckets: the buckets cut into at most `num_partitions` runs of
/// about the same number of suffixes.
fn cut_partitions(
    counts: &[usize],
    num_partitions: usize,
) -> (Vec<Range<usize>>, Vec<Range<usize>>) {
    let num_suffixes: usize = counts.iter().sum();
    let share = num_suffixes.div_ceil(num_partitions.max(1));
    let mut buckets = vec![];
    let mut ranges: Vec<Range<usize>> = vec![];
    let mut first_bucket = 0;
    let mut start = 0;
    let mut taken = 0;
    for (bucket, &count) in counts.iter().enumerate() {
        taken += count;

        // a partition closes once the suffixes taken so far
        // reach its cumulative share, so a bucket far over
        // its share is paid back by the partitions after it
        if taken > start && taken >= share * (ranges.len() + 1) {
            buckets.push(first_bucket..bucket + 1);
            ranges.push(start..taken);
            first_bucket = bucket + 1;
            start = taken;
        }
    }
    if taken > start {
        buckets.push(first_bucket..counts.len());
        ranges.push(start..taken);
    }
    (ranges, buckets)
}

/// Runs of partitions of at most `max_suffixes` suffixes each, or
/// one partition where that alone is more.
fn turns(ranges: &[Range<usize>], max_suffixes: usize) -> Vec<Range<usize>> {
    let mut turns = vec![];
    let mut first = 0;
    let mut taken = 0;
    for (p, range) in ranges.iter().enumerate() {
        if taken > 0 && taken + range.len() > max_suffixes {
            turns.push(first..p);
            first = p;
            taken = 0;
        }
        taken += range.len();
    }
    if first < ranges.len() {
        turns.push(first..ranges.len());
    }
    turns
}

/// `values` in the order of `key(i)` for each index `i`, every key
/// below `num_keys`, and where each key's run starts, with one more
/// entry for the end.
pub(super) fn counting_sort<T: Copy + Default>(
    values: &[T],
    num_keys: usize,
    key: impl Fn(usize) -> usize,
) -> (Vec<T>, Vec<usize>) {
    let mut starts = vec![0; num_keys + 1];
    for i in 0..values.len() {
        starts[key(i) + 1] += 1;
    }
    for k in 0..num_keys {
        starts[k + 1] += starts[k];
    }
    let mut next = starts.clone();
    let mut ordered = vec![T::default(); values.len()];
    for (i, &value) in values.iter().enumerate() {
        let slot = &mut next[key(i)];
        ordered[*slot] = value;
        *slot += 1;
    }
    (ordered, starts)
}

/// Sort the suffixes in `out` with their LCPs, using `work`, which
/// starts with the same suffixes, as scratch space.
#[allow(clippy::too_many_arguments)]
fn merge_sort<W: Int>(
    text: &[u8],
    work_sa: &mut [W],
    work_lcp: &mut [W],
    out_sa: &mut [W],
    out_lcp: &mut [W],
    mask: Option<&SeedMask>,
    max_query_len: usize,
    n_ranges: &[Range<usize>],
    lcp_cache: &ThreadLocal<RefCell<LcpCache>>,
) {
    let n = work_sa.len();
    if n == 1 {
        out_lcp[0] = W::default();
    } else {
        const NESTED_PAR_GRAIN_SIZE: usize = 1 << 13;
        let mid = n / 2;
        let (work_sa_l, work_sa_r) = work_sa.split_at_mut(mid);
        let (work_lcp_l, work_lcp_r) = work_lcp.split_at_mut(mid);
        let (out_sa_l, out_sa_r) = out_sa.split_at_mut(mid);
        let (out_lcp_l, out_lcp_r) = out_lcp.split_at_mut(mid);

        // each half is sorted into `work`, swapping the
        // roles, so the merge below reads `work`
        if mid < NESTED_PAR_GRAIN_SIZE {
            merge_sort(
                text,
                out_sa_l,
                out_lcp_l,
                work_sa_l,
                work_lcp_l,
                mask,
                max_query_len,
                n_ranges,
                lcp_cache,
            );
            merge_sort(
                text,
                out_sa_r,
                out_lcp_r,
                work_sa_r,
                work_lcp_r,
                mask,
                max_query_len,
                n_ranges,
                lcp_cache,
            );
        } else {
            rayon::join(
                || {
                    merge_sort(
                        text,
                        out_sa_l,
                        out_lcp_l,
                        work_sa_l,
                        work_lcp_l,
                        mask,
                        max_query_len,
                        n_ranges,
                        lcp_cache,
                    )
                },
                || {
                    merge_sort(
                        text,
                        out_sa_r,
                        out_lcp_r,
                        work_sa_r,
                        work_lcp_r,
                        mask,
                        max_query_len,
                        n_ranges,
                        lcp_cache,
                    )
                },
            );
        }

        merge(
            text,
            work_sa,
            work_lcp,
            mid,
            out_sa,
            out_lcp,
            mask,
            max_query_len,
            n_ranges,
            &mut lcp_cache.get_or_default().borrow_mut(),
        );
    }
}

/// Merge the sorted runs `from[..mid]` and `from[mid..]` into `into`.
#[allow(clippy::too_many_arguments)]
fn merge<W: Int>(
    text: &[u8],
    from_sa: &mut [W],
    from_lcp: &mut [W],
    mid: usize,
    target_sa: &mut [W],
    target_lcp: &mut [W],
    mask: Option<&SeedMask>,
    max_query_len: usize,
    n_ranges: &[Range<usize>],
    lcp_cache: &mut LcpCache,
) {
    let (mut x, mut y) = from_sa.split_at_mut(mid);
    let (mut lcp_x, mut lcp_y) = from_lcp.split_at_mut(mid);
    let mut len_x = x.len();
    let mut len_y = y.len();
    // the last LCP taken from x
    let mut m = W::default();
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
                let max_n = W::from_usize(text.len()) - shorter_suffix;

                let context = match mask {
                    Some(seed_mask) => W::from_usize(
                        seed_mask
                            .positions
                            .iter()
                            .filter(|&i| *i < max_n.to_usize())
                            .count(),
                    ),
                    None if max_query_len > 0 => {
                        min(W::from_usize(max_query_len), max_n)
                    }
                    None => max_n,
                };

                // LCP(X_i, Y_j)
                let (len_lcp, full_len_lcp) = if m < context {
                    let lcp = find_lcp(
                        text,
                        x[idx_x].to_usize(),
                        y[idx_y].to_usize(),
                        context - m,
                        m.to_usize(), // skip
                        mask,
                        max_query_len,
                        n_ranges,
                        Some(lcp_cache),
                    );
                    let full_lcp = match mask {
                        Some(seed_mask) => mask_offset(lcp.to_usize(), seed_mask),
                        None => lcp.to_usize(),
                    };
                    (lcp, W::from_usize(full_lcp))
                } else {
                    (context, context)
                };

                // an LCP reaching the context sorts the shorter suffix
                // first, and otherwise the byte after the LCP decides
                if len_lcp >= context {
                    target_sa[idx_target] = shorter_suffix;
                } else {
                    let cmp = text[(x[idx_x] + full_len_lcp).to_usize()]
                        .cmp(&text[(y[idx_y] + full_len_lcp).to_usize()]);

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

/// The LCP of the suffixes at `start1` and `start2`, which share their
/// first `skip` characters, comparing at most `len` characters.
#[allow(clippy::too_many_arguments)]
#[inline(always)]
fn find_lcp<W: Int>(
    text: &[u8],
    start1: usize,
    start2: usize,
    len: W,
    skip: usize,
    mask: Option<&SeedMask>,
    max_query_len: usize,
    n_ranges: &[Range<usize>],
    lcp_cache: Option<&mut LcpCache>,
) -> W {
    // with a mask, `len` is not consulted: the comparison
    // stops at the end of the mask or of the text
    if let Some(mask) = mask {
        return W::from_usize(
            skip + masked_lcp(
                &text[start1..],
                &text[start2..],
                &mask.positions[skip..],
            ),
        );
    }

    let mut skip = skip;
    if let (Some(end1), Some(end2)) =
        (find_n_run(n_ranges, start1), find_n_run(n_ranges, start2))
    {
        let (run1, run2) = (end1 - start1, end2 - start2);
        if run1 != run2 {
            // the byte after the shorter run is N against not N
            let run = min(run1, run2);
            return W::from_usize(if max_query_len > 0 {
                min(run, max_query_len)
            } else {
                run
            });
        }
        // equal runs end together, so compare on from there
        skip = max(skip, run1);
    }
    if max_query_len > 0 && skip >= max_query_len {
        return W::from_usize(max_query_len);
    }

    // `len` counts bytes past `skip`, so the cap leaves
    // max_query_len - skip of them
    let len = if max_query_len > 0 {
        max_query_len - skip
    } else {
        len.to_usize()
    };
    W::from_usize(suffix_lcp(text, start1, start2, skip, len, lcp_cache))
}

/// The end of the run of `N`s in `n_ranges` holding `suffix`, if it is
/// in one.
fn find_n_run(n_ranges: &[Range<usize>], suffix: usize) -> Option<usize> {
    n_ranges
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
        .map(|i| n_ranges[i].end)
}

#[cfg(test)]
mod test {
    use super::{
        counting_sort, cut_partitions, find_lcp, keep_table, read_sequence_file, turns,
        Alphabet, PrefixEncoder, SuffixArrayBuilder,
    };
    use crate::{format::SeedMask, SuffixArray};
    use anyhow::Result;
    use pretty_assertions::assert_eq;
    use std::path::Path;
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
    fn test_keep_table() {
        let kept = |alphabet| -> Vec<u8> {
            let keep = keep_table(&alphabet);
            (0..=255u8).filter(|&c| keep[c as usize]).collect()
        };
        let every_byte: Vec<u8> = (0..=255).collect();
        assert_eq!(
            kept(Alphabet::Dna {
                allow_ambiguity: false,
                ignore_softmask: false,
            }),
            b"$ACGT"
        );
        assert_eq!(
            kept(Alphabet::Dna {
                allow_ambiguity: true,
                ignore_softmask: false,
            }),
            every_byte
        );
        assert_eq!(
            kept(Alphabet::Amino {
                index_unknown: false
            }),
            every_byte
        );
        assert_eq!(
            kept(Alphabet::Amino {
                index_unknown: true
            }),
            every_byte
        );
        assert_eq!(kept(Alphabet::Bytes), every_byte);
    }

    #[test]
    fn test_write() -> Result<()> {
        let outfile = NamedTempFile::new()?;
        SuffixArrayBuilder::from_sequence_file("../data/inputs/3.fa", b'%')?
            .with_alphabet(Alphabet::Dna {
                allow_ambiguity: false,
                ignore_softmask: true,
            })
            .build_to_file(outfile.path())?;
        let suffix_array = SuffixArray::open(&outfile.path().to_string_lossy())?;
        let meta = suffix_array.metadata()?;
        assert_eq!(meta.text_len, 113);
        assert_eq!(meta.num_suffixes, 101);
        Ok(())
    }

    fn encoder(text: &[u8], max_query_len: usize) -> PrefixEncoder {
        let mut occupancy = [false; 256];
        for &b in text {
            occupancy[b as usize] = true;
        }
        PrefixEncoder::from_occupancy(occupancy, max_query_len)
    }

    /// `ranks` packed from the top of a u64 at `bits` each.
    fn packed(ranks: &[u64], bits: u32) -> u64 {
        ranks
            .iter()
            .enumerate()
            .map(|(j, &r)| r << (64 - bits * (j as u32 + 1)))
            .sum()
    }

    #[test]
    fn test_widths() {
        let dna = encoder(b"$%ACGNT", 0);
        assert_eq!((dna.bits, dna.key_chars, dna.bucket_chars), (3, 21, 5));
        assert_eq!(dna.num_buckets(), 32768);

        let protein = encoder(b"$%*-ABCDEFGHIKLMNPQRSTVWY", 0);
        assert_eq!(
            (protein.bits, protein.key_chars, protein.bucket_chars),
            (5, 12, 3)
        );

        let bytes: Vec<u8> = (0..=255).collect();
        let bytes = encoder(&bytes, 0);
        assert_eq!((bytes.bits, bytes.key_chars, bytes.bucket_chars), (8, 8, 2));
        assert_eq!(bytes.num_buckets(), 65536);

        let one = encoder(b"$", 0);
        assert_eq!((one.bits, one.key_chars, one.bucket_chars), (1, 64, 16));
    }

    #[test]
    fn test_encode() -> Result<()> {
        // ranks in byte order: $=0 A=1 C=2 G=3 N=4 T=5, 3 bits each
        let text = b"ACGTNA$";
        let encoder = encoder(text, 0);

        // A C G T N A $ and 14 characters past the end
        let key = encoder.encode(text, 0, None);
        assert_eq!(key, packed(&[1, 2, 3, 5, 4, 1, 0], 3));
        assert_eq!(
            encoder.bucket(key),
            packed(&[1, 2, 3, 5, 4], 3) as usize >> 49
        );

        // N A $ and the rest past the end
        assert_eq!(encoder.encode(text, 4, None), packed(&[4, 1, 0], 3));

        // care positions 0, 1 and 3, then past the mask: A C T
        let mask = SeedMask::new("1101")?;
        assert_eq!(encoder.encode(text, 0, Some(&mask)), packed(&[1, 2, 5], 3));
        // A $ and position 8 past the end
        assert_eq!(encoder.encode(text, 5, Some(&mask)), packed(&[1], 3));

        // a cap keeps that many characters
        let capped = super::test::encoder(text, 2);
        assert_eq!(
            capped.capped(capped.encode(text, 0, None)),
            packed(&[1, 2], 3)
        );
        Ok(())
    }

    #[test]
    fn test_roll_reaches_the_same_key_as_encode() {
        let text = b"ACGTNACGTTTGANNACGTACGTACGTAAA$";
        for max_query_len in [0, 4] {
            let encoder = encoder(text, max_query_len);
            let mut key = encoder.encode(text, 0, None);
            for i in 0..text.len() {
                assert_eq!(key, encoder.encode(text, i, None), "at {i}");
                key = encoder.roll(text, i, key);
            }
        }
    }

    #[test]
    fn test_shared_counts_equal_leading_characters() {
        let text = b"ACGTACGA$";
        let encoder = encoder(text, 0);
        let a = encoder.encode(text, 0, None);
        let b = encoder.encode(text, 4, None);
        // ACGTACGA$ and ACGA$ share ACG
        assert_eq!(encoder.shared(a, b), 3);
        // the lowest bit of a character deciding is still that character
        let c = packed(&[1, 2, 3, 4], 3);
        let d = packed(&[1, 2, 3, 5], 3);
        assert_eq!(encoder.shared(c, d), 3);
        assert_eq!(encoder.shared(a, encoder.encode(text, 8, None)), 0);
    }

    #[test]
    fn test_cut_at_cumulative_shares() {
        // 16 suffixes in 3 partitions: a share of 6, so the cuts
        // fall at 6 and 12
        let (ranges, buckets) = cut_partitions(&[3, 2, 1, 6, 2, 2], 3);
        assert_eq!(ranges, [0..6, 6..12, 12..16]);
        assert_eq!(buckets, [0..3, 3..4, 4..6]);

        // an empty bucket never opens or closes a partition
        let (ranges, buckets) = cut_partitions(&[0, 5, 0, 0, 5, 0], 2);
        assert_eq!(ranges, [0..5, 5..10]);
        assert_eq!(buckets, [0..2, 2..5]);

        // more partitions than suffixes
        let (ranges, _) = cut_partitions(&[1, 0, 1], 1024);
        assert_eq!(ranges, [0..1, 1..2]);
    }

    #[test]
    fn test_turns() {
        let (ranges, _) = cut_partitions(&[3, 2, 1, 6, 2, 2], 3);
        assert_eq!(turns(&ranges, 100), vec![0..3; 1]);
        assert_eq!(turns(&ranges, 6), [0..1, 1..2, 2..3]);
        assert_eq!(turns(&ranges, 10), [0..1, 1..3]);
        // a partition over the limit is a turn by itself
        assert_eq!(turns(&ranges, 1), [0..1, 1..2, 2..3]);
    }

    #[test]
    fn test_counting_sort() {
        let keys = [2, 0, 2, 1, 0];
        let (ordered, starts) = counting_sort(&[10u32, 11, 12, 13, 14], 4, |i| keys[i]);
        assert_eq!(ordered, [11, 14, 13, 10, 12]);
        assert_eq!(starts, [0, 2, 3, 5, 5]);
    }

    #[test]
    fn test_find_lcp_no_seed_mask() -> Result<()> {
        //           012345
        let text = b"TTTAGC";

        // 0: TTTAGC
        // 1:  TTAGC
        // 6: len of text
        assert_eq!(find_lcp::<u32>(text, 0, 1, 6, 0, None, 0, &[], None), 2);

        // 0: TTTAGC
        // 2:   TAGC
        // 6: len of text
        assert_eq!(find_lcp::<u32>(text, 0, 2, 6, 0, None, 0, &[], None), 1);

        // 0: TTTAGC
        // 1:  TTAGC
        // 1: max query len = 1
        assert_eq!(find_lcp::<u32>(text, 0, 1, 1, 0, None, 0, &[], None), 1);

        // 0: TTTAGC
        // 3:    AGC
        // 6: len of text
        assert_eq!(find_lcp::<u32>(text, 0, 3, 6, 0, None, 0, &[], None), 0);

        // 0: TTTTAG
        // 1:  TTTAG
        // skip 2 with max query len 2: nothing left to compare
        let text = b"TTTTAG";
        assert_eq!(find_lcp::<u32>(text, 0, 1, 6, 2, None, 2, &[], None), 2);

        // skip 1 with max query len 3: one byte left, W == W
        assert_eq!(find_lcp::<u32>(text, 0, 1, 6, 1, None, 3, &[], None), 3);

        // two runs of four N with the same remaining length: the
        // match goes on past them, GA against GAT
        let text = b"NNNNGA%NNNNGAT";
        let n_ranges = [0..4, 7..11];
        assert_eq!(
            find_lcp::<u32>(text, 0, 7, 14, 0, None, 0, &n_ranges, None),
            6
        );
        assert_eq!(
            find_lcp::<u32>(text, 0, 7, 14, 0, None, 5, &n_ranges, None),
            5
        );
        assert_eq!(
            find_lcp::<u32>(text, 0, 7, 14, 0, None, 3, &n_ranges, None),
            3
        );
        // runs of different remaining length end at the shorter
        assert_eq!(
            find_lcp::<u32>(text, 1, 7, 14, 0, None, 0, &n_ranges, None),
            3
        );

        // both suffixes inside a run of four N, max query len 2
        let text = b"NNNNA";
        // from_ref, since clippy reads the array literal [0..4]
        // as a typo for a range
        assert_eq!(
            find_lcp::<u32>(
                text,
                0,
                1,
                5,
                0,
                None,
                2,
                std::slice::from_ref(&(0..4)),
                None
            ),
            2
        );

        Ok(())
    }

    #[test]
    fn test_find_lcp_with_seed_mask() -> Result<()> {
        //           012345
        let text = b"TTTTTA";
        let mask = SeedMask::new("1101")?;
        let mask = Some(&mask);

        // 0: TTTTTA
        // 1:  TTTTA
        assert_eq!(find_lcp::<u32>(text, 0, 1, 3, 0, mask, 0, &[], None), 3);

        // 0: TTTTTA
        // 2:   TTTA
        assert_eq!(find_lcp::<u32>(text, 0, 2, 3, 0, mask, 0, &[], None), 2);

        // 0: TTTTTA
        // 5:      A
        assert_eq!(find_lcp::<u32>(text, 0, 5, 3, 0, mask, 0, &[], None), 0);

        Ok(())
    }
}
