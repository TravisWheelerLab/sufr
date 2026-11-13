use crate::{
    format::{Header, SeedMask, SuffixSortType, OUTFILE_VERSION},
    io::{as_bytes, write_all_at, Int},
    lcp::{masked_lcp, suffix_lcp, LcpCache},
};
use anyhow::{anyhow, bail, Result};
use log::info;
use needletail::parse_fastx_file;
use rayon::prelude::*;
use std::{
    cell::RefCell,
    cmp::{max, Ordering},
    fs::{self, File},
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
    /// Nucleotides, of which only suffixes starting with `A`, `C`, `G`
    /// or `T` are indexed.
    Dna {
        /// Whether to replace lowercase bytes with `N` rather than
        /// upper-case them.
        ignore_softmask: bool,
    },

    /// Amino acids.
    Amino,

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

    /// Sort the suffixes by the care positions of `mask` rather than by
    /// every character.
    pub fn with_seed_mask(mut self, mask: SeedMask) -> Self {
        self.seed_mask = Some(mask);
        self
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
        let packing = Packing::from_occupancy(occupancy);
        match &self.seed_mask {
            None => self.build::<W, _>(text, &keep, &UnmaskedKeys { packing }, path),
            Some(mask) => {
                let scheme =
                    MaskedKeys::new(packing, &mask.positions, mask.bytes.len())?;
                self.build::<W, _>(text, &keep, &scheme, path)
            }
        }
    }

    /// Build the suffix array of `text` with the keys of `scheme` and
    /// write it to `path`.
    fn build<W: Int, S: KeyScheme>(
        self,
        text: &[u8],
        keep: &[bool; 256],
        scheme: &S,
        path: &Path,
    ) -> Result<()> {
        let packing = scheme.packing();
        let bmi2 = has_bmi2();

        // ---
        // count the suffixes of each bucket, those whose first byte
        // `keep` marks, and cut the buckets into partitions

        let now = Instant::now();
        let num_buckets = packing.num_buckets();
        let counts = (0..text.len().div_ceil(CHUNK_SIZE))
            .into_par_iter()
            .fold(
                || vec![0; num_buckets],
                |mut counts, chunk| {
                    let start = chunk * CHUNK_SIZE;
                    let end = (start + CHUNK_SIZE).min(text.len());
                    #[cfg(target_arch = "x86_64")]
                    if bmi2 {
                        // SAFETY: bmi2 is the CPU check from has_bmi2
                        unsafe {
                            count_chunk_bmi2(
                                text,
                                keep,
                                scheme,
                                start..end,
                                &mut counts,
                            )
                        };
                        return counts;
                    }
                    count_chunk(text, keep, scheme, start..end, &mut counts);
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

        let (is_dna, ignore_softmask) = match self.alphabet {
            Alphabet::Dna { ignore_softmask } => (true, ignore_softmask),
            Alphabet::Amino | Alphabet::Bytes => (false, false),
        };
        let mut header = Header {
            version: OUTFILE_VERSION,
            is_dna,
            ignore_softmask,
            text_len: text.len(),
            num_suffixes,
            sort_type: match &self.seed_mask {
                Some(seed_mask) => SuffixSortType::Mask(seed_mask.clone()),
                None => SuffixSortType::MaxQueryLen(0),
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
                        let range = start..end;
                        let n = 'n: {
                            #[cfg(target_arch = "x86_64")]
                            if bmi2 {
                                // SAFETY: bmi2 is the CPU check from has_bmi2
                                break 'n unsafe {
                                    gather_chunk_bmi2(
                                        text, keep, scheme, range, lo, width, suffixes,
                                        keys,
                                    )
                                };
                            }
                            gather_chunk(
                                text, keep, scheme, range, lo, width, suffixes, keys,
                            )
                        };
                        let suffixes = &suffixes[..n];
                        let keys = &keys[..n];

                        let partition =
                            |i: usize| partition_of[packing.bucket(keys[i])] - first;
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
                    let (sa, lcp) =
                        sort_partition(text, scheme, suffixes, keys, &lcp_cache);

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
                    let (last, first) = (last_suffixes[i - 1], first_suffixes[i]);
                    W::from_usize(scheme.lcp(text, last, first))
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
    let now = Instant::now();
    let ignore_softmask = matches!(
        alphabet,
        Alphabet::Dna {
            ignore_softmask: true,
            ..
        }
    );
    let mut table = [0u8; 256];
    for (b, to) in table.iter_mut().enumerate() {
        let b = b as u8;
        *to = match (b.is_ascii_lowercase(), ignore_softmask) {
            (false, _) => b,
            (true, false) => b.to_ascii_uppercase(),
            (true, true) => b'N',
        };
    }
    let occupancy = text
        .par_chunks_mut(CHUNK_SIZE)
        .map(|chunk| {
            let mut occupancy = [false; 256];
            for b in chunk {
                *b = table[*b as usize];
                occupancy[*b as usize] = true;
            }
            occupancy
        })
        .reduce(
            || [false; 256],
            |mut total, occupancy| {
                for (total, seen) in total.iter_mut().zip(occupancy) {
                    *total |= seen;
                }
                total
            },
        );
    info!("Normalized the text in {:?}", now.elapsed());
    occupancy
}

/// The first bytes of the suffixes that are indexed: `ACGT` and the
/// sentinel for `Dna`, else every byte.
fn keep_table(alphabet: &Alphabet) -> [bool; 256] {
    match alphabet {
        Alphabet::Dna { .. } => {
            let mut keep = [false; 256];
            for &c in b"ACGT" {
                keep[c as usize] = true;
            }
            keep[SENTINEL_CHARACTER as usize] = true;
            keep
        }
        Alphabet::Amino | Alphabet::Bytes => [true; 256],
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
    let capacity = fs::metadata(path).map_or(0, |m| m.len() as usize);
    let mut seq: Vec<u8> = Vec::with_capacity(capacity);
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

/// The suffixes one task of a pass over the text takes.
const CHUNK_SIZE: usize = 1 << 19;

/// The number of bits in a bucket number.
//
// TODO: no commit or comment records why this is 16
const BUCKET_BITS: u32 = 16;

/// The maximum number of windows that a mask can contain.
const MAX_WINDOWS: usize = 4;

/// How the characters of one text become ranks, and how ranks pack into
/// a 64-bit key.
#[derive(Clone, Copy)]
pub(super) struct Packing {
    /// For each byte value, the position of that byte in increasing order
    /// among the bytes of the text.
    rank: [u8; 256],

    /// The number of bits for each character in a key.
    //
    //   ceil(log2(number of distinct bytes))
    //
    // at least 1, because a text of one record with no bases is only
    // the sentinel, and 0 bits would make chars_per_key divide by 0
    bits_per_char: u32,

    /// The number of characters in a key.
    //
    //   64 / bits_per_char, rounded down
    chars_per_key: usize,

    /// The number of characters in a bucket number. These characters are
    /// the first characters of the key.
    //
    //   BUCKET_BITS / bits_per_char, rounded down
    chars_per_bucket: usize,

    /// The number of low bits in a key that do not contain a character.
    //
    //   64 - bits_per_char * chars_per_key
    key_pad_bits: u32,
}

impl Packing {
    /// The packing for a text whose bytes `occupancy` marks.
    pub(super) fn from_occupancy(occupancy: [bool; 256]) -> Self {
        let mut rank = [0u8; 256];
        let mut distinct: usize = 0;
        for (c, &seen) in occupancy.iter().enumerate() {
            if seen {
                rank[c] = distinct as u8;
                distinct += 1;
            }
        }
        let bits_per_char = distinct.next_power_of_two().trailing_zeros().max(1);
        let chars_per_key = 64 / bits_per_char as usize;
        let chars_per_bucket = (BUCKET_BITS / bits_per_char) as usize;
        let key_pad_bits = 64 - bits_per_char * chars_per_key as u32;
        info!(
            "{distinct} distinct bytes at {bits_per_char} bits: \
             {chars_per_key}-character keys, {chars_per_bucket}-character buckets"
        );
        Self {
            rank,
            bits_per_char,
            chars_per_key,
            chars_per_bucket,
            key_pad_bits,
        }
    }

    /// The number of distinct bucket numbers.
    pub(super) fn num_buckets(&self) -> usize {
        1 << (self.bits_per_char * self.chars_per_bucket as u32)
    }

    /// The bucket number of a suffix with the key `key`. The bucket number
    /// is the first `chars_per_bucket` characters of the key.
    #[inline(always)]
    pub(super) fn bucket(&self, key: u64) -> usize {
        (key >> (64 - self.bits_per_char * self.chars_per_bucket as u32)) as usize
    }

    /// The number of characters at the start of the unequal keys `a` and
    /// `b` that are the same.
    #[inline(always)]
    pub(super) fn shared_chars(&self, a: u64, b: u64) -> usize {
        ((a ^ b).leading_zeros() / self.bits_per_char) as usize
    }

    /// The key of the suffix at `i`. The key contains the first
    /// `chars_per_key` characters of the suffix.
    pub(super) fn key_at(&self, text: &[u8], i: usize) -> u64 {
        // the ranks go in from the top bit down, so keys compare as the
        // characters do. for protein, with 5 bits per character and 12
        // characters per key:
        //
        //   bits  63-59  58-54  ...  8-4  3-0
        //         c0     c1     ...  c11  pad
        let mut key = 0u64;
        for j in 0..self.chars_per_key {
            // a character past the end of the text ranks 0, so the key
            // of a short suffix still has chars_per_key characters
            let rank = text.get(i + j).map_or(0, |&c| self.rank[c as usize] as u64);
            key = (key << self.bits_per_char) | rank;
        }
        key << self.key_pad_bits
    }

    /// The key of the suffix at `i + 1`, from the key `key` of the suffix
    /// at `i`. The first character moves out at the top. The character
    /// at `i + chars_per_key` comes in at the bottom.
    #[inline(always)]
    pub(super) fn roll(&self, text: &[u8], i: usize, key: u64) -> u64 {
        let incoming = text
            .get(i + self.chars_per_key)
            .map_or(0, |&c| self.rank[c as usize] as u64);
        (key << self.bits_per_char) | (incoming << self.key_pad_bits)
    }
}

/// How one build makes the key of each suffix.
pub(super) trait KeyScheme: Sync {
    /// What a scan holds at one suffix: enough to give the key of that
    /// suffix and to step to the next suffix cheaply.
    type RollingKeyState: Copy;

    /// The packing that the keys of this scheme use.
    fn packing(&self) -> &Packing;

    /// The rolling key state of the suffix at `i`.
    fn start(&self, text: &[u8], i: usize) -> Self::RollingKeyState;

    /// Step `state` from the suffix at `i` to the suffix at `i + 1`.
    fn advance(&self, text: &[u8], i: usize, state: &mut Self::RollingKeyState);

    /// The key of the suffix that `state` is at.
    fn key(&self, state: &Self::RollingKeyState) -> u64;

    /// Order a run of suffixes with equal keys, and give each suffix after
    /// the first its LCP with the one before it.
    fn order_ties<W: Int>(
        &self,
        text: &[u8],
        sa: &mut [W],
        lcp: &mut [W],
        lcp_cache: &ThreadLocal<RefCell<LcpCache>>,
    );

    /// The LCP of the suffixes at `a` and `b`.
    fn lcp(&self, text: &[u8], a: usize, b: usize) -> usize;
}

/// The keys of a build without a seed mask. The key of a suffix contains
/// its first `chars_per_key` characters.
pub(super) struct UnmaskedKeys {
    pub(super) packing: Packing,
}

impl KeyScheme for UnmaskedKeys {
    type RollingKeyState = u64;

    fn packing(&self) -> &Packing {
        &self.packing
    }

    #[inline(always)]
    fn start(&self, text: &[u8], i: usize) -> u64 {
        self.packing.key_at(text, i)
    }

    #[inline(always)]
    fn advance(&self, text: &[u8], i: usize, key: &mut u64) {
        *key = self.packing.roll(text, i, *key);
    }

    #[inline(always)]
    fn key(&self, key: &u64) -> u64 {
        *key
    }

    fn order_ties<W: Int>(
        &self,
        text: &[u8],
        sa: &mut [W],
        lcp: &mut [W],
        lcp_cache: &ThreadLocal<RefCell<LcpCache>>,
    ) {
        merge_sort(
            text,
            &mut sa.to_vec(),
            &mut lcp.to_vec(),
            sa,
            lcp,
            lcp_cache,
        );
    }

    fn lcp(&self, text: &[u8], a: usize, b: usize) -> usize {
        find_lcp(text, a, b, text.len(), 0, None)
    }
}

/// The keys of a build with a seed mask. The key of a suffix contains the
/// characters at the care positions of the mask.
pub(super) struct MaskedKeys {
    packing: Packing,

    /// The care positions of the mask, in increasing order.
    positions: Vec<usize>,

    /// The number of windows that contain the mask.
    //
    //   ceil((last care position + 1) / chars_per_key)
    windows: usize,

    /// For each window, a bit mask with a 1 at each bit of each care
    /// position.
    //
    // care position p is character p % chars_per_key of window
    // p / chars_per_key
    care_bits: [u64; MAX_WINDOWS],

    /// For each window, the number of care-position bits in the windows
    /// after it.
    //
    //   window_offset[w] = sum of care_bits[v].count_ones() for v > w
    //
    // masked_key moves the care positions of window w up by this number,
    // then moves the full key up by pad_bits. the care positions of
    // window 0 are at the top of the masked key
    window_offset: [u32; MAX_WINDOWS],

    /// The number of low bits in a masked key that do not contain a
    /// character.
    //
    //   64 - bits_per_char * (number of care positions)
    pad_bits: u32,

    /// Whether the CPU has BMI2, so that `key` can use pext.
    bmi2: bool,
}

impl MaskedKeys {
    /// The keys for the mask with the care positions `positions` over
    /// `span` characters. The result is an error if a key cannot hold the
    /// mask.
    pub(super) fn new(
        packing: Packing,
        positions: &[usize],
        span: usize,
    ) -> Result<Self> {
        // a key holds every care position of a mask, so the sort by key
        // is the whole sort under one; a longer mask is refused, as is
        // one spanning more windows than the scan rolls
        if positions.len() > packing.chars_per_key {
            bail!(
                "the seed mask has {} care positions, and a key holds {} characters of this alphabet",
                positions.len(),
                packing.chars_per_key
            );
        }
        if span > MAX_WINDOWS * packing.chars_per_key {
            bail!(
                "the seed mask spans {} characters, and the scan covers {} of this alphabet",
                span,
                MAX_WINDOWS * packing.chars_per_key
            );
        }

        let windows =
            (positions[positions.len() - 1] + 1).div_ceil(packing.chars_per_key);
        let char_ones = (1u64 << packing.bits_per_char) - 1;
        let mut care_bits = [0; MAX_WINDOWS];
        for &p in positions {
            let (window, j) = (p / packing.chars_per_key, p % packing.chars_per_key);
            care_bits[window] |=
                char_ones << (64 - packing.bits_per_char * (j as u32 + 1));
        }
        let mut window_offset = [0; MAX_WINDOWS];
        let mut offset = 0;
        for window in (0..windows).rev() {
            window_offset[window] = offset;
            offset += care_bits[window].count_ones();
        }
        Ok(Self {
            packing,
            positions: positions.to_vec(),
            windows,
            care_bits,
            window_offset,
            pad_bits: 64 - packing.bits_per_char * positions.len() as u32,
            bmi2: has_bmi2(),
        })
    }

    /// The masked key of the suffix with the windows `windows`. The masked
    /// key contains the characters at the care positions, in order, from
    /// the top bit.
    #[inline(always)]
    fn masked_key(&self, windows: &[u64; MAX_WINDOWS]) -> u64 {
        let bits_per_char = self.packing.bits_per_char;
        let mut key = 0;
        let per_window = self.care_bits.iter().zip(&self.window_offset);
        for (&window, (&care_bits, &offset)) in
            windows.iter().zip(per_window).take(self.windows)
        {
            let mut chars = 0;
            let mut care = care_bits;
            while care != 0 {
                let top_zeros = care.leading_zeros();
                let char_low = 64 - top_zeros - bits_per_char;
                chars = (chars << bits_per_char)
                    | ((window >> char_low) & ((1 << bits_per_char) - 1));
                care &= !(((1u64 << bits_per_char) - 1) << char_low);
            }
            key |= chars << offset;
        }
        key << self.pad_bits
    }

    /// `masked_key` by pext.
    #[cfg(target_arch = "x86_64")]
    #[target_feature(enable = "bmi2")]
    #[inline]
    fn masked_key_pext(&self, windows: &[u64; MAX_WINDOWS]) -> u64 {
        let mut key = 0;
        let per_window = self.care_bits.iter().zip(&self.window_offset);
        for (&window, (&care_bits, &offset)) in
            windows.iter().zip(per_window).take(self.windows)
        {
            key |= std::arch::x86_64::_pext_u64(window, care_bits) << offset;
        }
        key << self.pad_bits
    }
}

impl KeyScheme for MaskedKeys {
    /// The windows of the suffix. Window `w` is the key of the suffix that
    /// starts `w * chars_per_key` characters later.
    type RollingKeyState = [u64; MAX_WINDOWS];

    fn packing(&self) -> &Packing {
        &self.packing
    }

    #[inline(always)]
    fn start(&self, text: &[u8], i: usize) -> [u64; MAX_WINDOWS] {
        let mut windows = [0; MAX_WINDOWS];
        for (w, window) in windows.iter_mut().enumerate().take(self.windows) {
            *window = self
                .packing
                .key_at(text, i + w * self.packing.chars_per_key);
        }
        windows
    }

    #[inline(always)]
    fn advance(&self, text: &[u8], i: usize, windows: &mut [u64; MAX_WINDOWS]) {
        for (w, window) in windows.iter_mut().enumerate().take(self.windows) {
            *window =
                self.packing
                    .roll(text, i + w * self.packing.chars_per_key, *window);
        }
    }

    #[inline(always)]
    fn key(&self, windows: &[u64; MAX_WINDOWS]) -> u64 {
        #[cfg(target_arch = "x86_64")]
        if self.bmi2 {
            // SAFETY: new sets bmi2 from has_bmi2
            return unsafe { self.masked_key_pext(windows) };
        }
        self.masked_key(windows)
    }

    fn order_ties<W: Int>(
        &self,
        text: &[u8],
        sa: &mut [W],
        lcp: &mut [W],
        _lcp_cache: &ThreadLocal<RefCell<LcpCache>>,
    ) {
        // equal masked keys agree at every care position, so the run goes
        // shortest first, and the LCP is the care positions that fit in
        // the shorter suffix
        sa.sort_unstable_by(|a, b| b.cmp(a));
        let last = self.positions[self.positions.len() - 1];
        for k in 1..sa.len() {
            let shorter = sa[k - 1].to_usize();
            lcp[k] = W::from_usize(if shorter + last < text.len() {
                self.positions.len()
            } else {
                self.positions
                    .partition_point(|&p| shorter + p < text.len())
            });
        }
    }

    fn lcp(&self, text: &[u8], a: usize, b: usize) -> usize {
        masked_lcp(&text[a..], &text[b..], &self.positions)
    }
}

/// Whether the CPU has bmi2, for pext.
fn has_bmi2() -> bool {
    #[cfg(target_arch = "x86_64")]
    {
        is_x86_feature_detected!("bmi2")
    }
    #[cfg(not(target_arch = "x86_64"))]
    false
}

/// Count the suffixes of `positions` whose first byte `keep` marks into
/// `counts` by bucket.
#[inline(always)]
fn count_chunk<S: KeyScheme>(
    text: &[u8],
    keep: &[bool; 256],
    scheme: &S,
    positions: Range<usize>,
    counts: &mut [usize],
) {
    let packing = scheme.packing();
    let mut state = scheme.start(text, positions.start);
    for i in positions {
        if keep[text[i] as usize] {
            counts[packing.bucket(scheme.key(&state))] += 1;
        }
        scheme.advance(text, i, &mut state);
    }
}

/// `count_chunk`, compiled with BMI2.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "bmi2")]
fn count_chunk_bmi2<S: KeyScheme>(
    text: &[u8],
    keep: &[bool; 256],
    scheme: &S,
    positions: Range<usize>,
    counts: &mut [usize],
) {
    count_chunk(text, keep, scheme, positions, counts)
}

/// Gather the suffixes of `positions` whose first byte `keep` marks and
/// whose bucket is within `width` of `lo`, with their keys, into
/// `suffixes` and `keys`. The result is the number gathered.
#[allow(clippy::too_many_arguments)]
#[inline(always)]
fn gather_chunk<W: Int, S: KeyScheme>(
    text: &[u8],
    keep: &[bool; 256],
    scheme: &S,
    positions: Range<usize>,
    lo: usize,
    width: usize,
    suffixes: &mut [W],
    keys: &mut [u64],
) -> usize {
    let packing = scheme.packing();
    let mut state = scheme.start(text, positions.start);
    let mut n = 0;
    for i in positions {
        let key = scheme.key(&state);
        let bucket = packing.bucket(key);
        if keep[text[i] as usize] && bucket.wrapping_sub(lo) < width {
            suffixes[n] = W::from_usize(i);
            keys[n] = key;
            n += 1;
        }
        scheme.advance(text, i, &mut state);
    }
    n
}

/// `gather_chunk`, compiled with BMI2.
#[cfg(target_arch = "x86_64")]
#[allow(clippy::too_many_arguments)]
#[target_feature(enable = "bmi2")]
fn gather_chunk_bmi2<W: Int, S: KeyScheme>(
    text: &[u8],
    keep: &[bool; 256],
    scheme: &S,
    positions: Range<usize>,
    lo: usize,
    width: usize,
    suffixes: &mut [W],
    keys: &mut [u64],
) -> usize {
    gather_chunk(text, keep, scheme, positions, lo, width, suffixes, keys)
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

/// Sort the `suffixes` of a partition, each with its key, and give each
/// its LCP with the one before it: by key, then each run of equal keys
/// in the order of `scheme`.
fn sort_partition<W: Int, S: KeyScheme>(
    text: &[u8],
    scheme: &S,
    suffixes: Vec<W>,
    keys: Vec<u64>,
    lcp_cache: &ThreadLocal<RefCell<LcpCache>>,
) -> (Vec<W>, Vec<W>) {
    let len = suffixes.len();
    let mut keyed: Vec<(u64, W)> = keys.into_iter().zip(suffixes).collect();
    keyed.sort_unstable_by_key(|&(key, _)| key);
    let mut sa: Vec<W> = keyed.iter().map(|&(_, s)| s).collect();
    let mut lcp = vec![W::default(); len];
    let mut run = 0;
    for i in 1..=len {
        if i < len && keyed[i].0 == keyed[run].0 {
            continue;
        }
        if i - run > 1 {
            scheme.order_ties(text, &mut sa[run..i], &mut lcp[run..i], lcp_cache);
        }
        // the LCP across unequal keys is read off the keys
        if run > 0 {
            lcp[run] = W::from_usize(
                scheme
                    .packing()
                    .shared_chars(keyed[run - 1].0, keyed[run].0),
            );
        }
        run = i;
    }
    (sa, lcp)
}

/// Sort the suffixes in `out` with their LCPs, using `work`, which
/// starts with the same suffixes, as scratch space.
fn merge_sort<W: Int>(
    text: &[u8],
    work_sa: &mut [W],
    work_lcp: &mut [W],
    out_sa: &mut [W],
    out_lcp: &mut [W],
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
            merge_sort(text, out_sa_l, out_lcp_l, work_sa_l, work_lcp_l, lcp_cache);
            merge_sort(text, out_sa_r, out_lcp_r, work_sa_r, work_lcp_r, lcp_cache);
        } else {
            rayon::join(
                || {
                    merge_sort(
                        text, out_sa_l, out_lcp_l, work_sa_l, work_lcp_l, lcp_cache,
                    )
                },
                || {
                    merge_sort(
                        text, out_sa_r, out_lcp_r, work_sa_r, work_lcp_r, lcp_cache,
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
            &mut lcp_cache.get_or_default().borrow_mut(),
        );
    }
}

/// Merge the sorted runs `from[..mid]` and `from[mid..]` into `into`.
fn merge<W: Int>(
    text: &[u8],
    from_sa: &mut [W],
    from_lcp: &mut [W],
    mid: usize,
    target_sa: &mut [W],
    target_lcp: &mut [W],
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

                let context = max_n;

                // LCP(X_i, Y_j)
                let len_lcp = if m < context {
                    find_lcp(
                        text,
                        x[idx_x].to_usize(),
                        y[idx_y].to_usize(),
                        (context - m).to_usize(),
                        m.to_usize(), // skip
                        Some(lcp_cache),
                    )
                } else {
                    context.to_usize()
                };
                let len_lcp = W::from_usize(len_lcp);

                // an LCP reaching the context sorts the shorter suffix
                // first, and otherwise the byte after the LCP decides
                if len_lcp >= context {
                    target_sa[idx_target] = shorter_suffix;
                } else {
                    let cmp = text[(x[idx_x] + len_lcp).to_usize()]
                        .cmp(&text[(y[idx_y] + len_lcp).to_usize()]);

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
#[inline(always)]
fn find_lcp(
    text: &[u8],
    start1: usize,
    start2: usize,
    len: usize,
    skip: usize,
    lcp_cache: Option<&mut LcpCache>,
) -> usize {
    suffix_lcp(text, start1, start2, skip, len, lcp_cache)
}

#[cfg(test)]
mod test {
    use super::{
        counting_sort, cut_partitions, find_lcp, keep_table, read_sequence_file, turns,
        Alphabet, KeyScheme, MaskedKeys, Packing, SuffixArrayBuilder,
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
                ignore_softmask: false,
            }),
            b"$ACGT"
        );
        assert_eq!(kept(Alphabet::Amino), every_byte);
        assert_eq!(kept(Alphabet::Bytes), every_byte);
    }

    #[test]
    fn test_write() -> Result<()> {
        let outfile = NamedTempFile::new()?;
        SuffixArrayBuilder::from_sequence_file("../data/inputs/3.fa", b'%')?
            .with_alphabet(Alphabet::Dna {
                ignore_softmask: true,
            })
            .build_to_file(outfile.path())?;
        let suffix_array = SuffixArray::open(&outfile.path().to_string_lossy())?;
        let header = suffix_array.header();
        assert_eq!(header.text_len, 113);
        assert_eq!(header.num_suffixes, 101);
        Ok(())
    }

    fn packing(text: &[u8]) -> Packing {
        let mut occupancy = [false; 256];
        for &b in text {
            occupancy[b as usize] = true;
        }
        Packing::from_occupancy(occupancy)
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
        let dna = packing(b"$%ACGNT");
        assert_eq!(
            (dna.bits_per_char, dna.chars_per_key, dna.chars_per_bucket),
            (3, 21, 5)
        );
        assert_eq!(dna.num_buckets(), 32768);

        let protein = packing(b"$%*-ABCDEFGHIKLMNPQRSTVWY");
        assert_eq!(
            (
                protein.bits_per_char,
                protein.chars_per_key,
                protein.chars_per_bucket
            ),
            (5, 12, 3)
        );

        let bytes: Vec<u8> = (0..=255).collect();
        let bytes = packing(&bytes);
        assert_eq!(
            (
                bytes.bits_per_char,
                bytes.chars_per_key,
                bytes.chars_per_bucket
            ),
            (8, 8, 2)
        );
        assert_eq!(bytes.num_buckets(), 65536);

        let one = packing(b"$");
        assert_eq!(
            (one.bits_per_char, one.chars_per_key, one.chars_per_bucket),
            (1, 64, 16)
        );
    }

    #[test]
    fn test_key_at() -> Result<()> {
        // ranks in byte order: $=0 A=1 C=2 G=3 N=4 T=5, 3 bits each
        let text = b"ACGTNA$";
        let packing = packing(text);

        // A C G T N A $ and 14 characters past the end
        let key = packing.key_at(text, 0);
        assert_eq!(key, packed(&[1, 2, 3, 5, 4, 1, 0], 3));
        assert_eq!(
            packing.bucket(key),
            packed(&[1, 2, 3, 5, 4], 3) as usize >> 49
        );

        // N A $ and the rest past the end
        assert_eq!(packing.key_at(text, 4), packed(&[4, 1, 0], 3));

        // care positions 0, 1 and 3, then past the mask: A C T
        let mask = SeedMask::new("1101")?;
        let scheme = MaskedKeys::new(packing, &mask.positions, mask.bytes.len())?;
        let masked = |i| scheme.key(&scheme.start(text, i));
        assert_eq!(masked(0), packed(&[1, 2, 5], 3));
        // A $ and position 8 past the end
        assert_eq!(masked(5), packed(&[1], 3));
        Ok(())
    }

    /// The masked key by lookup: the characters at the care positions.
    fn lookup(packing: &Packing, text: &[u8], i: usize, positions: &[usize]) -> u64 {
        let mut key = 0u64;
        for &off in positions {
            let rank = text
                .get(i + off)
                .map_or(0, |&c| packing.rank[c as usize] as u64);
            key = (key << packing.bits_per_char) | rank;
        }
        key << (64 - packing.bits_per_char * positions.len() as u32)
    }

    #[test]
    fn masked_key_matches_the_lookup_at_every_position() -> Result<()> {
        let text = b"ACGTNACGTTTGANNACGTACGTACGTAAAGATTACAGATTACACCCGGGTTTAAANNNNACGT$";
        for mask in [
            "101",
            "1101",
            "11011011011011011",
            "111010010100110111",
            "11011011011011011011011011011",
            "1000000000000000000000000000000000000000001",
        ] {
            let mask = SeedMask::new(mask)?;
            let packing = packing(text);
            let scheme = MaskedKeys::new(packing, &mask.positions, mask.bytes.len())?;
            let mut windows = scheme.start(text, 0);
            for i in 0..text.len() {
                let expected = lookup(&packing, text, i, &mask.positions);
                assert_eq!(scheme.masked_key(&windows), expected, "{mask:?} at {i}");
                #[cfg(target_arch = "x86_64")]
                if super::has_bmi2() {
                    // SAFETY: bmi2 is the CPU check from has_bmi2
                    let pext = unsafe { scheme.masked_key_pext(&windows) };
                    assert_eq!(pext, expected, "{mask:?} at {i} by pext");
                }
                scheme.advance(text, i, &mut windows);
            }
        }
        Ok(())
    }

    #[test]
    fn test_roll_reaches_the_same_key_as_key_at() {
        let text = b"ACGTNACGTTTGANNACGTACGTACGTAAA$";
        let packing = packing(text);
        let mut key = packing.key_at(text, 0);
        for i in 0..text.len() {
            assert_eq!(key, packing.key_at(text, i), "at {i}");
            key = packing.roll(text, i, key);
        }
    }

    #[test]
    fn test_shared_chars_counts_equal_leading_characters() {
        let text = b"ACGTACGA$";
        let packing = packing(text);
        let a = packing.key_at(text, 0);
        let b = packing.key_at(text, 4);
        // ACGTACGA$ and ACGA$ share ACG
        assert_eq!(packing.shared_chars(a, b), 3);
        // the lowest bit of a character deciding is still that character
        let c = packed(&[1, 2, 3, 4], 3);
        let d = packed(&[1, 2, 3, 5], 3);
        assert_eq!(packing.shared_chars(c, d), 3);
        assert_eq!(packing.shared_chars(a, packing.key_at(text, 8)), 0);
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
    fn a_mask_past_the_key_is_refused() -> Result<()> {
        // 5 distinct bytes at 3 bits: 21 characters in a key
        let text = b"ACGTACGTAC$".to_vec();
        let outfile = NamedTempFile::new()?;
        let build = |mask: &str| {
            SuffixArrayBuilder::from_data(text.clone(), vec![(0, "s".to_string())])
                .with_seed_mask(SeedMask::new(mask)?)
                .build_to_file(outfile.path())
        };
        assert!(build(&format!("{}01", "1".repeat(20))).is_ok());
        assert!(build(&format!("{}01", "1".repeat(21))).is_err());
        Ok(())
    }

    #[test]
    fn test_find_lcp_no_seed_mask() -> Result<()> {
        //           012345
        let text = b"TTTAGC";

        // 0: TTTAGC
        // 1:  TTAGC
        // 6: len of text
        assert_eq!(find_lcp(text, 0, 1, 6, 0, None), 2);

        // 0: TTTAGC
        // 2:   TAGC
        // 6: len of text
        assert_eq!(find_lcp(text, 0, 2, 6, 0, None), 1);

        // 0: TTTAGC
        // 1:  TTAGC
        // 1: compare one byte
        assert_eq!(find_lcp(text, 0, 1, 1, 0, None), 1);

        // 0: TTTAGC
        // 3:    AGC
        // 6: len of text
        assert_eq!(find_lcp(text, 0, 3, 6, 0, None), 0);

        Ok(())
    }
}
