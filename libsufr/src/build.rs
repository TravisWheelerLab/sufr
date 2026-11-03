mod partition;
mod scratch;
mod sort;
mod write;

use crate::{
    format::{Header, SeedMask, SuffixSortType, OUTFILE_VERSION},
    index::{Index, Storage, Text},
    io::Int,
    SuffixArray,
};
use anyhow::Result;
use log::info;
use needletail::parse_fastx_file;
use partition::{group_buffers, partition, PrefixEncoder};
use scratch::ScratchKind;
use sort::{first_lcps, sort_partitions};
use std::{mem, ops::Range, path::Path, time::Instant};
use write::write;

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
    memory_scratch: bool,
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
            memory_scratch: false,
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

    /// Set whether to keep the partitions in memory during construction
    /// instead of in temporary files (default: `false`).
    pub fn with_memory_scratch(mut self, memory_scratch: bool) -> Self {
        self.memory_scratch = memory_scratch;
        self
    }

    // only Index::<u64> can read a small text written this
    // way, so this is for tests
    #[cfg(test)]
    pub(crate) fn with_force_64bit(mut self, force_64bit: bool) -> Self {
        self.force_64bit = force_64bit;
        self
    }

    /// Set whether to write the LCP array (default: `false`).
    pub fn with_write_lcp(mut self, write_lcp: bool) -> Self {
        self.write_lcp = write_lcp;
        self
    }

    /// Build the suffix array, write it to `path` and open the result.
    /// With `low_memory`, the suffixes stay on disk.
    pub fn build_to_file(
        self,
        path: impl AsRef<Path>,
        low_memory: bool,
    ) -> Result<SuffixArray> {
        let path = path.as_ref();
        let fits_u32 = (self.text.len() as u64) < u32::MAX as u64 && !self.force_64bit;
        if fits_u32 {
            Ok(SuffixArray::U32(construct::<u32>(
                path,
                self.text,
                self.sequences,
                self.alphabet,
                self.seed_mask,
                self.max_query_len,
                self.num_partitions,
                self.write_lcp,
                self.memory_scratch,
                low_memory,
            )?))
        } else {
            Ok(SuffixArray::U64(construct::<u64>(
                path,
                self.text,
                self.sequences,
                self.alphabet,
                self.seed_mask,
                self.max_query_len,
                self.num_partitions,
                self.write_lcp,
                self.memory_scratch,
                low_memory,
            )?))
        }
    }
}

/// Build the suffix array of `text` with `T` per suffix and LCP, write it
/// to `path` and open the result.
#[allow(clippy::too_many_arguments)]
fn construct<T: Int>(
    path: &Path,
    mut text: Vec<u8>,
    sequences: Vec<(usize, String)>,
    alphabet: Alphabet,
    mask: Option<SeedMask>,
    max_query_len: usize,
    num_partitions: usize,
    write_lcp: bool,
    memory_scratch: bool,
    low_memory: bool,
) -> Result<Index<T>> {
    let occupancy = normalize(&mut text, &alphabet);
    let keep = keep_table(&alphabet);
    let n_ranges = match alphabet {
        Alphabet::Dna {
            allow_ambiguity: true,
            ..
        } => find_n_ranges(&text),
        _ => vec![],
    };
    let encoder = PrefixEncoder::from_occupancy(occupancy);
    let mask = mask.as_ref();
    let scratch = if memory_scratch {
        ScratchKind::Memory
    } else {
        ScratchKind::Disk
    };

    let partitioned = partition::<T>(&text, &encoder, &keep, mask, scratch)?;
    let header = header(
        &alphabet,
        &sequences,
        text.len(),
        partitioned.num_suffixes,
        mask,
        max_query_len,
        write_lcp,
        mem::size_of::<T>(),
    );
    let groups = group_buffers(partitioned, num_partitions)?;
    let partitions =
        sort_partitions(&text, groups, mask, max_query_len, &n_ranges, write_lcp)?;
    let first_lcps = if write_lcp {
        first_lcps(&text, &partitions, mask, max_query_len, &n_ranges)
    } else {
        vec![]
    };
    write(path, &header, &text, partitions, &first_lcps)?;

    let suffixes = if low_memory {
        Storage::Disk
    } else {
        Storage::Memory
    };
    Index::new(
        &path.to_string_lossy(),
        header,
        Text::Memory(text),
        suffixes,
        None,
    )
}

/// The header of the output file, with `int_size` bytes per
/// suffix array value.
#[allow(clippy::too_many_arguments)]
fn header(
    alphabet: &Alphabet,
    sequences: &[(usize, String)],
    text_len: usize,
    num_suffixes: usize,
    mask: Option<&SeedMask>,
    max_query_len: usize,
    write_lcp: bool,
    int_size: usize,
) -> Header {
    let (is_dna, allow_ambiguity, ignore_softmask) = match *alphabet {
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
        text_len,
        num_suffixes,
        sort_type: match mask {
            Some(seed_mask) => SuffixSortType::Mask(seed_mask.clone()),
            None => SuffixSortType::MaxQueryLen(max_query_len),
        },
        sequence_starts: sequences.iter().map(|(start, _)| *start).collect(),
        sequence_names: sequences.iter().map(|(_, name)| name.clone()).collect(),
        int_size,
        text_pos: 0,
        suffix_array_pos: 0,
        lcp_pos: None,
    };
    header.set_positions(write_lcp);
    header
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

#[cfg(test)]
mod test {
    use super::{keep_table, read_sequence_file, Alphabet, SuffixArrayBuilder};
    use crate::SuffixArray;
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
    fn test_write() -> Result<()> {
        let outfile = NamedTempFile::new()?;
        SuffixArrayBuilder::from_sequence_file("../data/inputs/3.fa", b'%')?
            .with_alphabet(Alphabet::Dna {
                allow_ambiguity: false,
                ignore_softmask: true,
            })
            .build_to_file(outfile.path(), false)?;
        let suffix_array = SuffixArray::open(&outfile.path().to_string_lossy())?;
        let meta = suffix_array.metadata()?;
        assert_eq!(meta.text_len, 113);
        assert_eq!(meta.num_suffixes, 101);
        Ok(())
    }
}
