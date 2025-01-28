//! Common types

use anyhow::{bail, Result};
use chrono::{DateTime, Local};
use regex::Regex;
use std::{
    cmp::Ordering,
    fmt::{self, Debug, Display},
    hash::Hash,
    ops::Range,
    ops::{Add, AddAssign, Div, Sub},
};

// --------------------------------------------------
/// Serialization version
pub const OUTFILE_VERSION: u8 = 6;

/// The sentinel character placed at the end of the text
/// (and so must not occur in the given text)
pub const SENTINEL_CHARACTER: u8 = b'$';

// --------------------------------------------------
/// Describes the suffixes in a _.sufr_ file were sorted
#[derive(Debug, Clone, PartialEq)]
pub enum SuffixSortType {
    /// A maximum query length, which can be zero
    MaxQueryLen(usize),

    /// A seed mask of 1/0 for care/don't-care positions
    Mask(SeedMask),
}

// --------------------------------------------------
/// A struct describing a seed mask
#[derive(Debug, PartialEq, Clone)]
pub struct SeedMask {
    /// The given mask string of 1/0s
    pub mask: String,

    /// A vector of 1/0 `u8` valus representing the mask
    pub bytes: Vec<u8>,

    /// The offset positions of the "care" values (from the 1s)
    pub positions: Vec<usize>,

    /// A vector of "difference" values to add to the
    /// range 0..`weight` that will return the "care" positions
    pub differences: Vec<usize>,

    /// The number of "care" positions (popcount of 1)
    pub weight: usize,
}

// --------------------------------------------------
impl SeedMask {
    /// Create a new `SeedMask` from an input string.
    /// The string must:
    /// * be comprised entirely of 1 or 0
    /// * start and end with 1
    /// * contain at least one 0
    pub fn new(mask: &str) -> Result<Self> {
        if !Self::is_valid(mask) {
            bail!("Invalid seed mask '{mask}'")
        }

        let bytes = Self::parse(mask);
        let positions = Self::get_positions(&bytes);
        let differences = Self::get_differences(&positions);
        let weight = positions.len();

        Ok(Self {
            mask: mask.to_string(),
            bytes,
            positions,
            differences,
            weight,
        })
    }

    /// Instantiate a `SeedMask` from the byte representation
    /// from a `SufrFile`
    pub fn from_bytes(bytes: &[u8]) -> Result<Self> {
        let mask: String = bytes
            .iter()
            .filter_map(|b| match b {
                0 => Some('0'),
                1 => Some('1'),
                _ => None,
            })
            .collect();
        if mask.is_empty() {
            bail!("Bytes must be 1 or 0");
        }
        let positions = Self::get_positions(bytes);
        let differences = Self::get_differences(&positions);
        let weight = positions.len();

        Ok(Self {
            mask,
            bytes: bytes.to_vec(),
            positions,
            differences,
            weight,
        })
    }

    /// Determine if a seed mask is valid
    fn is_valid(mask: &str) -> bool {
        let seed_re = Regex::new("^1+0[01]*1$").unwrap();
        seed_re.is_match(mask)
    }

    /// Turn a valid string mask into a vector of `u8` values
    /// suitable for serialization.
    fn parse(mask: &str) -> Vec<u8> {
        mask.as_bytes()
            .iter()
            .flat_map(|b| match b {
                b'1' => Some(1),
                b'0' => Some(0),
                _ => None,
            })
            .collect()
    }

    /// Find the differences to add to each index to get the offsets
    /// of the "care" positions.
    fn get_differences(positions: &[usize]) -> Vec<usize> {
        // Mask: "1001101"
        // M: [0, 3, 4, 6]
        // U: [0, 1, 2, 3]
        // D: [0, 2, 2, 3]
        positions.iter().enumerate().map(|(i, &m)| m - i).collect()
    }

    /// Return the offsets of the "care" positions from the bytes
    fn get_positions(bytes: &[u8]) -> Vec<usize> {
        // [1, 0, 1] -> [0, 2]
        bytes
            .iter()
            .enumerate()
            .filter_map(|(i, &b)| (b == 1).then_some(i))
            .collect()
    }
}

impl Display for SeedMask {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.mask)
    }
}

// --------------------------------------------------
/// A struct for use in searching the suffix array
#[derive(Debug)]
pub struct SearchOptions {
    /// A vector of query strings
    pub queries: Vec<String>,

    /// A maximum query length to use.
    /// If the suffix array was sorted with a shorter MQL, that
    /// value will be used instead.
    pub max_query_len: Option<usize>,

    /// More memory will result in higher throughput/latency.
    /// When `true`, the suffix array will be placed
    /// into memory. When `false`, the suffix array will be
    /// read from disk. NB: initially reading the suffix array
    /// with `low_memory` will also cause the `text` to
    /// remain on disk.
    pub low_memory: bool,

    /// Whether or not to return location information or to simply count
    pub find_suffixes: bool,
}

// --------------------------------------------------
/// A struct representing the result of a suffix array search
#[derive(Debug, PartialEq)]
pub struct SearchResult<T>
where
    T: Int + FromUsize<T> + Sized + Send + Sync + serde::ser::Serialize,
{
    /// The ordinal position of the query
    pub query_num: usize,

    /// The query string
    pub query: String,

    /// A optional summary of the places where the query was found
    pub locations: Option<SearchResultLocations<T>>,
}

// --------------------------------------------------
/// The suffix locations matching a given query.
#[derive(Debug, PartialEq)]
pub struct SearchResultLocations<T> {
    /// The range of suffix ranks
    pub ranks: Range<usize>,

    /// A vector of the suffix locations
    pub suffixes: Vec<T>,
}

// --------------------------------------------------
/// A struct describing how a query compares to a suffix
#[derive(Debug)]
pub struct Comparison {
    /// Whether the suffix is greater, less, or equal
    pub cmp: Ordering,

    /// The length of the longest common prefix (LCP)
    /// between the query and the suffix
    pub lcp: usize,
}

// --------------------------------------------------
/// This struct is returned by `libsufr::utils::read_sequence_file`
/// for reading sequences from a FASTA/Q file.
#[derive(Debug)]
pub struct SequenceFileData {
    /// The sequence as a vector of bytes. Multiple sequences are
    /// separated by a user-supplied sequence delimiter.
    pub seq: Vec<u8>,

    /// The offsets where each sequence starts
    pub start_positions: Vec<usize>,

    /// The names of the sequences, will be the same length as `start_positions`
    pub sequence_names: Vec<String>,
}

// --------------------------------------------------
/// Trait to generically describe an "integer" of size `u8` (for text/bytes),
/// `u32` for suffix/LCP arrays over a text with a length < 2^32,
/// or `u64` for longer texts.
pub trait Int:
    Debug
    + AddAssign
    + Add<Output = Self>
    + Sub<Output = Self>
    + Div<Output = Self>
    + Copy
    + Default
    + Display
    + Ord
    + Hash
    + serde::ser::Serialize
{
    fn to_usize(&self) -> usize;
}

impl Int for u8 {
    fn to_usize(&self) -> usize {
        *self as usize
    }
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

pub trait FromUsize<T> {
    fn from_usize(val: usize) -> T;
}

impl FromUsize<u8> for u8 {
    fn from_usize(val: usize) -> u8 {
        val as u8
    }
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

// --------------------------------------------------
/// Options for counting the occurrences of suffixes
#[derive(Debug)]
pub struct CountOptions {
    /// Vector of query strings
    pub queries: Vec<String>,

    /// Maximum query length for search
    pub max_query_len: Option<usize>,

    /// Low memory, put suffix array into memory or not
    pub low_memory: bool,
}

// --------------------------------------------------
/// A struct representing the results of counting the occurrences of suffixes
#[derive(Debug, PartialEq)]
pub struct CountResult {
    /// The ordinal position of the original query
    pub query_num: usize,

    /// The query string
    pub query: String,

    /// Number of times a query was found
    pub count: usize,
}

// --------------------------------------------------
/// Options for extracting suffixes from a search
#[derive(Debug, Clone)]
pub struct ExtractOptions {
    /// Vector of query strings
    pub queries: Vec<String>,

    /// Maximum query length for search
    pub max_query_len: Option<usize>,

    /// Low memory options, put suffix array into memory
    pub low_memory: bool,

    /// Optional length of prefix to append before found suffixes (context)
    pub prefix_len: Option<usize>,

    /// Optional limit to the length of the suffix returned
    pub suffix_len: Option<usize>,
}

// --------------------------------------------------
/// The result of querying/extracting suffixes
#[derive(Debug, PartialEq)]
pub struct ExtractResult {
    /// The ordinal position of the original query
    pub query_num: usize,

    /// The query string
    pub query: String,

    /// A vector of the sequences containing the queries and locations
    pub sequences: Vec<ExtractSequence>,
}

// --------------------------------------------------
/// A struct describing a found query in the context of a sequence
/// This struct is used by the `sufr extract` command to print the
/// results of a search in FASTA format.
/// Normally, the query will be at the beginning of the result, so the
/// `sequence_start` will be 0.
/// When looking for alignment seeds, the user may request some prior context
/// via a `prefix_len`, in which case the sequence start will be whatever
/// length of prefix appended before the actual found hit. Note that the user
/// may request something like 100 characters of prefix but less than that is
/// appended because the hit was closer than 100 to the start of the sequence.
#[derive(Debug, PartialEq)]
pub struct ExtractSequence {
    /// The position of the suffix in the suffix array
    pub suffix: usize,

    /// The rank of the suffix in the suffix array
    pub rank: usize,

    /// The name of the sequence containing a query hit
    pub sequence_name: String,

    /// The start/offset of the containing sequence in the full `text`
    pub sequence_start: usize,

    /// The hit's relative start/stop range inside the sequence
    /// including the prefix/suffix lengths shown
    pub sequence_range: Range<usize>,

    /// The query hit's start position from the beginning of the shown context
    /// E.g., if the user requested a prefix of 10, then this value will be
    /// between 0-10, depending on the location of the hit inside the sequence.
    pub suffix_offset: usize,
}

// --------------------------------------------------
/// Arguments to sufr_file.list
pub struct ListOptions {
    /// Ranks of suffixes to show
    pub ranks: Vec<usize>,

    /// Show rank column
    pub show_rank: bool,

    /// Show suffix position column
    pub show_suffix: bool,

    /// Show LCP column
    pub show_lcp: bool,

    /// Low memory
    pub low_memory: bool,

    /// Length of suffixes to show
    pub len: Option<usize>,

    /// Number of suffixes to show
    pub number: Option<usize>,

    /// Output
    //pub output: &'a mut Box<dyn Write>,
    pub output: Option<String>,
}

// --------------------------------------------------
/// A struct for use in locating suffixes
#[derive(Debug)]
pub struct LocateOptions {
    /// A vector of query strings
    pub queries: Vec<String>,

    /// A maximum query length to use.
    /// If the suffix array was sorted with a shorter MQL, that
    /// value will be used instead.
    pub max_query_len: Option<usize>,

    /// More memory will result in higher throughput/latency.
    /// With `None`, the `text` and suffix array will be placed
    /// into memory. At low memory, the suffix array will be
    /// read from disk. At very low, the text will also be left
    /// on disk.
    pub low_memory: bool,
}

// --------------------------------------------------
/// A struct representing the results of a search that includes the
/// locations of the suffixes in their sequence context.
#[derive(Debug, PartialEq)]
pub struct LocateResult {
    /// The ordinal position of the original query
    pub query_num: usize,

    /// The query string
    pub query: String,

    /// A vector of positions where the query was found.
    /// This will be empty when the query was not present.
    pub positions: Vec<LocatePosition>,
}

// --------------------------------------------------
/// A struct representing the relative location of a query hit
/// in the context of a sequence.
#[derive(Debug, PartialEq)]
pub struct LocatePosition {
    /// The position of the suffix in the suffix array
    pub suffix: usize,

    /// The rank of the suffix in the suffix array
    pub rank: usize,

    /// The name of the sequence containing a query hit
    pub sequence_name: String,

    /// The start position of the hit in the sequence
    pub sequence_position: usize,
}

// --------------------------------------------------
/// The arguments for creating a `SufrBuilder` struct
#[derive(Clone, Debug)]
pub struct SufrBuilderArgs {
    /// Text as raw U8 bytes. cf. `libsufr::utils::read_sequence_file`
    /// Note: this value will be kept in memory during the build process.
    pub text: Vec<u8>,

    /// The path to the .sufr file that will be written.
    pub path: Option<String>,

    /// Use low memory when reading in resulting suffix array
    pub low_memory: bool,

    /// Maximum query length determines a prefix length of the suffixes.
    /// Without this value, suffixes will be fully sorted.
    pub max_query_len: Option<usize>,

    /// Indicates that the input is nucleotides, which has implications
    /// for ignoring ambiguity characters (not A, C, G, or T) and
    /// soft-masked/lowercase characters (usually indicating low-complexity
    /// regions).
    pub is_dna: bool,

    /// Whether or not to allow ambiguity characters (not A, C, G, or T)
    /// when handling nucleotides.
    pub allow_ambiguity: bool,

    /// Whether or not to ignore lowercased/softmasked nucleotide
    /// values (when `is_dna` is true).
    pub ignore_softmask: bool,

    /// When the `text` holds multiple sequences, this value contains
    /// the start positions of the sequences for locate queries.
    pub sequence_starts: Vec<usize>,

    /// When the `text` holds multiple sequences, this value contains
    /// the names of the sequences for locate queries.
    pub sequence_names: Vec<String>,

    /// The number of on-disk partitions to use when building the suffix array.
    /// Recommended to be at least the number of available CPUs, but
    /// a good number would place a 1-3 million suffixes into each partition,
    /// depending on the amount of available memory.
    /// The partitions are sorted independently and in parallel.
    /// Max memory usage will be determined by the average size of the
    /// partitions (which includes the number of suffixes in a partition
    /// and the integer size [`u32`, `u64`] to represent the suffixes)
    /// times the number of threads used to process concurrently.
    pub num_partitions: usize,

    /// An optional seed mask of 1/0 for care/don't-care positions,
    /// cf. `SeedMask`.
    pub seed_mask: Option<String>,

    /// A seed value for reproducibility when randomly choosing the
    /// suffixes for partitioning.
    pub random_seed: u64,
}

// --------------------------------------------------
/// A struct with metadata about the Sufr file
#[derive(Debug, PartialEq)]
pub struct SufrMetadata {
    /// Filename
    pub filename: String,

    /// Modified
    pub modified: DateTime<Local>,

    /// File size
    pub file_size: usize,

    /// File version
    pub file_version: usize,

    /// Nucleotides
    pub is_dna: bool,

    /// Allow ambiguity
    pub allow_ambiguity: bool,

    /// Ignore softmask
    pub ignore_softmask: bool,

    /// Text length
    pub text_len: usize,

    /// Number of suffixes
    pub len_suffixes: usize,

    /// Number of sequences
    pub num_sequences: usize,

    /// Start positions of sequences
    pub sequence_starts: Vec<usize>,

    /// Names of sequences
    pub sequence_names: Vec<String>,

    /// Sort type
    pub sort_type: SuffixSortType,
}

#[cfg(test)]
mod tests {
    use super::SeedMask;
    use anyhow::Result;
    use pretty_assertions::assert_eq;

    #[test]
    fn test_create_seed_mask() -> Result<()> {
        // No 0 at beginning
        let res = SeedMask::new("0101");
        assert!(res.is_err());

        // No 0 at end
        let res = SeedMask::new("1010");
        assert!(res.is_err());

        // Nothing other than 0/1
        let res = SeedMask::new("1021");
        assert!(res.is_err());

        // Must have a 0
        let res = SeedMask::new("1111");
        assert!(res.is_err());

        // Valid
        let res = SeedMask::new("101");
        assert!(res.is_ok());

        let expected = SeedMask {
            mask: "101".to_string(),
            bytes: vec![1, 0, 1],
            positions: vec![0, 2],
            differences: vec![0, 1],
            weight: 2,
        };
        assert_eq!(res.unwrap(), expected);

        let res = SeedMask::new("11101101101000011");
        assert!(res.is_ok());

        let expected = SeedMask {
            mask: "11101101101000011".to_string(),
            bytes: vec![1, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 0, 0, 0, 0, 1, 1],
            positions: vec![0, 1, 2, 4, 5, 7, 8, 10, 15, 16],
            differences: vec![0, 0, 0, 1, 1, 2, 2, 3, 7, 7],
            weight: 10,
        };
        assert_eq!(res.unwrap(), expected);

        Ok(())
    }

    #[test]
    fn test_display_seed_mask() -> Result<()> {
        for val in &["101", "11101101101000011"] {
            let mask = SeedMask::new(val).unwrap();
            assert_eq!(format!("{mask}"), val.to_string());
        }
        Ok(())
    }

    #[test]
    fn test_valid_seed_mask() -> Result<()> {
        let valid = ["101", "1001", "1101", "10101", "1110110110100001"];
        for pattern in valid {
            assert!(SeedMask::is_valid(pattern));
        }

        let invalid = [
            "", "abc", "1", "11", "111", "0", "00", "0111", "11100", "1a01",
        ];
        for pattern in invalid {
            assert!(!SeedMask::is_valid(pattern));
        }

        Ok(())
    }

    #[test]
    fn test_seed_mask_positions() -> Result<()> {
        assert_eq!(SeedMask::get_positions(&[1, 0, 1]), [0, 2]);
        assert_eq!(SeedMask::get_positions(&[1, 1, 0, 1, 1]), [0, 1, 3, 4]);
        assert_eq!(SeedMask::get_positions(&[1, 0, 1, 1, 0, 1]), [0, 2, 3, 5]);
        Ok(())
    }

    #[test]
    fn test_seed_mask_difference() -> Result<()> {
        // Empty is not a failure
        assert_eq!(SeedMask::get_differences(&[]), []);

        // "11011" -> [0, 1, 3, 4]
        //           - 0  1  2  3
        //            ------------
        //            [0, 0, 1, 1]
        assert_eq!(SeedMask::get_differences(&[0, 1, 3, 4]), [0, 0, 1, 1]);

        // "100001" -> [0, 5]
        //            - 0  1
        //            --------
        //             [0, 4]
        assert_eq!(SeedMask::get_differences(&[0, 5]), [0, 4]);

        // "1001101" -> [0, 3, 4, 6]
        //             - 0  1  2  3
        //             -------------
        //              [0, 2, 2, 3]
        assert_eq!(SeedMask::get_differences(&[0, 3, 4, 6]), [0, 2, 2, 3]);
        Ok(())
    }
}
