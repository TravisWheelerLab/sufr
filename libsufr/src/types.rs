//! Common types

use anyhow::{bail, Result};
use chrono::{DateTime, Local};
use regex::Regex;
use std::{
    cmp::Ordering,
    fmt::{self, Debug, Display},
    hash::Hash,
    ops::Range,
    ops::{Add, AddAssign, Sub},
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
    ///
    /// ```
    /// use anyhow::Result;
    /// use libsufr::types::SeedMask;
    ///
    /// fn main() -> Result<()> {
    ///     let mask = SeedMask::new("110110101")?;
    ///     let expected = SeedMask {
    ///         mask: "110110101".to_string(),
    ///         bytes: vec![1u8, 1u8, 0u8, 1u8, 1u8, 0u8, 1u8, 0u8, 1u8],
    ///         positions: vec![0, 1, 3, 4, 6, 8],
    ///         differences: vec![0, 0, 1, 1, 2, 3],
    ///         weight: 6,
    ///     };
    ///     assert_eq!(mask, expected);
    ///     Ok(())
    /// }
    /// ```
    ///
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
    ///
    /// ```
    /// use anyhow::Result;
    /// use libsufr::types::SeedMask;
    ///
    /// fn main() -> Result<()> {
    ///     let bytes: Vec<u8> = vec![1, 1, 0, 1, 1, 0, 1, 0, 1];
    ///     let mask = SeedMask::from_bytes(&bytes)?;
    ///     let expected = SeedMask {
    ///         mask: "110110101".to_string(),
    ///         bytes,
    ///         positions: vec![0, 1, 3, 4, 6, 8],
    ///         differences: vec![0, 0, 1, 1, 2, 3],
    ///         weight: 6,
    ///     };
    ///     assert_eq!(mask, expected);
    ///     Ok(())
    /// }
    /// ```
    ///
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
    ///
    /// ```
    /// use anyhow::Result;
    /// use libsufr::types::SeedMask;
    ///
    /// fn main() -> Result<()> {
    ///     assert!(SeedMask::new("").is_err());
    ///     assert!(SeedMask::new("0").is_err());
    ///     assert!(SeedMask::new("01").is_err());
    ///     assert!(SeedMask::new("10").is_err());
    ///     assert!(SeedMask::new("11").is_err());
    ///     assert!(SeedMask::new("101").is_ok());
    ///     Ok(())
    /// }
    /// ```
    ///
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

mod sealed {
    pub trait Sealed {}
    impl Sealed for u8 {}
    impl Sealed for u32 {}
    impl Sealed for u64 {}
}

// --------------------------------------------------
/// Trait to generically describe an "integer" of size `u8` (for text/bytes),
/// `u32` for suffix/LCP arrays over a text with a length < 2^32,
/// or `u64` for longer texts.
///
/// This trait is `Sealed` and cannot be implemented by downstream crates,
/// due to guarantees and assumptions in sufr's `unsafe` code.
pub trait Int:
    Debug
    + AddAssign
    + Add<Output = Self>
    + Sub<Output = Self>
    + Copy
    + Default
    + Display
    + Ord
    + Hash
    + Send
    + Sync
    + sealed::Sealed
{
    /// Convert an `Int` to a `usize`
    fn to_usize(&self) -> usize;

    /// Convert a `usize` to an `Int`
    fn from_usize(val: usize) -> Self;
}

impl Int for u8 {
    fn to_usize(&self) -> usize {
        *self as usize
    }

    fn from_usize(val: usize) -> Self {
        val as Self
    }
}

impl Int for u32 {
    fn to_usize(&self) -> usize {
        *self as usize
    }

    fn from_usize(val: usize) -> Self {
        val as Self
    }
}

impl Int for u64 {
    fn to_usize(&self) -> usize {
        *self as usize
    }

    fn from_usize(val: usize) -> Self {
        val as Self
    }
}

// --------------------------------------------------
/// Where a part of a _.sufr_ file is read from during searches
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub enum Load {
    /// Read from the file as needed
    #[default]
    Disk,

    /// Read into memory when the file is opened
    Memory,
}

// --------------------------------------------------
/// Options for opening a _.sufr_ file for searching
#[derive(Debug, Clone, Default)]
pub struct OpenOptions {
    /// Where the text is read from
    pub text: Load,

    /// Where the suffix array is read from
    pub suffix_array: Load,

    /// The number of leading query bytes a suffix must match to be a hit.
    /// `None` uses the length the file was sorted to, and a longer
    /// value is lowered to it.
    pub max_query_len: Option<usize>,
}

// --------------------------------------------------
/// How much of the surrounding sequence to include with an extracted hit
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct Context {
    /// The number of bytes to include before the hit
    pub before: usize,

    /// The number of bytes to include from the start of the hit.
    /// `None` runs to the end of the sequence.
    pub after: Option<usize>,
}

// --------------------------------------------------
/// A query hit in the context of the sequence containing it
#[derive(Debug, PartialEq)]
pub struct Location<'a> {
    /// The rank of the suffix in the suffix array
    pub rank: usize,

    /// The position of the suffix in the text
    pub suffix: usize,

    /// The name of the sequence containing the hit
    pub sequence_name: &'a str,

    /// The start of the containing sequence in the text
    pub sequence_start: usize,

    /// The position of the hit in the sequence
    pub sequence_position: usize,
}

// --------------------------------------------------
/// A query hit with a window of its sequence around it.
/// The window may be shorter than requested where it meets
/// the ends of the sequence.
#[derive(Debug, PartialEq)]
pub struct Extraction<'a> {
    /// The hit
    pub location: Location<'a>,

    /// The window, relative to the start of the sequence
    pub sequence_range: Range<usize>,

    /// The position of the hit in the window
    pub suffix_offset: usize,
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

    /// Whether to write the LCP array. Searching does not use it;
    /// listing LCP values does, and query-time max query lengths shorter
    /// than the build-time value are faster with it.
    pub write_lcp: bool,
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

    /// Whether or not the file has an on-disk LCP array
    pub has_lcp: bool,
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
