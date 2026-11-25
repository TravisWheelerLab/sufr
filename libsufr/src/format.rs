//! The _.sufr_ file layout

use crate::{
    error::{InvalidSeedMask, OpenError},
    io::write_all_at,
};
use std::{
    fmt::{self, Display},
    fs::File,
    io::{self, Read, Seek, SeekFrom},
    mem,
    path::Path,
    str::FromStr,
};

/// The serialization version.
pub(crate) const OUTFILE_VERSION: u8 = 7;

// the file is laid out as below, where u64 fields are
// little-endian, and the sequence starts, suffix array
// and LCP array hold `int_size` bytes per value
// in native byte order
//
//               offset  field
//                    0  version (u8)
//                    1  soft mask (u8, 0 fold, 1 preserve, 2 replace)
//                    2  soft mask replacement (u8, 0 unless replace)
//                    3  int_size (u8, 4 or 8)
//                    4  text_len (u64)
//                   12  text_pos (u64)
//                   20  suffix_array_pos (u64)
//                   28  lcp_pos (u64, 0 when there is no LCP array)
//                   36  num_suffixes (u64)
//                   44  num_sequences (u64)
//                   52  sequence starts (num_sequences values)
//                    .  seed mask length (u64)
//                    .  seed mask (u8 per position, 0/1)
//                    .  starts (32 bytes, bit b % 8 of byte b / 8 set
//                       when suffixes starting with byte b are indexed)
//             text_pos  text (text_len bytes)
//     suffix_array_pos  suffix array (num_suffixes values)
//              lcp_pos  LCP array (num_suffixes values), optional
//                    .  sequence names (bincode Vec<String>)
const VERSION: usize = 0;
const SOFT_MASK: usize = 1;
const SOFT_MASK_BYTE: usize = 2;
const INT_SIZE: usize = 3;
const TEXT_LEN: usize = 4;
const TEXT_POS: usize = 12;
const SUFFIX_ARRAY_POS: usize = 20;
const LCP_POS: usize = 28;
const NUM_SUFFIXES: usize = 36;
const NUM_SEQUENCES: usize = 44;
const SEQUENCE_STARTS: usize = 52;
const STARTS_LEN: usize = 32;

/// How the build writes lowercase (soft-masked) bytes into the text.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
#[non_exhaustive]
pub enum SoftMask {
    /// Upper-case them.
    #[default]
    Fold,

    /// Keep them as they are.
    Preserve,

    /// Write this byte over each of them.
    Replace(u8),
}

/// The header of a _.sufr_ file.
#[derive(Debug, PartialEq)]
pub struct Header {
    /// The serialization version.
    pub version: u8,

    /// How lowercase bytes were written into the text.
    pub soft_mask: SoftMask,

    /// The first bytes of the indexed suffixes besides the sentinel, or
    /// `None` when every suffix is indexed.
    pub starts: Option<Vec<u8>>,

    /// The length of the text.
    pub text_len: usize,

    /// The number of suffixes in the suffix array.
    pub num_suffixes: usize,

    /// The seed mask the suffixes were sorted by, if any.
    pub(crate) seed_mask: Option<SeedMask>,

    /// The start of each sequence in the text.
    pub sequence_starts: Vec<usize>,

    /// The name of each sequence.
    pub sequence_names: Vec<String>,

    /// The bytes per sequence start, suffix, and LCP value: 4 or 8.
    pub(crate) int_size: usize,

    /// The byte position of the text.
    pub(crate) text_pos: usize,

    /// The byte position of the suffix array.
    pub(crate) suffix_array_pos: usize,

    /// The byte position of the LCP array, if the file has one.
    pub(crate) lcp_pos: Option<usize>,
}

impl Header {
    /// Read the header of `path`.
    pub(crate) fn read(path: impl AsRef<Path>) -> Result<Header, OpenError> {
        let path = path.as_ref();
        let invalid = |reason: String| OpenError::Invalid {
            path: path.to_path_buf(),
            reason,
        };
        let read_error = |source| OpenError::io(path, source);
        let mut file = File::open(path).map_err(read_error)?;

        let mut fixed = [0u8; SEQUENCE_STARTS];
        file.read_exact(&mut fixed).map_err(read_error)?;
        let version = fixed[VERSION];
        if version != OUTFILE_VERSION {
            return Err(OpenError::Version {
                path: path.to_path_buf(),
                found: version,
            });
        }
        let int_size = fixed[INT_SIZE] as usize;
        if int_size != 4 && int_size != 8 {
            return Err(invalid(format!("integer size {int_size} is not 4 or 8")));
        }
        let u64_at =
            |pos: usize| u64::from_le_bytes(fixed[pos..pos + 8].try_into().unwrap());

        let text_len = u64_at(TEXT_LEN);
        if mem::size_of::<usize>() == 4 && text_len > usize::MAX as u64 {
            return Err(invalid(format!(
                "a text of {text_len} bytes is too long for this platform"
            )));
        }

        let num_sequences = u64_at(NUM_SEQUENCES) as usize;
        let mut starts = vec![0u8; num_sequences * int_size];
        file.read_exact(&mut starts).map_err(read_error)?;
        let sequence_starts = starts
            .chunks_exact(int_size)
            .map(|bytes| match int_size {
                4 => u32::from_ne_bytes(bytes.try_into().unwrap()) as usize,
                _ => u64::from_ne_bytes(bytes.try_into().unwrap()) as usize,
            })
            .collect();

        let mut buffer = [0u8; 8];
        file.read_exact(&mut buffer).map_err(read_error)?;
        let mut seed_mask = vec![0u8; u64::from_le_bytes(buffer) as usize];
        file.read_exact(&mut seed_mask).map_err(read_error)?;
        let seed_mask = if seed_mask.is_empty() {
            None
        } else {
            let mask = SeedMask::from_bytes(&seed_mask);
            Some(mask.ok_or_else(|| {
                invalid(format!("invalid seed mask bytes {seed_mask:?}"))
            })?)
        };

        let mut bitmap = [0u8; STARTS_LEN];
        file.read_exact(&mut bitmap).map_err(read_error)?;
        let starts: Vec<u8> = (0..=255u8)
            .filter(|&b| bitmap[b as usize / 8] & (1 << (b % 8)) != 0)
            .collect();
        let starts = (starts.len() < 256).then_some(starts);

        let soft_mask = match fixed[SOFT_MASK] {
            0 => SoftMask::Fold,
            1 => SoftMask::Preserve,
            2 => SoftMask::Replace(fixed[SOFT_MASK_BYTE]),
            other => return Err(invalid(format!("unknown soft mask {other}"))),
        };

        let lcp_pos = u64_at(LCP_POS) as usize;
        let mut header = Header {
            version: fixed[VERSION],
            soft_mask,
            starts,
            text_len: text_len as usize,
            num_suffixes: u64_at(NUM_SUFFIXES) as usize,
            seed_mask,
            sequence_starts,
            sequence_names: vec![],
            int_size,
            text_pos: u64_at(TEXT_POS) as usize,
            suffix_array_pos: u64_at(SUFFIX_ARRAY_POS) as usize,
            lcp_pos: (lcp_pos > 0).then_some(lcp_pos),
        };

        file.seek(SeekFrom::Start(header.names_pos() as u64))
            .map_err(read_error)?;
        let mut names = vec![];
        file.read_to_end(&mut names).map_err(read_error)?;
        header.sequence_names = bincode::deserialize(&names)
            .map_err(|e| invalid(format!("unreadable sequence names: {e}")))?;

        Ok(header)
    }

    /// Set the text, suffix array, and LCP positions from the
    /// header's contents, placing an LCP array when `write_lcp`.
    pub(crate) fn set_positions(&mut self, write_lcp: bool) {
        let seed_mask_len = self.seed_mask.as_ref().map_or(0, |m| m.len());
        self.text_pos = SEQUENCE_STARTS
            + self.sequence_starts.len() * self.int_size
            + mem::size_of::<u64>()
            + seed_mask_len
            + STARTS_LEN;
        self.suffix_array_pos = self.text_pos + self.text_len;
        let suffix_array_end =
            self.suffix_array_pos + self.num_suffixes * self.int_size;
        self.lcp_pos = write_lcp.then_some(suffix_array_end);
    }

    /// Whether the file has an LCP array.
    pub fn has_lcp(&self) -> bool {
        self.lcp_pos.is_some()
    }

    /// The seed mask the suffixes were sorted by, or `None` when they were
    /// sorted on whole suffixes.
    pub fn seed_mask(&self) -> Option<&SeedMask> {
        self.seed_mask.as_ref()
    }

    /// The number of characters the suffixes were sorted by: the text
    /// length for a full sort, the weight of a seed mask.
    pub(crate) fn sorted_len(&self) -> usize {
        self.seed_mask
            .as_ref()
            .map_or(self.text_len, |m| m.weight())
    }

    /// The byte position of the sequence names, after the arrays.
    pub(crate) fn names_pos(&self) -> usize {
        let num_arrays = if self.lcp_pos.is_some() { 2 } else { 1 };
        self.suffix_array_pos + num_arrays * self.num_suffixes * self.int_size
    }

    /// Write the header and the sequence names into `file`, leaving
    /// the text and arrays to the caller. Returns the file's length.
    pub(crate) fn write(&self, file: &File) -> io::Result<usize> {
        let seed_mask = self.seed_mask.as_ref().map_or(vec![], |m| m.to_bytes());

        let (soft_mask, soft_mask_byte) = match self.soft_mask {
            SoftMask::Fold => (0, 0),
            SoftMask::Preserve => (1, 0),
            SoftMask::Replace(b) => (2, b),
        };
        let mut bytes =
            vec![self.version, soft_mask, soft_mask_byte, self.int_size as u8];
        bytes.extend((self.text_len as u64).to_le_bytes());
        bytes.extend((self.text_pos as u64).to_le_bytes());
        bytes.extend((self.suffix_array_pos as u64).to_le_bytes());
        bytes.extend((self.lcp_pos.unwrap_or(0) as u64).to_le_bytes());
        bytes.extend((self.num_suffixes as u64).to_le_bytes());
        bytes.extend((self.sequence_starts.len() as u64).to_le_bytes());
        for &start in &self.sequence_starts {
            match self.int_size {
                4 => bytes.extend((start as u32).to_ne_bytes()),
                _ => bytes.extend((start as u64).to_ne_bytes()),
            }
        }
        bytes.extend((seed_mask.len() as u64).to_le_bytes());
        bytes.extend_from_slice(&seed_mask);
        let mut bitmap = [0u8; STARTS_LEN];
        for b in 0..=255u8 {
            if self
                .starts
                .as_ref()
                .is_none_or(|starts| starts.contains(&b))
            {
                bitmap[b as usize / 8] |= 1 << (b % 8);
            }
        }
        bytes.extend(bitmap);
        debug_assert_eq!(bytes.len(), self.text_pos);
        write_all_at(file, &bytes, 0)?;

        let names =
            bincode::serialize(&self.sequence_names).map_err(io::Error::other)?;
        write_all_at(file, &names, self.names_pos() as u64)?;

        Ok(self.names_pos() + names.len())
    }
}

/// A seed mask of care (1) and don't-care (0) positions.
///
/// The first and last positions must be a 1, and at
/// least one position must be a 0.
#[derive(Debug, PartialEq, Eq, Clone)]
pub struct SeedMask {
    /// The offsets of the care positions, ascending.
    positions: Vec<usize>,
}

impl SeedMask {
    /// The mask stored in a _.sufr_ header as one 1 or 0 byte per
    /// position.
    pub(crate) fn from_bytes(bytes: &[u8]) -> Option<Self> {
        let cares = bytes.iter().map(|&b| match b {
            1 => Some(true),
            0 => Some(false),
            _ => None,
        });
        Self::from_cares(cares)
    }

    fn from_cares(cares: impl Iterator<Item = Option<bool>>) -> Option<Self> {
        let cares: Vec<bool> = cares.collect::<Option<_>>()?;
        let valid = cares.first() == Some(&true)
            && cares.last() == Some(&true)
            && cares.contains(&false);
        valid.then(|| Self {
            positions: (0..cares.len()).filter(|&i| cares[i]).collect(),
        })
    }

    /// The offsets of the care positions, ascending.
    pub fn positions(&self) -> &[usize] {
        &self.positions
    }

    /// The number of care positions.
    pub fn weight(&self) -> usize {
        self.positions.len()
    }

    /// The number of positions, care and don't-care.
    pub(crate) fn len(&self) -> usize {
        self.positions[self.positions.len() - 1] + 1
    }

    /// The mask as one 1 or 0 byte per position.
    pub(crate) fn to_bytes(&self) -> Vec<u8> {
        let mut bytes = vec![0; self.len()];
        for &i in &self.positions {
            bytes[i] = 1;
        }
        bytes
    }
}

impl FromStr for SeedMask {
    type Err = InvalidSeedMask;

    fn from_str(mask: &str) -> Result<Self, InvalidSeedMask> {
        let cares = mask.bytes().map(|b| match b {
            b'1' => Some(true),
            b'0' => Some(false),
            _ => None,
        });
        Self::from_cares(cares).ok_or(InvalidSeedMask)
    }
}

impl Display for SoftMask {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            SoftMask::Fold => write!(f, "fold"),
            SoftMask::Preserve => write!(f, "preserve"),
            SoftMask::Replace(b) => write!(f, "replace with {}", *b as char),
        }
    }
}

impl Display for SeedMask {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for b in self.to_bytes() {
            write!(f, "{b}")?;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::{Header, SeedMask, SoftMask, OUTFILE_VERSION};
    use crate::error::OpenError;
    use pretty_assertions::assert_eq;
    use tempfile::NamedTempFile;

    type Result<T, E = Box<dyn std::error::Error>> = std::result::Result<T, E>;

    #[test]
    fn test_header_round_trip() -> Result<()> {
        let cases = [
            (
                4,
                Some("1101".parse()?),
                true,
                SoftMask::Replace(b'N'),
                None,
            ),
            (8, None, false, SoftMask::Preserve, Some(b"$ACGT".to_vec())),
            (4, None, false, SoftMask::Fold, Some(vec![])),
        ];
        for (int_size, seed_mask, write_lcp, soft_mask, starts) in cases {
            let mut header = Header {
                version: OUTFILE_VERSION,
                soft_mask,
                starts,
                text_len: 20,
                num_suffixes: 17,
                seed_mask,
                sequence_starts: vec![0, 7, 13],
                sequence_names: vec!["a".to_string(), "b".into(), "c".into()],
                int_size,
                text_pos: 0,
                suffix_array_pos: 0,
                lcp_pos: None,
            };
            header.set_positions(write_lcp);

            let outfile = NamedTempFile::new()?;
            let len = header.write(outfile.as_file())?;
            assert_eq!(len as u64, outfile.as_file().metadata()?.len());

            let path = outfile.path().to_string_lossy().to_string();
            assert_eq!(Header::read(&path)?, header);

            // the same file with a width that is neither 4 nor 8
            crate::io::write_all_at(outfile.as_file(), &[2], 3)?;
            let err = Header::read(&path).unwrap_err();
            assert!(
                matches!(&err, OpenError::Invalid { reason, .. } if reason == "integer size 2 is not 4 or 8"),
                "{err}"
            );

            // the same file with the next version in its first byte
            crate::io::write_all_at(outfile.as_file(), &[OUTFILE_VERSION + 1], 0)?;
            let err = Header::read(&path).unwrap_err();
            assert!(
                matches!(err, OpenError::Version { found, .. } if found == OUTFILE_VERSION + 1),
                "{err}"
            );

            // the original version and width, cut short inside the
            // sequence starts at 52
            crate::io::write_all_at(outfile.as_file(), &[OUTFILE_VERSION], 0)?;
            crate::io::write_all_at(outfile.as_file(), &[int_size as u8], 3)?;
            outfile.as_file().set_len(60)?;
            let err = Header::read(&path).unwrap_err();
            assert!(
                matches!(&err, OpenError::Invalid { reason, .. } if reason == "the file ends early"),
                "{err}"
            );
        }
        Ok(())
    }

    #[test]
    fn test_seed_mask() -> Result<()> {
        let cases: [(&str, &[usize]); 4] = [
            ("101", &[0, 2]),
            ("11011", &[0, 1, 3, 4]),
            ("110110101", &[0, 1, 3, 4, 6, 8]),
            ("11101101101000011", &[0, 1, 2, 4, 5, 7, 8, 10, 15, 16]),
        ];
        for (text, positions) in cases {
            let mask: SeedMask = text.parse()?;
            assert_eq!(mask.positions(), positions);
            assert_eq!(mask.weight(), positions.len());
            assert_eq!(mask.len(), text.len());
            assert_eq!(mask.to_string(), text);

            let bytes = mask.to_bytes();
            let expected: Vec<u8> = text.bytes().map(|b| b - b'0').collect();
            assert_eq!(bytes, expected);
            assert_eq!(SeedMask::from_bytes(&bytes), Some(mask));
        }
        Ok(())
    }

    #[test]
    fn test_invalid_seed_mask() {
        let invalid = [
            "", "abc", "1", "11", "111", "1111", "0", "00", "0101", "0111", "1010",
            "11100", "1021", "1a01",
        ];
        for text in invalid {
            assert!(text.parse::<SeedMask>().is_err(), "{text}");
        }

        let invalid: [&[u8]; 6] = [&[], &[1], &[1, 1], &[0, 1], &[1, 0], &[1, 2, 1]];
        for bytes in invalid {
            assert!(SeedMask::from_bytes(bytes).is_none(), "{bytes:?}");
        }
    }
}
