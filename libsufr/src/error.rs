//! The errors of building, opening and querying a suffix array

use crate::format::OUTFILE_VERSION;
use std::{
    error::Error,
    fmt::{self, Display},
    io,
    path::{Path, PathBuf},
};

/// An error building a suffix array from a sequence file.
#[derive(Debug)]
#[non_exhaustive]
pub enum BuildError {
    /// The sequence file could not be read or parsed as FASTA/FASTQ.
    Input {
        /// The sequence file.
        path: PathBuf,
        /// The parser's error.
        //
        // boxed, so the parser's error type stays out of the API
        source: Box<dyn Error + Send + Sync>,
    },

    /// The output file could not be created or written.
    Io {
        /// The output file.
        path: PathBuf,
        /// The I/O error.
        source: io::Error,
    },

    /// The seed mask has more care positions or spans more characters than
    /// a key of this text holds.
    MaskLength {
        /// The most care positions a mask may have.
        max_weight: usize,
        /// The most characters a mask may span.
        max_span: usize,
    },
}

impl Display for BuildError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            BuildError::Input { path, source } => {
                write!(f, "{}: {source}", path.display())
            }
            BuildError::Io { path, source } => {
                write!(f, "{}: {source}", path.display())
            }
            BuildError::MaskLength {
                max_weight,
                max_span,
            } => write!(
                f,
                "a seed mask for this text may have at most {max_weight} care \
                 positions and span at most {max_span} characters"
            ),
        }
    }
}

impl Error for BuildError {
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        match self {
            BuildError::Input { source, .. } => Some(source.as_ref()),
            BuildError::Io { source, .. } => Some(source),
            BuildError::MaskLength { .. } => None,
        }
    }
}

/// An error opening a _.sufr_ file.
#[derive(Debug)]
#[non_exhaustive]
pub enum OpenError {
    /// The file, or its subsample cache, could not be read or written.
    Io {
        /// The _.sufr_ file or the cache file.
        path: PathBuf,
        /// The I/O error.
        source: io::Error,
    },

    /// The file was written in a format version this build does not read.
    Version {
        /// The _.sufr_ file.
        path: PathBuf,
        /// The version in the file.
        found: u8,
    },

    /// The file is not a valid _.sufr_ file of this version.
    Invalid {
        /// The _.sufr_ file.
        path: PathBuf,
        /// What is wrong with the file.
        reason: String,
    },
}

impl OpenError {
    /// `source` as an error opening `path`: `Invalid` when the file ends
    /// before a read, else `Io`.
    pub(crate) fn io(path: &Path, source: io::Error) -> Self {
        match source.kind() {
            io::ErrorKind::UnexpectedEof => OpenError::Invalid {
                path: path.to_path_buf(),
                reason: "the file ends early".to_string(),
            },
            _ => OpenError::Io {
                path: path.to_path_buf(),
                source,
            },
        }
    }
}

impl Display for OpenError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            OpenError::Io { path, source } => write!(f, "{}: {source}", path.display()),
            OpenError::Version { path, found } => write!(
                f,
                "{}: format version {found}, and this build reads only version \
                 {OUTFILE_VERSION}",
                path.display()
            ),
            OpenError::Invalid { path, reason } => {
                write!(f, "{}: {reason}", path.display())
            }
        }
    }
}

impl Error for OpenError {
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        match self {
            OpenError::Io { source, .. } => Some(source),
            OpenError::Version { .. } | OpenError::Invalid { .. } => None,
        }
    }
}

/// An error reading from an open suffix array.
#[derive(Debug)]
#[non_exhaustive]
pub enum QueryError {
    /// The file could not be read.
    Io {
        /// The _.sufr_ file.
        path: PathBuf,
        /// The I/O error.
        source: io::Error,
    },

    /// The ranks or text positions asked for are outside the suffix array
    /// or the text.
    OutOfRange {
        /// The number of suffixes, or the length of the text.
        len: usize,
    },

    /// The LCPs were asked for, and the suffix array has none.
    NoLcp,
}

impl Display for QueryError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            QueryError::Io { path, source } => {
                write!(f, "{}: {source}", path.display())
            }
            QueryError::OutOfRange { len } => {
                write!(f, "the range asked for is outside 0..{len}")
            }
            QueryError::NoLcp => write!(f, "the suffix array has no LCP array"),
        }
    }
}

impl Error for QueryError {
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        match self {
            QueryError::Io { source, .. } => Some(source),
            QueryError::OutOfRange { .. } | QueryError::NoLcp => None,
        }
    }
}

/// A seed mask string that is not 1s and 0s, starting and ending with 1,
/// with at least one 0.
#[derive(Debug, PartialEq, Eq)]
#[non_exhaustive]
pub struct InvalidSeedMask;

impl Display for InvalidSeedMask {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "a seed mask must be 1s and 0s, start and end with 1, and hold \
             at least one 0"
        )
    }
}

impl Error for InvalidSeedMask {}
