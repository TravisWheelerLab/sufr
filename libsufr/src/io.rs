//! Positioned file reads and writes, and the byte casts of the arrays

use anyhow::{bail, Result};
use rayon::prelude::*;
use std::{
    fmt::Debug,
    fs::File,
    io,
    ops::{Add, Sub},
    slice,
};

mod sealed {
    pub trait Sealed {}
    impl Sealed for u8 {}
    impl Sealed for u32 {}
    impl Sealed for u64 {}
}

/// An element of a _.sufr_ array: `u8` for the text, and `u32` or `u64`
/// for the suffix and LCP arrays.
///
/// The trait is sealed because the byte casts in this module rely on
/// `u8`, `u32` and `u64` having no padding and no invalid bit patterns.
pub trait Int:
    Copy
    + Default
    + Ord
    + Debug
    + Send
    + Sync
    + 'static
    + Add<Output = Self>
    + Sub<Output = Self>
    + Into<u64>
    + sealed::Sealed
{
    /// The value as a `usize`.
    fn to_usize(self) -> usize {
        self.into() as usize
    }

    /// The value of `val`, truncated to `Self`.
    fn from_usize(val: usize) -> Self;
}

impl Int for u8 {
    fn from_usize(val: usize) -> Self {
        val as Self
    }
}

impl Int for u32 {
    fn from_usize(val: usize) -> Self {
        val as Self
    }
}

impl Int for u64 {
    fn from_usize(val: usize) -> Self {
        val as Self
    }
}

/// The bytes of `buffer`, in native byte order.
pub(crate) fn as_bytes_mut<T: Int>(buffer: &mut [T]) -> &mut [u8] {
    // TODO: rework - endianness-dependent
    // SAFETY: T is u8, u32 or u64, which have no padding and
    //         accept every bit pattern, and the byte length is
    //         the length of `buffer`
    unsafe {
        slice::from_raw_parts_mut(
            buffer.as_mut_ptr() as *mut _,
            std::mem::size_of_val(buffer),
        )
    }
}

/// The bytes of `vec`, in native byte order.
pub(crate) fn as_bytes<T: Int>(vec: &[T]) -> &[u8] {
    // TODO: rework - endianness-dependent
    // SAFETY: T is u8, u32 or u64, which have no padding, and
    //         the byte length is the length of `vec`
    unsafe {
        slice::from_raw_parts(vec.as_ptr() as *const _, std::mem::size_of_val(vec))
    }
}

// the positioned reads and writes below take an absolute offset
// and leave the file cursor alone, so several threads may use
// one file at once

/// Write from `buf` at `offset`, returning the number of bytes written,
/// which may be fewer than requested.
#[cfg(unix)]
fn write_at(file: &File, buf: &[u8], offset: u64) -> io::Result<usize> {
    use std::os::unix::fs::FileExt;

    file.write_at(buf, offset)
}

#[cfg(windows)]
fn write_at(file: &File, buf: &[u8], offset: u64) -> io::Result<usize> {
    use std::os::windows::fs::FileExt;

    file.seek_write(buf, offset)
}

#[cfg(unix)]
pub(crate) fn read_exact_at(
    file: &File,
    buf: &mut [u8],
    offset: u64,
) -> io::Result<()> {
    use std::os::unix::fs::FileExt;

    file.read_exact_at(buf, offset)
}

// the loop of std's unix read_exact_at, over seek_read
#[cfg(windows)]
pub(crate) fn read_exact_at(
    file: &File,
    mut buf: &mut [u8],
    mut offset: u64,
) -> io::Result<()> {
    use std::os::windows::fs::FileExt;

    while !buf.is_empty() {
        match file.seek_read(buf, offset) {
            Ok(0) => break,
            Ok(n) => {
                buf = &mut buf[n..];
                offset += n as u64;
            }
            Err(ref e) if e.kind() == io::ErrorKind::Interrupted => {}
            Err(e) => return Err(e),
        }
    }
    if !buf.is_empty() {
        Err(io::Error::new(
            io::ErrorKind::UnexpectedEof,
            "failed to fill whole buffer",
        ))
    } else {
        Ok(())
    }
}

/// Write all of `buf` at `offset`, in parallel 64 MB chunks.
pub(crate) fn write_all_at(file: &File, buf: &[u8], offset: u64) -> Result<()> {
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

#[cfg(test)]
mod tests {
    use super::{as_bytes, as_bytes_mut};
    use anyhow::Result;
    use pretty_assertions::assert_eq;

    #[test]
    fn test_as_bytes_mut() -> Result<()> {
        let mut buf: Vec<u32> = vec![1];
        as_bytes_mut(&mut buf).copy_from_slice(&[0, 0, 0, 0]);
        assert_eq!(buf, &[0u32]);

        let mut buf: Vec<u64> = vec![1];
        as_bytes_mut(&mut buf).copy_from_slice(&[0, 0, 0, 0, 0, 0, 0, 0]);
        assert_eq!(buf, &[0u64]);

        let mut buf: Vec<u32> = vec![1];
        as_bytes_mut(&mut buf).copy_from_slice(&[1, 0, 0, 0]);
        assert_eq!(buf, &[1u32]);

        let mut buf: Vec<u64> = vec![1];
        as_bytes_mut(&mut buf).copy_from_slice(&[1, 0, 0, 0, 0, 0, 0, 0]);
        assert_eq!(buf, &[1u64]);

        let mut buf: Vec<u32> = vec![1];
        as_bytes_mut(&mut buf).copy_from_slice(&[255, 255, 255, 255]);
        assert_eq!(buf, &[u32::MAX]);

        let mut buf: Vec<u64> = vec![1];
        as_bytes_mut(&mut buf)
            .copy_from_slice(&[255, 255, 255, 255, 255, 255, 255, 255]);
        assert_eq!(buf, &[u64::MAX]);

        let mut buf: Vec<u32> = vec![1, 1];
        as_bytes_mut(&mut buf).copy_from_slice(&[0, 0, 0, 0, 255, 255, 255, 255]);
        assert_eq!(buf, &[0u32, u32::MAX]);

        let mut buf: Vec<u64> = vec![1, 1];
        as_bytes_mut(&mut buf).copy_from_slice(&[
            0, 0, 0, 0, 0, 0, 0, 0, 255, 255, 255, 255, 255, 255, 255, 255,
        ]);
        assert_eq!(buf, &[0u64, u64::MAX]);

        Ok(())
    }

    #[test]
    fn test_as_bytes() -> Result<()> {
        let res = as_bytes(&[0u32]);
        assert_eq!(res, &[0, 0, 0, 0]);

        let res = as_bytes(&[0u64]);
        assert_eq!(res, &[0, 0, 0, 0, 0, 0, 0, 0]);

        let res = as_bytes(&[1u32]);
        assert_eq!(res, &[1, 0, 0, 0]);

        let res = as_bytes(&[1u64]);
        assert_eq!(res, &[1, 0, 0, 0, 0, 0, 0, 0]);

        let res = as_bytes(&[u32::MAX]);
        assert_eq!(res, &[255, 255, 255, 255]);

        let res = as_bytes(&[u64::MAX]);
        assert_eq!(res, &[255, 255, 255, 255, 255, 255, 255, 255]);

        let res = as_bytes(&[0u32, u32::MAX]);
        assert_eq!(res, &[0, 0, 0, 0, 255, 255, 255, 255]);

        let res = as_bytes(&[0u64, u64::MAX]);
        assert_eq!(
            res,
            &[0, 0, 0, 0, 0, 0, 0, 0, 255, 255, 255, 255, 255, 255, 255, 255]
        );

        Ok(())
    }
}
