//! # Low memory access to Sufr's on-disk arrays (text/SA/LCP)

use crate::{types::Int, util::slice_int_to_slice_u8};
use anyhow::{bail, Result};
use std::{cmp::min, fs::File, io, marker::PhantomData, mem, ops::Range};

// --------------------------------------------------
/// Struct to mediate file access to on-disk arrays of text, suffix/LCP arrays
#[derive(Debug)]
pub struct FileAccess<T: Int> {
    /// A read-only filehandle to the _.sufr_ file
    file: File,

    /// The size in bytes for the entire text/SA/LCP
    pub size: usize,

    /// The starting byte position of the structure being read (text/SA/LCP)
    start_position: u64,

    /// The final byte position of the structure being read (text/SA/LCP)
    end_position: u64,

    _marker: PhantomData<fn() -> T>,
}

impl<T: Int> FileAccess<T> {
    /// Create a read-only file access to a portion of a _.sufr_ file
    /// representing the text, suffix array, or LCP array.
    /// This struct must be initialized using an `Int` of `u8` for the `text`
    /// or `u32`/`u64` for the SA/LCP.
    /// The metadata needed to create this can be found in the `SufrFile`.
    ///
    /// Args:
    /// * `filename`: the _.sufr_ filename
    /// * `start`: the byte position in the file of the array
    /// * `num_elements`: the length of the text/SA/LCP
    pub fn new(filename: &str, start: u64, num_elements: usize) -> Result<Self> {
        let file = File::open(filename)?;
        let size = num_elements * mem::size_of::<T>();
        Ok(FileAccess {
            file,
            size,
            start_position: start,
            end_position: start + size as u64,
            _marker: PhantomData,
        })
    }

    /// Create a `FileAccessIter` iterator.
    pub fn iter(&self) -> FileAccessIter<'_, T> {
        FileAccessIter {
            file_access: self,
            buffer: vec![],
            buffer_pos: 0,
            current_position: self.start_position,
            exhausted: false,
        }
    }

    // --------------------------------------------------
    /// Return a value (`u8`/character from text or a SA/LCP value)
    ///
    /// Args:
    /// * `pos`: position in the array
    //
    // TODO: Ignoring lots of Results to return Option
    pub fn get(&self, pos: usize) -> Option<T> {
        // Don't bother looking for something beyond the end
        let seek = self.start_position + (pos * mem::size_of::<T>()) as u64;
        if seek < self.end_position {
            let mut val = T::from_usize(0);
            read_exact_at(
                &self.file,
                slice_int_to_slice_u8(std::slice::from_mut(&mut val)),
                seek,
            )
            .unwrap();
            Some(val)
        } else {
            None
        }
    }

    // --------------------------------------------------
    /// Return a range of values (`u8`/characters from text or a SA/LCP values)
    ///
    /// Args:
    /// * `range`: start/stop positions in the array
    pub fn get_range(&self, range: Range<usize>) -> Result<Vec<T>> {
        assert!(range.start <= range.end);
        let start = self.start_position as usize + (range.start * mem::size_of::<T>());
        let end = self.start_position as usize + (range.end * mem::size_of::<T>());
        let valid = self.start_position as usize..self.end_position as usize + 1;
        if valid.contains(&start) && valid.contains(&end) {
            let num_vals = range.len();
            let mut buffer: Vec<T> = vec![T::from_usize(0); num_vals];
            read_exact_at(
                &self.file,
                slice_int_to_slice_u8(&mut buffer),
                start as u64,
            )?;
            Ok(buffer)
        } else {
            bail!("Invalid range: {range:?}")
        }
    }
}

// --------------------------------------------------
/// An iterator over the values from a `FileAccess`
#[derive(Debug)]
pub struct FileAccessIter<'a, T: Int> {
    file_access: &'a FileAccess<T>,
    /// Internal buffer for reading a portion of the file
    buffer: Vec<T>,
    /// The current position when reading the buffer
    buffer_pos: usize,
    /// The current position after reading a portion of the structure
    /// from disk and into the `buffer`
    current_position: u64,
    /// Whether or not the last read off disk made it to the end of the array
    exhausted: bool,
}

impl<'a, T: Int> FileAccessIter<'a, T> {
    /// The maximum size in bytes of the buffer (currently 2^24, 16 MiB)
    const BUFFER_SIZE: usize = 2usize.pow(24);
}

impl<T: Int> Iterator for FileAccessIter<'_, T> {
    type Item = T;

    fn next(&mut self) -> Option<Self::Item> {
        if self.exhausted {
            None
        } else {
            // Fill the buffer
            if self.buffer.is_empty() || self.buffer_pos == self.buffer.len() {
                if self.current_position >= self.file_access.end_position {
                    self.exhausted = true;
                    return None;
                }

                // Read whole elements, at most BUFFER_SIZE bytes, straight
                // into the reused typed buffer
                let bytes_wanted = min(
                    Self::BUFFER_SIZE,
                    (self.file_access.end_position - self.current_position) as usize,
                );
                let num_vals = bytes_wanted / mem::size_of::<T>();
                self.buffer.resize(num_vals, T::from_usize(0));
                read_exact_at(
                    &self.file_access.file,
                    slice_int_to_slice_u8(&mut self.buffer),
                    self.current_position,
                )
                .unwrap();

                self.current_position += (num_vals * mem::size_of::<T>()) as u64;
                self.buffer_pos = 0;
            }

            let val = self.buffer.get(self.buffer_pos).copied();

            self.buffer_pos += 1;
            val
        }
    }
}

// Position-independent reads and writes for both Unix and Windows

/// Write from `buf` at an absolute `offset` without using the file cursor,
/// so several threads may write disjoint regions of one file at once.
/// Returns the number of bytes written, which may be fewer than requested.
#[cfg(unix)]
pub(crate) fn write_at(file: &File, buf: &[u8], offset: u64) -> io::Result<usize> {
    use std::os::unix::fs::FileExt;

    file.write_at(buf, offset)
}

// Windows version of `write_at`; `seek_write` also takes an absolute
// offset per call, so concurrent calls place their data correctly.
#[cfg(windows)]
pub(crate) fn write_at(file: &File, buf: &[u8], offset: u64) -> io::Result<usize> {
    use std::os::windows::fs::FileExt;

    file.seek_write(buf, offset)
}

#[cfg(unix)]
fn read_exact_at(file: &File, buf: &mut [u8], offset: u64) -> io::Result<()> {
    use std::os::unix::fs::FileExt;

    file.read_exact_at(buf, offset)
}

// Essentially the same as the Unix read_exact_at implementation, but using seek_read
#[cfg(windows)]
fn read_exact_at(file: &File, mut buf: &mut [u8], mut offset: u64) -> io::Result<()> {
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
