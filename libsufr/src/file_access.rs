use anyhow::{bail, Result};
use crate::{
    util::slice_u8_to_vec,
    types::{FromUsize, Int}
};
use std::{
    mem,
    cmp::min,
    fs::File,
    io::{Read, Seek, SeekFrom},
    ops::Range,
};

// --------------------------------------------------
#[derive(Debug)]
pub struct FileAccess<T>
where
    T: Int + FromUsize<T> + Sized + Send + Sync + serde::ser::Serialize,
{
    file: File,
    buffer: Vec<T>,
    buffer_size: usize,
    buffer_pos: usize,
    pub size: usize,
    start_position: u64,
    current_position: u64,
    end_position: u64,
    exhausted: bool,
}

impl<T> FileAccess<T>
where
    T: Int + FromUsize<T> + Sized + Send + Sync + serde::ser::Serialize,
{
    pub fn new(filename: &str, start: u64, num_elements: usize) -> Result<Self> {
        let file = File::open(filename)?;
        let size = num_elements * mem::size_of::<T>();
        Ok(FileAccess {
            file,
            buffer: vec![],
            buffer_size: 2usize.pow(30),
            buffer_pos: 0,
            size,
            start_position: start,
            current_position: start,
            end_position: start + size as u64,
            exhausted: false,
        })
    }

    pub fn reset(&mut self) {
        self.buffer = vec![];
        self.buffer_pos = 0;
        self.current_position = self.start_position;
        self.exhausted = false;
    }

    pub fn iter(&mut self) -> FileAccessIter<T> {
        FileAccessIter { file_access: self }
    }

    // --------------------------------------------------
    // TODO: Ignoring lots of Results to return Option
    pub fn get(&mut self, pos: usize) -> Option<T> {
        // Don't bother looking for something beyond the end
        let seek = self.start_position + (pos * mem::size_of::<T>()) as u64;
        if seek < self.end_position {
            let _ = self.file.seek(SeekFrom::Start(seek));
            let mut buffer: Vec<u8> = vec![0; mem::size_of::<T>()];
            let bytes_read = self.file.read(&mut buffer).unwrap();
            (bytes_read == mem::size_of::<T>()).then(|| {
                let res = unsafe {
                    std::slice::from_raw_parts(buffer.as_ptr() as *const _, 1)
                };
                res[0]
            })
        } else {
            None
        }
    }

    // --------------------------------------------------
    pub fn get_range(&mut self, range: Range<usize>) -> Result<Vec<T>> {
        let start = self.start_position as usize + (range.start * mem::size_of::<T>());
        let end = self.start_position as usize + (range.end * mem::size_of::<T>());
        let valid = self.start_position as usize..self.end_position as usize + 1;
        if valid.contains(&start) && valid.contains(&end) {
            self.file.seek(SeekFrom::Start(start as u64))?;
            let mut buffer: Vec<u8> = vec![0; end - start];
            let bytes_read = self.file.read(&mut buffer)?;
            let num_vals = bytes_read / mem::size_of::<T>();
            Ok(slice_u8_to_vec(&buffer, num_vals))
        } else {
            bail!("Invalid range: {range:?}")
        }
    }
}

// --------------------------------------------------
#[derive(Debug)]
pub struct FileAccessIter<'a, T>
where
    T: Int + FromUsize<T> + Sized + Send + Sync + serde::ser::Serialize,
{
    file_access: &'a mut FileAccess<T>,
}

impl<T> Iterator for FileAccessIter<'_, T>
where
    T: Int + FromUsize<T> + Sized + Send + Sync + serde::ser::Serialize,
{
    type Item = T;

    fn next(&mut self) -> Option<Self::Item> {
        if self.file_access.exhausted {
            None
        } else {
            // Fill the buffer
            if self.file_access.buffer.is_empty()
                || self.file_access.buffer_pos == self.file_access.buffer.len()
            {
                if self.file_access.current_position >= self.file_access.end_position {
                    self.file_access.exhausted = true;
                    return None;
                }

                self.file_access
                    .file
                    .seek(SeekFrom::Start(self.file_access.current_position))
                    .unwrap();

                let bytes_wanted = min(
                    self.file_access.buffer_size * mem::size_of::<T>(),
                    (self.file_access.end_position - self.file_access.current_position)
                        as usize,
                );

                let mut buffer: Vec<u8> = vec![0; bytes_wanted];
                self.file_access.file.read_exact(&mut buffer).unwrap();
                self.file_access.current_position =
                    self.file_access.file.stream_position().unwrap();

                let num_vals = bytes_wanted / mem::size_of::<T>();
                self.file_access.buffer = slice_u8_to_vec(&buffer, num_vals);
                self.file_access.buffer_pos = 0;
            }

            let val = self
                .file_access
                .buffer
                .get(self.file_access.buffer_pos)
                .copied();

            self.file_access.buffer_pos += 1;
            val
        }
    }
}
