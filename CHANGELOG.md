# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## Unreleased

### Added

- Start `CHANGELOG.md` for changes made after `0.7.12`
- Add `libsufr::SuffixArrayBuilder`, a builder-style API for creating a
  suffix array from a sequence file (`from_sequence_file`) or from in-memory
  text (`from_data`)
- Add `libsufr::Alphabet` (`Dna`, `Amino`, `Bytes`) and
  `SuffixArrayBuilder::with_alphabet` to choose how the text is read
- Add `SuffixArrayBuilder::with_memory_budget` and `sufr create
  --memory-budget`, the most bytes of suffixes held in memory while
  sorting. The default is the length of the text. A text whose suffixes
  need more is scanned once per budget.
- Add `SuffixArray::open` and `SuffixArray::open_with`. `OpenOptions` and
  `SourceKind` choose whether the text and the suffix array are held in
  memory or read from disk.
- Add `SuffixArray::search`, `search_all` and `par_search_all`, which
  return `Hits` for each query. `Hits` gives the count, the ranks, the
  suffixes, `locate` and `extract`.
- Add `Header` and `SuffixArray::header`, `suffixes`, `lcps` and `text`,
  which read a file's header and its arrays by range. The byte positions
  of the arrays in the header are not public.
- Add the cargo feature `typed`, which exports `Int`, `TypedSuffixArray`
  and the storage types `Source`, `MemorySource` and `DiskSource`
- Add the example `libsufr/examples/search.rs`
- Re-export `SuffixArray` at the crate root as `libsufr::SuffixArray`
- Implement vectorized LCP calculation (AVX2 on x86-64, NEON on AArch64), significantly
  improving build time on some inputs
- Implement an LCP cache, significantly improving build time on some inputs with very long
  duplications or tandem repeats
- Add `SuffixArrayBuilder::with_write_lcp` and `sufr create --write-lcp`.
  The LCP array is now optional and off by default. `list --show-lcp`
  fails on a file without one. The header's LCP offset is 0 when there is
  no LCP array, and `Header::has_lcp` and `summarize` report whether one
  is present. A search with a query-time max query
  length reads the LCP array, when the file has one, to subsample the
  suffix array.
- Subsample the suffix array for a query-time max query length without an
  LCP array by streaming the suffix array in blocks of four million
  elements and comparing each suffix with its predecessor in the text

### Changed

- (breaking) Remove `SufrMetadata`. `SuffixArray::header` gives the facts
  about the file that it held. `sufr summarize` reads the file's size and
  time itself.
- (breaking) `count`, `extract` and `locate` read queries from a file
  only with `-q,--query-file <FILE>`, and `-q -` reads them from stdin.
  A positional argument is always a query, even when a file by that name
  exists. A missing query file, or no queries at all, is an error.
- (breaking) `-t,--threads`, `--log` and `--log-file` follow the
  subcommand, as in `sufr create -t 16`, and every subcommand takes them.
  `--log` no longer has the short form `-l`.
- `sufr summarize` no longer prints a "Max query len" row. A file built
  with a seed mask still shows its "Seed mask" row.
- (breaking) Simplify the `Int` trait definition and bounds, subsuming
  `FromUsize`. `Int` is exported only under the `typed` feature.
- (breaking) Make `Int` unable to implement from other crates, due to guarantees required by sufr's `unsafe` code
- (breaking) Change the partitioning strategy to be based on a prefix of
  each suffix, instead of using random sampling. A first pass counts the
  suffixes of each prefix, so the size of every partition and its place
  in the output are known before sorting.
- Sort each partition by a 64-bit key holding the first characters of
  each suffix, and merge sort only the suffixes whose keys are equal.
  Under a seed mask the key holds the care positions.
- Change the default for `sufr create --num-partitions` from 16 to 1024.
  It caps the number of sort tasks. `--memory-budget` sets the memory a
  build uses.
- Construction writes no temporary files. It holds the text and up to the
  memory budget of suffixes in memory. Peak memory is higher: on 1 GiB
  of DNA at 8 threads it went from 1.4 GiB to 2.4 GiB, measured when the
  change was made.
- Construction is faster. In our measurements on 1 GiB of DNA at 8
  threads, a build went from 57.7 s to 8.8 s, and inputs with long
  repeats went from minutes to seconds.
- Reduce some unnecessary internal allocations
- During build, the merge sort makes better use of the available threads
  with work-stealing via `rayon::join`
- Write the output file with positioned writes (`pwrite` on Unix,
  `seek_write` on Windows) in parallel 64 MiB chunks. The header and text
  are written first, and each sorted partition is written, with its LCPs
  if any, straight to its place in the file.
- `sufr create` refuses a seed mask with more care positions than one
  64-bit key holds, or one spanning more than four keys of characters.
  A key typically holds 21 characters of DNA.
- (breaking) `sufr create -i` requires `--dna`. It no longer converts
  lowercase protein to `X`.
- (breaking) Raise the minimum supported Rust version from 1.80 to 1.88
- Byte 2 of the header records whether suffixes starting with a byte
  other than A, C, G or T were indexed. It is 1 for all but DNA input.
  A non-DNA file differs from before in that byte. Readers ignore it.
- `sufr locate` orders the hits by the order of the sequences in the
  input, not by their names
- An `extract` window holds only its sequence's own bytes. It no longer
  includes the sequence delimiter or the sentinel, and no longer runs into
  the next sequence.

### Removed

- (breaking) Remove the `sufr_builder` module, `SufrBuilder`,
  `SufrBuilderArgs`, `SuffixArray::new`, and `SuffixArray::write`. Use
  `libsufr::SuffixArrayBuilder` instead.
- (breaking) Remove `SuffixArray::read`, `count`, `extract`, `locate`,
  `list` and `string_at`, with their option and result structs
  (`CountOptions`, `CountResult`, `ExtractOptions`, `ExtractResult`,
  `ExtractSequence`, `LocateOptions`, `LocateResult`, `LocatePosition`,
  `ListOptions`, `SearchOptions`, `SearchResult`,
  `SearchResultLocations`). Use `SuffixArray::search` and `Hits` instead.
- (breaking) Remove `SufrFile`. Use `SuffixArray::open_with` instead.
- (breaking) Remove the public modules. `SeedMask`, `SuffixSortType` and
  `OUTFILE_VERSION` are at the crate root.
  `SENTINEL_CHARACTER`, `Comparison`, `SequenceFileData`,
  `read_sequence_file`, `read_text_length` and `find_lcp_full_offset` are
  private or gone.
- (breaking) Remove the build-time max query length, `sufr create -m`.
  Opening a file built with one is an error. The query-time `-m` of
  `count`, `locate` and `extract` stays.
- Remove `sufr create -r,--random-seed`,
  which is not used for the prefix-based partitioning
- (breaking) Remove `sufr create -a,--allow-ambiguity` and the "Allow
  Ambiguity" row of `sufr summarize`. The LCP cache made
  them unnecessary. A DNA build indexes only suffixes starting with A,
  C, G or T.
- Remove and make private the `unsafe` byte-reading utilities
  `slice_u8_to_vec`, `usize_to_bytes`, `vec_to_slice_u8`.

### Fixed

- `sufr` exits 0 with nothing on stderr when the reader of its output
  stops early, as `head` does. `summarize` panicked there, and `list`,
  `locate` and `extract` reported "Broken pipe" and exited 1.
- Disk I/O is updated to reduce intermediate allocations and remove unaligned memory reads/writes
- Searches (count, extract, or locate) with a `max_query_len`/`--max-query-len` could create cache
  files containing extraneous garbage memory; they no longer do so. The extra data was never read,
  so existing cache files should still be readable/writable across versions.
- A query that sorts before every suffix no longer subtracts below rank
  0. Debug builds panicked. Release builds found no hits by chance.
- A query that matches a whole suffix, sentinel included, and has bytes
  left no longer panics. It finds nothing.
- Under a query-time max query length, a hit over several groups of the
  subsampled suffix array no longer loses most of its last group. A hit
  on the last group alone no longer panics. On `smol.fa`, `count A` gave
  93 without a cap and 71 with `-m 2`.
- `sufr list -o` writes to the named file
- `count`, `locate` and `extract` stop with an error when a read fails,
  instead of skipping the query
