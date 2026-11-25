# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## Unreleased

### Added

- Start `CHANGELOG.md` for changes made after `0.7.12`.
- Add `libsufr::SuffixArrayBuilder`, a builder-style API for creating a
  suffix array from a sequence file (`from_sequence_file`).
- Add `SoftMask` and `SuffixArrayBuilder::with_soft_mask`, which set how
  the build writes lowercase bytes into the text. `Header::soft_mask`
  records the setting, and `SoftMask` is `#[non_exhaustive]`.
- Add `SuffixArrayBuilder::with_starts`, which indexes only the suffix
  at the sentinel and those starting with the given bytes.
  `Header::starts` records the bytes.
- Add `sufr create --soft-mask` and `--starts`. `--dna` is short for
  `--starts ACGT`, and `-i` for `--soft-mask N`.
- Add `SuffixArrayBuilder::with_memory_budget` and
  `sufr create --memory-budget` (in MiB), which cap the memory for
  suffixes while sorting. The default is the length of the text; with a
  smaller budget, the build takes more passes over the text.
- Add `SuffixArray::open` and `SuffixArray::open_with`. `OpenOptions`
  sets whether the text and the suffix array are held in memory or read
  from disk (`SourceKind`), and the max query length.
- Add `SuffixArray::search`, `search_all` and `par_search_all`, which
  return `Hits` for each query. `Hits` gives the count, the ranks, the
  suffixes, `locate` and `extract`. `Location` and `Extraction` are
  `#[non_exhaustive]`.
- Add `Header` and `SuffixArray::header`, `suffixes`, `lcps` and `text`,
  which return a file's header, and its suffixes, LCPs and text by
  range.
- Add the cargo feature `unstable` for low-level API that may break in
  any release: `Int`, `TypedSuffixArray` and its storage types,
  `SuffixArrayBuilder::from_data` and `with_force_64bit`, and
  `From<TypedSuffixArray>` for `SuffixArray`.
- Add the example `libsufr/examples/search.rs`.
- Add `SuffixArrayBuilder::with_write_lcp` and
  `sufr create --write-lcp`, which write the LCP array.
  `Header::has_lcp` and `summarize` report whether a file has one.

### Changed

- (breaking) `count`, `extract` and `locate` read queries from a file
  only with `-q,--query-file <FILE>`, and `-q -` reads them from stdin.
  A positional argument is always a query, even when a file by that name
  exists. A missing query file, or no queries at all, is an error.
- (breaking) `-t,--threads`, `--log` and `--log-file` follow the
  subcommand, as in `sufr create -t 16`, and every subcommand takes
  them. `--log` no longer has the short form `-l`.
- (breaking) `sufr create` and `SuffixArrayBuilder` write the LCP array
  only when asked with `--write-lcp` or `with_write_lcp`.
  `list --show-lcp` fails on a file without one; searches with a max
  query length still work.
- (breaking) `Int` now includes `FromUsize`, is sealed, and is exported
  only under `unstable`.
- Partition suffixes by their leading characters instead of random
  sampling, sort them by a 64-bit key of those characters, and merge
  sort only ties.
- Change the default for `sufr create --num-partitions` from 16 to 1024.
  It caps the number of sort tasks.
- Construction is several times faster on DNA in our tests, and can be
  orders of magnitude faster on inputs with long repeats, helped by an
  LCP cache and vectorized LCP calculation on x86-64 and AArch64. It
  writes the output in parallel, writes no temporary files, and uses
  significantly less peak memory on DNA in our tests.
- (breaking) `sufr create` and `SuffixArrayBuilder::build_to_file`
  refuse a seed mask too long for the 64-bit sort key: for typical DNA,
  more than 21 care positions or a span over 84 characters. We expect
  few masks to reach the limit.
- (breaking) `sufr create -i` requires `--dna`. In 0.7.12, `-i` also
  wrote N over lowercase protein.
- (breaking) Raise the minimum supported Rust version from 1.80 to 1.88.
- (breaking) `sufr summarize` prints the rows `Soft Mask` and `Starts`
  in place of `DNA` and `Ignore Softmask`, and no longer prints
  `Max query len`.
- (breaking) Raise the _.sufr_ format version from 6 to 7. The header
  records the soft mask, the indexed starts and the integer width, and
  no longer holds the build-time max query length. This release does not
  open a file written by 0.7.12 or earlier; rebuild it with
  `sufr create`.
- `sufr locate` lists hits in the order of the input sequences rather
  than by name.
- An `extract` window no longer runs past its sequence into the
  delimiter, the sentinel or the next sequence.
- (breaking) libsufr returns its own error types, `BuildError`,
  `OpenError`, `QueryError` and `InvalidSeedMask`, in place of
  `anyhow::Error`. An error about a file carries its path, and the enums
  are `#[non_exhaustive]`.
- `SuffixArrayBuilder::from_sequence_file` replaces bytes in a record ID
  that are not valid UTF-8 with U+FFFD. It used to fail the build.

### Removed

- (breaking) Remove `SufrMetadata` and `SuffixArray::metadata`.
  `SuffixArray::header` replaces it, and `sufr summarize` reads the
  file's size and modification time itself.
- (breaking) Remove the `sufr_builder` module, `SufrBuilder`,
  `SufrBuilderArgs`, `SuffixArray::new`, and `SuffixArray::write`. Use
  `libsufr::SuffixArrayBuilder` instead.
- (breaking) Remove `SuffixArray::read`, `count`, `extract`, `locate`,
  `list` and `string_at`, with their option and result structs
  (`CountOptions`, `CountResult`, `ExtractOptions`, `ExtractResult`,
  `ExtractSequence`, `LocateOptions`, `LocateResult`, `LocatePosition`,
  `ListOptions`, `SearchOptions`, `SearchResult`,
  `SearchResultLocations`). Use `SuffixArray::search` and `Hits`
  instead.
- (breaking) Remove `SufrFile`. Use `SuffixArray::open_with` instead.
- (breaking) Remove the public modules. `SuffixArray` and `SeedMask` are
  at the crate root. `OUTFILE_VERSION`, `SENTINEL_CHARACTER`,
  `Comparison`, `SequenceFileData`, `read_sequence_file`,
  `read_text_length` and `find_lcp_full_offset` are private or gone.
- (breaking) Remove `SuffixSortType`. `Header::seed_mask` returns the
  seed mask a file was sorted by, or `None` for a file sorted on whole
  suffixes.
- (breaking) Remove the build-time max query length, `sufr create -m`.
  The query-time `-m` of `count`, `locate` and `extract` stays.
- (breaking) Remove `sufr create -r,--random-seed`, which the new
  partitioning does not use.
- (breaking) Remove `sufr create -a,--allow-ambiguity` and the "Allow
  Ambiguity" row of `sufr summarize`. A DNA build indexes only suffixes
  starting with A, C, G or T. With `-a`, 0.7.12 could misorder suffixes
  in runs of 1000 or more N and store LCPs there that were too short.
- (breaking) Remove the byte-conversion utilities `slice_u8_to_vec`,
  `usize_to_bytes` and `vec_to_slice_u8`.
- (breaking) Make the fields of `SeedMask` private, and remove
  `differences` and `SeedMask::from_bytes`. `SeedMask::positions` and
  `SeedMask::weight` return the care positions and their number.
- (breaking) Replace `SeedMask::new` with `impl FromStr for SeedMask`.
  Parse a mask with `"101".parse::<SeedMask>()`.

### Fixed

- `count`, `locate` and `extract` exit with an error instead of
  panicking when they need the subsample cache for `-m` and there is no
  home directory.
- The error for a file of another format version names the file.
- `sufr` and `SuffixArray::open` refuse a file whose seed mask would not
  parse as a `SeedMask`. 0.7.12 accepted malformed masks.
- `sufr` exits 0 with nothing on stderr when the reader of its output
  stops early, as `head` does.
- Reading the suffix and LCP arrays no longer risks unaligned integer
  reads, which are undefined behaviour in Rust.
- Search cache files for `-m` no longer contain stray bytes after the
  data. The extra bytes were never read.
- A query that sorts before every suffix no longer panics in debug
  builds.
- A query with more bytes after a `$` no longer panics; it finds
  nothing.
- Searches with a max query length no longer miss some hits, and no
  longer panic on some queries.
- `sufr list -o` writes to the named file.
- `count`, `locate` and `extract` stop with an error when a read fails,
  instead of skipping the query.
