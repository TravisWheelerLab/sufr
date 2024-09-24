# Parallel Construction of Suffix Arrays in Rust

This code is inspired by [Cache-friendly, Parallel, and Samplesort-based Constructor for Suffix Arrays and LCP Arrays](doi.org/10.4230/LIPIcs.WABI.2023.16).
I copied many ideas from the original C++ implementation [CaPS-SA](https://github.com/jamshed/CaPS-SA.git), most notably the mergesort that constructs the longest common prefix (LCP).

The basic ideas are as follow:

* Read the input file as `u8` (unsigned 8-bit integer values). 
* Select the suffixes, which are normally 0 to the length of the text but there is the option to skip suffixes starting with _N_. Note: The original C++ uses 32-bit integers if the input length is less than 2^32 and 64-bit integers, otherwise. This Rust implementation currently only uses `usize` values, which defaults to 64-bit on 64-bit architectures. Copying this idea would certainly use less disk space for the resulting SA and may result in greater performance and less memory usage.
* Split the input into subarrays and sort, also producing arrays for the LCP and ordered candidate suffixes for global pivots.
* Merge the candidate pivots and downsample to select suffixes for global pivots.
* Use the global pivots to partition the original subarrays into sub-subarrays containing the suffixes that fall into given ranges. E.g., all the values less than the first/lowest pivot suffix, then all the values greater than or equal to the first pivot and less than the second, and so on.
* The resulting sub-subarrays are sorted that are merged. Because the values fall into nonoverlapping ranges, these subarrays can be appended in order to produce the final SA.
* Produce a binary-encoded output file of `usize` values containing the length of the SA in the first position and the sorted SA values following.

Some advantages to this algorithm:

* The various partitioned subarrays can be processed independently by separate threads, and no thread will ever have to merge the entire input.
* Suffix comparisons are made faster by caching LCPs.
* Keeping the input text as an array of 8-bit integers vs `char` (UTF-8) results in lower memory usage.

## Setup

* [Install Rust](https://www.rust-lang.org/tools/install)
* Execute `cargo run` to build a debug version of the program from source

```
$ cargo run
Usage: sufr [COMMAND]

Commands:
  check   Create suffix array
  create  Create suffix array
  read    Read suffix array and extract sequences
  help    Print this message or the help of the given subcommand(s)

Options:
  -h, --help  Print help
```

* Execute `cargo build --release` to build an optimized binary in _./target/release/sufr_, which you can execute directly and copy into your `$PATH`:

```
$ cargo build --release
   Compiling sufr v0.1.0 (/Users/kyclark/wheelerlab/sufr)
    Finished `release` profile [optimized] target(s) in 10.57s

$ ./target/release/sufr -h
Usage: sufr [COMMAND]

Commands:
  create   Create suffix array
  check    Check correctness of suffix array/LCP
  extract  Read suffix array and extract sequences
  help     Print this message or the help of the given subcommand(s)

Options:
  -h, --help  Print help
```

## Code Overview

The code is organized into a Cargo workspace (https://doc.rust-lang.org/book/ch14-03-cargo-workspaces.html):

* _libsufr/src/lib.rs_: Core functionality to sort a suffix array and create LCP
* _sufr/src/main.rs_: The main entry point for the `sufr` CLI
* _sufr/src/lib.rs_: A library that implements the CLI functions

As the CLI usage shows, `sufr` currently supports three actions, create, check, and extract.

### Create a suffix array

To begin, you must create a suffix array (SA) using the `create` action:

```
$ cargo run -- create -h
Create suffix array

Usage: sufr create [OPTIONS] <INPUT>

Arguments:
  <INPUT>  Input file

Options:
  -n, --num-partitions <NUM_PARTS>  Subproblem count [default: 16]
  -m, --max-context <CONTEXT>       Max context
  -t, --threads <THREADS>           Number of threads
  -o, --output <OUTPUT>             Output file
  -d, --dna                         Input is DNA, ignore sequences starting with 'N'
  -c, --check                       Verify order
  -l, --log <LOG>                   Log level [possible values: info, debug]
      --log-file <LOG_FILE>         Log file
  -h, --help                        Print help
  -V, --version                     Print version
```

The input file is a required positional argument and should be a FASTA/Q-formatted file with one or more sequences.

```
$ cat sufr/tests/inputs/2.fa
>ABC
acgtacgt
>DEF
acgtacgt
```

The algorithm works as follows.
First, read the input sequences into a string of `u8` bytes and uppercase the characters.
Each sequence is separated by a dollar sign (`$`), and a final hash sign (`#`) is appended to the end.
For instance, the preceding sequence has the following bytes and suffix positions:

```
seq:    [ A,  C,  G,  T,  A,  C,  G,  T,  $,  A,  C,  G,  T,  A,  C,  G,  T,  #]
bytes:  [65, 67, 71, 84, 65, 67, 71, 84, 36, 65, 67, 71, 84, 65, 67, 71, 84, 35]
suffix: [ 0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17]
```

NOTE: If the `--dna` flag is present, suffixes are skipped if they begin with any character other than _A_, _C_, _G_, or _T_.

Next, we partition the suffixes into `N` partitions by randomly selecting `--num-partitions` - 1 pivot suffixes, sorting them, and using the pivots to place each suffix into the highest bounded partition.
The partitions are sorted using a merge sort algorithm that also generates an LCP (longest common prefix) array.
The sorted suffix/LCP arrays are then concatenated to produce the final output.

The preceding example is sorted into the following order/LCP:

```
   Pos  LCP   Suffix
    13    0   ACGT#
     4    4   ACGT$ACGTACGT#
     9    4   ACGTACGT#
     0    8   ACGTACGT$ACGTACGT#
    14    0   CGT#
     5    3   CGT$ACGTACGT#
    10    3   CGTACGT#
     1    7   CGTACGT$ACGTACGT#
    15    0   GT#
     6    2   GT$ACGTACGT#
    11    2   GTACGT#
     2    6   GTACGT$ACGTACGT#
    16    0   T#
     7    1   T$ACGTACGT#
    12    1   TACGT#
     3    5   TACGT$ACGTACGT#
```

### Read a suffix array

```
$ cargo run -- read -h
Read suffix array and extract sequences

Usage: sufr read [OPTIONS] --array <SA> --sequence <SEQ>

Options:
  -a, --array <SA>         Suffix array file
  -s, --sequence <SEQ>     Sequence file
  -m, --max-len <MAX>      Maximum length of sequence [default: 0]
  -e, --extract <EXTRACT>  Extract positions [default: 1]
  -n, --number             Number output
  -o, --output <OUTPUT>    Output
  -h, --help               Print help
  -V, --version            Print version
```

### Check a suffix array

```
$ cargo run -- check -h
Create suffix array

Usage: sufr check --array <SA> --sequence <SEQ>

Options:
  -a, --array <SA>      Suffix array file
  -s, --sequence <SEQ>  Sequence file
  -h, --help            Print help
  -V, --version         Print version
```

## Testing

Run **`cargo test`**.

## Author

Ken Youens-Clark <kyclark@arizona.edu>
