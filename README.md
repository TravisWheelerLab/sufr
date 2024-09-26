# Parallel Construction of Suffix Arrays in Rust

This code is inspired by [Cache-friendly, Parallel, and Samplesort-based Constructor for Suffix Arrays and LCP Arrays](doi.org/10.4230/LIPIcs.WABI.2023.16).
We copied many ideas from the original C++ implementation [CaPS-SA](https://github.com/jamshed/CaPS-SA.git), most notably the mergesort that constructs the longest common prefix (LCP).

## Setup

* [Install Rust](https://www.rust-lang.org/tools/install)
* Run **`cargo install sufr`** to install the CLI
* Alternately, execute `cargo run` in the source code directory to build a debug version of the program from source

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

The `sufr` CLI will create an output file containing a binary-encoded representation of the sorted suffix/LCP arrays along with the original sequence data and other metadata used to generate the arrays.
For instance, with the _1.fa_ file, the default output file will be _1.sufr_:

```
$ cargo run -- create --log debug sufr/tests/inputs/1.fa -n 2 --dna --check
    Finished `dev` profile [unoptimized + debuginfo] target(s) in 0.03s
     Running `target/debug/sufr create --log debug sufr/tests/inputs/1.fa -n 2 --dna --check`
[2024-09-24T21:04:14Z INFO  sufr] Using 8 threads
[2024-09-24T21:04:14Z INFO  sufr] Read raw input of len 11 in 6.081625ms
[2024-09-24T21:04:14Z DEBUG sufr] Raw input '[65, 67, 71, 84, 78, 78, 65, 67, 71, 84, 35]'
[2024-09-24T21:04:14Z INFO  libsufr] Created unsorted suffix array of len 8 in 25.209µs
[2024-09-24T21:04:14Z INFO  libsufr] Selected/sorted 1 pivots in 151.708µs
[2024-09-24T21:04:14Z INFO  libsufr] Split into 2 partitions (avg 4) in 309µs
[2024-09-24T21:04:14Z INFO  libsufr] Sorted partitions in 162.75µs
[2024-09-24T21:04:14Z INFO  libsufr] Concatenated partitions in 17.041µs
[2024-09-24T21:04:14Z INFO  libsufr] Fixed LCP boundaries in 6.708µs
[2024-09-24T21:04:14Z INFO  libsufr] Total time to create suffix array: 890.208µs
[2024-09-24T21:04:14Z DEBUG sufr] Sorted = [6, 0, 7, 1, 8, 2, 9, 3]
[2024-09-24T21:04:14Z DEBUG sufr] Suffixes = [
        "ACGT#",
        "ACGTNNACGT#",
        "CGT#",
        "CGTNNACGT#",
        "GT#",
        "GTNNACGT#",
        "T#",
        "TNNACGT#",
    ]
[2024-09-24T21:04:14Z DEBUG sufr] LCP = [
        0,
        4,
        0,
        3,
        0,
        2,
        0,
        1,
    ]
[2024-09-24T21:04:14Z INFO  sufr] Checked order, found 0 errors in 6.584µs
[2024-09-24T21:04:14Z INFO  sufr] Checked LCP, found 0 errors in 25.042µs
[2024-09-24T21:04:14Z INFO  sufr] Wrote 130 bytes to '1.sufr' in 3.253ms
```

### Extracting suffix from a sufr file

You can use the `extract` action to view the sorted arrays:

```
$ cargo run -- extract -h
    Finished `dev` profile [unoptimized + debuginfo] target(s) in 0.03s
     Running `target/debug/sufr extract -h`
Read suffix array and extract sequences

Usage: sufr extract [OPTIONS] <SUFR>

Arguments:
  <SUFR>  Sufr file

Options:
  -m, --max-len <MAX>      Maximum length of sequence
  -e, --extract <EXTRACT>  Extract positions [default: 1]
  -n, --number             Number output
  -o, --output <OUTPUT>    Output
  -h, --help               Print help
  -V, --version            Print version
```

For example, to view the first 10 suffixes from the _1.sufr_ file:

```
$ cargo run -- extract 1.sufr -e 1-10 -n
  0: ACGT#
  1: ACGTNNACGT#
  2: CGT#
  3: CGTNNACGT#
  4: GT#
  5: GTNNACGT#
  6: T#
  7: TNNACGT#
```

### Check a suffix array

Use the `check` action to verify the order of the suffix array:

```
$ cargo run -- check -h
Check correctness of suffix array/LCP

Usage: sufr check [OPTIONS] <SUFR>

Arguments:
  <SUFR>  Sufr file

Options:
  -v, --verbose  List errors
  -h, --help     Print help
  -V, --version  Print version
```

For instance:

```
$ cargo run -- check 1.sufr
Found 0 errors in suffix array.
Found 0 errors in LCP
Finished checking in 264.709µs
```

## Testing

Run **`cargo test`**.

## Authors

* Ken Youens-Clark <kyclark@arizona.edu>
* Jack Roddy <jroddy@pharmacy.arizona.edu>
* Travis Wheeler <twheeler@arizona.edu>
