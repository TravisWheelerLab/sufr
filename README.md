# Parallel Construction of Suffix Arrays in Rust

This code is inspired by [Cache-friendly, Parallel, and Samplesort-based Constructor for Suffix Arrays and LCP Arrays](doi.org/10.4230/LIPIcs.WABI.2023.16).
I copied many ideas from the original C++ implementation [CaPS-SA](https://github.com/jamshed/CaPS-SA.git), most notably the mergesort that constructs the longest common prefix (LCP).

The basic ideas are as follow:

* Read the input file as `u8` (unsigned 8-bit integer values). Note: The original C++ uses 32-bit integers if the number of suffixes falls is less than 2^32 and 64-bit integers, otherwise. This Rust implementation currently only uses `usize` values, which defaults to 64-bit on 64-bit architectures. Copying this idea would certainly use less disk space for the resulting SA and may result in greater performance and less memory usage when creating.
* Select the suffixes, which are normally 0 to the length of the text but there is the option to skip suffixes starting with _N_.
* Split the input into subarrays and sort, also producing LCPs and candidate suffixes for global pivots
* Select and order global pivots.
* Use the global pivots to partition the original subarrays into new subarrays where all the values fall into given ranges. E.g., all the values less than the first/lowest pivot suffix, then all the values greater than or equal to the first pivot and less than the second, and so on.
* The resulting sub-subarrays will have sorted selections from the original subarrays that can be merged into one array. Because the values fall into nonoverlapping ranges, these subarrays can be appended in order to produce the final SA.
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

* Execute `cargo build --release` to build an optimized binary in _./target/release/sufr_, whic you can execute directly and copy into your `$PATH`:

```
$ cargo build --release
   Compiling sufr v0.1.0 (/Users/kyclark/wheelerlab/sufr)
    Finished `release` profile [optimized] target(s) in 10.57s

$ ./target/release/sufr -h
Usage: sufr [COMMAND]

Commands:
  check   Create suffix array
  create  Create suffix array
  read    Read suffix array and extract sequences
  help    Print this message or the help of the given subcommand(s)

Options:
  -h, --help  Print help
```

## Code Overview

The code is organized like standard Rust programs in the _src_ directory:

* _src/main.rs_: The main entry point for the `sufr` CLI
* _src/lib.rs_: A library that implements the CLI functions
* _src/suffix_array.rs_: A crate for interacting with suffix arrays

As the CLI usage shows, `sufr` currently supports three actions:

### Create a suffix array

To begin, you must create a suffix array (SA) using the `create` action:

```
$ cargo run -- create -h
Create suffix array

Usage: sufr create [OPTIONS] <INPUT>

Arguments:
  <INPUT>  Input file

Options:
  -s, --subproblem-count <SUBPROBLEMS>  Subproblem count [default: 16]
  -m, --max-context <CONTEXT>           Max context
  -t, --threads <THREADS>               Number of threads [default: 16]
  -o, --output <OUTPUT>                 Output file [default: sufr.sa]
  -i, --ignore-start-n                  Ignore sequences starting with N
  -c, --check                           Verify order
  -l, --log <LOG>                       Log level [possible values: info, debug]
      --log-file <LOG_FILE>             Log file
  -h, --help                            Print help
  -V, --version                         Print version
```

The input file is a required positional argument.
Currently, the file should be a plain text file with a single genomic sequence on one line.
The program will likely read FASTA/Q files in the future.
E.g.:

```
$ cat tests/inputs/seq1.txt
CTCACC
```

The algorithm works in the following steps.
First, read the input file as a string of `u8` bytes (note this could therefore fail if provided UTF-8) and uppercase the characters. 
At present, no filtering is applied, i.e., screening for a particular genomic alphabet (DNA/RNA/AA) or ambiguity codes (IUPAC), which is why the sequence should include no line breaks (at present).
To illustrate, consider the following input:

```
$ cat sample.txt
CTNNCACC
```

The raw input to the program looks like the following:

```
[67, 84, 78, 78, 67, 65, 67, 67, 36]
  C   A   N   N   C   A   C   C   $
```

Next, determine the start positions of the suffixes starting from the first character and going to the end, which is noted with an appended `$`.
The suffix array for the preceding sequence is as follows:

```
[ 0, 1, 2, 3, 4, 5, 6, 7, 8 ]
  C  A  N  N  C  A  C  C  $
```

The suffixes are as follows:

```
0 CTNNCACC$
1 TNNCACC$
2 NNCACC$
3 NCACC$
4 CACC$
5 ACC$
6 CC$
7 C$
8 $
```

If `-i|--ignore-start-n` is true, then suffixes beginning with `N` are skipped.
Using the same example sequence, the suffixes would be as follows:

```
0 CTNNCACC$
1 TNNCACC$
4 CACC$
5 ACC$
6 CC$
7 C$
8 $
```

Split the input sequence into `-s|--subproblem-count` arrays of equal size and sort, optionally using the `-m|--max-context` argument to limit string comparisons to a length less than the entire string.
For example, if you are aligning very short sequences, you might care to set a context of 16 or 20 base pairs.
This will speed up the creation of the SA, but the resulting SA will be only partially sorted.

The `sort_subarrays` function returns the sorted SA, the LCP, and an evenly selected sample of suffixes that will be used to find pivots later in the algorithm.
For example the previous sequence split into 3 subproblems will produce the following:

```
[
    (
        [      // This is the sorted SA
            0, // CTNNCACC$ sorts before
            1, // TNNCACC$
        ],
        [      // This is the LCP for the sorted SA
            0, // The first element in the LCP is always 0
            0, // TNNCACC$ has an LCP of 0 to CTNNCACC$
        ],
        [      // This is a sample from the SA for candidate pivots
            1, // The suffix "TNNCACC$" was selected
        ],
    ),
    (
        [      
            5, // ACC$ sorts before
            4, // CACC$
        ],
        [
            0, // First LCP always 0
            0, // No common prefix from CACC$ to preceding value ACC$
        ],
        [
            4, // CACC$ was selected for pivot
        ],
    ),
    (
        [
            8, // $ sorts before
            7, // C$ sorts before
            6, // CC$
        ],
        [
            0, // First LCP always 0
            0, // No LCP from C$ to $
            1, // LCP of 1 from CC$ to C$
        ],
        [
            6, // CC$ was selected for pivot
        ],
    ),
]
```

Next, select pivot suffixes from the subarray pivots.
Using the preceding example, this would happen by first constructing an LCP array for the sampled suffixes to use in merging them all into a single sorted SA:

```
[
    4, // CACC$
    6, // CC$
    1, // TNNCACC$
],
```

Then select the final global pivots:

```
[
    4, // CACC$
    6, // CC$
],
```

Next, search the original subarrays for ranges of values delimited by the pivots.
For instance:

```
[
    [             // For the subarray [ 0, 1 ] => [ CTNNCACC$, TNNCACC$ ]
        None,     // Values <= CACC$
        None,     // Values <= CC$
        Some(     // Values >  CC$
            0..2,
        ),
    ],
    [             // For the subarray [ 5, 4 ] => [ ACC$, CACC$ ]
        Some(     // Values <= CACC$
            0..2, 
        ),
        None,     // Values <= CC$
        None,     // Values >  CC$
    ],
    [             // For the subarray [ 8, 7, 6 ] => [ $, C$, CC$ ]
        Some(     // Values <= CACC$
            0..2,
        ),
        Some(     // Values <= CC$
            2..3,
        ),
        None,     // Values >  CC$
    ],
]
```

Use these ranges to partition the original subarrays:

```
[
    [             // For the subarray [ 0, 1 ] => [ CTNNCACC$, TNNCACC$ ]
        None,     // Values <= CACC$
        None,     // Values <= CC$
        Some(     // Values >  CC$
            [
                0,
                1,
            ],
        ),
    ],
    [             // For the subarray [ 5, 4 ] => [ ACC$, CACC$ ]
        Some(     // Values <= CACC$
            [
                5,
                4,
            ],
        ),
        None,     // Values <= CC$
        None,     // Values >  CC$
    ],
    [             // For the subarray [ 8, 7, 6 ] => [ $, C$, CC$ ]
        Some(     // Values <= CACC$
            [
                8,
                7,
            ],
        ),
        Some(    // Values <= CC$
            [
                6,
            ],
        ),
        None,    // Values >  CC$
    ],
]
```

Transpose these values to create sub-subarrays of relatively sorted suffixes:

```
[
    [           // Values <= CACC$
        [
            5,  // ACC$
            4,  // CACC$
        ],
        [
            8,  // $
            7,  // C$
        ],
    ],
    [           // Values <= CC$
        [
            6,  // CC$
        ],
    ],
    [           // Values >  CC$
        [
            0,  // CTNNCACC$
            1,  // TNNCACC$
        ],
    ],
]
```

Merge the sub-subarrays:

```
[
    [      // Values <= CACC$
        8, // $
        5, // ACC$
        7, // C$
        4, // CACC$
    ],
    [      // Values <= CC$
        6, // CC$
    ],
    [      // Values >  CC$
        0, // CTNNCACC$
        1, // TNNCACC$
    ],
]
```

Finally, concatenate the suffixes into the final sorted SA:

```
[
    8, // $
    5, // ACC$
    7, // C$
    4, // CACC$
    6, // CC$
    0, // CTNNCACC$
    1, // TNNCACC$
]
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
