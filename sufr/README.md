# Parallel Construction of Suffix Arrays in Rust

This code is inspired by [Cache-friendly, Parallel, and Samplesort-based Constructor for Suffix Arrays and LCP Arrays](doi.org/10.4230/LIPIcs.WABI.2023.16).
We copied many ideas from the original C++ implementation [CaPS-SA](https://github.com/jamshed/CaPS-SA.git), most notably the mergesort that constructs the longest common prefix (LCP).

## Setup

* [Install Rust](https://www.rust-lang.org/tools/install)
* Run **`cargo install sufr`** to install the CLI
* Alternately, execute `cargo run` in the source code directory to build a debug version of the program from source

```
$ cargo run
Parallel Construction of Suffix Arrays in Rust

Usage: sufr [OPTIONS] [COMMAND]

Commands:
  create     Create sufr file
  check      Check sufr file for correctness
  extract    Extract suffixes from a sufr file
  list       List the suffix array from a sufr file
  count      Count occurrences of sequences in a sufr file
  locate     Locate sequences in a sufr file
  summarize  Summarize sufr file
  help       Print this message or the help of the given subcommand(s)

Options:
  -t, --threads <THREADS>    Number of threads
  -l, --log <LOG>            Log level [possible values: info, debug]
      --log-file <LOG_FILE>  Log file
  -h, --help                 Print help
  -V, --version              Print version
```

* Execute `cargo build --release` to build an optimized binary in _./target/release/sufr_, which you can execute directly and copy into your `$PATH`:

```
$ cargo build --release

$ ./target/release/sufr -h
Parallel Construction of Suffix Arrays in Rust

Usage: sufr [OPTIONS] [COMMAND]

Commands:
  create     Create sufr file
  check      Check sufr file for correctness
  extract    Extract suffixes from a sufr file
  list       List the suffix array from a sufr file
  count      Count occurrences of sequences in a sufr file
  locate     Locate sequences in a sufr file
  summarize  Summarize sufr file
  help       Print this message or the help of the given subcommand(s)

Options:
  -t, --threads <THREADS>    Number of threads
  -l, --log <LOG>            Log level [possible values: info, debug]
      --log-file <LOG_FILE>  Log file
  -h, --help                 Print help
  -V, --version              Print version
```

Some of the commands may produce debug/info messages.
Use the `-l|--log` option to view these on STDOUT or optionally write to a given `--log-file`.

### Create a suffix array

To begin, you must create a _.sufr_ using the `create` (`cr`) action:

```
$ sufr create -h
Create sufr file

Usage: sufr create [OPTIONS] <INPUT>

Arguments:
  <INPUT>  Input file

Options:
  -n, --num-partitions <NUM_PARTS>  Subproblem count [default: 16]
  -m, --max-query-len <CONTEXT>     Max context
  -o, --output <OUTPUT>             Output file
  -d, --dna                         Input is DNA
  -a, --allow-ambiguity             Allow suffixes starting with ambiguity codes
  -i, --ignore-softmask             Ignore suffixes in soft-mask/lowercase regions
  -D, --sequence-delimiter <DELIM>  Character to separate sequences [default: %]
  -h, --help                        Print help
```

The resulting output file will contain:

* metadata about the input
* the entire input text encoded as `u8` (bytes)
* a fully sorted suffix array (SA)
* an array of the LCP (longest common prefix) for the SA

The input file is a required positional argument and should be a FASTA/Q-formatted file with one or more sequences.

```
$ cat sufr/tests/inputs/2.fa
>ABC
acgtacgt
>DEF
acgtacgt
```

First, read the input sequences into a string of `u8` bytes and uppercase the characters.
Each sequence is separated by a specified character (`%` by default).
A final dollar sign (`$`) is appended to the end of the input text.
If the `--dna` flag is present, suffixes are skipped if they begin with any character other than _A_, _C_, _G_, or _T_ unless the `--allow-ambiguity` flag is present.
Additionally, soft-masked (lowercase) input is converted to uppercase unless the `--ignore-softmask` flag is present, in which case it is converted to `N` and ignored.

Next, we partition the suffixes into some number `N` partitions by randomly selecting `--num-partitions` - 1 pivot suffixes, sorting them, and using the pivots to place each suffix into the highest bounded partition.
The partitions are sorted using a merge sort algorithm that also generates an LCP (longest common prefix) array.
The sorted suffix/LCP arrays are then concatenated to produce the final output.

The `sufr` CLI will create an output file containing a binary-encoded representation of the sorted suffix/LCP arrays along with the original sequence data and other metadata used to generate the arrays.
For instance, with the _1.fa_ file, the default output file will be _1.sufr_:

```
$ sufr --log debug create sufr/tests/inputs/1.fa -n 2 --dna
[2024-11-22T22:28:45Z INFO  sufr] Using 8 threads
[2024-11-22T22:28:45Z INFO  sufr] Read input of len 11 in 1.0535ms
[2024-11-22T22:28:45Z INFO  libsufr::sufr_builder] Selected 1 pivot in 59.25µs
[2024-11-22T22:28:45Z INFO  libsufr::sufr_builder] Wrote 9 unsorted suffixes to partition in 356.833µs
[2024-11-22T22:28:45Z INFO  libsufr::sufr_builder] Sorted 9 suffixes in 2 partitions (avg 4) in 834.542µs
[2024-11-22T22:28:45Z INFO  sufr] Wrote 164 bytes to '1.sufr' in 593.875µs
[2024-11-22T22:28:45Z INFO  sufr] Total time: 4.518709ms
```

### Summarize a sufr file

Use the `summarize` (`su`) action to view metadata about a _.sufr_ file:

```
$ sufr su -h
Summarize sufr file

Usage: sufr summarize <SUFR>

Arguments:
  <SUFR>  Sufr file

Options:
  -h, --help  Print help
```

For instance:

```
$ sufr su 1.sufr
+-----------------+------------------+
| Filename        | 1.sufr           |
+-----------------+------------------+
| Modified        | 2024-11-15 09:21 |
+-----------------+------------------+
| File Size       | 164 bytes        |
+-----------------+------------------+
| File Version    | 4                |
+-----------------+------------------+
| DNA             | true             |
+-----------------+------------------+
| Allow Ambiguity | false            |
+-----------------+------------------+
| Ignore Softmask | false            |
+-----------------+------------------+
| Text Length     | 11               |
+-----------------+------------------+
| Num Suffixes    | 9                |
+-----------------+------------------+
| Max query len   | 0                |
+-----------------+------------------+
| Num sequences   | 1                |
+-----------------+------------------+
| Sequence starts | 0                |
+-----------------+------------------+
| Headers         | 1                |
+-----------------+------------------+
```

### Check a suffix array

Use the `check` (`ch`) action to verify the order of the suffix array:

```
$ sufr check -h
Check sufr file for correctness

Usage: sufr check [OPTIONS] <SUFR>

Arguments:
  <SUFR>  Sufr file

Options:
  -v, --verbose  List errors
  -h, --help     Print help
```

For instance:

```
$ sufr check 1.sufr
Checked 9 suffixes, found 0 errors in suffix array.
Finished checking in 1.112417ms.
```

### Listing suffixes in a sufr file

You can use the `list` (`ls`) action to view the sorted arrays by their _rank_:

```
$ sufr ls -h
List the suffix array from a sufr file

Usage: sufr list [OPTIONS] <SUFR> [RANK]...

Arguments:
  <SUFR>     Sufr file
  [RANK]...  Ranks of suffixes to show

Options:
  -l, --len <LEN>     Length of suffixes to show
  -o, --output <OUT>  Output
  -h, --help          Print help
```

For example, to view the _1.sufr_ file:

```
$ sufr ls 1.sufr
 R  S  L
 0 10  0: $
 1  6  0: ACGT$
 2  0  4: ACGTNNACGT$
 3  7  0: CGT$
 4  1  3: CGTNNACGT$
 5  8  0: GT$
 6  2  2: GTNNACGT$
 7  9  0: T$
 8  3  1: TNNACGT$
```

An optional positional argument for the _ranks_ allows you to select only a portion:

```
$ sufr ls 1.sufr 0-5
 R  S  L
 0 10  0: $
 1  6  0: ACGT$
 2  0  4: ACGTNNACGT$
 3  7  0: CGT$
 4  1  3: CGTNNACGT$
 5  8  0: GT$
```

As suffixes can get quite long, use the `-l|--len` option to restrict the length of the shown suffix:

```
$ sufr ls 1.sufr 2-3 -l 3
 R  S  L
 2  0  4: ACG
 3  7  0: CGT
```

### Count occurrences of suffixes

Use the `count` (`co`) command to find the number of occurrences of suffixes:

```
$ sufr co -h
Count occurrences of sequences in a sufr file

Usage: sufr count [OPTIONS] <SUFR> <QUERY>...

Arguments:
  <SUFR>      Sufr file
  <QUERY>...  Query

Options:
  -m, --max-query-len <LEN>  Maximum query length
  -o, --output <OUT>         Output
  -l, --low-memory           Low-memory
  -h, --help                 Print help
```

The `-l|--low-memory` option will force the suffixes to be read from disk rather than loaded into memory. 
The time to load a large index (e.g., human genome) into memory may take longer than the actual search, so you may find this is faster for only a few queries but slower for a large number.

For example:

```
$ sufr co 1.sufr AC GT X
AC 2
GT 2
X not found
```

Suffixes that do not occur are printed to `STDERR` such as the "X" in the preceding output.

### Locate suffixes

Use the `locate` (`lo`) command to find the positions of a given suffix:

```
$ sufr lo -h
Locate sequences in a sufr file

Usage: sufr locate [OPTIONS] <SUFR> <QUERY>...

Arguments:
  <SUFR>      Sufr file
  <QUERY>...  Query

Options:
  -m, --max-query-len <LEN>  Maximum query length
  -o, --output <OUT>         Output
  -l, --low-memory           Low-memory
  -a, --abs                  Show absolute position in text
  -h, --help                 Print help
```

For instance, given the _3.fa_ input:

```
$ cat sufr/tests/inputs/3.fa
>1
AGCTTTTCATTCTGACTGCAACGG
GCAATANNNNNNNNNNTGTCTC
>2
TGTGTGGATTAAAAAAAGAGTGTC
>3
TGATAGCAGCTTCTGAACTGGTTA
CCTGCCGTGAGTAAAT
```

The suffix "ACT" is found at position 14 in sequence _1_ and 16 in sequence _3_, "GT" is found at position 20 in sequence _3_, and "X" is not found and is printed to `STDERR`:

```
$ sufr lo 3.sufr ACT GTT
ACT
1 14
3 16
//
GTT
3 20
//
X not found
```

Use the `-a|--absolute` flag if you prefer to see the raw suffix position:

```
$ sufr lo 3.sufr ACT GTT -a
ACT 14 88
GTT 92
```

By default, the entire suffix array will be loaded into memory.
If your machine lacks suffixient resources, you can use the `-m|--max-query-len` option to create a down-sampled suffix array.
Or use the `-l|--low-memory` option to search the suffix array on disk, which is slower but requires only enough memory to hold the original text/sequences.

### Extract suffixes

Use the `extract` (`ex`) to get suffixes by their position in the text:

```
$ sufr ex -h
Extract suffixes from a sufr file

Usage: sufr extract [OPTIONS] <SUFR> [SUFFIX]...

Arguments:
  <SUFR>       Sufr file
  [SUFFIX]...  Suffixes to extract

Options:
  -p, --prefix-len <PREFIX_LEN>  Prefix length
  -s, --suffix-len <SUFFIX_LEN>  Suffix length
  -l, --lcp                      Show LCP
  -o, --output <OUT>             Output
  -h, --help                     Print help
```

Using the preceding `locate` results, I can extract the suffixes at those positions.
As the suffixes can get quite long, I will use the `-s|--suffix-len` option to limit them to 15bp:

```
$ sufr ex 3.sufr 14 88 92 -s 15
ACTGCAACGGGCAAT
ACTGGTTACCTGCCG
GTTACCTGCCGTGAG
```

Combine this with the `-p|--prefix-len` to control the amount of preceding text, which might be useful when identifying alignment seeds:

```
$ sufr ex 3.sufr 14 88 92 -s 15 -p 5
TTCTGACTGCAACGGGCAAT
TCTGAACTGGTTACCTGCCG
AACTGGTTACCTGCCGTGAG
```

## Testing

Run **`cargo test`**.

## Authors

* Ken Youens-Clark <kyclark@arizona.edu>
* Jack Roddy <jroddy@pharmacy.arizona.edu>
* Travis Wheeler <twheeler@arizona.edu>
