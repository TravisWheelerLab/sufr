# Parallel Construction of Suffix Arrays in Rust

This code is inspired by [Cache-friendly, Parallel, and Samplesort-based Constructor for Suffix Arrays and LCP Arrays](http://doi.org/10.4230/LIPIcs.WABI.2023.16).
We copied many ideas from the original C++ implementation [CaPS-SA](https://github.com/jamshed/CaPS-SA.git), most notably the mergesort that constructs the longest common prefix (LCP).

## Setup

* [Install Rust](https://www.rust-lang.org/tools/install)   (this is easy)
* Run **`cargo install sufr`** to install the CLI
* Alternately, execute `cargo run` in the source code directory to build a debug version of the program from source

```
$ cargo run
Parallel Construction of Suffix Arrays in Rust

Usage: sufr [OPTIONS] [COMMAND]

Commands:
  create     Create sufr file
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
$ sufr cr -h
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
  -s, --seed-mask <MASK>            Spaced seeds mask
  -r, --random-seed <RANDSEED>      Random seed [default: 42]
  -h, --help                        Print help
```

The resulting binary-encoded output file will contain:

* metadata about the input
* the entire input text encoded as `u8` (bytes)
* a fully sorted suffix array (SA)
* an array of the LCP (longest common prefix) for the SA

The input file is a required positional argument and should be a FASTA/Q-formatted file with one or more sequences, e.g.:

```
$ cat data/inputs/2.fa
>ABC
acgtacgt
>DEF
acgtacgt
```

To briefly describe our algorigthm, the input sequences is read into a vector of `u8` bytes.
Unless we are told to ignore softmasked input, lowercase values are converted to uppercase; in the case that we do ignore softmask data, lowercase are converted to an ambiguity character (_N_ for nucleotides and _X_ otherwise, e.g., protein).
Multiple sequence are separated by a specified character (`%` by default).
A sentinel character is appended to the end (default `$`) is appended to the end of the input text.
If the `--dna` flag is present, suffixes are skipped if they begin with any character other than _A_, _C_, _G_, or _T_ unless the `--allow-ambiguity` flag is present.

Next, we partition the suffixes into some number partitions by randomly selecting `--num-partitions` - 1 pivot suffixes, sorting them, and using the pivots to place each suffix into the highest bounded partition.
The partitions are sorted using a merge sort algorithm that also generates an LCP (longest common prefix) array.
The sorted suffix/LCP arrays are then concatenated to produce the final output.

The `sufr` CLI will create an output file containing a binary-encoded representation of the sorted suffix/LCP arrays along with the original sequence data and other metadata used to generate the arrays.
For instance, with the _1.fa_ file, the default output file will be _1.sufr_:

```
$ sufr --log debug create data/inputs/1.fa -n 2 --dna
[2025-01-29T18:56:00Z INFO  sufr] Using 8 threads
[2025-01-29T18:56:00Z INFO  sufr] Read input of len 11 in 980.5µs
[2025-01-29T18:56:00Z INFO  libsufr::sufr_builder] Selected 1 pivot in 51.917µs
[2025-01-29T18:56:00Z INFO  libsufr::sufr_builder] Wrote 9 unsorted suffixes to partition in 282.208µs
[2025-01-29T18:56:00Z INFO  libsufr::sufr_builder] Sorted 9 suffixes in 2 partitions (avg 4) in 530.292µs
[2025-01-29T18:56:00Z INFO  sufr] Wrote 172 bytes to '1.sufr' in 1.822333ms
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
| Modified        | 2025-01-29 11:56 |
+-----------------+------------------+
| File Size       | 172 bytes        |
+-----------------+------------------+
| File Version    | 6                |
+-----------------+------------------+
| DNA             | true             |
+-----------------+------------------+
| Allow Ambiguity | false            |
+-----------------+------------------+
| Ignore Softmask | false            |
+-----------------+------------------+
| Text Length     | 11               |
+-----------------+------------------+
| Len Suffixes    | 9                |
+-----------------+------------------+
| Max query len   | 0                |
+-----------------+------------------+
| Num sequences   | 1                |
+-----------------+------------------+
| Sequence starts | 0                |
+-----------------+------------------+
| Sequence names  | 1                |
+-----------------+------------------+
```

### Listing suffixes in a sufr file

You can use the `list` (`ls`) action to view the sorted arrays by their _rank_:

```
$ sufr ls -h
List the suffix array from a sufr file

Usage: sufr list [OPTIONS] <FILE> [RANK]...

Arguments:
  <FILE>     Sufr file
  [RANK]...  Ranks of suffixes to show

Options:
  -r, --show-rank        Show rank column
  -s, --show-suffix      Show suffix position column
  -p, --show-lcp         Show LCP column
  -v, --very-low-memory  Very low memory
      --len <LEN>        Length of suffixes to show
  -n, --number <LEN>     Number of suffixes to show
  -o, --output <OUT>     Output
  -h, --help             Print help
```

For example, to view the _1.sufr_ file including rank, suffix, and LCP:

```
$ sufr ls -rsp 1.sufr
 0 10  0 $
 1  6  0 ACGT$
 2  0  4 ACGTNNACGT$
 3  7  0 CGT$
 4  1  3 CGTNNACGT$
 5  8  0 GT$
 6  2  2 GTNNACGT$
 7  9  0 T$
 8  3  1 TNNACGT$
```

An optional positional argument for the _ranks_ allows you to select only a portion:

```
$ sufr ls 1.sufr 0-5
$
ACGT$
ACGTNNACGT$
CGT$
CGTNNACGT$
GT$
```

As suffixes can get quite long, use the `-l|--len` option to restrict the length of the shown suffix:

```
$ sufr ls 1.sufr 2-3 --len 3
ACG
CGT
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
  -l, --low-memory           Low memory
  -v, --very-low-memory      Very low memory
  -h, --help                 Print help
```

.Low/very low memory usage
****
The `-l|--low-memory` option will force the suffixes to be read from disk rather than loaded into memory. 
The `-v|--very-low-memory` option will also force the text to be read from disk rather than loaded into memory. 
The time to load a large index (e.g., human genome) into memory may take longer than the actual search, so you may find this is faster for only a few queries but slower for a large number.
Also consider the `-m|--max-query-len` option to create a down-sampled suffix array.
****

For example:

```
$ sufr co 1.sufr AC GT X
AC 2
GT 2
X 0
```

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
  -o, --output <OUT>         Output
  -m, --max-query-len <LEN>  Maximum query length
  -l, --low-memory           Low memory
  -v, --very-low-memory      Very low memory
  -a, --abs                  Show absolute position in text
  -h, --help                 Print help
```

For instance, given the _3.fa_ input:

```
$ cat data/inputs/3.fa
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
$ sufr lo data/expected/3.sufr ACT GTT X
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

### Extract suffixes

Use the `extract` (`ex`) to print suffixes in FASTA format:

```
$ sufr ex -h
Extract suffixes from a sufr file

Usage: sufr extract [OPTIONS] <SUFR> <QUERY>...

Arguments:
  <SUFR>      Sufr file
  <QUERY>...  Query

Options:
  -m, --max-query-len <LEN>      Maximum query length
  -l, --low-memory               Low memory
  -v, --very-low-memory          Very low memory
  -p, --prefix-len <PREFIX_LEN>  Prefix length
  -s, --suffix-len <SUFFIX_LEN>  Suffix length
  -o, --output <OUT>             Output
  -h, --help                     Print help
```

Using the preceding `locate` results, I can extract the suffixes at those positions.
As the suffixes can get quite long, I will use the `-s|--suffix-len` option to limit them to 15bp:

```
$ sufr ex data/expected/3.sufr ACT GTT -s 15
>1:14-29 ACT 0
ACTGCAACGGGCAAT
>3:16-31 ACT 0
ACTGGTTACCTGCCG
>3:20-35 GTT 0
GTTACCTGCCGTGAG
```

Combine this with the `-p|--prefix-len` to control the amount of preceding text, which might be useful when identifying alignment seeds:

```
$ sufr ex data/expected/3.sufr ACT GTT -s 15 -p 3
>1:11-29 ACT 3
CTGACTGCAACGGGCAAT
>3:13-31 ACT 3
TGAACTGGTTACCTGCCG
>3:17-35 GTT 3
CTGGTTACCTGCCGTGAG
```

The FASTA header contains the matching sequence ID, a colon, the start/stop position of the extracted sequence, followed by the query, and finally the location of the query in the extracted sequence, which is only relevant if you included a prefix.

## Testing

Run **`cargo test`**.

## Authors

* Ken Youens-Clark <kyclark@arizona.edu>
* Jack Roddy <jroddy@arizona.edu>
* Travis Wheeler <twheeler@arizona.edu>
