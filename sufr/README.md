# Parallel Construction of Suffix Arrays in Rust

This is the CLI tool for creating a suffix array using `libsufr`.

## Usage

Use `-h|--help` to view documentation:

```
$ sufr -h
Usage: sufr [COMMAND]

Commands:
  create   Create suffix array
  check    Check correctness of suffix array/LCP
  extract  Read suffix array and extract sequences
  help     Print this message or the help of the given subcommand(s)

Options:
  -h, --help  Print help
```

### Create a Suffix Array

Use the `create` action to generate the sorted suffix/LCP arrays from FASTA/Q input:

```
$ sufr create -h
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

For example:

```
$ sufr create sufr/tests/inputs/1.fa --log info
[2024-09-26T17:16:30Z INFO  sufr] Using 8 threads
[2024-09-26T17:16:30Z INFO  sufr] Read raw input of len 11 in 1.835292ms
[2024-09-26T17:16:30Z INFO  libsufr] Created unsorted suffix array of len 11 in 500ns
[2024-09-26T17:16:30Z INFO  libsufr] Selected/sorted 11 pivots in 42.625µs
[2024-09-26T17:16:30Z INFO  libsufr] Split into 16 partitions (avg 0) in 25.208µs
[2024-09-26T17:16:30Z INFO  libsufr] Sorted partitions in 20.042µs
[2024-09-26T17:16:30Z INFO  libsufr] Concatenated partitions in 3.292µs
[2024-09-26T17:16:30Z INFO  libsufr] Fixed LCP boundaries in 292ns
[2024-09-26T17:16:30Z INFO  libsufr] Total time to create suffix array: 115.25µs
[2024-09-26T17:16:30Z INFO  sufr] Wrote 154 bytes to '1.sufr' in 751.708µs
```

### Check Output

The resulting _.sufr_ file is a binary-encoded representation of the suffix/LCP arrays and the original sequence and other metadata from the `create` action.
Use the `check` action to verify that the suffix array is correctly sorted and that the LCP values are accurate:

```
$ sufr check -h
Check correctness of suffix array/LCP

Usage: sufr check [OPTIONS] <SUFR>

Arguments:
  <SUFR>  Sufr file

Options:
  -v, --verbose  List errors
  -h, --help     Print help
  -V, --version  Print version
```

For example:

```
$ sufr check 1.sufr
Found 0 errors in suffix array.
Found 0 errors in LCP
Finished checking in 32.958µs.
```

### Extract Suffixes

The `extract` action will display a range of suffixes from the _.sufr_ file:

```
$ sufr extract -h
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

For example:

```
$ sufr extract 1.sufr -e 1-10
#
ACGT#
ACGTNNACGT#
CGT#
CGTNNACGT#
GT#
GTNNACGT#
NACGT#
NNACGT#
T#
```

## Authors

* Ken Youens-Clark <kyclark@arizona.edu>
* Jack Roddy <jroddy@pharmacy.arizona.edu>
* Travis Wheeler <twheeler@arizona.edu>
