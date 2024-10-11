# Parallel Construction of Suffix Arrays in Rust

This is the CLI tool for creating a suffix array using `libsufr`.

## Usage

Use `-h|--help` to view documentation:

```
$ sufr --help
Parallel Construction of Suffix Arrays in Rust

Usage: sufr [OPTIONS] [COMMAND]

Commands:
  create   Create suffix array
  check    Check correctness of suffix array/LCP
  extract  Read suffix array and extract sequences
  search   Search a suffix array
  help     Print this message or the help of the given subcommand(s)

Options:
  -l, --log <LOG>            Log level [possible values: info, debug]
      --log-file <LOG_FILE>  Log file
  -h, --help                 Print help
  -V, --version              Print version
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
  -h, --help                        Print help
  -V, --version                     Print version
```

For example, using the human T2T chr1 (248M bp):

```
/usr/bin/time -l sufr --log info create --dna -n 64 chr1.fa
[2024-10-11T20:24:20Z INFO  sufr] Using 8 threads
[2024-10-11T20:24:20Z INFO  sufr] Read raw input of len 248,387,329 in 258.480792ms
[2024-10-11T20:24:20Z INFO  libsufr] Selected 639 pivots in 535.25µs
[2024-10-11T20:24:25Z INFO  libsufr] Wrote unsorted partitions in 5.393286208s
[2024-10-11T20:24:41Z INFO  libsufr] Sorted 64 partitions (avg 3881052) in 15.999000708s
[2024-10-11T20:24:43Z INFO  sufr] Wrote 2,235,486,019 bytes to 'chr1.sufr' in 1.407683834s
[2024-10-11T20:24:43Z INFO  sufr] Total time: 23.106772459s
       23.32 real       143.65 user        10.12 sys
          1177272320  maximum resident set size
                   0  average shared memory size
                   0  average unshared data size
                   0  average unshared stack size
              245897  page reclaims
                   4  page faults
                   0  swaps
                   0  block input operations
                   0  block output operations
                   0  messages sent
                   0  messages received
                   0  signals received
                4911  voluntary context switches
             1386865  involuntary context switches
        794098220487  instructions retired
        497221392750  cycles elapsed
           966626944  peak memory footprint
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
$ sufr check chr1.sufr
Checked 248,387,329 suffixes, found 0 errors in suffix array.
Finished checking in 4.4269505s.
```

### Extract Suffixes

The `extract` action will display a range of suffixes from the _.sufr_ file:

```
$ sufr extract --help
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

For example, to select the one-millionth suffix from chr1 with a maximum length of 30 bp:

```
$ sufr extract -n -e 1000000 -m 30 chr1.sufr
  16003693: AAAAAGAAAACACATCATGTAACTGACTGT
```

### Search

The `search` action will look for a query string in the _.sufr_ file:

```
$ sufr search --help
Search a suffix array

Usage: sufr search <QUERY> <SUFR>

Arguments:
  <QUERY>  Query
  <SUFR>   Sufr file

Options:
  -h, --help     Print help
  -V, --version  Print version
```

For instance:

```
$ sufr search AAAAAGAAAACACATCATGTAACTGACTGT chr1.sufr
Query 'AAAAAGAAAACACATCATGTAACTGACTGT' found in range 1000000..1000001 in 413.042µs
```

## Authors

* Ken Youens-Clark <kyclark@arizona.edu>
* Jack Roddy <jroddy@pharmacy.arizona.edu>
* Travis Wheeler <twheeler@arizona.edu>
