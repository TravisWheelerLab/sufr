# Parallel Construction of Suffix Arrays in Rust

See https://docs.rs/libsufr/latest/libsufr/ for API documentation.

This is a Rust library crate for creating suffix arrays.

The basic ideas are as follow:

* Read the input file as `u8` (unsigned 8-bit integer values). 
* Select the suffixes, which are normally 0 to the length of the text but there is the option to skip suffixes that don't start with A/C/G/T if the input is DNA. Note: Following the C++ implementation, we use 32-bit integers if the input text length is less than 2^32 and 64-bit integers, otherwise.
* Create partitions by randomly choosing suffixes and copying the suffixes to the highest possible partition bounded by any pivot.
* Sort the partitions in parallel.
* Merge the partitions. Because the values fall into nonoverlapping ranges, these subarrays can be appended in order to produce the final SA.
* Produce a binary-encoded output file with the suffix/LCP arrays and other metadata.

Some advantages to this algorithm:

* The various partitioned subarrays can be processed independently by separate threads, and no thread will ever have to merge the entire input.
* Suffix comparisons are made faster by caching LCPs.
* Using `u8` for the input text and 32-bits (when possible) for SA/LCP results in lower memory usage.

See the repository for documentation: https://github.com/TravisWheelerLab/sufr

## See Also

Use **`cargo install sufr`** for a CLI.

## Authors

* Ken Youens-Clark <kyclark@arizona.edu>
* Jack Roddy <jroddy@pharmacy.arizona.edu>
* Travis Wheeler <twheeler@arizona.edu>
