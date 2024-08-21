use anyhow::Result;
use clap::Parser;
use std::{
    cmp::min,
    fs::{self, File},
    io::Read,
    mem,
    path::PathBuf,
    time::Instant,
};
use substring::Substring;

/// Read suffix array
#[derive(Debug, Parser)]
#[command(author, version, about)]
pub struct Args {
    /// Suffix array filea
    #[arg(short('a'), long, value_name = "SA")]
    pub sa: PathBuf,

    /// Sequence file
    #[arg(short, long, value_name = "SEQ")]
    pub sequence: PathBuf,

    /// Extract positions
    #[arg(
        short, 
        long, 
        value_name = "EXTRACT", 
        value_parser = clap::value_parser!(u64).range(0..),
        default_value = "2",
        num_args(0..)
    )]
    pub extract: Vec<u64>,
}

// --------------------------------------------------
fn main() {
    if let Err(e) = run(Args::parse()) {
        eprintln!("Error: {e}");
        std::process::exit(1);
    }
}

// --------------------------------------------------
pub fn run(args: Args) -> Result<()> {
    let mut sa_file = File::open(args.sa)?;

    // The first 64-bits of the file contain the size of the SA
    let mut buffer = [0; 8];
    sa_file.read_exact(&mut buffer)?;
    //dbg!(&buffer);

    // Convert the Vec<u8> to a usize
    let sa_len = usize::from_ne_bytes(buffer);
    println!("Length of SA: {sa_len}");

    let start = Instant::now();
    // Allocate a buffer to hold the data
    let mut buffer2 = vec![0u8; sa_len * mem::size_of::<u32>()];
    //let mut buffer2 = vec![0u8; sa_len * mem::size_of::<usize>()];

    // Read the bytes into the buffer
    sa_file.read_exact(&mut buffer2)?;

    // Convert the buffer into a slice of i32 integers
    let integers: &[u32] = unsafe {
        std::slice::from_raw_parts(buffer2.as_ptr() as *const u32, sa_len)
    };

    //let integers: &[usize] = unsafe {
    //    std::slice::from_raw_parts(buffer2.as_ptr() as *const usize, sa_len)
    //};

    println!("SA construction finished in {:?}", start.elapsed());
    println!("{integers:?}");

    // Print the integers
    //for &integer in integers {
    //    println!("{}", integer);
    //}

    //let seq: Vec<u8> = fs::read_to_string(args.sequence)?.bytes().collect();
    let start = Instant::now();
    let seq = fs::read_to_string(args.sequence)?.to_uppercase();
    let seq = format!("{}$", seq.trim());
    println!("Read sequence finished in {:?}", start.elapsed());

    //for start in [0, 100_000, 2_000_000, 3_000_000] {
    for &start in &args.extract {
        let start = start as usize;
        if let Some(&pos) = integers.get(start) {
            let pos = pos as usize;
            //let stop = pos + (sa_len - pos);
            //let len = min(sa_len - 1 - start, 8);
            //println!("start = {pos}, len = {len}");
            //let stop = pos + len;
            //println!("start = {pos}, stop = {stop}");
            //println!(
            //    "({start}) {pos} = {}",
            //    String::from_utf8_lossy(&seq[pos..stop])
            //);
            println!("({start}) {pos} = {}", seq.substring(pos, pos + 8));
        }
    }

    //if let Some(later) = integers.get(2_000_000) {
    //    let start = *later as usize;
    //    println!(
    //        "later = {}",
    //        String::from_utf8_lossy(&seq[start..start + 8])
    //    );
    //}

    // Read the rest of the file into memory
    //let mut contents: Vec<u8> = vec![];
    //sa_file.read_to_end(&mut contents)?;
    //dbg!(&contents);

    // Determine how many bits to take
    //let mult = if (n as u64) < u32::MAX as u64 { 4 } else { 8 };
    //dbg!(&mult);

    //let raw_sa: Vec<_> = contents.drain(0..(n * mult)).collect();
    //dbg!(&raw_sa[..20]);

    //let chunked: Vec<_> = raw_sa.chunks(mult).collect();
    //dbg!(&chunked[..20]);

    //let converted: Vec<_> = chunked
    //    .iter()
    //    .map(|&v| usize::from_ne_bytes(v.into()))
    //    .collect();

    //.chunk_by(|chunk| usize::from_ne_bytes(chunk))
    //.collect();

    ////dbg!(&seq);

    //for pos in suffix_array {
    //    let p = pos as usize;
    //    println!("{pos} = {}", seq[p]);
    //    break;
    //}

    //let mut lcp: Vec<u8> = vec![0; *n as usize];
    //file.read_exact(&mut lcp)?;
    //dbg!(&lcp);

    Ok(())
}
