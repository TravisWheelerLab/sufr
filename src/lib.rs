mod suffix_array;

use anyhow::{anyhow, bail, Result};
use clap::{builder::PossibleValue, Parser, ValueEnum};
use format_num::NumberFormat;
use log::{debug, info};
use needletail::parse_fastx_file;
use regex::Regex;
use std::{
    ffi::OsStr,
    fmt::Debug,
    fs::File,
    io::{BufWriter, Read},
    mem,
    num::NonZeroUsize,
    ops::Range,
    path::PathBuf,
    time::Instant,
};
//use substring::Substring;
use suffix_array::{FromUsize, Int, SuffixArrayBuilder};
//use u4::{AsNibbles, U4x2, U4};

// --------------------------------------------------
#[derive(Parser, Debug)]
#[command(arg_required_else_help = true)]
pub struct Cli {
    #[command(subcommand)]
    pub command: Option<Command>,
}

#[derive(Parser, Debug)]
pub enum Command {
    /// Create suffix array
    #[clap(alias = "ch")]
    Check(CheckArgs),

    /// Create suffix array
    #[clap(alias = "cr")]
    Create(CreateArgs),

    /// Read suffix array and extract sequences
    #[clap(alias = "re")]
    Read(ReadArgs),
}

#[derive(Debug, Parser)]
#[command(author, version, about)]
pub struct CheckArgs {
    /// Suffix array file
    #[arg(short, long, value_name = "SA")]
    pub array: String,

    /// Sequence file
    #[arg(short, long, value_name = "SEQ")]
    pub sequence: String,
}

#[derive(Debug, Parser)]
#[command(author, version, about)]
pub struct CreateArgs {
    /// Input file
    #[arg(value_name = "INPUT")]
    pub input: String,

    /// Subproblem count
    #[arg(short, long, value_name = "NUM_PARTS", default_value = "16")]
    pub num_partitions: usize,

    /// Max context
    #[arg(short, long, value_name = "CONTEXT")]
    pub max_context: Option<usize>,

    /// Number of threads
    #[arg(short, long, value_name = "THREADS", default_value = "16")]
    pub threads: Option<usize>,

    /// Output file
    #[arg(short, long, value_name = "OUTPUT")]
    pub output: Option<String>,

    /// Input is DNA, ignore sequences starting with 'N'
    #[arg(short('d'), long("dna"))]
    pub is_dna: bool,

    /// Verify order
    #[arg(short, long)]
    pub check: bool,

    /// Log level
    #[arg(short, long)]
    pub log: Option<LogLevel>,

    /// Log file
    #[arg(long)]
    pub log_file: Option<String>,
}

/// Read suffix array
#[derive(Debug, Parser)]
#[command(author, version, about)]
pub struct ReadArgs {
    /// Sufr file
    #[arg(value_name = "FILE")]
    pub filename: String,

    /// Maximum length of sequence
    #[arg(short, long, value_name = "MAX", default_value = "0")]
    pub max_len: usize,

    /// Extract positions
    #[arg(short, long, value_name = "EXTRACT", default_value = "1")]
    pub extract: String,

    /// Number output
    #[arg(short, long)]
    pub number: bool,

    /// Output
    #[arg(short, long)]
    pub output: Option<String>,
}

#[derive(Debug, Clone)]
pub enum LogLevel {
    Info,
    Debug,
}

impl ValueEnum for LogLevel {
    fn value_variants<'a>() -> &'a [Self] {
        &[LogLevel::Info, LogLevel::Debug]
    }

    fn to_possible_value<'a>(&self) -> Option<PossibleValue> {
        Some(match self {
            LogLevel::Info => PossibleValue::new("info"),
            LogLevel::Debug => PossibleValue::new("debug"),
        })
    }
}

type PositionList = Vec<Range<usize>>;

#[derive(Debug)]
struct SuffixArray {
    version: u8,
    is_dna: u8,
    len: usize,
    suffix_array: Vec<u32>,
    lcp: Vec<u32>,
}

//#[derive(Debug)]
//pub struct PackedSeq {
//    data: Vec<u8>,
//}

//impl PackedSeq {
//    pub fn new(text: Vec<u8>) -> Self {
//        let u4s: Vec<_> = text.iter().filter_map(Self::byte_to_u4).collect();
//        let mut data: Vec<u8> = vec![0u8; text.len() / 2];
//        AsNibbles(&mut data).pack_from_slice(&u4s);
//        Self { data }
//    }

//    pub fn values(&self) -> AsNibbles<&Vec<u8>> {
//        AsNibbles(&self.data)
//    }

//    pub fn get(&self, index: usize) -> U4 {
//        let idx = index / 2;
//        let rem = index % 2;
//        let val = U4x2::from_byte(self.data[idx]);
//        if rem == 0 {
//            val.left()
//        } else {
//            val.right()
//        }
//    }

//    fn byte_to_u4(byte: &u8) -> Option<U4> {
//        match byte {
//            36 => Some(U4::Dec00), // $
//            65 => Some(U4::Dec01), // A
//            66 => Some(U4::Dec02), // B
//            67 => Some(U4::Dec03), // C
//            68 => Some(U4::Dec04), // D
//            71 => Some(U4::Dec05), // G
//            72 => Some(U4::Dec06), // H
//            75 => Some(U4::Dec07), // K
//            77 => Some(U4::Dec08), // M
//            78 => Some(U4::Dec09), // N
//            82 => Some(U4::Dec10), // R
//            83 => Some(U4::Dec11), // S
//            84 => Some(U4::Dec12), // T
//            86 => Some(U4::Dec13), // V
//            87 => Some(U4::Dec14), // W
//            89 => Some(U4::Dec15), // Y
//            _ => None,             // Anything else
//        }
//    }

//    fn u4_to_byte(val: U4) -> u8 {
//        match val {
//            U4::Dec00 => 36,
//            U4::Dec01 => 65,
//            U4::Dec02 => 66,
//            U4::Dec03 => 67,
//            U4::Dec04 => 68,
//            U4::Dec05 => 71,
//            U4::Dec06 => 72,
//            U4::Dec07 => 75,
//            U4::Dec08 => 77,
//            U4::Dec09 => 78,
//            U4::Dec10 => 82,
//            U4::Dec11 => 83,
//            U4::Dec12 => 84,
//            U4::Dec13 => 86,
//            U4::Dec14 => 87,
//            U4::Dec15 => 89,
//        }
//    }
//}

// --------------------------------------------------
pub fn check(args: &CheckArgs) -> Result<()> {
    let len = read_suffix_array_len(&args.array)?;
    dbg!(len);
    //let sa = if len < u32::MAX as u64 {
    //    let sa = read_suffix_array<u32>(&args.array)?;
    //    sa
    //} else {
    //    let sa = read_suffix_array<u64>(&args.array)?;
    //    sa
    //};
    //dbg!(&sa);
    //let sa_len = sa.len();
    //let seq = read_input(open(&args.sequence)?);

    //if sa_len != seq.len() {
    //    bail!("SA len {sa_len} does not match sequence len {}", seq.len());
    //}

    //let start = Instant::now();
    //let mut previous: Option<usize> = None;
    //let mut num_errors = 0;
    //for (i, &cur) in sa.iter().enumerate() {
    //    if let Some(p) = previous {
    //        if !is_less(&seq[p..sa_len], &seq[cur..sa_len]) {
    //            num_errors += 1;
    //            println!("POS {}", i + 1);
    //        }
    //    }
    //    previous = Some(cur);
    //}

    //println!(
    //    "Found {num_errors} error{} in {:?}.",
    //    if num_errors == 1 { "" } else { "s" },
    //    start.elapsed()
    //);

    Ok(())
}

// --------------------------------------------------
pub fn read_sequence_input(filename: &str) -> Result<(Vec<u8>, Vec<String>)> {
    let mut reader = parse_fastx_file(filename)?;
    let mut seqs = vec![];
    let mut headers = vec![];
    let mut i = 0;
    while let Some(rec) = reader.next() {
        let rec = rec?;
        if i > 0 {
            // Sequence delimiter
            seqs.push(b'?');
        }
        // Uppercase (mask w/32)
        let mut seq: Vec<u8> =
            rec.seq().iter().map(|b| b & 0b1011111).collect();
        seqs.append(&mut seq);
        i += 1;

        headers.push(String::from_utf8(rec.id().to_vec())?);
    }

    // File delimiter
    seqs.push(b'#');

    Ok((seqs, headers))
}

// --------------------------------------------------
//pub fn read_input_u4(input: impl BufRead) -> PackedSeq {
//    //text.push(b'$');
//    let mut text: Vec<u8> = input
//        .bytes()
//        .map_while(Result::ok)
//        .filter(|&b| (0b1000001..=0b1011010).contains(&b)) // A-Z
//        .map(|b| b & 0b1011111) // Uppercase (mask w/32)
//        .collect();
//    text.push(b'$');

//    // The len of text must be even b/c we pack 2 u4's into each each u8
//    if text.len() % 2 == 1 {
//        text.push(b'$');
//    }
//    PackedSeq::new(text)
//}

// --------------------------------------------------
pub fn create(args: &CreateArgs) -> Result<()> {
    // Read sequence input
    let now = Instant::now();
    let (text, _headers) = read_sequence_input(&args.input)?;
    let len = text.len() as u64;
    let num_fmt = NumberFormat::new();
    info!(
        "Read raw input of len {} in {:?}",
        num_fmt.format(",.0", len as f64),
        now.elapsed()
    );
    debug!("Raw input '{:?}'", text);

    //println!("{:?}", String::from_utf8(text.clone()));
    //println!("{}", headers.join(", "));
    //let packed = read_input_u4(&mut BufReader::new(File::open(&args.input)?));
    //println!("Input as U4 {:?}", &packed);
    //for (i, val) in packed.values().iter().enumerate() {
    //    println!(
    //        "{i}: {val:04b} => {} => {}",
    //        PackedSeq::u4_to_byte(val),
    //        packed.get(i),
    //    );
    //}

    if len < u32::MAX as u64 {
        let sa: SuffixArrayBuilder<u32> = SuffixArrayBuilder::new(
            text,
            len as u32,
            args.max_context.map(|val| val as u32),
            args.is_dna,
        );
        _create(sa, args)?;
    } else {
        let sa: SuffixArrayBuilder<u64> = SuffixArrayBuilder::new(
            text,
            len,
            args.max_context.map(|val| val as u64),
            args.is_dna,
        );
        _create(sa, args)?;
    }

    Ok(())
}

// --------------------------------------------------
pub fn _create<T>(sa: SuffixArrayBuilder<T>, args: &CreateArgs) -> Result<()>
where
    T: Int + FromUsize<T> + Sized + Send + Sync + Debug,
{
    env_logger::Builder::new()
        .filter_level(match args.log {
            Some(LogLevel::Debug) => log::LevelFilter::Debug,
            Some(LogLevel::Info) => log::LevelFilter::Info,
            _ => log::LevelFilter::Off,
        })
        .target(match args.log_file {
            // Optional log file, default to STDOUT
            Some(ref filename) => env_logger::Target::Pipe(Box::new(
                BufWriter::new(File::create(filename)?),
            )),
            _ => env_logger::Target::Stdout,
        })
        .init();

    if let Some(num) = args.threads {
        info!("Using {num} thread{}", if num == 1 { "" } else { "s" });
        rayon::ThreadPoolBuilder::new()
            .num_threads(num)
            .build_global()
            .unwrap();
    }

    let total_start = Instant::now();

    // Sort
    let (sorted_sa, lcp) = sa.sort(args.num_partitions);
    info!("Suffix generated in {:?}s", total_start.elapsed());

    // Check
    if args.check {
        debug!("Sorted = {:?}", &sorted_sa);
        debug!(
            "Suffixes = {:#?}",
            sorted_sa
                .iter()
                .map(|v| sa.string_at(v.to_usize()))
                .collect::<Vec<_>>()
        );
        debug!("LCP = {lcp:#?}");
        let now = Instant::now();
        let num_errors = sa.check_order(&sorted_sa).len();
        info!(
            "Checked order, found {num_errors} error{} in {:?}",
            if num_errors == 1 { "" } else { "s" },
            now.elapsed()
        );
        let lcp_errors = sa.check_lcp(&sorted_sa, &lcp);
        let num_errors = lcp_errors.len();
        info!(
            "Checked LCP, found {num_errors} error{} in {:?}",
            if num_errors == 1 { "" } else { "s" },
            now.elapsed()
        );
        if num_errors > 0 {
            debug!("{lcp_errors:?}");
        }
    }

    // Write out suffix array length
    let now = Instant::now();
    let outfile = &args.output.clone().unwrap_or(format!(
        "{}.sufr",
        PathBuf::from(&args.input)
            .file_stem()
            .unwrap_or(OsStr::new("out"))
            .to_string_lossy()
    ));
    sa.write(&sorted_sa, &lcp, outfile)?;
    info!("Wrote output file '{outfile}' in {:?}", now.elapsed());

    Ok(())
}

// --------------------------------------------------
// Parse an index from a string representation of an integer.
// Ensures the number is non-zero.
// Ensures the number does not start with '+'.
// Returns an index, which is a non-negative integer that is
// one less than the number represented by the original input.
#[allow(dead_code)]
fn parse_index(input: &str) -> Result<usize> {
    let value_error = || anyhow!(r#"illegal list value: "{input}""#);
    input
        .starts_with('+')
        .then(|| Err(value_error()))
        .unwrap_or_else(|| {
            input
                .parse::<NonZeroUsize>()
                .map(|n| usize::from(n) - 1)
                .map_err(|_| value_error())
        })
}

// --------------------------------------------------
#[allow(dead_code)]
fn parse_pos(range: &str) -> Result<PositionList> {
    let range_re = Regex::new(r"^(\d+)-(\d+)$").unwrap();
    range
        .split(',')
        .map(|val| {
            parse_index(val).map(|n| n..n + 1).or_else(|e| {
                range_re.captures(val).ok_or(e).and_then(|captures| {
                    let n1 = parse_index(&captures[1])?;
                    let n2 = parse_index(&captures[2])?;
                    if n1 >= n2 {
                        bail!(
                            "First number in range ({}) \
                            must be lower than second number ({})",
                            n1 + 1,
                            n2 + 1
                        );
                    }
                    Ok(n1..n2 + 1)
                })
            })
        })
        .collect::<Result<_, _>>()
        .map_err(From::from)
}

// --------------------------------------------------
pub fn read(args: &ReadArgs) -> Result<()> {
    let now = Instant::now();
    let sufr_file = _read(&args.filename)?;
    dbg!(sufr_file);

    //let sa = read_suffix_array(&args.array)?;
    //let sa_len = sa.len();
    //info!("Read SA in {:?}", start.elapsed());

    //let (seq, _headers) = read_input(&args.sequence)?;

    //if seq.len() != sa_len {
    //    bail!("SA len {sa_len} does not match sequence len {}", seq.len());
    //}

    //let positions: Vec<_> = parse_pos(&args.extract)?
    //    .into_iter()
    //    .flat_map(|r| r.collect::<Vec<_>>())
    //    .collect();

    //let mut output: Box<dyn Write> = match &args.output {
    //    Some(out_name) => Box::new(File::create(out_name)?),
    //    _ => Box::new(io::stdout()),
    //};

    //for start in positions {
    //    if let Some(&pos) = sa.get(start) {
    //        //let pos = pos as usize;
    //        let end = if args.max_len > 0 {
    //            pos + args.max_len
    //        } else {
    //            sa_len
    //        };

    //        if args.number {
    //            writeln!(
    //                output,
    //                "{:3}: {}",
    //                start + 1,
    //                seq.substring(pos, end)
    //            )?;
    //        } else {
    //            writeln!(output, "{}", seq.substring(pos, end))?;
    //        }
    //    }
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

// --------------------------------------------------
fn read_suffix_array_len(filename: &str) -> Result<u64> {
    // The first 64-bits of the file contain the size of the SA
    let mut sa_file =
        File::open(filename).map_err(|e| anyhow!("{filename}: {e}"))?;
    let mut buffer = [0; 8];
    sa_file.read_exact(&mut buffer)?;
    Ok(u64::from_ne_bytes(buffer))
}

// --------------------------------------------------
fn _read(filename: &str) -> Result<SuffixArray> {
    //fn read_suffix_array<T>(filename: &str) -> Result<Vec<T>>
    //where
    //    T: Int + FromUsize<T> + Sized + Send + Sync + Debug,
    //{
    let mut file =
        File::open(filename).map_err(|e| anyhow!("{filename}: {e}"))?;

    let mut buffer = [0; 2];
    file.read_exact(&mut buffer)?;
    let version = buffer[0];
    let is_dna = buffer[1];

    // The next 64-bits of the file contain the size of the SA
    let mut buffer = [0; 8];
    file.read_exact(&mut buffer)?;

    // Convert the Vec<u8> to a usize
    let len = usize::from_ne_bytes(buffer);

    // Allocate a buffer to hold the data
    //let suffix_array = if sa_len < u32::MAX as usize {
    //    let mut buffer = vec![0u8; sa_len * mem::size_of::<u32>()];
    //    sa_file.read_exact(&mut buffer)?;
    //    let sa: &[u32] = unsafe {
    //        std::slice::from_raw_parts(buffer.as_ptr() as *const u32, sa_len)
    //    };
    //    sa
    //} else {
    //    let mut buffer = vec![0u8; sa_len * mem::size_of::<u64>()];
    //    sa_file.read_exact(&mut buffer)?;
    //    let sa: &[u64] = unsafe {
    //        std::slice::from_raw_parts(buffer.as_ptr() as *const u64, sa_len)
    //    };
    //    sa
    //};

    // Read the bytes into the buffer
    let mut buffer = vec![0u8; len * mem::size_of::<u32>()];
    file.read_exact(&mut buffer)?;
    let suffix_array: Vec<u32> = unsafe {
        std::slice::from_raw_parts(buffer.as_ptr() as *const u32, len)
            .to_vec()
    };

    let mut buffer = vec![0u8; len * mem::size_of::<u32>()];
    file.read_exact(&mut buffer)?;
    let lcp: Vec<u32> = unsafe {
        std::slice::from_raw_parts(buffer.as_ptr() as *const u32, len)
            .to_vec()
    };

    // How can I have either 32 or 64 vectors?
    // Convert the buffer into a slice of i32 integers
    //let suffix_array: &[u32] = unsafe {
    //    std::slice::from_raw_parts(buffer.as_ptr() as *const u32, sa_len)
    //};

    //let suffix_array: &[usize] = unsafe {
    //    std::slice::from_raw_parts(buffer.as_ptr() as *const usize, sa_len)
    //};

    Ok(SuffixArray {
        version,
        is_dna,
        len,
        suffix_array,
        lcp,
    })
}
