use anyhow::{anyhow, Result};
use clap::{builder::PossibleValue, Parser, ValueEnum};
use log::{debug, info};
use std::{
    cmp::{max, min},
    //collections::HashMap,
    fs::{self, File},
    //io::{BufRead, BufWriter, Write},
    io::{BufRead, BufWriter},
    mem,
    path::PathBuf,
    //slice,
    time::Instant,
};
use substring::Substring;

// --------------------------------------------------
/// Suffix array
#[derive(Debug, Parser)]
#[command(author, version, about)]
pub struct Args {
    /// Input file
    #[arg(value_name = "INPUT")]
    pub input: PathBuf,

    /// Search pattern
    #[arg(value_name = "PATTERN")]
    pub pattern: Option<String>,

    /// Subproblem count
    #[arg(short, long, value_name = "SUBPROBLEMS", default_value = "8192")]
    pub subproblem_count: usize,

    /// Max context
    #[arg(short, long, value_name = "CONTEXT")]
    pub max_context: Option<usize>,

    /// Output file
    #[arg(short, long, value_name = "OUTPUT", default_value = "sufr.sa")]
    pub output: PathBuf,

    /// Log level
    #[arg(short, long)]
    pub log: Option<LogLevel>,

    /// Log file
    #[arg(long)]
    pub log_file: Option<String>,
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

struct SuffixArray {
    pub text: Vec<u8>,
    pub len: usize,
    pub max_context: usize,
}

impl SuffixArray {
    pub fn new(
        input: &PathBuf,
        max_context: Option<usize>,
    ) -> Result<SuffixArray> {
        let mut text: Vec<u8> = fs::read_to_string(input)
            .map_err(|e| anyhow!("{}: {e}", input.display()))?
            .trim()
            .to_uppercase()
            .bytes()
            .collect();
        text.push(b'$');
        let len = text.len();

        Ok(SuffixArray {
            text,
            len,
            max_context: max_context.unwrap_or(len),
        })
    }

    pub fn merge_sort(
        self: &Self,
        x: &mut [usize],
        y: &mut [usize],
        n: usize,
        lcp: &mut [usize],
        lcp_w: &mut [usize],
    ) {
        debug!("x {x:?} y {y:?} n {n} lcp {lcp:?} lcp_w {lcp_w:?}");
        if n == 1 {
            lcp[0] = 0;
        } else {
            let mid = n / 2;
            self.merge_sort(
                &mut y[..mid],
                &mut x[..mid],
                mid,
                &mut lcp_w[..mid],
                &mut lcp[..mid],
            );

            self.merge_sort(
                &mut y[mid..],
                &mut x[mid..],
                n - mid,
                &mut lcp_w[mid..],
                &mut lcp[mid..],
            );

            self.merge(x, mid, lcp_w, y, lcp);
        }
        dbg!(x);
        dbg!(y);
        dbg!(lcp_w);
    }

    //) -> (Vec<usize>, Vec<usize>, Vec<usize>, Vec<usize>) {
    fn merge(
        self: &Self,
        suffix_array: &mut [usize],
        mid: usize,
        lcp_w: &mut [usize],
        z: &mut [usize],
        lcp_z: &mut [usize],
    ) {
        debug!("\n>>> MERGE {suffix_array:?} mid {mid} <<<");
        let mut x = suffix_array[..mid].to_vec();
        let mut y = suffix_array[mid..].to_vec();
        let mut lcp_x = lcp_w[..mid].to_vec();
        let mut lcp_y = lcp_w[mid..].to_vec();
        debug!("LCP_x {lcp_x:?}");
        debug!("LCP_y {lcp_y:?}");
        let mut len_x = x.len();
        let mut len_y = y.len();

        // Bookkeeping
        let mut m = 0;
        let mut i = 0; // Index into x
        let mut j = 0; // Index into y
        let mut k = 0; // Index into z

        // LCP of last compared pair
        // In the original code, l_x is intially set to null
        //#[allow(unused_assignments)]
        //let mut l_x = 0; // LCP(X_i, X_{i-1})

        //debug!("MERGE x {x:?} y {y:?} z {z:?} len_x {len_x} len_y {len_y}");
        //debug!("MERGE lcp_x {lcp_x:?} lcp_y {lcp_y:?} lcp_z {lcp_z:?}");

        //let mut num_swaps = 0;
        while i < len_x && j < len_y {
            debug!("\n>>> x {} y {} <<<", x[i], y[j]);
            let l_x = lcp_x[i];

            if l_x > m {
                debug!("BRANCH 1 i '{i}' l_x '{l_x}' > m '{m}'");
                debug!("sa    {suffix_array:?}");
                debug!("lcp_w {lcp_w:?}");
                debug!("lcp_x {lcp_x:?}");
                debug!("lcp_y {lcp_y:?}");
                z[k] = x[i];
                debug!(">> k1 = '{k}' TW <<");
                lcp_z[k] = l_x;
            } else if l_x < m {
                debug!("BRANCH 2 l_x '{l_x}' < m '{m}'");
                z[k] = y[j];
                debug!(">> k2 = '{k}' TW <<");
                lcp_z[k] = m;
                m = l_x;
            } else {
                debug!("BRANCH 3");
                // Length of shorter suffix
                let max_n = self.len - max(x[i], y[j]);

                // Prefix-context length for the suffixes
                let context = min(self.max_context, max_n);

                // LCP(X_i, Y_j)
                let this_lcp =
                    lcp(&self.text[x[i]..], &self.text[y[j]..], context - m);
                //debug!("CALLING LCP x[i] = {} y[j] = {}", x[i], y[j]);
                let n = m + lcp(
                    &self.text[(x[i] + m)..],
                    &self.text[(y[j] + m)..],
                    context - m,
                );
                debug!("m       '{m}'");
                debug!("lcp     '{this_lcp}'");
                debug!("n       '{n}'");
                debug!("max_n   '{max_n}'");
                debug!("context '{context}'");
                debug!(
                    "left  {} {:?}",
                    x[i],
                    String::from_utf8(self.text[x[i]..].to_vec())
                );
                debug!(
                    "right {} {:?}",
                    y[j],
                    String::from_utf8(self.text[y[j]..].to_vec())
                );
                //debug!("LCP = {n}, max_n = {max_n}, k = {k}, n = {n}");

                // If the len of the LCP is the entire shorter
                // sequence, take that.
                // Else, look at the next char after the LCP
                // to determine order.
                //z[k] = if n == max_n {
                //    max(x[i], y[j])
                //} else if self.text[x[i] + n] < self.text[y[j] + n] {
                //    x[i]
                //} else {
                //    y[j]
                //};

                if n == max_n {
                    z[k] = max(x[i], y[j])
                } else if self.text[x[i] + n] < self.text[y[j] + n] {
                    debug!(
                        "Checking next chars '{}' and '{}'",
                        self.text[x[i] + n] as char,
                        self.text[y[j] + n] as char
                    );
                    z[k] = x[i]
                } else {
                    z[k] = y[j]
                }

                //lcp_z[k] = if z[k] == x[i] { l_x } else { m };
                if z[k] == x[i] {
                    lcp_z[k] = l_x;
                    debug!(">> k3 = '{k}' TW <<");
                } else {
                    debug!(">> m = '{m}' TW <<");
                    debug!(">> n = '{n}' TW <<");
                    debug!(">> k3 = '{k}' TW <<");
                    lcp_z[k] = m
                }
                m = n;
            }
            debug!("x {} y {} => z {}", x[i], y[j], z[k]);

            if z[k] == x[i] {
                i += 1;
            } else {
                j += 1;
                //num_swaps += 1;
                debug!("SWAP!");
                mem::swap(&mut x, &mut y);
                mem::swap(&mut len_x, &mut len_y);
                mem::swap(&mut lcp_x, &mut lcp_y);
                mem::swap(&mut i, &mut j);
            }

            k += 1;
        }

        // Copy rest of the data from X to Z.
        while i < len_x {
            z[k] = x[i];
            lcp_z[k] = lcp_x[i];
            i += 1;
            k += 1;
        }

        // Copy rest of the data from Y to Z.
        //dbg!(&len_y);
        //dbg!(&j);
        if j < len_y {
            debug!("FINISHING Y");
            debug!("m '{m}'");
            debug!("k '{k}'");
            z[k] = y[j];
            lcp_z[k] = m;
            j += 1;
            k += 1;
            while j < len_y {
                //dbg!(&k);
                //dbg!(&z);
                //dbg!(&y);
                debug!("..k '{k}' j '{j}' y '{}'", y[j]);
                z[k] = y[j];
                lcp_z[k] = lcp_y[j];
                j += 1;
                k += 1;
            }
        }

        debug!("z     {z:?}");
        debug!("LCP_z {lcp_z:?}");
    }
}

// --------------------------------------------------
pub fn run(args: Args) -> Result<()> {
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

    info!("Reading input {}", &args.input.display());

    //let file = BufReader::new(
    //    File::open(&args.input)
    //        .map_err(|e| anyhow!("{}: {e}", args.input.display()))?,
    //);

    let start = Instant::now();
    let suf_arr = SuffixArray::new(&args.input, args.max_context)?;
    info!("Reading/LCP finished in {:?}", start.elapsed());

    info!("Input length '{}'", suf_arr.len);
    debug!("Raw input '{:?}'", suf_arr.text);
    let text = String::from_utf8(suf_arr.text.clone())?;
    debug!("Input '{text}'");
    let mut sa: Vec<usize> = (0..suf_arr.len).collect();
    debug!("Initial SA '{:?}'", sa);

    let mut sa_w = sa.clone();
    let mut lcp = vec![0; suf_arr.len];
    let mut lcp_w = vec![0; suf_arr.len];

    let start = Instant::now();
    suf_arr.merge_sort(&mut sa_w, &mut sa, suf_arr.len, &mut lcp, &mut lcp_w);
    info!("merge_sort finished in {:?}", start.elapsed());
    //info!("Sorted SA = {:?}", sa);

    debug!("Suffixes of '{text}'");
    for &pos in &sa {
        debug!("pos {pos:4} suffix {}", text.substring(pos, suf_arr.len));
    }

    //debug!("SA = {:?}", sa);
    //debug!("SA_w = {:?}", sa_w);
    //debug!("LCP = {:?}", lcp);
    //debug!("LCP_w = {:?}", lcp_w);

    //let mut out = File::create(args.output)?;
    //let slice_u8: &[u8] = unsafe {
    //    slice::from_raw_parts(
    //        sa_w.as_ptr() as *const _,
    //        suf_arr.len * mem::size_of::<usize>(),
    //    )
    //};
    //out.write(&usize_to_bytes(suf_arr.len))?;
    //out.write_all(slice_u8)?;

    Ok(())
}

// --------------------------------------------------
fn lcp(s1: &[u8], s2: &[u8], len: usize) -> usize {
    for i in 0..len {
        if s1[i] != s2[i] {
            return i;
        }
    }
    len
}

// --------------------------------------------------
#[test]
fn test_lcp() -> Result<()> {
    assert_eq!(lcp(&[b'A'], &[b'C'], 1), 0);

    assert_eq!(lcp(&[b'A'], &[b'A'], 1), 1);

    assert_eq!(lcp(&[b'A', b'A'], &[b'A', b'A', b'C'], 2), 2);

    assert_eq!(lcp(&[b'A', b'C'], &[b'A', b'A', b'C'], 2), 1);

    assert_eq!(lcp(&[b'A', b'A'], &[b'A', b'A', b'C'], 1), 1);
    Ok(())
}

// --------------------------------------------------
//fn usize_to_bytes(value: usize) -> Vec<u8> {
//    // Determine the size of usize in bytes
//    let size = std::mem::size_of::<usize>();

//    // Create a vector to hold the bytes
//    let mut bytes = Vec::with_capacity(size);

//    // Convert usize to bytes
//    for i in 0..size {
//        bytes.push((value >> (i * 8)) as u8);
//    }

//    bytes
//}

// --------------------------------------------------
//#[inline]
//fn byte_to_base(byte: &u8) -> Option<u8> {
//    match byte {
//        65 | 97 => Some(0),  // A
//        67 | 99 => Some(1),  // C
//        71 | 103 => Some(3), // G
//        84 | 116 => Some(2), // T
//        _ => None,
//    }
//}

// --------------------------------------------------
/// Reads input file into an array of 2-bit encoded values
/// E.g., A = 0, C = 1, G = 2, T = 3
fn _read_input(mut input: impl BufRead) -> Result<Vec<u8>> {
    let mut ret = vec![];
    let start = Instant::now();
    let buffer_size = 1024 * 128;

    // Mask with 32 to uppercase ASCII values
    //let ret: Vec<_> = input
    //    .bytes()
    //    .map_while(Result::ok)
    //    .map(|b| b & 0b1011111)
    //    .collect();

    // 1.05537875s
    //let mut buffer = vec![0; buffer_size];
    //loop {
    //    let bytes_read = input.read(&mut buffer)?;
    //    if bytes_read == 0 {
    //        break;
    //    }
    //    ret.extend(
    //        buffer[..bytes_read].iter().filter_map(|b| byte_to_base(b)),
    //    );
    //}

    // 1.699109792s in debug, 1.031828417s in release
    let mut buffer = vec![0; buffer_size];
    loop {
        let bytes_read = input.read(&mut buffer)?;
        if bytes_read == 0 {
            break;
        }
        for byte in buffer[..bytes_read].iter() {
            if let Some(checked) = match byte {
                65 | 97 => Some(0),
                67 | 99 => Some(1),
                71 | 103 => Some(3),
                84 | 116 => Some(2),
                _ => None,
            } {
                ret.push(checked);
            }
        }
    }

    // 53.436009s in debug, 1.731694375s in release
    //let lookup = HashMap::<u8, u8>::from([
    //    (65, 0),  // A
    //    (97, 0),  // a
    //    (67, 1),  // C
    //    (99, 1),  // c
    //    (71, 3),  // G
    //    (103, 3), // g
    //    (84, 2),  // T
    //    (116, 2), // t
    //]);
    //let ret: Vec<u8> = input
    //    .bytes()
    //    .filter_map(|b| b.ok())
    //    .filter_map(|b| lookup.get(&b))
    //    .cloned()
    //    .collect();

    // 35.509148416s in debug, 1.198498916s in release
    //let lookup = HashMap::from([
    //    (65, 0),  // A
    //    (97, 0),  // a
    //    (67, 1),  // C
    //    (99, 1),  // c
    //    (71, 3),  // G
    //    (103, 3), // g
    //    (84, 2),  // T
    //    (116, 2), // t
    //]);

    //let mut buffer = vec![0; 65536];
    //loop {
    //    let bytes_read = input.read(&mut buffer)?;
    //    if bytes_read == 0 {
    //        break;
    //    }
    //    for byte in buffer[..bytes_read].iter() {
    //        if let Some(val) = lookup.get(byte) {
    //            ret.push(*val);
    //        }
    //    }
    //}

    // 10.935326709s in debug, 1.262571875s in release
    //for byte in input.bytes() {
    //    if let Some(val) = match byte.unwrap() {
    //        65 | 97 => Some(0),
    //        67 | 99 => Some(1),
    //        71 | 103 => Some(3),
    //        84 | 116 => Some(2),
    //        _ => None,
    //    } {
    //        ret.push(val);
    //    }
    //}

    // 10.899836542s in debug, 675.702666ms in release
    //for byte in input.bytes() {
    //    let base = BASE_LOOKUP[byte.unwrap() as usize];
    //    if base != 255 {
    //        ret.push(base);
    //    }
    //}

    // 18.955611916s in debug, 857.590666ms in release
    //let ret: Vec<u8> = input
    //    .bytes()
    //    .filter_map(|b| b.ok())
    //    .map(|b| BASE_LOOKUP[b as usize])
    //    .filter(|&b| b < 255)
    //    .collect();

    debug!("Read input in {:?}", start.elapsed());

    Ok(ret)
}

// --------------------------------------------------
#[cfg(test)]
mod tests {
    //use super::{caps_sa, read_input};
    //use super::_read_input;
    //use pretty_assertions::assert_eq;
    //use std::io::Cursor;

    //#[test]
    //fn test_read_input() {
    //    let cursor = Cursor::new("");
    //    let res = read_input(cursor);
    //    assert!(res.is_ok());
    //    assert_eq!(res.unwrap(), &[]);

    //    let cursor = Cursor::new("ACGT");
    //    let res = read_input(cursor);
    //    assert!(res.is_ok());
    //    assert_eq!(res.unwrap(), &[0, 1, 2, 3]);

    //    let cursor = Cursor::new("TGCA");
    //    let res = read_input(cursor);
    //    assert!(res.is_ok());
    //    assert_eq!(res.unwrap(), &[3, 2, 1, 0]);

    //    let cursor = Cursor::new("acgt");
    //    let res = read_input(cursor);
    //    assert!(res.is_ok());
    //    assert_eq!(res.unwrap(), &[0, 1, 2, 3]);

    //    let cursor = Cursor::new("tgca");
    //    let res = read_input(cursor);
    //    assert!(res.is_ok());
    //    assert_eq!(res.unwrap(), &[3, 2, 1, 0]);

    //    let cursor = Cursor::new("NNN");
    //    let res = read_input(cursor);
    //    assert!(res.is_ok());
    //    assert_eq!(res.unwrap(), &[]);
    //}

    //#[test]
    //fn test_caps_sa() {
    //    let (sa, lcp) = gen("AACTGCGGAT$"); // "CCATGGAG$"
    //    assert_eq!(sa, &[10, 0, 1, 8, 5, 2, 7, 4, 6, 9, 3,]);
    //    assert_eq!(lcp, &[0, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1,]);

    //    let (sa, lcp) = gen("CCATGGAG$");
    //    assert_eq!(sa, &[8, 6, 2, 1, 0, 7, 5, 4, 3]);
    //    assert_eq!(lcp, &[0, 0, 1, 0, 1, 0, 1, 1, 0]);
    //}
}
