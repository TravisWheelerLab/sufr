use anyhow::Result;
use clap::{builder::PossibleValue, Parser, ValueEnum};
use log::{debug, info};
use std::{
    cmp::{max, min},
    //collections::HashMap,
    fs::{self, File},
    //fs::{self, File},
    //io::{BufWriter, Read},
    io::{BufRead, BufWriter, Write},
    mem,
    //mem::swap,
    path::PathBuf,
    slice,
    time::Instant,
};
use substring::Substring;

//const BASE_LOOKUP: [u8; 256] = [
//    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
//    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
//    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
//    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
//    255, 255, 255, 255, 255, 255, 255, 255, 255, 0, 255, 1, 255, 255, 255, 2,
//    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 3, 255, 255,
//    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 0, 255, 1, 255, 255,
//    255, 2, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 3,
//    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
//    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
//    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
//    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
//    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
//    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
//    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
//    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
//    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
//    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
//];

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
    pub text: Vec<char>,
    pub len: usize,
    //sa: Vec<usize>,
    //lcp: Vec<usize>,
    pub max_context: usize,
}

impl SuffixArray {
    pub fn new(input: &PathBuf, max_context: usize) -> Result<SuffixArray> {
        let mut text: Vec<char> = fs::read_to_string(input)?
            .trim()
            .to_uppercase()
            .chars()
            .collect();
        text.push('$');
        let len = text.len();

        Ok(SuffixArray {
            text,
            len,
            //sa: (0..len).collect(),
            //lcp: vec![0; len],
            max_context,
        })
    }

    pub fn merge_sort2(
        self: &Self,
        x: &mut [usize],
        y: &mut [usize],
        n: usize,
        lcp: &mut [usize],
        lcp_w: &mut [usize],
    ) {
        println!("x {x:?} y {y:?} n {n} lcp {lcp:?} lcp_w {lcp_w:?}");
        if n == 1 {
            lcp[0] = 0;
        } else {
            let mid = n / 2;
            self.merge_sort2(
                &mut y[..mid],
                &mut x[..mid],
                mid,
                &mut lcp_w[..mid],
                &mut lcp[..mid],
            );

            self.merge_sort2(
                &mut y[mid..],
                &mut x[mid..],
                n - mid,
                &mut lcp_w[mid..],
                &mut lcp[mid..],
            );

            self.merge2(x, mid, lcp_w, y, lcp)
        }
    }

    fn merge2(
        self: &Self,
        x_prime: &mut [usize],
        mid: usize,
        lcp_w: &mut [usize],
        z: &mut [usize],
        lcp_z: &mut [usize],
    ) {
        let mut x = x_prime[..mid].to_vec();
        let mut y = x_prime[mid..].to_vec();
        let mut len_x = x.len();
        let mut len_y = y.len();
        let mut lcp_x = lcp_w[..mid].to_vec();
        let mut lcp_y = lcp_w[mid..].to_vec();
        let mut m = 0; // LCP of last compared pair
                       // In the original code, l_x is intially set to null
        #[allow(unused_assignments)]
        let mut l_x = 0; // LCP(X_i, X_{i-1})
        let mut i = 0; // Index into x
        let mut j = 0; // Index into y
        let mut k = 0; // Index into z
        println!("MERGE x = {x:?} y = {y:?} len_x = {len_x} len_y = {len_y}");
        println!(
            "MERGE lcp_x = {lcp_x:?} lcp_y = {lcp_y:?} lcp_z = {lcp_z:?}"
        );

        while i < len_x && j < len_y {
            l_x = lcp_x[i];
            println!("i = {i} l_x = {l_x} m = {m}");

            if l_x > m {
                println!("BRANCH 1");
                z[k] = x[i];
                lcp_z[k] = l_x;
            } else if l_x < m {
                println!("BRANCH 2");
                z[k] = y[j];
                lcp_z[k] = m;
                m = l_x;
            } else {
                println!("BRANCH 3");
                // Length of shorter suffix
                let max_n = self.len - max(x[i], y[j]);

                // Prefix-context length for the suffixes
                //let context = min(self.max_context, max_n);

                // LCP(X_i, Y_j)
                println!("CALLING LCP x[i] = {} y[j] = {}", x[i], y[j]);
                let n = m + lcp(
                    &self.text[x[i]..],
                    &self.text[y[j]..],
                    max_n, //context - m,
                );
                println!("LCP = {n}");

                // If the len of the LCP is the entire shorter
                // sequence, take that.
                // Else, look at the next char after the LCP
                // to determine order.
                z[k] = if n == max_n {
                    max(x[i], y[j])
                } else if &self.text[x[i] + n] < &self.text[y[j] + n] {
                    x[i]
                } else {
                    y[j]
                };

                lcp_z[k] = if z[k] == x[i] { l_x } else { m };
                dbg!(&z);
                dbg!(&lcp_z);
                m = n;
            }

            if z[k] == x[i] {
                i += 1;
            } else {
                j += 1;
                mem::swap(&mut x, &mut y);
                mem::swap(&mut len_x, &mut len_y);
                mem::swap(&mut lcp_x, &mut lcp_y);
                mem::swap(&mut i, &mut j);
            }

            k += 1;
        }

        // Copy rest of the data from X to Z.
        while i < len_x {
            i += 1;
            k += 1;
            z[k] = x[i];
            lcp_z[k] = lcp_x[i];
        }

        // Copy rest of the data from Y to Z.
        dbg!(&len_y);
        dbg!(&j);
        if j < len_y {
            z[k] = y[j];
            lcp_z[k] = m;
            while j < len_y {
                dbg!(&k);
                dbg!(&z);
                dbg!(&y);
                z[k] = y[j];
                lcp_z[k] = lcp_y[j];
                j += 1;
                k += 1;
            }
        }

        // Copy back to the original values?
        x.append(&mut y);
        for i in 0..x.len() {
            x_prime[i] = x[i];
        }

        lcp_x.append(&mut lcp_y);
        for i in 0..lcp_x.len() {
            lcp_w[i] = lcp_x[i];
        }
    }

    pub fn merge_sort(self: &Self, arr: &mut [usize]) {
        if arr.len() > 1 {
            let mid = arr.len() / 2;

            // Sort the left half recursively.
            self.merge_sort(&mut arr[..mid]);

            // Sort the right half recursively.
            self.merge_sort(&mut arr[mid..]);

            // Combine the two halves.
            self.merge(arr, mid);
        }
    }

    fn merge(self: &Self, arr: &mut [usize], mid: usize) {
        // Create temporary vectors to support the merge.
        let left_half = arr[..mid].to_vec();
        let right_half = arr[mid..].to_vec();

        // Indexes to track the positions while merging.
        let mut l = 0;
        let mut r = 0;

        for v in arr {
            // Choose either the smaller element, or from
            // whichever vec is not exhausted.
            if r == right_half.len()
                || (l < left_half.len()
                    && self.text[left_half[l]..] < self.text[right_half[r]..])
            {
                *v = left_half[l];
                l += 1;
            } else {
                *v = right_half[r];
                r += 1;
            }
        }
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

    let suf_arr =
        SuffixArray::new(&args.input, args.max_context.unwrap_or(0))?;
    info!("Input length {}", suf_arr.len);
    debug!("Input {}", suf_arr.text.iter().cloned().collect::<String>());
    //info!("SA {:?}", sa.sa);
    //info!("LCP {:?}", sa.lcp);

    let start = Instant::now();
    let mut sa: Vec<usize> = (0..suf_arr.len).collect();
    let mut sa_w = sa.clone();
    let mut lcp = vec![0; suf_arr.len];
    let mut lcp_w = vec![0; suf_arr.len];
    suf_arr.merge_sort2(
        &mut sa_w,
        &mut sa,
        suf_arr.len,
        &mut lcp,
        &mut lcp_w,
    );
    //info!("SA_w = {:?}", sa_w);
    //sa.merge_sort(&mut sa_w);
    info!("merge_sort finished in {:?}", start.elapsed());
    debug!("SA = {:?}", sa);
    debug!("SA_w = {:?}", sa_w);
    debug!("LCP = {:?}", lcp);
    debug!("LCP_w = {:?}", lcp_w);

    let text: String = suf_arr.text.iter().cloned().collect();
    for &pos in &sa_w {
        println!("pos {pos:4} suffix {}", text.substring(pos, suf_arr.len));
    }

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
fn lcp(s1: &[char], s2: &[char], len: usize) -> usize {
    dbg!(&s1);
    dbg!(&s2);
    for i in 0..len {
        if s1[i] != s2[i] {
            return i;
        }
    }
    len
}

#[test]
fn test_lcp() -> Result<()> {
    assert_eq!(lcp(&['A'], &['C'], 1), 0);
    assert_eq!(lcp(&['A'], &['A'], 1), 1);
    assert_eq!(lcp(&['A', 'A'], &['A', 'A', 'C'], 2), 2);
    assert_eq!(lcp(&['A', 'C'], &['A', 'A', 'C'], 2), 1);
    Ok(())
}

// --------------------------------------------------
fn usize_to_bytes(value: usize) -> Vec<u8> {
    // Determine the size of usize in bytes
    let size = std::mem::size_of::<usize>();

    // Create a vector to hold the bytes
    let mut bytes = Vec::with_capacity(size);

    // Convert usize to bytes
    for i in 0..size {
        bytes.push((value >> (i * 8)) as u8);
    }

    bytes
}

// --------------------------------------------------
//fn merge(
//    x: &[usize],
//    len_x: usize,
//    y: &[usize],
//    len_y: usize,
//    lcp_x: &[usize],
//    lcp_y: &[usize],
//    x: &[usize],
//    lcp_z: &[usize],
//    max_context: usize,
//) {
//    let m = 0; // LCP of the last compared pair.
//    let l_x = 0; // LCP(X_i, X_{i - 1}).
//    let i = 0; // Index into `X`.
//    let j = 0; // Index into `Y`.
//    let k = 0; // Index into `Z`.

//    while i < len_x && j < len_y {
//        l_x = lcp_x[i];

//        if l_x > m {
//            z[k] = x[i];
//            lcp_z[k] = l_x;
//            m = m;
//        } else if l_x < m {
//            z[k] = y[j];
//            lcp_z[k] = m;
//            m = l_x;
//        } else
//        // Compute LCP of X_i and Y_j through linear scan.
//        {
//            let max_n = n_ - max(x[i], y[j]); // Length of the s  horter suffix.
//            let context = min(max_context, max_n); // Prefix-cont  ext length for the suffixes.
//                                                   // LCP(X_i, Y_j)
//            let n = m + lcp_opt_avx(
//                T_ + (X[i] + m),
//                T_ + (Y[j] + m),
//                context - m,
//            );

//            // Whether the shorter suffix is a prefix of the longer one.
//            z[k] = if (n == max_n) {
//                max(x[i], y[j])
//            } else {
//                if (T_[x[i] + n] < T_[y[j] + n]) {
//                    x[i]
//                } else {
//                    y[j]
//                }
//            };
//            lcp_z[k] = if z[k] == x[i] { l_x } else { m };
//            m = n;
//        }

//        if (z[k] == x[i]) {
//            i += 1;
//        } else {
//            // Swap X and Y (and associated data structures)
//            // when Y_j gets   pulled into Z.
//            j += 1;
//            swap(x, y);
//            swap(len_x, len_y);
//            swap(lcp_x, lcp_y);
//            swap(i, j);
//        }

//        k += 1;
//    }

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
//pub fn caps_sa(
//    _text: &str,
//    n: usize,
//    subproblem_count: usize,
//    _max_context: usize,
//) -> (SuffixArray, LCP) {
//    let sa: Vec<_> = (0..n).collect();
//    let sa_w = sa.clone();
//    let mut lcp: Vec<_> = (0..n).collect();
//    let lcp_w = lcp.clone();
//    let m = n / subproblem_count;
//    for i in 0..m {
//        let b = i * m;
//        let e = (i + 1) * m;
//        merge_sort(
//            &sa_w[b..e],
//            &sa[b..e],
//            m + (if i < (subproblem_count - 1) {
//                0
//            } else {
//                n % subproblem_count
//            }),
//            &mut lcp[b..e],
//            &lcp_w[b..e],
//        );
//    }

//    (sa, lcp)
//}

// --------------------------------------------------
//pub fn merge_sort(
//    x: &[usize],
//    y: &[usize],
//    n: usize,
//    lcp: &mut [usize],
//    lcp_w: &[usize],
//) {
//    if n == 1 {
//        lcp[0] = 0;
//    } else {
//        let m = n / 2;
//        merge_sort(y, x, m, &mut lcp_w, lcp);
//        merge_sort(&y[0..m], &x[0..m], n - m, &mut lcp_w[0..m], &lcp[0..m]);

//        //merge(x, m, x + m, n - m, lcp_w, lcp_w + m, y, lcp);
//    }
//}

// --------------------------------------------------
//pub fn naive_gen(text: &str) -> (SuffixArray, LCP) {
//    let mut suffixes: Vec<_> = (0..text.len())
//        .map(|i| (text[i..].to_string(), i))
//        .collect();
//    suffixes.sort();

//    let sa: Vec<_> = suffixes.iter().map(|(_, i)| *i).collect();
//    //dbg!(&sa);

//    let suffixes: Vec<_> = suffixes.iter().map(|(t, _)| t).collect();
//    dbg!(&suffixes);

//    let mut lcp: Vec<usize> = Vec::with_capacity(sa.len());
//    lcp.push(0);
//    //dbg!(&lcp);

//    for i in 1..suffixes.len() {
//        let cur = suffixes[i];
//        let prev = suffixes[i - 1];
//        let same = cur
//            .chars()
//            .zip(prev.chars())
//            .map_while(|(c1, c2)| if c1 == c2 { Some(true) } else { None })
//            .count();
//        lcp.push(same);
//    }

//    (sa, lcp)
//}

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
