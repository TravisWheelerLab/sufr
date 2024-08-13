use clap::Parser;
use sufr::{run, Args};

// --------------------------------------------------
fn main() {
    if let Err(e) = run(Args::parse()) {
        eprintln!("Error: {e}");
        std::process::exit(1);
    }
}

//fn main() {
//    for text in ["AACTGCGGAT$", "CCATGGAG$"] {
//        let enum_text: Vec<_> = text
//            .chars()
//            .enumerate()
//            .map(|(i, c)| format!("{i:2}: {c}"))
//            .collect();
//        println!("{}", enum_text.join("\n"));
//        let (sa, lcp) = sa::gen(text);
//        dbg!(&sa);
//        dbg!(&lcp);
//    }
//}
