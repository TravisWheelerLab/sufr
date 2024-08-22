use anyhow::Result;
use clap::Parser;
use sufr::{Cli, Command};

// --------------------------------------------------
fn main() {
    if let Err(e) = run(Cli::parse()) {
        eprintln!("Error: {e}");
        std::process::exit(1);
    }
}

// --------------------------------------------------
fn run(args: Cli) -> Result<()> {
    match &args.command {
        Some(Command::Create(args)) => {
            sufr::create(args)?;
            Ok(())
        }
        Some(Command::Read(args)) => {
            sufr::read(args)?;
            Ok(())
        }
        _ => unreachable!(),
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
