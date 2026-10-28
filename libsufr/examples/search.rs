use anyhow::Result;
use libsufr::{suffix_array::SuffixArray, types::Context};

fn main() -> Result<()> {
    // text "ACGTNNACGT$"
    let suffix_array =
        SuffixArray::open(concat!(env!("CARGO_MANIFEST_DIR"), "/../data/inputs/1.sufr"))?;

    let hits = suffix_array.search(b"ACG")?;
    println!("count: {}", hits.count());
    println!("locate: {:?}", hits.locate()?);
    println!("extract: {:?}", hits.extract(Context { before: 1, after: Some(3) })?);

    for hits in suffix_array.search_many(&["AC", "GT", "XX"]) {
        println!("count: {}", hits?.count());
    }

    Ok(())
}
