use anyhow::Result;
use libsufr::{Context, SuffixArray};

fn main() -> Result<()> {
    // the text is ACGTNNACGT$
    let path = concat!(env!("CARGO_MANIFEST_DIR"), "/../data/inputs/1.sufr");
    let suffix_array = SuffixArray::open(path)?;

    let hits = suffix_array.search(b"ACG")?;
    println!("count: {}", hits.count());
    println!("locate: {:?}", hits.locate()?);

    let context = Context {
        before: 1,
        after: Some(3),
    };
    for extraction in hits.extract(context)? {
        let start = extraction.location.sequence_start;
        let window = extraction.sequence_range;
        let text = suffix_array.text(start + window.start..start + window.end)?;
        println!(
            "extract: {} with the hit at {}",
            String::from_utf8_lossy(&text),
            extraction.suffix_offset
        );
    }

    for hits in suffix_array.search_all(&["AC", "GT", "XX"]) {
        println!("count: {}", hits?.count());
    }

    Ok(())
}
