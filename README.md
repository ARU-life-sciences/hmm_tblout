# `hmm_tblout`

Simple parsing of tabular output from `HMMER::nhmmer --tblout ...`.

## Example

Run this example using `cargo run --example print_coordinates ./data/test.tbl`.

```rust
use hmm_tblout;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // get the command line args, only parse the
    // first one which should be a fasta file
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 2 {
        println!("Usage: print_coordinates <tblout_file>");
        std::process::exit(1);
    }

    // the Reader::from_path is probably the most convenient API.
    let reader = hmm_tblout::Reader::from_path(args[1].clone())?;

    // I'm iterating over the owned records here, but you can borrow
    // too using `.records()`.
    for record in reader.into_records() {
        let r = record?;
        // there is a function to retrieve each field.
        let tname = r.target_name();
        let strand = r.strand();
        let alifrom = r.ali_from();
        let alito = r.ali_to();

        println!("{}\t{}\t{}\t{}", tname, strand, alifrom, alito);
    }

    Ok(())
}
```


## Yet to implement

May handle these in the future. Or feel free to contribute!

- Does not handle the description column, as this may contain spaces. 
- Metadata ignored.
