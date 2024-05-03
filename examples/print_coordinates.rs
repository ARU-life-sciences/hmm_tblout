extern crate hmm_tblout;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // get the command line args, only parse the
    // first one which should be a fasta file
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 2 {
        println!("Usage: print_coordinates <tblout_file>");
        std::process::exit(1);
    }

    let reader = hmm_tblout::Reader::from_path(args[1].clone())?;

    for record in reader.into_records() {
        let r = record?;
        let tname = r.target_name();
        let strand = r.strand();
        let alifrom = r.ali_from();
        let alito = r.ali_to();

        println!("{}\t{}\t{}\t{}", tname, strand, alifrom, alito);
    }

    Ok(())
}
