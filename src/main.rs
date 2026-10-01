use chloroplast_ir::{
    arrange_from_contigs, arrange_single_contig, classify, fasta, find_inverted_repeats,
    ArrangedGenome, IrFinderConfig,
};
use std::env;
use std::fs;
use std::process::ExitCode;

/*
Gaurav Sablok
gsablok@proton.me
*/

fn usage() -> String {
    r#"chloroplast-ir — arrange an assembled chloroplast genome FASTA into the
canonical end-to-end LSC -> IRb -> SSC -> IRa order, filling N only where an
assembly gap genuinely leaves the junction unknown.

Gaurav Sablok
gsablok@proton.me

USAGE:
    chloroplast-ir --input <FASTA> --output <FASTA> [OPTIONS]

INPUT SHAPES:
    * One FASTA record  -> treated as a single already-circularized contig.
      The genome is rotated (never base-altered) to the canonical start.
    * 2+ FASTA records   -> treated as separate contigs to classify and
      stitch together (see --lsc/--ir/--ir2/--ssc to override auto-detection).

OPTIONS:
    --input <FILE>        Input FASTA (required)
    --output <FILE>        Output FASTA (required)
    --id <STRING>          Sequence id for the output record (default: derived)
    --min-len <N>           Minimum IR length to accept, bp (default: 1000)
    --kmer <N>              Seed k-mer size, <=31 (default: 21)
    --gap-n <N>             N's inserted at each inter-contig junction in
                             multi-contig mode (default: 0)
    --wrap <N>              FASTA line width, 0 = no wrap (default: 70)
    --lsc <SEQ_ID>          Force which input record is the LSC (multi-contig)
    --ir <SEQ_ID>           Force which input record is IR copy 1 (multi-contig)
    --ir2 <SEQ_ID>          Force which input record is IR copy 2 (multi-contig)
    --ssc <SEQ_ID>          Force which input record is the SSC (multi-contig)
    -h, --help              Show this help

EXIT CODES:
    0 success, 1 usage/argument error, 2 no inverted repeat could be found.
"#
    .to_string()
}

struct Args {
    input: String,
    output: String,
    id: Option<String>,
    min_len: usize,
    kmer: usize,
    gap_n: usize,
    wrap: usize,
    force_lsc: Option<String>,
    force_ir: Option<String>,
    force_ir2: Option<String>,
    force_ssc: Option<String>,
}

fn parse_args() -> Result<Args, String> {
    let mut a = Args {
        input: String::new(),
        output: String::new(),
        id: None,
        min_len: 1000,
        kmer: 21,
        gap_n: 0,
        wrap: 70,
        force_lsc: None,
        force_ir: None,
        force_ir2: None,
        force_ssc: None,
    };
    let mut it = env::args().skip(1);
    let mut have_input = false;
    let mut have_output = false;
    while let Some(flag) = it.next() {
        macro_rules! val {
            () => {
                it.next()
                    .ok_or_else(|| format!("{flag} requires a value"))?
            };
        }
        match flag.as_str() {
            "-h" | "--help" => {
                print!("{}", usage());
                std::process::exit(0);
            }
            "--input" => {
                a.input = val!();
                have_input = true;
            }
            "--output" => {
                a.output = val!();
                have_output = true;
            }
            "--id" => a.id = Some(val!()),
            "--min-len" => {
                a.min_len = val!().parse().map_err(|_| "--min-len must be an integer")?
            }
            "--kmer" => a.kmer = val!().parse().map_err(|_| "--kmer must be an integer")?,
            "--gap-n" => a.gap_n = val!().parse().map_err(|_| "--gap-n must be an integer")?,
            "--wrap" => a.wrap = val!().parse().map_err(|_| "--wrap must be an integer")?,
            "--lsc" => a.force_lsc = Some(val!()),
            "--ir" => a.force_ir = Some(val!()),
            "--ir2" => a.force_ir2 = Some(val!()),
            "--ssc" => a.force_ssc = Some(val!()),
            other => return Err(format!("unknown argument: {other}")),
        }
    }
    if !have_input || !have_output {
        return Err("--input and --output are required".to_string());
    }
    if a.kmer == 0 || a.kmer > 31 {
        return Err("--kmer must be between 1 and 31".to_string());
    }
    Ok(a)
}

fn run() -> Result<(), (String, u8)> {
    let args = match parse_args() {
        Ok(a) => a,
        Err(e) => {
            eprintln!("error: {e}\n\n{}", usage());
            return Err((String::new(), 1));
        }
    };

    let records = fasta::read_file(&args.input)
        .map_err(|e| (format!("failed to read {}: {e}", args.input), 1))?;

    let cfg = IrFinderConfig {
        kmer: args.kmer,
        min_len: args.min_len,
        ..Default::default()
    };

    let (arranged, out_id, report): (ArrangedGenome, String, String) = if records.len() == 1 {
        let rec = &records[0];
        let hits = find_inverted_repeats(&rec.seq, &cfg);
        let best = hits.first().ok_or_else(|| {
            (
                format!(
                    "no inverted repeat >= {} bp found in '{}'. Try lowering --min-len, or if this \
                     genome is split across multiple contigs, pass all of them in one FASTA.",
                    args.min_len, rec.id
                ),
                2,
            )
        })?;
        let arranged = arrange_single_contig(&rec.seq, best);
        let r = &arranged.regions;
        let report = format!(
            "input: single contig '{}' ({} bp)\nIR match: {} bp, {:.2}% identity ({} matches / {} mismatches)\nLSC {}..{} ({} bp)\nIRb {}..{} ({} bp)\nSSC {}..{} ({} bp)\nIRa {}..{} ({} bp)\n",
            rec.id, rec.seq.len(),
            best.region1_len(), best.identity() * 100.0, best.matches, best.mismatches,
            r.lsc.0, r.lsc.1, r.lsc.1 - r.lsc.0,
            r.ir_b.0, r.ir_b.1, r.ir_b.1 - r.ir_b.0,
            r.ssc.0, r.ssc.1, r.ssc.1 - r.ssc.0,
            r.ir_a.0, r.ir_a.1, r.ir_a.1 - r.ir_a.0,
        );
        let id = args
            .id
            .clone()
            .unwrap_or_else(|| format!("{}_arranged", rec.id));
        (arranged, id, report)
    } else {
        let by_id = |wanted: &str| {
            records
                .iter()
                .find(|r| r.id == wanted)
                .map(|r| r.seq.as_slice())
                .ok_or_else(|| (format!("no input record with id '{wanted}'"), 1))
        };

        let (lsc, ir, ssc, ir2, used_heuristic, lsc_id, ir_id, ssc_id, ir2_id): (
            &[u8],
            &[u8],
            &[u8],
            Option<Vec<u8>>,
            bool,
            String,
            String,
            String,
            Option<String>,
        ) = if args.force_lsc.is_some() || args.force_ir.is_some() || args.force_ssc.is_some() {
            let lsc_id = args.force_lsc.clone().ok_or((
                "--ir/--ssc were given but --lsc was not; all of --lsc/--ir/--ssc are required together"
                    .to_string(),
                1,
            ))?;
            let ir_id = args.force_ir.clone().ok_or((
                "--lsc/--ssc were given but --ir was not; all of --lsc/--ir/--ssc are required together"
                    .to_string(),
                1,
            ))?;
            let ssc_id = args.force_ssc.clone().ok_or((
                "--lsc/--ir were given but --ssc was not; all of --lsc/--ir/--ssc are required together"
                    .to_string(),
                1,
            ))?;
            let lsc = by_id(&lsc_id)?;
            let ir = by_id(&ir_id)?;
            let ssc = by_id(&ssc_id)?;
            let ir2 = match &args.force_ir2 {
                Some(id) => Some(by_id(id)?.to_vec()),
                None => None,
            };
            (
                lsc,
                ir,
                ssc,
                ir2,
                false,
                lsc_id,
                ir_id,
                ssc_id,
                args.force_ir2.clone(),
            )
        } else {
            let named: Vec<(&str, &[u8])> = records
                .iter()
                .map(|r| (r.id.as_str(), r.seq.as_slice()))
                .collect();
            let c = classify(&named, &cfg).map_err(|e| (e, 2))?;
            (
                c.lsc,
                c.ir,
                c.ssc,
                c.ir_second.clone(),
                c.used_length_heuristic,
                c.lsc_id.to_string(),
                c.ir_id.to_string(),
                c.ssc_id.to_string(),
                c.ir_second_id.map(|s| s.to_string()),
            )
        };

        let arranged = arrange_from_contigs(lsc, ir, ssc, ir2.as_deref(), args.gap_n);
        let r = &arranged.regions;
        let mut report = format!(
            "input: {} contigs\nLSC = '{}' ({} bp)\nIR  = '{}' ({} bp){}\nSSC = '{}' ({} bp)\n",
            records.len(),
            lsc_id,
            lsc.len(),
            ir_id,
            ir.len(),
            match &ir2_id {
                Some(id) => format!(
                    "  [second IR copy: '{id}', {} bp]",
                    ir2.as_ref().unwrap().len()
                ),
                None => "  [second IR copy not found in input; synthesized as reverse complement]"
                    .to_string(),
            },
            ssc_id,
            ssc.len(),
        );
        if used_heuristic {
            report.push_str(
                "NOTE: role assignment used a length-based heuristic (longest=LSC, \
                 middle=IR, shortest=SSC) because no matching inverted-repeat contig pair was \
                 found. Verify this is correct for your data, or override with --lsc/--ir/--ssc.\n",
            );
        }
        if args.gap_n > 0 {
            report.push_str(&format!(
                "Inserted {} N's at each of the inter-contig junctions (unknown exact junction sequence).\n",
                args.gap_n
            ));
        }
        report.push_str(&format!(
            "Output layout: LSC {}..{}  IRb {}..{}  SSC {}..{}  IRa {}..{}\n",
            r.lsc.0, r.lsc.1, r.ir_b.0, r.ir_b.1, r.ssc.0, r.ssc.1, r.ir_a.0, r.ir_a.1,
        ));
        let id = args
            .id
            .clone()
            .unwrap_or_else(|| "chloroplast_arranged".to_string());
        (arranged, id, report)
    };

    let n_count = arranged
        .seq
        .iter()
        .filter(|&&b| b.to_ascii_uppercase() == b'N')
        .count();
    let total = arranged.seq.len();
    let full_report = format!(
        "{report}Total length: {total} bp ({n_count} bp N, {:.4}%)\n",
        100.0 * n_count as f64 / total.max(1) as f64
    );

    let desc = format!(
        "arranged_chloroplast_genome length={total} lsc={}-{} irb={}-{} ssc={}-{} ira={}-{}",
        arranged.regions.lsc.0 + 1,
        arranged.regions.lsc.1,
        arranged.regions.ir_b.0 + 1,
        arranged.regions.ir_b.1,
        arranged.regions.ssc.0 + 1,
        arranged.regions.ssc.1,
        arranged.regions.ir_a.0 + 1,
        arranged.regions.ir_a.1,
    );
    let out_record = fasta::FastaRecord {
        id: out_id,
        description: desc,
        seq: arranged.seq,
    };
    let out_text = fasta::write_string(&[out_record], args.wrap);
    fs::write(&args.output, out_text)
        .map_err(|e| (format!("failed to write {}: {e}", args.output), 1))?;

    eprint!("{full_report}");
    eprintln!("Wrote {}", args.output);
    Ok(())
}

fn main() -> ExitCode {
    match run() {
        Ok(()) => ExitCode::SUCCESS,
        Err((msg, code)) => {
            if !msg.is_empty() {
                eprintln!("error: {msg}");
            }
            ExitCode::from(code)
        }
    }
}
