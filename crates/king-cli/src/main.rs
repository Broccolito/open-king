//! `king` — drop-in CLI matching KING 2.3.2's command line and console output.
//!
//! The startup sequence below is the reference binary's own order, established by
//! probing it; each step names the observation that fixes its position.
//!
//!  1. banner and the "parameters in effect" block;
//!  2. the WARNING block, if the command line had problems — parsing never aborts;
//!  3. FATAL ERROR if no `-b` fileset was named, and nothing else is printed;
//!  4. either the "please specify an analysis" notice (no analysis requested) or the
//!     "will run separately" line (more than one);
//!  5. `KING starts at <time>`;
//!  6. the analysis-parameter checks, in this order: `--maxP`, `--sexchr`,
//!     `--seglength`, `--minConc`, `--risk`/`--model`;
//!  7. FATAL ERROR if the `.bed` cannot be opened;
//!  8. the analyses themselves.
//!
//! Everything goes to stdout; the reference writes nothing to stderr, and neither do we.
//! Every exit here is status 1, which is what the reference returns for all four
//! captured command lines.

#![forbid(unsafe_code)]

use std::io::Write;

use king_cli::{cli, console, load};

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let parsed = cli::parse(&args);
    let opts = &parsed.options;

    let mut out = console::startup(opts);
    out.push_str(&console::warning_block(&parsed.warnings));
    emit(&out);

    if opts.bed.is_empty() {
        fatal(console::GENOTYPE_FILES_REQUIRED);
    }

    if !opts.any_analysis() {
        emit(&console::no_analysis_block());
    } else {
        let separate = opts.separate_analyses();
        if separate.len() > 1 {
            emit(&console::separate_analyses_block(&separate));
        }
    }

    emit(&console::king_starts_at(console::now_local()));

    // Analysis-parameter validation, in the reference's own order.
    if opts.was_given(cli::Opt::MaxP) {
        // Both tails go through the inverse normal, and the first bad one is fatal.
        let p = opts.double(cli::Opt::MaxP);
        for tail in [p / 2.0, 1.0 - p / 2.0] {
            if tail <= 0.0 {
                fatal(&console::p_value_out_of_range(tail));
            }
        }
    }

    let sexchr = opts.int(cli::Opt::Sexchr);
    if sexchr < console::MIN_SEX_CHROMOSOME {
        fatal(&console::sex_chromosome_out_of_range(sexchr));
    }
    if sexchr != console::HUMAN_SEX_CHROMOSOME {
        emit(&console::non_human_notice(sexchr));
    }

    if opts.was_given(cli::Opt::Seglength) {
        let seglength = opts.double(cli::Opt::Seglength);
        if seglength > console::SEGLENGTH_MIN && seglength < console::SEGLENGTH_MAX {
            emit(&console::segment_length_notice(seglength));
        } else {
            emit(console::SEGLENGTH_OUT_OF_RANGE);
        }
    }

    let min_conc = opts.double(cli::Opt::MinConc);
    if !(0.0..=1.0).contains(&min_conc) {
        emit(console::MINCONC_OUT_OF_RANGE);
    }

    if opts.flag(cli::Opt::Risk) && opts.string(cli::Opt::Model).is_empty() {
        fatal(console::RISK_MODEL_REQUIRED);
    }

    // The fileset load owns its own console sequence and its own `.bed` probe — including
    // the reference's quirk of skipping that probe on the `--risk --model` path.
    let loaded = match load::load(opts, &mut std::io::stdout()) {
        Ok(loaded) => loaded,
        Err(message) => fatal(&message),
    };

    // ------------------------------------------------------------------
    // SEAM: the analysis engines plug in here.
    //
    // TODO(engine): dispatch the analyses named by `opts.analyses_in_effect()` into
    // `king_core`. The reference reprints `loaded.preamble()` — the
    // "Autosome genotypes stored in …" line and the blank line under it — ahead of
    // *each* pass, then that pass's `console::options_in_effect`, then its body; a
    // `--kinship --ibs` run shows the pair twice. `console::king_ends_at` closes the run.
    //
    // The single-analysis preamble below is what the loader is obliged to print and no
    // more; the per-analysis loop replaces it.
    // ------------------------------------------------------------------

    emit(&loaded.preamble());

    std::process::exit(0);
}

/// Write to stdout and flush, so ordering survives the `exit` calls below.
fn emit(s: &str) {
    let mut stdout = std::io::stdout();
    let _ = stdout.write_all(s.as_bytes());
    let _ = stdout.flush();
}

/// Print a FATAL ERROR block and leave with the reference's exit status of 1.
fn fatal(message: &str) -> ! {
    emit(&console::fatal_block(message));
    std::process::exit(1);
}
