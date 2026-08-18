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

use king_cli::{analysis, cli, console, load};

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
    // Analysis dispatch. Each pass emits, in order: `loaded.preamble()` when
    // `load::prints_preamble` says that analysis opens with it, then that pass's
    // `console::options_in_effect`, then its body — so a `--kinship --ibs` run will
    // show the preamble twice and a `--build` run not at all.
    //
    // SEAM: `--kinship`, `--ibs` and `--duplicate` are implemented. Every other analysis
    // still falls through to the bare preamble, which is as far as the loader's own
    // obligation reaches.
    // ------------------------------------------------------------------

    let mut ran = false;
    // `--related` leads the dispatch because it leads `SEPARATE_ANALYSES`: a
    // `--related --ibdseg` run prints the `--related` pass first.
    if opts.flag(cli::Opt::Related) {
        if analysis::related::downgrades_to_kinship(loaded.fileset.samples.len()) {
            let plain = opts.without_degree();
            emit(&analysis::related::small_sample_notice());
            emit(&loaded.preamble());
            emit(&console::options_in_effect(&["--kinship".to_string()]));
            analysis::kinship::run(&plain, &loaded, &mut std::io::stdout());
            ran = true;
        } else {
            emit(&loaded.preamble());
            emit(&console::options_in_effect(
                &analysis::related::options_in_effect(opts),
            ));
            analysis::related::run(opts, &loaded, &mut std::io::stdout());
            ran = true;
        }
    }
    if opts.flag(cli::Opt::Kinship) {
        emit(&loaded.preamble());
        emit(&console::options_in_effect(
            &analysis::kinship::options_in_effect(opts),
        ));
        analysis::kinship::run(opts, &loaded, &mut std::io::stdout());
        ran = true;
    }
    // `--autoQC` follows `--kinship` in the reference's own pass order, and opens with a
    // line of the loader's preamble's shape but not its arithmetic: its word count is
    // `ceil(m / 16)`, the denser packing the QC pipeline uses.
    if opts.flag(cli::Opt::AutoQc) {
        emit(&analysis::autoqc::preamble(&loaded));
        emit(&console::options_in_effect(&analysis::options_in_effect(
            opts,
            cli::Opt::AutoQc,
        )));
        analysis::autoqc::run(opts, &loaded, &mut std::io::stdout());
        ran = true;
    }
    for (opt, body) in [
        (
            cli::Opt::Duplicate,
            analysis::duplicate::run as fn(&cli::Options, &load::Loaded, &mut dyn Write),
        ),
        (cli::Opt::Ibs, analysis::ibs::run),
    ] {
        if !opts.flag(opt) {
            continue;
        }
        if load::prints_preamble(opt) {
            emit(&loaded.preamble());
        }
        emit(&console::options_in_effect(&analysis::options_in_effect(
            opts, opt,
        )));
        body(opts, &loaded, &mut std::io::stdout());
        ran = true;
    }
    // The five passes that open with no preamble: `--unrelated`, `--build` and
    // `--cluster` run their own family clustering, and the two QC reports write
    // `allsegs.txt` silently before their body.
    for (opt, body) in [
        (
            cli::Opt::Unrelated,
            analysis::unrelated::run as fn(&cli::Options, &load::Loaded, &mut dyn Write),
        ),
        (cli::Opt::Build, analysis::build::run),
        (cli::Opt::Bysample, analysis::qc::run_bysample),
        (cli::Opt::BySnp, analysis::qc::run_bysnp),
        (cli::Opt::Cluster, analysis::cluster::run),
    ] {
        if !opts.flag(opt) {
            continue;
        }
        // The two QC reports write `allsegs.txt` silently before their body; the three
        // clustering passes emit it themselves, inside their own console block.
        if matches!(opt, cli::Opt::Bysample | cli::Opt::BySnp) {
            let _ = analysis::ibdseg::segment_prepass(opts, &loaded);
        }
        // None of the three opens with the preamble, so the blank line that always
        // precedes `Options in effect:` has to come from here rather than from it.
        emit("\n");
        emit(&console::options_in_effect(&analysis::options_in_effect(
            opts, opt,
        )));
        body(opts, &loaded, &mut std::io::stdout());
        ran = true;
    }
    if opts.flag(cli::Opt::Ibdseg) {
        // Under five samples the reference quietly becomes a `--kinship` run: it says so,
        // prints the preamble `--ibdseg` never prints, and takes the kinship path whole.
        if analysis::ibdseg::downgrades_to_kinship(loaded.fileset.samples.len()) {
            let plain = opts.without_degree();
            emit(&analysis::ibdseg::small_sample_notice());
            emit(&loaded.preamble());
            emit(&console::options_in_effect(&["--kinship".to_string()]));
            analysis::kinship::run(&plain, &loaded, &mut std::io::stdout());
        } else {
            // The splitped line is the one console line that precedes `Options in
            // effect:` rather than following it.
            if analysis::splitped::is_generated(&loaded.fileset.samples) {
                emit(&analysis::ibdseg::splitped_notice(
                    opts.string(cli::Opt::Prefix),
                ));
            } else {
                emit("\n");
            }
            emit(&console::options_in_effect(
                &analysis::ibdseg::options_in_effect(opts),
            ));
            analysis::ibdseg::run(opts, &loaded, &mut std::io::stdout());
        }
        ran = true;
    }
    if !ran
        && cli::all()
            .any(|o| o.kind() == cli::Kind::Flag && opts.flag(o) && load::prints_preamble(o))
    {
        emit(&loaded.preamble());
    }

    // The reference closes a completed run with the timestamp and a blank line; a fatal
    // exit never reaches it.
    if ran {
        emit(&console::king_ends_at(console::now_local()));
        emit("\n");
    }

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
