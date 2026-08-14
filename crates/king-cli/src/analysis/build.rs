//! `--build`: family clustering followed by pedigree reconstruction.
//!
//! The pass opens with the same clustering prologue as `--cluster` and `--unrelated`
//! ([`unrelated::clustering_prologue`]) and then reconstructs pedigrees inside the
//! families that clustering built:
//!
//! ```text
//! Pedigree reconstruction starts at <t>
//! Reconstructing pedigree...
//! Age information not provided.
//! Total length of <k> chromosomal segments usable for IBD segment analysis is <d> Mb.
//!   Information of these chromosomal segments can be found in file <p>allsegs.txt
//!
//! <the reconstruction log, echoed to the console as well as to the file>
//!
//! Details of pedigree reconstruction are available in log file <p>build.log
//! Update-ID information is saved in file <p>updateids.txt
//! No pedigrees can be reconstructed.
//! Pedigree reconstruction ends at <t>
//! ```
//!
//! Note the segment total is printed **twice** — once by the clustering prologue and
//! again here — and that the two blank lines above `Details of …` are the segment
//! block's own blank plus the one that closes an empty log.
//!
//! # What gets reconstructed
//!
//! Only the **newly clustered** families do. Every capture in which clustering joined
//! nothing — nine of the thirteen, `threegen`'s single twelve-member three-generation
//! family included — produces a **zero-byte** `<p>build.log`, a **zero-byte**
//! `<p>updateparents.txt`, no `<p>updateids.txt` at all despite the console line saying
//! otherwise, and the closing `No pedigrees can be reconstructed.`. `bigish`, the one
//! fileset large enough for the hundred-sample clustering gate to fire, is also the one
//! that logs rules and inferences and writes all three files.
//!
//! **The reconstruction itself is unimplemented**; this module reproduces the
//! no-reconstruction outcome and stops. See `docs/PARITY.md` §11.

use std::io::Write;
use std::path::Path;

use crate::analysis::{ibdseg, out_path, unrelated};
use crate::cli::Options;
use crate::console;
use crate::load::Loaded;

/// Run the pass. The caller has already printed `Options in effect:`.
pub fn run(opts: &Options, loaded: &Loaded, out: &mut dyn Write) {
    let clustering = unrelated::clustering_prologue(opts, loaded, out);
    if clustering.tiny {
        // Under ten samples the reference stops at the disabled notice: no
        // reconstruction block, and not one file.
        return;
    }

    let _ = out.write_all(
        format!(
            "Pedigree reconstruction starts at {}\n",
            console::ctime(console::now_local())
        )
        .as_bytes(),
    );
    let _ = out.write_all(b"Reconstructing pedigree...\n");
    // No `.fam` carries ages, and the corpus never provides one, so this notice is
    // unconditional here.
    let _ = out.write_all(b"Age information not provided.\n");
    let _ = out.write_all(ibdseg::segment_prepass(opts, loaded).as_bytes());

    // The log, echoed to stdout and written to the file. Empty unless clustering merged
    // families, which is the branch this module does not implement.
    let log = String::new();
    let _ = out.write_all(log.as_bytes());
    let _ = out.write_all(b"\n");

    let log_path = out_path(opts, "build.log");
    let parents_path = out_path(opts, "updateparents.txt");
    let _ = std::fs::write(Path::new(&log_path), log.as_bytes());
    let _ = std::fs::write(Path::new(&parents_path), b"");

    let _ = out.write_all(
        format!("Details of pedigree reconstruction are available in log file {log_path}\n")
            .as_bytes(),
    );
    // Announced whether or not the file is written: on an unmerged run the reference
    // prints this line and leaves no `updateids.txt` behind.
    let _ = out.write_all(
        format!(
            "Update-ID information is saved in file {}\n",
            out_path(opts, "updateids.txt")
        )
        .as_bytes(),
    );
    if clustering.any_merged {
        let _ = out.write_all(
            format!("Update-parent information is saved in file {parents_path}\n").as_bytes(),
        );
    } else {
        let _ = out.write_all(b"No pedigrees can be reconstructed.\n");
    }
    let _ = out.write_all(
        format!(
            "Pedigree reconstruction ends at {}\n",
            console::ctime(console::now_local())
        )
        .as_bytes(),
    );
}
