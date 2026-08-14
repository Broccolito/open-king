//! `.bed` — the PLINK 1 binary genotype matrix, and the fileset entry point.
//!
//! # File layout
//!
//! Three header bytes, then the genotype block:
//!
//! ```text
//! 0x6c 0x1b   magic
//! 0x01        mode: 1 = SNP-major, 0 = individual-major
//! ...         ceil(n_samples / 4) bytes per variant, one row per variant
//! ```
//!
//! Within a row, four samples are packed per byte and the **low-order** bit pair holds the
//! **first** sample of the four. The two-bit codes and the bit planes they become are
//! given by [`planes_for_code`]; the table there was verified byte for byte against a
//! fileset built by PLINK 1.9 from a `.ped` with known genotypes (see the
//! `plink_written_genotypes_round_trip` test, which carries that hexdump inline).
//!
//! # Transposition
//!
//! The reader transposes as it goes: the file is variant-major, the bit planes are
//! sample-major, with bit `m` of sample `s`'s words holding the `m`-th **retained**
//! variant. Filtered-out rows are read and discarded, so `m` counts retained variants
//! only and the planes never carry a variant the caller asked to drop. Words are
//! zero-initialised and only ever OR-ed into, which is what keeps the tail of the final
//! word zero as [`Genotypes`] requires.

use std::cmp::min;
use std::fs::File;
use std::io::{BufReader, Read};
use std::path::{Path, PathBuf};

use crate::{bim, fam, Fileset, Genotypes, IoError, Result, Variant, VariantFilter};

/// The two magic bytes every PLINK 1 `.bed` starts with.
pub const MAGIC: [u8; 2] = [0x6c, 0x1b];
/// Mode byte for individual-major order, which KING does not accept.
pub const MODE_INDIVIDUAL_MAJOR: u8 = 0x00;
/// Mode byte for SNP-major order — the only layout PLINK has written since 1.00.
pub const MODE_SNP_MAJOR: u8 = 0x01;
/// Length of the magic + mode header.
pub const HEADER_LEN: usize = 3;

/// Bytes one variant's row occupies: four samples per byte, rounded up.
pub fn bytes_per_variant(n_samples: usize) -> usize {
    n_samples.div_ceil(4)
}

/// Total size a well-formed `.bed` must have for this many samples and variants.
pub fn expected_bed_len(n_samples: usize, n_variants: usize) -> u64 {
    HEADER_LEN as u64 + n_variants as u64 * bytes_per_variant(n_samples) as u64
}

/// PLINK's two-bit genotype code to the `(plane0, plane1)` bit pair of [`Genotypes`].
///
/// | code | genotype           | `plane0` | `plane1` |
/// |------|--------------------|----------|----------|
/// | `00` | homozygous A1/A1   | 1        | 1        |
/// | `01` | missing            | 0        | 0        |
/// | `10` | heterozygous       | 0        | 1        |
/// | `11` | homozygous A2/A2   | 1        | 0        |
///
/// Note that PLINK's `01` is *missing*, not a genotype — the codes are not a count of A2
/// alleles, and reading them as one silently turns every no-call into a homozygote.
pub const fn planes_for_code(code: u8) -> (u64, u64) {
    match code & 0b11 {
        0b00 => (1, 1), // hom A1/A1
        0b10 => (0, 1), // het
        0b11 => (1, 0), // hom A2/A2
        // 0b01 is missing; the `&` above makes any other value impossible.
        _ => (0, 0),
    }
}

/// Whether a variant survives `filter`.
pub fn is_kept(variant: &Variant, filter: VariantFilter) -> bool {
    match filter {
        VariantFilter::All => true,
        VariantFilter::Autosomes => variant.is_autosome(),
        VariantFilter::Chromosome(c) => variant.chrom_code() == Some(c),
    }
}

/// The `.bed`, `.bim` and `.fam` paths implied by a prefix or by an explicit `.bed` path.
///
/// A trailing `.bed` (any case) is stripped to get the prefix; otherwise the argument *is*
/// the prefix. The suffix is appended rather than substituted, so a prefix that already
/// contains a dot — `study.v2` — resolves to `study.v2.bed`, not `study.bed`.
pub fn paths_for(prefix_or_bed: &Path) -> (PathBuf, PathBuf, PathBuf) {
    let is_bed = prefix_or_bed
        .extension()
        .is_some_and(|e| e.eq_ignore_ascii_case("bed"));
    if is_bed {
        let prefix = prefix_or_bed.with_extension("");
        (
            prefix_or_bed.to_path_buf(),
            with_suffix(&prefix, "bim"),
            with_suffix(&prefix, "fam"),
        )
    } else {
        (
            with_suffix(prefix_or_bed, "bed"),
            with_suffix(prefix_or_bed, "bim"),
            with_suffix(prefix_or_bed, "fam"),
        )
    }
}

fn with_suffix(prefix: &Path, ext: &str) -> PathBuf {
    let mut s = prefix.as_os_str().to_os_string();
    s.push(".");
    s.push(ext);
    PathBuf::from(s)
}

/// Load a `.bed`/`.bim`/`.fam` triple.
///
/// `prefix_or_bed` may be either the shared prefix or the `.bed` path itself.
/// `fam_override` / `bim_override` correspond to the reference binary's `--fam` / `--bim`
/// options and replace the derived paths when given.
///
/// Only variants passing `filter` are transposed into [`Fileset::genotypes`];
/// [`Fileset::variants`] still holds every row of the `.bim`, and [`Fileset::kept`] says
/// which of them the bit planes carry, in order.
pub fn read_fileset(
    prefix_or_bed: &Path,
    filter: VariantFilter,
    fam_override: Option<&Path>,
    bim_override: Option<&Path>,
) -> Result<Fileset> {
    let (bed_path, default_bim, default_fam) = paths_for(prefix_or_bed);
    let fam_path = fam_override.map_or(default_fam, Path::to_path_buf);
    let bim_path = bim_override.map_or(default_bim, Path::to_path_buf);

    let samples = fam::read_fam(&fam_path)?;
    fam::check_duplicates(&samples)?;
    let variants = bim::read_bim(&bim_path)?;
    let (genotypes, kept) = read_bed(&bed_path, samples.len(), &variants, filter)?;

    Ok(Fileset {
        samples,
        variants,
        genotypes,
        kept,
    })
}

/// Read a `.bed` whose sample and variant tables have already been loaded.
///
/// Returns the bit planes and the indices into `variants` they carry.
pub fn read_bed(
    path: &Path,
    n_samples: usize,
    variants: &[Variant],
    filter: VariantFilter,
) -> Result<(Genotypes, Vec<usize>)> {
    let file = File::open(path).map_err(|source| IoError::Open {
        path: path.to_path_buf(),
        source,
    })?;
    let found = file
        .metadata()
        .map_err(|source| IoError::Read {
            path: path.to_path_buf(),
            source,
        })?
        .len();

    let mut reader = BufReader::new(file);
    let mut header = [0u8; HEADER_LEN];
    if found < HEADER_LEN as u64 {
        return Err(IoError::SizeMismatch {
            path: path.to_path_buf(),
            expected: expected_bed_len(n_samples, variants.len()),
            found,
        });
    }
    reader
        .read_exact(&mut header)
        .map_err(|source| IoError::Read {
            path: path.to_path_buf(),
            source,
        })?;
    check_header_and_len(&header, found, path, n_samples, variants.len())?;

    decode_snp_major(&mut reader, path, n_samples, variants, filter)
}

/// Decode an in-memory `.bed` image, header included.
///
/// Same contract as [`read_bed`]; `path` is used only to label errors.
pub fn read_bed_bytes(
    bytes: &[u8],
    path: &Path,
    n_samples: usize,
    variants: &[Variant],
    filter: VariantFilter,
) -> Result<(Genotypes, Vec<usize>)> {
    let found = bytes.len() as u64;
    if bytes.len() < HEADER_LEN {
        return Err(IoError::SizeMismatch {
            path: path.to_path_buf(),
            expected: expected_bed_len(n_samples, variants.len()),
            found,
        });
    }
    let (header, mut body) = bytes.split_at(HEADER_LEN);
    check_header_and_len(header, found, path, n_samples, variants.len())?;
    decode_snp_major(&mut body, path, n_samples, variants, filter)
}

/// Validate magic, mode and total length, in that order.
///
/// Order matters: a file that is not a `.bed` at all should be reported as such rather
/// than as a size mismatch, and an individual-major `.bed` is a *format* problem even
/// though its length happens to be checkable too.
fn check_header_and_len(
    header: &[u8],
    found: u64,
    path: &Path,
    n_samples: usize,
    n_variants: usize,
) -> Result<()> {
    if header[0] != MAGIC[0] || header[1] != MAGIC[1] {
        return Err(IoError::BadMagic {
            path: path.to_path_buf(),
        });
    }
    match header[2] {
        MODE_SNP_MAJOR => {}
        MODE_INDIVIDUAL_MAJOR => {
            return Err(IoError::IndividualMajor {
                path: path.to_path_buf(),
            })
        }
        mode => {
            return Err(IoError::BadMode {
                path: path.to_path_buf(),
                mode,
            })
        }
    }
    let expected = expected_bed_len(n_samples, n_variants);
    if found != expected {
        return Err(IoError::SizeMismatch {
            path: path.to_path_buf(),
            expected,
            found,
        });
    }
    Ok(())
}

/// Stream the genotype block, transposing retained rows into the bit planes.
fn decode_snp_major<R: Read>(
    reader: &mut R,
    path: &Path,
    n_samples: usize,
    variants: &[Variant],
    filter: VariantFilter,
) -> Result<(Genotypes, Vec<usize>)> {
    let row_len = bytes_per_variant(n_samples);
    let n_kept = variants.iter().filter(|v| is_kept(v, filter)).count();
    let words = n_kept.div_ceil(64);

    let mut plane0 = vec![vec![0u64; words]; n_samples];
    let mut plane1 = vec![vec![0u64; words]; n_samples];
    let mut kept = Vec::with_capacity(n_kept);
    let mut row = vec![0u8; row_len];

    for (index, variant) in variants.iter().enumerate() {
        reader
            .read_exact(&mut row)
            .map_err(|source| IoError::Read {
                path: path.to_path_buf(),
                source,
            })?;
        if !is_kept(variant, filter) {
            continue;
        }
        let m = kept.len();
        let (word, bit) = (m / 64, m % 64);

        for (b, &byte) in row.iter().enumerate() {
            // 0x55 is four missing calls, which set nothing in either plane.
            if byte == 0x55 {
                continue;
            }
            let base = b * 4;
            let lanes = min(4, n_samples - base);
            for lane in 0..lanes {
                // The LOW-order bit pair is the FIRST sample in the byte.
                let code = (byte >> (2 * lane)) & 0b11;
                let (p0, p1) = planes_for_code(code);
                plane0[base + lane][word] |= p0 << bit;
                plane1[base + lane][word] |= p1 << bit;
            }
        }
        kept.push(index);
    }

    Ok((
        Genotypes {
            plane0,
            plane1,
            n_samples,
            n_variants: n_kept,
        },
        kept,
    ))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;
    use std::sync::atomic::{AtomicUsize, Ordering};

    /// The genotype a `(plane0, plane1)` bit pair denotes, for readable assertions.
    #[derive(Debug, Clone, Copy, PartialEq, Eq)]
    enum Call {
        HomA1,
        Het,
        HomA2,
        Missing,
    }

    fn call_at(g: &Genotypes, sample: usize, m: usize) -> Call {
        let (word, bit) = (m / 64, m % 64);
        let p0 = (g.plane0[sample][word] >> bit) & 1;
        let p1 = (g.plane1[sample][word] >> bit) & 1;
        match (p0, p1) {
            (1, 1) => Call::HomA1,
            (0, 1) => Call::Het,
            (1, 0) => Call::HomA2,
            _ => Call::Missing,
        }
    }

    fn hex(s: &str) -> Vec<u8> {
        let digits: Vec<u8> = s
            .bytes()
            .filter(|b| !b.is_ascii_whitespace())
            .map(|b| match b {
                b'0'..=b'9' => b - b'0',
                b'a'..=b'f' => b - b'a' + 10,
                b'A'..=b'F' => b - b'A' + 10,
                other => panic!("not a hex digit: {:?}", other as char),
            })
            .collect();
        assert_eq!(digits.len() % 2, 0, "odd number of hex digits");
        digits.chunks(2).map(|c| (c[0] << 4) | c[1]).collect()
    }

    fn variant(chrom: &str, id: &str, bp: i64, a1: &str, a2: &str) -> Variant {
        Variant {
            chrom: chrom.to_string(),
            id: id.to_string(),
            cm: 0.0,
            bp,
            a1: a1.to_string(),
            a2: a2.to_string(),
        }
    }

    fn tmp_dir(tag: &str) -> PathBuf {
        static N: AtomicUsize = AtomicUsize::new(0);
        let n = N.fetch_add(1, Ordering::Relaxed);
        let dir =
            std::env::temp_dir().join(format!("open-king-io-{}-{}-{}", std::process::id(), tag, n));
        fs::create_dir_all(&dir).unwrap();
        dir
    }

    // ---------------------------------------------------------------------------------
    // The PLINK-written fixture.
    //
    // Built with PLINK 1.9 from this .ped (6 samples, all female so chrX stays diploid):
    //
    //   fam1 s1 0 0 2 -9   A A   C C   0 0   G G   A A
    //   fam1 s2 0 0 2 -9   A G   C C   0 0   G T   A A
    //   fam2 s3 0 0 2 -9   G G   C T   0 0   T T   A C
    //   fam2 s4 0 0 2 -9   0 0   T T   0 0   G G   C C
    //   fam3 s5 0 0 2 -9   G G   0 0   0 0   G T   A C
    //   fam3 s6 0 0 2 -9   G G   C T   0 0   G T   0 0
    //
    // and this .map (chr1 x3, chr2 x1, chrX x1). PLINK chose A1 = the minor allele, so
    // the .bim it wrote is:
    //
    //   1  v1  0  1000  A  G
    //   1  v2  0  2000  T  C
    //   1  v3  0  3000  0  0     (all calls missing -> both alleles "0")
    //   2  v4  0  1000  T  G
    //   23 vx  0  1000  C  A
    //
    // `xxd x.bed` gave exactly the 13 bytes below.
    // ---------------------------------------------------------------------------------
    const TINY_BED_HEX: &str = "6c1b01 780f 2f09 5505 cb0a 2f06";

    fn tiny_variants() -> Vec<Variant> {
        vec![
            variant("1", "v1", 1000, "A", "G"),
            variant("1", "v2", 2000, "T", "C"),
            variant("1", "v3", 3000, "0", "0"),
            variant("2", "v4", 1000, "T", "G"),
            variant("23", "vx", 1000, "C", "A"),
        ]
    }

    const TINY_FAM: &str = "\
fam1 s1 0 0 2 -9
fam1 s2 0 0 2 -9
fam2 s3 0 0 2 -9
fam2 s4 0 0 2 -9
fam3 s5 0 0 2 -9
fam3 s6 0 0 2 -9
";

    const TINY_BIM: &str = "\
1\tv1\t0\t1000\tA\tG
1\tv2\t0\t2000\tT\tC
1\tv3\t0\t3000\t0\t0
2\tv4\t0\t1000\tT\tG
23\tvx\t0\t1000\tC\tA
";

    /// The genotype matrix as written into the `.ped`, re-expressed against the A1/A2 that
    /// PLINK assigned. `expected[s][v]`.
    fn tiny_expected() -> [[Call; 5]; 6] {
        use Call::*;
        [
            //  v1(A/G)  v2(T/C)  v3(0/0)  v4(T/G)  vx(C/A)
            [HomA1, HomA2, Missing, HomA2, HomA2], // s1: AA CC 00 GG AA
            [Het, HomA2, Missing, Het, HomA2],     // s2: AG CC 00 GT AA
            [HomA2, Het, Missing, HomA1, Het],     // s3: GG CT 00 TT AC
            [Missing, HomA1, Missing, HomA2, HomA1], // s4: 00 TT 00 GG CC
            [HomA2, Missing, Missing, Het, Het],   // s5: GG 00 00 GT AC
            [HomA2, Het, Missing, Het, Missing],   // s6: GG CT 00 GT 00
        ]
    }

    #[test]
    fn plink_written_genotypes_round_trip() {
        let bytes = hex(TINY_BED_HEX);
        assert_eq!(bytes.len(), 13);
        // header, then 5 rows of ceil(6/4) = 2 bytes
        assert_eq!(bytes.len() as u64, expected_bed_len(6, 5));
        assert_eq!(&bytes[..3], &[0x6c, 0x1b, 0x01]);

        let variants = tiny_variants();
        let (g, kept) =
            read_bed_bytes(&bytes, Path::new("x.bed"), 6, &variants, VariantFilter::All).unwrap();
        assert_eq!(kept, vec![0, 1, 2, 3, 4]);
        assert_eq!((g.n_samples, g.n_variants), (6, 5));

        let expected = tiny_expected();
        for (s, row) in expected.iter().enumerate() {
            for (v, want) in row.iter().enumerate() {
                assert_eq!(
                    call_at(&g, s, v),
                    *want,
                    "sample {s} variant {} ({})",
                    v,
                    variants[v].id
                );
            }
        }
    }

    #[test]
    fn low_order_bit_pair_is_the_first_sample_in_the_byte() {
        // Row 1 of the fixture is 0x78 0x0f. 0x78 = 0b01_11_10_00, so reading the byte
        // from the low end gives s1..s4 = 00, 10, 11, 01 = HomA1, Het, HomA2, Missing —
        // which is exactly what the .ped says. Reading it from the high end would give
        // the reverse and is what this test exists to rule out.
        let bytes = hex(TINY_BED_HEX);
        assert_eq!(bytes[3], 0x78);
        let (g, _) = read_bed_bytes(
            &bytes,
            Path::new("x.bed"),
            6,
            &tiny_variants(),
            VariantFilter::All,
        )
        .unwrap();
        assert_eq!(call_at(&g, 0, 0), Call::HomA1);
        assert_eq!(call_at(&g, 1, 0), Call::Het);
        assert_eq!(call_at(&g, 2, 0), Call::HomA2);
        assert_eq!(call_at(&g, 3, 0), Call::Missing);
        // second byte of the row holds samples 5 and 6, again from the low end
        assert_eq!(bytes[4], 0x0f);
        assert_eq!(call_at(&g, 4, 0), Call::HomA2);
        assert_eq!(call_at(&g, 5, 0), Call::HomA2);
    }

    #[test]
    fn code_to_plane_table_is_the_documented_one() {
        assert_eq!(planes_for_code(0b00), (1, 1)); // hom A1/A1
        assert_eq!(planes_for_code(0b01), (0, 0)); // MISSING, not a genotype
        assert_eq!(planes_for_code(0b10), (0, 1)); // het
        assert_eq!(planes_for_code(0b11), (1, 0)); // hom A2/A2
                                                   // the two derived masks the kernels rely on
        for code in 0..4u8 {
            let (p0, p1) = planes_for_code(code);
            let is_hom = code == 0b00 || code == 0b11;
            let non_missing = code != 0b01;
            assert_eq!(p0 == 1, is_hom, "plane0 is the is-homozygous mask");
            assert_eq!(
                (p0 | p1) == 1,
                non_missing,
                "plane0|plane1 is the called mask"
            );
        }
    }

    #[test]
    fn missing_calls_are_zero_in_both_planes() {
        let bytes = hex(TINY_BED_HEX);
        let (g, _) = read_bed_bytes(
            &bytes,
            Path::new("x.bed"),
            6,
            &tiny_variants(),
            VariantFilter::All,
        )
        .unwrap();
        // v3 is missing for everyone: its bit is clear in both planes for every sample
        for s in 0..6 {
            assert_eq!(g.plane0[s][0] >> 2 & 1, 0, "sample {s} plane0 at v3");
            assert_eq!(g.plane1[s][0] >> 2 & 1, 0, "sample {s} plane1 at v3");
        }
        // the individually missing calls too: s4/v1, s5/v2, s6/vx
        for (s, v) in [(3, 0), (4, 1), (5, 4)] {
            assert_eq!(call_at(&g, s, v), Call::Missing);
        }
        // and a missing call is distinguishable from every real genotype
        assert_ne!(call_at(&g, 3, 1), Call::Missing);
    }

    #[test]
    fn variant_filter_selects_which_rows_reach_the_planes() {
        let bytes = hex(TINY_BED_HEX);
        let variants = tiny_variants();
        let read = |f| read_bed_bytes(&bytes, Path::new("x.bed"), 6, &variants, f).unwrap();

        let (all, kept) = read(VariantFilter::All);
        assert_eq!(kept, vec![0, 1, 2, 3, 4]);
        assert_eq!(all.n_variants, 5);

        // autosomes: chr1 x3 + chr2 x1; the chrX variant is dropped
        let (auto, kept) = read(VariantFilter::Autosomes);
        assert_eq!(kept, vec![0, 1, 2, 3]);
        assert_eq!(auto.n_variants, 4);

        let (chr1, kept) = read(VariantFilter::Chromosome(1));
        assert_eq!(kept, vec![0, 1, 2]);
        assert_eq!(chr1.n_variants, 3);

        let (chr23, kept) = read(VariantFilter::Chromosome(23));
        assert_eq!(kept, vec![4]);
        assert_eq!(chr23.n_variants, 1);

        let (none, kept) = read(VariantFilter::Chromosome(5));
        assert!(kept.is_empty());
        assert_eq!(none.n_variants, 0);
        assert_eq!(none.words_per_sample(), 0);
        assert!(none.plane0[0].is_empty());

        // Retained variants are renumbered densely: with chr2 only, the single kept row
        // sits at bit 0, not at its .bim position 3.
        let (chr2, kept) = read(VariantFilter::Chromosome(2));
        assert_eq!(kept, vec![3]);
        for (s, row) in tiny_expected().iter().enumerate() {
            assert_eq!(call_at(&chr2, s, 0), row[3], "sample {s}");
        }
        // ... and with chrX only, s3's het at vx lands at bit 0
        assert_eq!(call_at(&chr23, 2, 0), Call::Het);
    }

    #[test]
    fn filtering_does_not_disturb_the_surviving_rows() {
        let bytes = hex(TINY_BED_HEX);
        let variants = tiny_variants();
        let (auto, kept) = read_bed_bytes(
            &bytes,
            Path::new("x.bed"),
            6,
            &variants,
            VariantFilter::Autosomes,
        )
        .unwrap();
        for (m, &v) in kept.iter().enumerate() {
            for (s, row) in tiny_expected().iter().enumerate() {
                assert_eq!(call_at(&auto, s, m), row[v], "sample {s} variant {v}");
            }
        }
    }

    #[test]
    fn trailing_bits_of_the_final_word_are_zero() {
        // 70 variants -> 2 words, so bits 70..128 must be clear. Every call is hom A1/A1
        // (code 00, byte 0x00), which sets a bit in BOTH planes, so any stray tail bit
        // would show up in both.
        let n_samples = 3;
        let n_variants = 70;
        let mut bytes = vec![0x6c, 0x1b, 0x01];
        bytes.extend(std::iter::repeat(0x00).take(n_variants * bytes_per_variant(n_samples)));
        let variants: Vec<Variant> = (0..n_variants)
            .map(|i| variant("1", &format!("v{i}"), i as i64 + 1, "A", "G"))
            .collect();

        let (g, _) = read_bed_bytes(
            &bytes,
            Path::new("t.bed"),
            n_samples,
            &variants,
            VariantFilter::All,
        )
        .unwrap();
        assert_eq!(g.n_variants, 70);
        assert_eq!(g.words_per_sample(), 2);
        for s in 0..n_samples {
            assert_eq!(g.plane0[s].len(), 2);
            assert_eq!(g.plane0[s][0], u64::MAX);
            assert_eq!(g.plane1[s][0], u64::MAX);
            // bits 0..6 of word 1 set, 6..64 clear
            assert_eq!(g.plane0[s][1], (1u64 << 6) - 1, "sample {s} plane0 tail");
            assert_eq!(g.plane1[s][1], (1u64 << 6) - 1, "sample {s} plane1 tail");
        }

        // Same again, but reaching 70 retained variants by filtering 100 rows, so the
        // tail is determined by the retained count and not by the file's row count.
        let mut mixed = vec![0x6c, 0x1b, 0x01];
        mixed.extend(std::iter::repeat(0x00).take(100 * bytes_per_variant(n_samples)));
        let mixed_variants: Vec<Variant> = (0..100)
            .map(|i| {
                let chrom = if i < 70 { "1" } else { "23" };
                variant(chrom, &format!("v{i}"), i as i64 + 1, "A", "G")
            })
            .collect();
        let (g, kept) = read_bed_bytes(
            &mixed,
            Path::new("t.bed"),
            n_samples,
            &mixed_variants,
            VariantFilter::Autosomes,
        )
        .unwrap();
        assert_eq!(kept.len(), 70);
        assert_eq!(g.words_per_sample(), 2);
        for s in 0..n_samples {
            assert_eq!(g.plane0[s][1], (1u64 << 6) - 1);
            assert_eq!(g.plane1[s][1], (1u64 << 6) - 1);
        }
    }

    #[test]
    fn padding_bits_in_the_last_byte_of_a_row_are_ignored() {
        // 6 samples -> 2 bytes per row, of which the top 4 bits of byte 2 are padding.
        // PLINK writes them as 0 (= code 00 = hom A1/A1); a reader that decoded them
        // would invent samples 7 and 8.
        let bytes = hex(TINY_BED_HEX);
        let (g, _) = read_bed_bytes(
            &bytes,
            Path::new("x.bed"),
            6,
            &tiny_variants(),
            VariantFilter::All,
        )
        .unwrap();
        assert_eq!(g.plane0.len(), 6);
        assert_eq!(g.plane1.len(), 6);
        assert_eq!(g.n_samples, 6);
    }

    #[test]
    fn individual_major_files_are_rejected() {
        let mut bytes = hex(TINY_BED_HEX);
        bytes[2] = MODE_INDIVIDUAL_MAJOR;
        match read_bed_bytes(
            &bytes,
            Path::new("x.bed"),
            6,
            &tiny_variants(),
            VariantFilter::All,
        ) {
            Err(IoError::IndividualMajor { path }) => assert_eq!(path, Path::new("x.bed")),
            other => panic!("expected IndividualMajor, got {other:?}"),
        }
    }

    #[test]
    fn bad_magic_and_bad_mode_are_rejected() {
        let mut wrong_magic = hex(TINY_BED_HEX);
        wrong_magic[1] = 0x1c;
        assert!(matches!(
            read_bed_bytes(
                &wrong_magic,
                Path::new("x.bed"),
                6,
                &tiny_variants(),
                VariantFilter::All
            ),
            Err(IoError::BadMagic { .. })
        ));

        let mut wrong_mode = hex(TINY_BED_HEX);
        wrong_mode[2] = 0x02;
        match read_bed_bytes(
            &wrong_mode,
            Path::new("x.bed"),
            6,
            &tiny_variants(),
            VariantFilter::All,
        ) {
            Err(IoError::BadMode { mode, .. }) => assert_eq!(mode, 0x02),
            other => panic!("expected BadMode, got {other:?}"),
        }

        // magic is checked before length, so a non-.bed file of the wrong size still
        // reports the real problem
        assert!(matches!(
            read_bed_bytes(
                b"not a bed file at all",
                Path::new("x.bed"),
                6,
                &tiny_variants(),
                VariantFilter::All
            ),
            Err(IoError::BadMagic { .. })
        ));
    }

    #[test]
    fn truncated_and_overlong_files_are_rejected() {
        let full = hex(TINY_BED_HEX);
        // one byte short of the 13 a 6-sample, 5-variant fileset needs
        match read_bed_bytes(
            &full[..12],
            Path::new("x.bed"),
            6,
            &tiny_variants(),
            VariantFilter::All,
        ) {
            Err(IoError::SizeMismatch {
                expected, found, ..
            }) => {
                assert_eq!((expected, found), (13, 12));
            }
            other => panic!("expected SizeMismatch, got {other:?}"),
        }
        // a whole row short
        assert!(matches!(
            read_bed_bytes(
                &full[..11],
                Path::new("x.bed"),
                6,
                &tiny_variants(),
                VariantFilter::All
            ),
            Err(IoError::SizeMismatch { .. })
        ));
        // header only
        match read_bed_bytes(
            &full[..3],
            Path::new("x.bed"),
            6,
            &tiny_variants(),
            VariantFilter::All,
        ) {
            Err(IoError::SizeMismatch { found, .. }) => assert_eq!(found, 3),
            other => panic!("expected SizeMismatch, got {other:?}"),
        }
        // shorter than the header
        match read_bed_bytes(
            &full[..2],
            Path::new("x.bed"),
            6,
            &tiny_variants(),
            VariantFilter::All,
        ) {
            Err(IoError::SizeMismatch { found, .. }) => assert_eq!(found, 2),
            other => panic!("expected SizeMismatch, got {other:?}"),
        }
        // trailing junk
        let mut long = full.clone();
        long.push(0);
        match read_bed_bytes(
            &long,
            Path::new("x.bed"),
            6,
            &tiny_variants(),
            VariantFilter::All,
        ) {
            Err(IoError::SizeMismatch {
                expected, found, ..
            }) => {
                assert_eq!((expected, found), (13, 14));
            }
            other => panic!("expected SizeMismatch, got {other:?}"),
        }
        // a .bed whose .fam has the wrong sample count fails the same way: 7 samples
        // would need 2 bytes per row as well, but 8 samples need 2 and 9 need 3
        assert!(matches!(
            read_bed_bytes(
                &full,
                Path::new("x.bed"),
                9,
                &tiny_variants(),
                VariantFilter::All
            ),
            Err(IoError::SizeMismatch {
                expected: 18,
                found: 13,
                ..
            })
        ));
    }

    #[test]
    fn size_helpers_round_up_to_whole_bytes() {
        assert_eq!(bytes_per_variant(0), 0);
        assert_eq!(bytes_per_variant(1), 1);
        assert_eq!(bytes_per_variant(4), 1);
        assert_eq!(bytes_per_variant(5), 2);
        assert_eq!(bytes_per_variant(8), 2);
        assert_eq!(bytes_per_variant(9), 3);
        assert_eq!(expected_bed_len(6, 5), 13);
        assert_eq!(expected_bed_len(8, 320), 643);
        assert_eq!(expected_bed_len(0, 10), 3);
    }

    #[test]
    fn paths_are_derived_from_a_prefix_or_a_bed_path() {
        let from_prefix = paths_for(Path::new("/data/study"));
        assert_eq!(from_prefix.0, Path::new("/data/study.bed"));
        assert_eq!(from_prefix.1, Path::new("/data/study.bim"));
        assert_eq!(from_prefix.2, Path::new("/data/study.fam"));

        let from_bed = paths_for(Path::new("/data/study.bed"));
        assert_eq!(from_bed, from_prefix);
        // uppercase extension is still a .bed path
        assert_eq!(
            paths_for(Path::new("/data/study.BED")).1,
            Path::new("/data/study.bim")
        );

        // a dotted prefix gets the suffix appended, not substituted
        let dotted = paths_for(Path::new("/data/study.v2"));
        assert_eq!(dotted.0, Path::new("/data/study.v2.bed"));
        assert_eq!(dotted.2, Path::new("/data/study.v2.fam"));
        // ... and the .bed form of that same fileset resolves identically
        assert_eq!(paths_for(Path::new("/data/study.v2.bed")), dotted);
    }

    #[test]
    fn read_fileset_loads_a_triple_from_disk() {
        let dir = tmp_dir("triple");
        fs::write(dir.join("x.bed"), hex(TINY_BED_HEX)).unwrap();
        fs::write(dir.join("x.bim"), TINY_BIM).unwrap();
        fs::write(dir.join("x.fam"), TINY_FAM).unwrap();

        // by prefix and by .bed path, identically
        for arg in [dir.join("x"), dir.join("x.bed")] {
            let fs_ = read_fileset(&arg, VariantFilter::Autosomes, None, None).unwrap();
            assert_eq!(fs_.samples.len(), 6);
            assert_eq!(fs_.samples[0].fid, "fam1");
            assert_eq!(fs_.samples[5].iid, "s6");
            // every .bim row is retained in `variants` ...
            assert_eq!(fs_.variants.len(), 5);
            // ... but only the autosomal ones reach the planes
            assert_eq!(fs_.kept, vec![0, 1, 2, 3]);
            assert_eq!(fs_.genotypes.n_variants, 4);
            assert_eq!(fs_.genotypes.words_per_sample(), 1);
            for (m, &v) in fs_.kept.iter().enumerate() {
                for (s, row) in tiny_expected().iter().enumerate() {
                    assert_eq!(
                        call_at(&fs_.genotypes, s, m),
                        row[v],
                        "sample {s} variant {v}"
                    );
                }
            }
        }
        fs::remove_dir_all(&dir).unwrap();
    }

    #[test]
    fn fam_and_bim_overrides_replace_the_derived_paths() {
        let dir = tmp_dir("override");
        fs::write(dir.join("x.bed"), hex(TINY_BED_HEX)).unwrap();
        fs::write(dir.join("x.bim"), TINY_BIM).unwrap();
        fs::write(dir.join("x.fam"), TINY_FAM).unwrap();
        // alternates: same shape, different content
        let alt_fam = TINY_FAM.replace("fam1", "ALT").replace("fam2", "ALT");
        fs::write(dir.join("alt.fam"), &alt_fam).unwrap();
        let alt_bim = TINY_BIM.replace("\n1\tv3", "\n23\tv3");
        fs::write(dir.join("alt.bim"), &alt_bim).unwrap();

        let with_fam = read_fileset(
            &dir.join("x"),
            VariantFilter::Autosomes,
            Some(&dir.join("alt.fam")),
            None,
        )
        .unwrap();
        assert_eq!(with_fam.samples[0].fid, "ALT");
        assert_eq!(with_fam.kept, vec![0, 1, 2, 3]);

        // moving v3 to chrX in the override changes what the autosome filter keeps
        let with_bim = read_fileset(
            &dir.join("x"),
            VariantFilter::Autosomes,
            None,
            Some(&dir.join("alt.bim")),
        )
        .unwrap();
        assert_eq!(with_bim.samples[0].fid, "fam1");
        assert_eq!(with_bim.kept, vec![0, 1, 3]);
        assert_eq!(with_bim.genotypes.n_variants, 3);

        // a missing override is an Open error against the override's own path
        match read_fileset(
            &dir.join("x"),
            VariantFilter::All,
            Some(&dir.join("nope.fam")),
            None,
        ) {
            Err(IoError::Open { path, .. }) => assert_eq!(path, dir.join("nope.fam")),
            other => panic!("expected Open error, got {other:?}"),
        }
        fs::remove_dir_all(&dir).unwrap();
    }

    #[test]
    fn read_fileset_rejects_duplicate_samples() {
        // The reference binary makes this fatal: "Family fam1: Person s1 is duplicated".
        let dir = tmp_dir("dup");
        fs::write(dir.join("x.bed"), hex(TINY_BED_HEX)).unwrap();
        fs::write(dir.join("x.bim"), TINY_BIM).unwrap();
        fs::write(dir.join("x.fam"), TINY_FAM.replace("s2", "s1")).unwrap();
        match read_fileset(&dir.join("x"), VariantFilter::All, None, None) {
            Err(IoError::DuplicateSample { fid, iid }) => {
                assert_eq!((fid.as_str(), iid.as_str()), ("fam1", "s1"));
            }
            other => panic!("expected DuplicateSample, got {other:?}"),
        }
        fs::remove_dir_all(&dir).unwrap();
    }

    #[test]
    fn streamed_and_in_memory_readers_agree() {
        let dir = tmp_dir("stream");
        let bytes = hex(TINY_BED_HEX);
        let path = dir.join("x.bed");
        fs::write(&path, &bytes).unwrap();
        let variants = tiny_variants();

        let (from_file, kept_file) = read_bed(&path, 6, &variants, VariantFilter::All).unwrap();
        let (from_mem, kept_mem) =
            read_bed_bytes(&bytes, &path, 6, &variants, VariantFilter::All).unwrap();
        assert_eq!(kept_file, kept_mem);
        assert_eq!(from_file.plane0, from_mem.plane0);
        assert_eq!(from_file.plane1, from_mem.plane1);

        // and the on-disk reader reports the same failures
        fs::write(&path, &bytes[..12]).unwrap();
        assert!(matches!(
            read_bed(&path, 6, &variants, VariantFilter::All),
            Err(IoError::SizeMismatch {
                expected: 13,
                found: 12,
                ..
            })
        ));
        match read_bed(&dir.join("absent.bed"), 6, &variants, VariantFilter::All) {
            Err(IoError::Open { .. }) => {}
            other => panic!("expected Open error, got {other:?}"),
        }
        fs::remove_dir_all(&dir).unwrap();
    }

    // ---------------------------------------------------------------------------------
    // Cross-check against KING 2.3.2.
    //
    // 8 samples x 320 variants (200 on chr1, 100 on chr2, 20 on chrX), built by PLINK 1.9
    // from a .ped with ~5% missing calls. The reference binary reported
    //
    //   Genotype data consist of 300 autosome SNPs, 20 X-chromosome SNPs
    //   Autosome genotypes stored in 5 words for each of 8 individuals
    //
    // confirming both that its default relatedness path is autosomes-only and that its
    // word count is ceil(300/64) = 5, i.e. `Genotypes::words_per_sample`.
    //
    // `king -b big.bed --ibs` then printed, for all 28 pairs, the N_SNP / N_IBS0 /
    // NHetHet / NHomHom / N_Het1 / N_Het2 counts tabulated below. Recomputing them here
    // from the decoded bit planes pins the code->plane mapping to the reference's own
    // reading of the same bytes.
    // ---------------------------------------------------------------------------------
    const BIG_BED_HEX: &str = "
         6c1b01fabeecbef7ff2fbfbffcce8baefeffda2beeaaee2beffaeff3fafaffbf
         9b82fffe3fafbafebab4f8f73ebbcea3baabeefefeb9fbceeef9afeec2aaeae5
         aaf7bfe23ebfefebaaeac7a3beafabbbe8a2bdffacbb9efaeab726afabcbbffe
         eefffbf6bbefbe3eaaf2ef378bebbe8fabfeaeba5eeefa5efaafbe36fbffbeeb
         acedd63aeffe38ee2b786bff9aaefaffabaeffb7bbadbfbfef3fafbbebbbbb3a
         eaabab3bbbfae39b7efebfaafbfaebcadeefbfbe1ffe7bfbef973ebefeeef27e
         aa2e4ee3fb1b8be0faba9c9befaeeaef7fbabbebcfaa9bcbccfa8d7aaabbdd60
         fefbf6efbbf1fefebeafdff66eeaefdceebba0fbeabefabbbcffbbbffe8fbe3d
         bc7abb32aba7faf2b6faec33eefbe19abffab3ffafaefbbbb9effffcaafbe3af
         affbbfceca87bfeeeaeabaefabeba8f33afeabeec2fee2faaf8ffff3baccb33e
         0bfbfebacabeeaafdbbcacfbffbf8bdf8baeaffffbc3abfe9fb7afeffb8bb3eb
         3e52ece2bfef638f3bfbabba6bfeefab3ee3effcffffee8ffbaf28abd1fafa2e
         feaabffef6fe32ae2cffebfebaaebae0b7bc9eab2a78f2bffbbb6bfbef0bbbae
         fbfeeae3aab8fbfbeb2fa2ef3a72f4e2acecbae78c6e3fab8ce3afeaaeb2befc
         aaffbeabeac3d2daeeab9ffff0fefbe3fef3afb2b4ffaecdefaebeeab8cff3f3
         beffffadcebfdf45e7fbfaff624eaf4f23bfaffefeef6bfbffbffbfabaeb7aab
         6bbae238fbbeefa3d3aefbffee2abefabbbf73dfbef4cfaeeaeaeaeefb8d6bea
         3aa6fea3af9fcffffe2b72e98a3b2eefff8bbb2f3bbc6bcacbba3ea38ff4affe
         8aff3be7a3bbff3b2afb82efffaaeaaefdaffcaff3fbbfaf2affdcf7aafe6ffe
         bbeb6ebaeaa3fffe7b02fb8ada1fbf2eafcabe3ac9eeaa3effff2efee33feeef
         f8eeff";

    /// `(N_SNP, N_IBS0, NHetHet, NHomHom, N_Het1, N_Het2)`, in the reference's own column
    /// order.
    type Counts = (u32, u32, u32, u32, u32, u32);
    /// One `king.ibs` row: the two sample indices, then the counts.
    type IbsRow = (usize, usize, Counts);

    /// Every row `king -b big.bed --ibs` wrote to `king.ibs`, transcribed.
    const KING_IBS: [IbsRow; 28] = [
        (0, 1, (269, 28, 42, 99, 110, 102)),
        (0, 2, (276, 28, 42, 103, 110, 105)),
        (0, 3, (276, 25, 37, 102, 106, 105)),
        (0, 4, (281, 23, 58, 109, 111, 119)),
        (0, 5, (280, 26, 50, 103, 109, 118)),
        (0, 6, (268, 23, 41, 96, 108, 105)),
        (0, 7, (271, 25, 40, 87, 109, 115)),
        (1, 2, (268, 27, 34, 97, 100, 105)),
        (1, 3, (267, 29, 40, 105, 101, 101)),
        (1, 4, (272, 23, 42, 97, 100, 117)),
        (1, 5, (271, 33, 46, 99, 103, 115)),
        (1, 6, (258, 31, 44, 102, 97, 103)),
        (1, 7, (262, 29, 41, 96, 97, 110)),
        (2, 3, (277, 26, 49, 114, 108, 104)),
        (2, 4, (281, 35, 52, 108, 107, 118)),
        (2, 5, (278, 23, 46, 98, 107, 119)),
        (2, 6, (267, 21, 38, 96, 103, 106)),
        (2, 7, (271, 24, 44, 96, 103, 116)),
        (3, 4, (280, 24, 37, 94, 104, 119)),
        (3, 5, (279, 28, 47, 102, 107, 117)),
        (3, 6, (266, 25, 36, 104, 95, 103)),
        (3, 7, (270, 19, 41, 94, 101, 116)),
        (4, 5, (283, 28, 51, 94, 119, 121)),
        (4, 6, (274, 22, 43, 92, 116, 109)),
        (4, 7, (275, 19, 43, 84, 117, 117)),
        (5, 6, (270, 22, 48, 93, 117, 108)),
        (5, 7, (272, 26, 51, 93, 113, 117)),
        (6, 7, (262, 28, 42, 92, 100, 112)),
    ];

    fn big_variants() -> Vec<Variant> {
        let mut v = Vec::with_capacity(320);
        for i in 0..200 {
            v.push(variant(
                "1",
                &format!("snp{i}"),
                (i + 1) as i64 * 1000,
                "A",
                "G",
            ));
        }
        for i in 200..300 {
            v.push(variant(
                "2",
                &format!("snp{i}"),
                (i - 199) as i64 * 1000,
                "A",
                "G",
            ));
        }
        for i in 0..20 {
            v.push(variant(
                "23",
                &format!("xsnp{i}"),
                (i + 1) as i64 * 1000,
                "C",
                "T",
            ));
        }
        v
    }

    /// [`Counts`] straight off the bit planes, using only the semantics documented on
    /// `Genotypes`.
    fn pair_counts(g: &Genotypes, i: usize, j: usize) -> Counts {
        let (mut n_snp, mut ibs0, mut het_het, mut hom_hom, mut het_i, mut het_j) =
            (0, 0, 0, 0, 0, 0);
        for w in 0..g.words_per_sample() {
            let (a0, a1) = (g.plane0[i][w], g.plane1[i][w]);
            let (b0, b1) = (g.plane0[j][w], g.plane1[j][w]);
            let called = (a0 | a1) & (b0 | b1);
            // heterozygous is plane0 = 0, plane1 = 1
            let ha = !a0 & a1;
            let hb = !b0 & b1;
            n_snp += called.count_ones();
            het_i += (ha & called).count_ones();
            het_j += (hb & called).count_ones();
            het_het += (ha & hb).count_ones();
            hom_hom += (a0 & b0).count_ones();
            // opposite homozygotes: both homozygous (plane0 set) but different plane1
            ibs0 += (a0 & b0 & (a1 ^ b1)).count_ones();
        }
        (n_snp, ibs0, het_het, hom_hom, het_i, het_j)
    }

    #[test]
    fn pairwise_counts_match_king_2_3_2() {
        let bytes = hex(BIG_BED_HEX);
        assert_eq!(bytes.len() as u64, expected_bed_len(8, 320));
        let variants = big_variants();
        let (g, kept) = read_bed_bytes(
            &bytes,
            Path::new("big.bed"),
            8,
            &variants,
            VariantFilter::Autosomes,
        )
        .unwrap();
        // matches "Genotype data consist of 300 autosome SNPs" and
        // "Autosome genotypes stored in 5 words for each of 8 individuals"
        assert_eq!(kept.len(), 300);
        assert_eq!(g.n_variants, 300);
        assert_eq!(g.words_per_sample(), 5);

        for &(i, j, counts) in KING_IBS.iter() {
            assert_eq!(
                pair_counts(&g, i, j),
                counts,
                "pair id{}/id{}",
                i + 1,
                j + 1
            );
        }
    }

    #[test]
    fn king_cross_check_is_sensitive_to_the_plane_mapping() {
        // Guards the test above: if the codes were swapped (11 read as hom A1/A1 and 00 as
        // hom A2/A2, the classic off-by-one in this format) the IBS0 counts would change.
        // Swapping plane1 does exactly that, since plane1 is what distinguishes the two
        // homozygote classes.
        let bytes = hex(BIG_BED_HEX);
        let (mut g, _) = read_bed_bytes(
            &bytes,
            Path::new("big.bed"),
            8,
            &big_variants(),
            VariantFilter::Autosomes,
        )
        .unwrap();
        let truth = pair_counts(&g, 0, 1);
        assert_eq!((truth.0, truth.1), (269, 28));
        // flip het <-> hom A2/A2 for sample 0 by clearing its plane0
        for w in 0..g.words_per_sample() {
            g.plane0[0][w] = 0;
        }
        assert_ne!(pair_counts(&g, 0, 1), truth);
    }
}
