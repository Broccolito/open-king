Supplementary runs backing the analysis in ../INDEX.md. Same capture format as the
core matrix (cmd.txt / stdout.txt / stderr.txt / exitcode.txt + output files).

  bigish__*__cpus1 / cpus4   --cpus determinism check (INDEX.md section 3)
  *__prefixZZZ               --prefix filename mapping (INDEX.md section 2)
  multifam__related_degree1__noscreen   evidence for quirk Q3

SPACE NOTE: in the bigish cpus1/cpus4 directories the output FILES were deleted after
verifying they are byte-identical to the corresponding core-matrix run
(../<same-name-without-__cpusN>/). MD5SUMS.txt in each directory records the hashes, so
the identity claim stays checkable without storing three copies of a 2.3 MB king.ibs0.
stdout.txt is kept in full — it is the only thing that differs between --cpus values.
