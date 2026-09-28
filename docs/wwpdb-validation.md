# Phase A: independent wwPDB Ramachandran validation

RamplotR now supports residue-matched comparison against **official wwPDB
validation XML**, separately from its existing independent Bio3D angle
comparison. Unlike the synthetic unit-test fixtures, the public reports are
independent reference assessments.

## Important differences between methods

RamplotR's existing density grids and four groups are General, GLY, PRO and
preProline. Its four regions use cumulative-density contours of 85%, 98%
and 99.95%. The legacy dataset descends from Lovell et al. (2003);
other bundled datasets have separate provenance.

wwPDB's MolProbity-derived validation uses three categories (Favored, Allowed,
OUTLIER) and six backbone classes: general, Ile/Val, Gly, pre-Pro, trans-Pro
and cis-Pro. MolProbity uses different underlying reference populations and
density thresholds: favored about 98%, allowed/outlier approximately 99.95%.
See https://doi.org/10.1107/S0907444909042073 and
https://pmc.ncbi.nlm.nih.gov/articles/PMC5734394/.

For a *descriptive* 3-way contingency table only, collapse RamplotR Favoured
and Allowed into favored, Generously allowed into allowed, and Not allowed
into outlier. This label crosswalk does NOT make the underlying methods
scientifically interchangeable. Report disagreements by residue group, and
do not assert MolProbity equivalence based on agreement percentages.

## Public, reproducible corpus

The committed `validation/manifest.csv` identifies 1CRN and 1UBQ
(X-ray crystallography), 6VXX (cryo-EM viral spike complex), 2DQ4
(a challenging crystallographic example with reported Ramachandran outliers),
and 1D3Z (solution-NMR ensemble, model one). The dedicated
CI job downloads each original RCSB mmCIF and the separately generated
wwPDB validation XML. Each run preserves both source files, their SHA256
hashes, exact download URLs, git commit and R package versions.

For every real structure the analysis produces:

- `residue_comparison.csv`: all residues, matched identifiers including chain
  and insertion code, independent phi/psi comparisons and both class labels.
- `group_contingency.csv`: all 3×3 class combinations per RamplotR group.
- `summary.csv`: finite-angle coverage, maximum circular angular differences
  and observed class agreement/disagreement counts.
- `provenance.csv` and `r-session.txt`: exact input hashes and R environment.

The numerical gate requires at least 90% matching of RamplotR's finite
phi/psi pairs and at most 1.5° absolute circular deviation. Genuine label
differences across independent distribution families are RECORDED, not
treated as failed tests. If source coordinates or reports are re-released,
investigate provenance before changing tolerance thresholds.

Tests also include a plainly marked synthetic XML fixture with insertion
codes, alternate conformations, ±180° angles, class disagreements and
deliberate negative cases. This fixture tests parsing only and is not
independent scientific validation.

## Reproduce

Run from the repository root with R installed:

~~~bash
Rscript -e 'install.packages(c("bio3d", "xml2", "digest"))'
Rscript tests/wwpdb.R
Rscript benchmarks/compare_wwpdb.R 1CRN original benchmarks/output/wwpdb/1CRN
~~~

For exact replay, supply the downloaded compressed XML and mmCIF as fourth
and fifth CLI arguments. The independent wwPDB GitHub Actions workflow
runs all three samples. It stores full source and output artifacts for
90 days; freeze the hashes and outputs in a versioned release or DOI-backed
archive for publications.

## Scope

Only model 1 and unambiguous residue identities are compared. Alternate
conformations, aliases, modified residues and missing atoms reduce
comparable coverage rather than being silently declared matched.
Comparative results are specific to the original RamplotR references.
The AlphaFold and ESMFold confidence workflows belong to Phase B and
are not substitutes for independently measured experimental coordinates.
