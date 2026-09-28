# Phase A: observed independent wwPDB validation results

**Date:** September 28, 2026. **Version:** original bundled reference
distributions, residue-aware classification, model 1. All samples were
processed in the same five-way independent workflow:
https://github.com/BiKC/RamplotR/actions/runs/36475593015

The pinned output metrics are in [baseline-2026-09-28.csv](../validation/baseline-2026-09-28.csv).
The source URLs and both official-file SHA256 checksums per sample are in
[baseline-source-hashes-2026-09-28.csv](../validation/baseline-source-hashes-2026-09-28.csv).
See [the method and limitations](wwpdb-validation.md) before quoting
classification agreement.

| Structure | Method | Independently matched finite angle pairs | Maximum circular difference | Agreement after descriptive three-way label mapping |
|---|---|---:|---:|---:|
| 1CRN | X-ray | 44/44 | 0.052° | 43/44 |
| 1UBQ | X-ray | 74/74 | 0.054° | 74/74 |
| 6VXX | cryo-EM | 2,844/2,844 | 0.055° | 2,796/2,844 |
| 2DQ4 | X-ray, outlier-rich | 682/682 | 0.055° | 651/682 |
| 1D3Z | NMR, model 1 | 74/74 | 0.055° | 74/74 |

All 3,718 eligible phi/psi pairs match independently reported wwPDB
angles to within 0.055° (the report prints its angles to one decimal
place). This is numerical validation of backbone geometry, not
equivalence of Ramachandran classification rules.

The classification comparisons differ at **80 of 3,718 residues** under
the explicitly documented, explanatory four-to-three-label crosswalk.
Of particular interest, **2DQ4 has seven official wwPDB outliers**, but
none of those same seven residues is outlier-labelled by RamplotR's
original density references. RamplotR independently flags one other
residue which wwPDB does not call an outlier. We preserve the seven
independent outlier identities in
[2dq4-wwpdb-outlier-examples.csv](../validation/2dq4-wwpdb-outlier-examples.csv)
for future method development. This is evidence of materially different
outlier-detection behaviour, **not** a reason to change the benchmark
gate to require perfect classification agreement.

The independent source XML and coordinate CIF are archived with every
CI artifact, along with all joined residues, reference-group confusion
matrices, source checksums and R session metadata. All input hashes and
quantitative snapshots are additionally pinned in the repository;
[tests/phase-a-baseline.R](../tests/phase-a-baseline.R) fails if an
official source changes without review or if the same pinned sources
produce different numbers.

## What these findings mean

The original RamplotR class labels cannot be cited or presented as
MolProbity/wwPDB-equivalent. The methods have different reference
populations, residue-type treatments and contour semantics. A future
MolProbity-compatible mode should source its validated six-class
reference grids, respect cis/trans proline and Ile/Val and compare its
residue-wise classifications against these pinned independent data,
including the seven 2DQ4 outliers.

For existing users and publications, keep the current mode available
for reproducibility, label its scientific provenance, and refer to
official wwPDB/MolProbity quality validation when making outlier
claims. The historical published release is separately preserved as
v0.1.0-legacy.
