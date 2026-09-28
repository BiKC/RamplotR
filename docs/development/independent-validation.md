# Independent experimental validation: implementation record

Base: main at a09d3eed4800bbe5c3aea37869fb181ea8fc8a6f. Historical version v0.1.0-legacy and reference grids are unchanged.

- [x] Compare phi/psi to the public wwPDB residue-level validation reports and Bio3D across curated structures.
- [x] Report coverage, circular angle differences, group-specific Ramachandran-class agreement and an explicit confusion matrix. Use independent source classifications, not regenerated labels from RamplotR.
- [x] Validate selected structures spanning crystallography and cryo-EM; capture report checksums, source URLs, R versions and analysis settings in reproducible artifacts.
- [x] Test XML parsing, stable residue matching, alt locations, insertion codes, missing angles and unequal classification vocabularies.
- [x] Add a public manifest and an easy-to-run verification workflow; fail CI on coverage/angle regressions, not on genuine differences between reference models.
- [x] Document why the four-region, four-reference-group RamplotR original method is not directly MolProbity-equivalent to its three regions and six groups.
- [x] Pass Ubuntu/Windows scientific regression, independent five-structure wwPDB validation, structure-validation benchmark and full live-browser CI, and merge PR #13 into main. AlphaFold and ESMFold prediction confidence are documented separately.

## Independent findings

[Five-structure results](../validation-results.md) and
[pinned SHA256 source manifests](../../validation/baseline-source-hashes-2026-09-28.csv)
were generated from official wwPDB validation reports. 2DQ4's seven official
outliers are not classified as outliers by the original RamplotR references.
The independent numerical geometry and group-specific class comparisons are
tested separately. Do not claim MolProbity-equivalent categories. Prediction-confidence analysis for AlphaFold and ESMFold is documented
in the separate [prediction-confidence guide](../prediction-confidence.md).

Completed integration: https://github.com/BiKC/RamplotR/pull/13
Reviewed independent validation workflow: https://github.com/BiKC/RamplotR/actions/runs/36476627652
