# Phase A: independent structural validation

Base: main at a09d3eed4800bbe5c3aea37869fb181ea8fc8a6f. Historical version v0.1.0-legacy and reference grids are unchanged.

- [ ] Compare phi/psi to the public wwPDB residue-level validation reports and Bio3D across curated structures.
- [ ] Report coverage, circular angle differences, group-specific Ramachandran-class agreement and an explicit confusion matrix. Use independent source classifications, not regenerated labels from RamplotR.
- [ ] Validate selected structures spanning crystallography and cryo-EM; capture report checksums, source URLs, R versions and analysis settings in reproducible artifacts.
- [ ] Test XML parsing, stable residue matching, alt locations, insertion codes, missing angles and unequal classification vocabularies.
- [ ] Add a public manifest and an easy-to-run verification workflow; fail CI on coverage/angle regressions, not on genuine differences between reference models.
- [ ] Document why the four-region, four-reference-group RamplotR original method is not directly MolProbity-equivalent to its three regions and six groups.
- [ ] Pass CI and merge the independent validation workflow into main. ESMFold and AlphaFold are explicitly Phase B.
