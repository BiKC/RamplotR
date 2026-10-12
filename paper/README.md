# RamplotR preprint draft

This directory contains the manuscript draft for RamplotR `main` at commit `6d1e4e95cdf5df4ae63bdbeccb4a6eeb7af41a1a`. Update this pin whenever the manuscript is resynchronized to a newer implementation revision.

## Scope

The paper describes the functionality used by the current application:

- native residue-aware RamplotR density regions for exploratory visualization;
- separate six-class Rama8000 standard validation (General, Gly, cis-Pro, trans-Pro, pre-Pro and Ile/Val);
- linked Ramachandran plot, true PDB residue numbering, sequence navigator, residue table and NGL 3D inspection;
- AlphaFold 2/ColabFold, AlphaFold 3 and ESMFold confidence handling;
- pairwise Conformational Change Explorer with wrapped delta-phi/delta-psi, alignment identity/coverage, linked dual-3D inspection, local hetero-residue context and paired pLDDT/ΔpLDDT when available;
- prediction ensembles and multi-model structural ensembles using circular backbone statistics, including AlphaFold 3 sample sets and residue-level disagreement returned to the shared inspector;
- group-level conformation comparison for biologically defined structure sets;
- sequence-based discovery of candidate experimental PDB counterparts for predicted models;
- optional official wwPDB validation XML as independent deposited-structure evidence;
- batch/reporting workflows and the Shinylive/webR browser deployment.

The manuscript treats these evidence layers separately. Native RamplotR `Not allowed`, Rama8000/wwPDB `Outlier`, prediction confidence, ensemble disagreement and external validation are not collapsed into a single quality score.

## Validation interpretation

The pinned five-structure validation set tests the numerical backbone and standard-validation implementation:

- all **3,718/3,718** eligible phi/psi pairs agree with independently reported wwPDB angles to within **0.055 degrees**;
- all **3,718/3,718** comparable residues receive exactly the same Rama8000 Favored/Allowed/Outlier category as the corresponding official wwPDB validation report;
- all seven official 2DQ4 Ramachandran outliers are independently recovered as Rama8000 Outlier.

The native RamplotR density regions remain an exploratory layer and are not remapped into standard outlier categories.

## Build

With a standard TeX installation:

```bash
cd paper
latexmk -pdf main.tex
```

The manuscript currently uses the real browser-test screenshot at `../docs/screenshots/overview.png`.

## Before submission

- Review the author list and contribution statement.
- Decide whether that revision should receive a publication release tag.
- Create an immutable software archive/DOI if desired.
- Add funding and acknowledgements where applicable.
- Consider expanding the independent validation set beyond the current five functional test structures.
- Consider adding a focused workflow figure for pairwise/prediction-ensemble analysis rather than additional UI screenshots.
