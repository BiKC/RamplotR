# Phase C: structural verification, ensembles and batch analysis

Phase C extends the established RamplotR Ramachandran analysis without
changing its reference densities, four-region labels, or the historical
`v0.1.0-legacy` release.

## Extended native geometry

Once a structure loads, expand **Extended structure verification** below the
main plots. RamplotR calculates the peptide dihedral **omega** (CA–C–N–CA)
only across geometrically connected peptide bonds, and **chi1** (N–CA–CB–X1)
for residues with an appropriate first side-chain atom. Cβ measurements
include the observed **CA–CB distance** (Å) and signed N–CA–C–CB tetrahedral
volume (Å³) for residues with all four atoms. These are descriptive
measurements, not independently validated Cβ-deviation or chirality-outlier
classifications. All angles are in degrees. Missing atoms, chain breaks and terminal residues have undefined
measurements; they are not scored as outliers.

For exploration, omega within 30° of 0° is labelled *cis*, omega within 30°
of ±180° is labelled *trans*, and other measured values are marked *twisted*.
These are descriptive flags, **not an independent MolProbity assessment**.
Chi1 is a measurement, not a rotamer outlier prediction.

## Independent experimental validation

For a deposited structure, obtain the official wwPDB validation XML or
`XML.gz` for the **same deposited entry and structural model**. Attach it in
the expandable verification panel and confirm that it belongs to the
loaded experimental structure. The app rejects official reports for declared
predicted-model inputs, including AlphaFold DB. For local or accession-loaded
experimental structures, you must still verify accession/model provenance:
matching sequence numbering alone cannot establish that two deposits are the
same experiment. The report is parsed locally and matched by **model,
chain, residue number, insertion code and residue type**, with unlabelled
alternate conformations preferred over alternate A.

When available, the inspector and detailed CSV show independent wwPDB
Ramachandran, side-chain rotamer, local clash, symmetry clash, bond-length
outlier, bond-angle outlier, RSCC and RSRZ annotations. A missing independent
record remains missing. The app reports exact matching coverage and preserves
the source file's MD5 checksum in its HTML report.

The original RamplotR contour interpretation and the wwPDB/MolProbity
reference systems are **not equivalent**. The independent outlier-rich
[Phase A results](phase-a-results.md) explicitly show genuine discrepancies.
Importing a report does not rewrite the original region classification.
Do not attach experimental wwPDB validation reports to an unrelated
AlphaFold/ESMFold prediction, even if their sequences are similar.

Official validation information:
https://www.wwpdb.org/validation/validation-reports

## Optional cryo-EM map overlay

Open **Local cryo-EM density map** beneath the NGL viewer, choose your own
CCP4/MRC file, and select *Show map*. The map is read by NGL directly in the
browser, with no external map-fitting service. It is rendered as a translucent
teal isosurface and can be adjusted from 0.5 to 5 sigma. *Remove* clears
the map; loading a different structure also clears any previous overlay.
Only files up to 64 MB are supported to protect ordinary laptop sessions.
NGL can be slow with large maps; the overlay is a qualitative visual aid,
**not an RSCC/Q-score or local map-model-fit measurement**. Use the
corresponding official report or a validated external map-fit tool when
quantitative claims are required.

NGL stage and volume API:
https://nglviewer.org/ngl/api/class/src/stage/stage.js~Stage.html

## NMR / multi-model ensemble

For an input containing multiple atom-compatible models, expand **Ensemble
analysis** on the Summary tab, then select *Analyse ensemble*. It analyses
up to the first 30 models on demand, using the currently chosen reference
dataset, classification mode and plotting background. The tool matches
residues by chain, residue number, insertion code and amino-acid identity,
not by row number. It reports circular means and standard deviations of
phi/psi angles, observed-model counts and the proportion of models agreeing
on a class. Single-angle observations have undefined variability.

A changed reference invalidates the previous analysis until recalculated.
Clicking an ensemble residue selects it in the shared inspector, plot and 3D
viewer. Use **Export ensemble CSV** to keep the per-residue results. A
prediction ensemble is an ensemble of output conformations, not proof of
experimental flexibility or pLDDT uncertainty.

## Offline batch mode

From the repository root, with R and Bio3D installed, run:

```bash
Rscript scripts/ramplotr-batch.R --input structures/ --output results/ \
  --reference original --mode residue --model 1 --report
```

The command accepts one local PDB/mmCIF file or a nonrecursive directory of
up to 1000 supported structure files. By default it produces residue-level
CSV and JSON, a batch-summary CSV, and optionally a standalone SVG plus HTML
report. Each file's output is named from its source filename, and existing
files are protected unless `--overwrite` is specified. The full CSV
includes omega/chi1 and available prediction-confidence or independent
validation annotations.

For a consistent NMR structure, request an ensemble table:

```bash
Rscript scripts/ramplotr-batch.R --input 1D3Z.pdb --output results/ \
  --model 1 --ensemble-models 20 --report
```

For **one** experimental structure, add its official XML:

```bash
Rscript scripts/ramplotr-batch.R --input 1CRN.pdb --output results/ \
  --validation-xml 1crn_validation.xml.gz --report
```

For an ESMFold prediction, explicitly set
`--prediction-source esmfold`; never request this for an experimental
structure. For AlphaFold 2/3, the CLI can read declared pLDDT from compatible
B-factor files; AF3 confidence sidecar parsing remains available in the
interactive Phase B uploader.

Use `--no-json` to avoid the optional jsonlite dependency; `--report`
requires htmltools and an SVG-capable R graphics device, while wwPDB XML
requires xml2. The CLI never downloads structures or sends coordinates to
external servers. Nonzero exit status indicates at least one failed input,
and `batch-summary.csv` records any per-file errors.

## Reproducibility and methodological limits

The outputs record the current model, file MD5 checksum, reference file MD5,
R/Bio3D versions, chosen scientific mode and prediction provenance where
declared. The HTML report distinguishes computed geometry from imported
independent evidence and prints model-ensemble statistics when available.
Use the pinned Phase A corpus and reference-validation reports before
making formal validation-performance claims. Avoid presenting any visual
density overlay or geometric heuristic as an official experimental
fit or MolProbity-equivalent score.
