# RamplotR

RamplotR is an R Shiny application for interactive Ramachandran plots of protein structures. It includes residue and chain filters, multiple reference-density datasets, regional statistics, an optional residue-aware classification and a 3D molecular viewer.

## Historical published version

The tagged version [v0.1.0-legacy](https://github.com/BiKC/RamplotR/releases/tag/v0.1.0-legacy) preserves the implementation available before the 2026 scientific corrections. Use this tag to reproduce analysis performed with the earlier application. New scientific behaviour and reference classifications may differ and should be recorded in subsequent analyses.

## Usage

The original app has been hosted at https://bioit.shinyapps.io/RamplotR/. The hosted deployment may not match the current repository.

To run this repository locally, install a current R 4.x and the following CRAN packages:

```r
install.packages(c("shiny", "shinyWidgets", "colourpicker", "bio3d", "NGLVieweR", "DT", "jsonlite"))
shiny::runApp("shinyRam")
```

The input panel lets you select a four-character PDB identifier or upload a local PDB or mmCIF structure (.pdb, .ent, .cif, .mcif, .mmcif). Uploads retain insertion codes and atom alternate locations for backbone processing. Structures with consistent multiple models offer a model selector. Both PDB and mmCIF uploads are supported; see the [inspection guide](docs/inspection-user-guide.md) for the current model and comparison limitations.

## Scientific interpretation

Reference density grids are in `shinyRam/static/`. They include the original distributions derived from the protein selection discussed by [Lovell et al. (2003)](https://pubmed.ncbi.nlm.nih.gov/12557186/), plus additional datasets. The chosen background controls the visual plot; residue-aware classification uses corresponding General, GLY, PRO and preProline grids from the selected dataset. The explicitly labelled legacy mode reproduces classification against the displayed background.

Reference distributions and density-percentile thresholds in RamplotR must not be described as equivalent to MolProbity quality metrics without independent validation. Save the selected reference dataset, scientific mode, threshold settings and application version alongside published results.

## Interface

The primary plot and 3D viewer share a selected-residue inspector, which remains visible when you switch tabs. A compact sequence overview directly beneath these views shows **every selected chain at once**; expand it to browse individual residues without leaving the plot. The loaded interface also offers a condensed laptop layout and a reversible focus mode for hiding analysis settings. Click a point, row, sequence letter or residue in the molecular structure to show the corresponding residue across all views. The inspector includes **Show in plot**, **Clear**, and outlier-review navigation. Changing residue filters, density reference, classification mode or palette updates the result without refetching the structure.

The residue table has readable angles, combined scientific/review filters and CSV export. The 3D viewer supports cartoon, ribbon, sticks, ball-and-stick and surface representations; ligand, DNA, RNA, spin and rock controls remain visible beneath the viewer as modern switches. RamplotR's own ordered publication palette is selected by default, with legacy colour schemes and custom colours still available.

The optional comparison tab aligns a chain from each of two structures by sequence and displays angular and classification differences, including insertions and deletions. Multi-model structures can be inspected model by model. The Summary tab exports vector SVG and 300-dpi PNG plots plus a self-contained report with reproducibility settings.

### AlphaFold and ESMFold confidence (Phase B)

AlphaFold DB accession lookup and explicit AlphaFold 2/3, ColabFold and ESMFold upload provenance can add pLDDT tracks beneath every chain. AlphaFold PAE and AF3 confidence JSON are optional, strictly matched to the selected coordinates, and shown as a linked, collapsible heatmap. ESMFold local PDBs expose pLDDT from their B-factor fields; standard ESMFold does not provide PAE. Experimental B-factors are never interpreted as prediction confidence. See the [prediction guide](docs/phase-b-guide.md).

For screenshots, limitations and the complete workflow see the [inspection and publication guide](docs/inspection-user-guide.md). The interface layout and browser verification history are documented in [interface-refresh.md](docs/interface-refresh.md).

## Independent Phase A validation

The [independent validation protocol](docs/wwpdb-validation.md) uses official
wwPDB residue-level validation reports to compare RamplotR's angles and
reference-dependent classifications. A [versioned public corpus manifest](validation/manifest.csv)
covers five experimental structures: 1CRN, 1UBQ, 6VXX, 2DQ4 and 1D3Z. The companion
[GitHub Actions workflow](.github/workflows/wwpdb-validation.yml) archives
source XML/mmCIF files, SHA256 hashes, exact software versions, per-residue
results and reference-group contingency matrices. Differences in categories
are explicitly reported instead of falsely claiming the four-group RamplotR
method and six-group MolProbity are identical. The [first five-structure
results](docs/phase-a-results.md), [pinned numerical results](validation/baseline-2026-09-28.csv)
and [pinned official source SHA256 hashes](validation/baseline-source-hashes-2026-09-28.csv)
are available for independent review. In particular, the outlier-rich 2DQ4
case exposes important differences between RamplotR and wwPDB outlier calls.

AlphaFold and ESMFold predicted-structure ingestion and linked confidence
assessment are planned for Phase B; they are deliberately kept out of the
independent experimental-validation benchmark.

## Regression tests

Pure-R tests do not require Shiny or Bio3D. From the repository root:

```sh
Rscript tests/scientific.R
Rscript tests/peptide.R
Rscript tests/input.R
Rscript tests/classification.R
Rscript tests/inspection.R
```

GitHub Actions runs source parsing, scientific regression tests and inspection/alignment tests on Ubuntu and Windows. Its full-browser job runs the Shiny app with real PDB input, checks table alignment and linked views, and exercises figure/report export and an actual multi-model NMR fixture. The separate structure validation workflow checks real structures and reference classifications.

## Project files

- `shinyRam/app.R`: user interface and Shiny server
- `shinyRam/R/`: structure loading, backbone analysis and classifications
- `shinyRam/static/`: reference-density matrices
- `shinyRam/www/`: plotting JavaScript and bundled viewer assets
- `tests/`: isolated scientific and input tests
- `docs/upgrade-worklog.md`: original scientific upgrade plan and verification status
- `docs/upgrade-inspection-worklog.md`: current inspection-upgrade checklist
- `docs/inspection-user-guide.md`: user workflow, colour scheme, export and comparison guidance

MIT license; see LICENSE.
