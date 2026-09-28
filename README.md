# RamplotR

RamplotR is an R Shiny application for interactive Ramachandran plots of protein structures. It includes residue and chain filters, multiple reference-density datasets, regional statistics, an optional residue-aware classification and a 3D molecular viewer.

## Historical published version

The tagged version [v0.1.0-legacy](https://github.com/BiKC/RamplotR/releases/tag/v0.1.0-legacy) preserves the implementation available before the 2026 scientific corrections. Use this tag to reproduce analysis performed with the earlier application. New scientific behaviour and reference classifications may differ and should be recorded in subsequent analyses.

## Usage

The original app has been hosted at https://bioit.shinyapps.io/RamplotR/. The hosted deployment may not match the current repository.

To run this repository locally, install a current R 4.x and the following CRAN packages:

```r
install.packages(c("shiny", "shinyWidgets", "colourpicker", "bio3d", "NGLVieweR", "DT"))
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

For screenshots, limitations and the complete workflow see the [inspection and publication guide](docs/inspection-user-guide.md). The interface layout and browser verification history are documented in [interface-refresh.md](docs/interface-refresh.md).

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
