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

The input panel lets you select a four-character PDB identifier or upload a local PDB or mmCIF structure (.pdb, .ent, .cif, .mcif, .mmcif). Uploads retain insertion codes and atom alternate locations for backbone processing. For now the app reads the first structural model.

## Scientific interpretation

Reference density grids are in `shinyRam/static/`. They include the original distributions derived from the protein selection discussed by [Lovell et al. (2003)](https://pubmed.ncbi.nlm.nih.gov/12557186/), plus additional datasets. The chosen background controls the visual plot; residue-aware classification uses corresponding General, GLY, PRO and preProline grids from the selected dataset. The explicitly labelled legacy mode reproduces classification against the displayed background.

Reference distributions and density-percentile thresholds in RamplotR must not be described as equivalent to MolProbity quality metrics without independent validation. Save the selected reference dataset, scientific mode, threshold settings and application version alongside published results.

## Interface

Once a structure is loaded, changing filters, the reference dataset, the
classification mode or plot/chain colours updates the plot and statistics
automatically. Click any residue on the Ramachandran plot, in the 3D structure,
or in the searchable residue list to highlight it across all three views.
Use **Clear selection** to dismiss the highlight. Press **Analyze structure**
only when loading another PDB accession or uploaded structure.


The current version has a responsive scientific workspace with the Ramachandran plot and 3D viewer alongside one another on wide screens. Reference controls and residue filters are grouped in the sidebar. The full UI changes and manual browser checks are documented in [interface-refresh.md](docs/interface-refresh.md).

## Regression tests

Pure-R tests do not require Shiny or Bio3D. From the repository root:

```sh
Rscript tests/scientific.R
Rscript tests/peptide.R
Rscript tests/input.R
Rscript tests/classification.R
```

GitHub Actions runs these tests and checks R source syntax on Linux and Windows. Integration checks using actual PDB/mmCIF structures and an independent structural validator are still planned.

## Project files

- `shinyRam/app.R`: user interface and Shiny server
- `shinyRam/R/`: structure loading, backbone analysis and classifications
- `shinyRam/static/`: reference-density matrices
- `shinyRam/www/`: plotting JavaScript and bundled viewer assets
- `tests/`: isolated scientific and input tests
- `docs/upgrade-worklog.md`: upgrade plan and verification status

MIT license; see LICENSE.
