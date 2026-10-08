# RamplotR

**Interactive protein-structure inspection and Ramachandran analysis, from experimental structures to AlphaFold and ESMFold predictions.**

RamplotR is an open-source R Shiny application that brings backbone geometry, a 3D molecular viewer, residue-level validation and prediction confidence into one linked workspace. Click a residue in the plot, sequence, table or structure to inspect it everywhere. Compare conformations, examine multi-model ensembles and export figures or reproducible reports.

![RamplotR showing a Ramachandran plot, molecular structure and linked sequence navigator](docs/screenshots/overview.png)

*RamplotR's default publication palette. Example: PDB [1CRN](https://www.rcsb.org/structure/1CRN). [More screenshots](#screenshots).*

[**Try the browser app**](https://bikc.be/RamplotR/) · [**Get started**](#get-started) · [**Explore the features**](#what-you-can-do) · [**Scientific interpretation**](#scientific-interpretation) · [**Documentation**](#documentation)

## What you can do

| Workflow | Capabilities |
| --- | --- |
| **Explore a structure** | Interactive φ/ψ plots with several reference-density datasets, native RamplotR density regions and a parallel six-class Rama8000 standard validation that matches the current cctbx/Phenix categories. |
| **Inspect residues in context** | Synchronized Ramachandran plot, searchable residue table, all-chain sequence navigator and NGL 3D viewer. For PDB/AlphaFold DB inputs, available canonical UniProt coordinates are shown alongside the original structure numbering. The issue queue distinguishes Rama8000 outliers from native RamplotR `Not allowed` regions and missing angles. |
| **Work with predicted models** | Import AlphaFold DB models or upload AlphaFold 2/3, ColabFold and ESMFold structures. Examine pLDDT/PAE, analyse AF2/ColabFold, AF3 or ESMFold model/seed ensembles with residue-level backbone-state and confidence agreement, and search experimental PDB counterparts for direct comparison. |
| **Examine structural geometry** | Explore peptide ω, side-chain χ1 and descriptive Cβ measurements. Optionally attach the matching deposited structure's official wwPDB validation report for independent rotamer, clash and geometry annotations. |
| **Inspect experimental evidence** | Overlay a local CCP4/MRC cryo-EM map in the 3D viewer as a qualitative aid, without uploading the map to a separate service. |
| **Discover experimental structures** | Use the new **Atlas** tab to find experimental PDB polymer entities by UniProt accession, inspect experimental method, resolution and reported sequence coverage, then open a candidate in Compare. The Atlas currently provides an initial capped inventory; automatic conformational-state clustering is on the roadmap. |
| **Compare models** | Sequence-align two structures with the Conformational Change Explorer using the same Ramachandran density background as the main plot; inspect wrapped φ/ψ shifts and broad backbone-state transitions, swap primary/comparison roles instantly, or compare biologically defined structure sets (for example apo/holo or WT/mutant) using residue-level circular φ/ψ means, dispersion and between-group backbone shifts. |
| **Publish or automate** | Export SVG and high-resolution PNG figures, filtered CSV tables and standalone HTML reports. Run the offline R command-line tool on individual files or a directory of structures. |

Advanced analysis stays in collapsible panels or dedicated comparison/summary views, keeping the everyday 2D/3D inspection screen uncluttered.

## Screenshots

All images below are taken from [real browser tests](docs/screenshots/README.md). The prediction-confidence example uses **synthetic test data** to demonstrate the interface, not a biological result.

| Linked residue inspection | All-chain sequence navigator |
| :---: | :---: |
| <img src="docs/screenshots/residue-inspection.png" alt="A selected residue highlighted on the Ramachandran plot and in 3D" width="520"> | <img src="docs/screenshots/all-chains.png" alt="Expandable sequence tracks for four chains of 1BBB" width="420"> |
| **Prediction confidence and PAE** | **Multi-model ensemble analysis** |
| <img src="docs/screenshots/prediction-pae.png" alt="Synthetic AlphaFold-style PAE heatmap linked to the structure viewer" width="520"> | <img src="docs/screenshots/ensemble.png" alt="Circular angle and classification consistency analysis of NMR models" width="520"> |

## Get started

### Run the interactive app

Install [R](https://www.r-project.org/) (a current R 4.x release), clone the repository and start the Shiny application:

```bash
git clone https://github.com/BiKC/RamplotR.git
cd RamplotR
```

```r
install.packages(c(
  "shiny", "shinyWidgets", "colourpicker", "bio3d",
  "NGLVieweR", "DT", "jsonlite", "htmltools", "xml2"
))
shiny::runApp("shinyRam")
```

Enter a four-character **PDB ID** (for example, `1CRN`), **upload** your own PDB/mmCIF structure or choose **AlphaFold DB** to retrieve an available model by UniProt accession. Local uploads support `.pdb`, `.ent`, `.cif`, `.mcif` and `.mmcif`.

The project has also been hosted at [bioit.shinyapps.io/RamplotR](https://bioit.shinyapps.io/RamplotR/), but that deployment may not reflect the latest GitHub version. Running locally is the most reliable way to use the current implementation. Public-accession retrieval requires an internet connection; uploaded coordinates and local map files can be inspected without an external folding service.

### Use the browser version

[**Open RamplotR at bikc.be/RamplotR**](https://bikc.be/RamplotR/). The public version runs through **Shinylive** on static one.com hosting. R runs in your browser through webR, so no R installation is required. A first visit downloads and starts webR and its R packages; subsequent visits can reuse cached assets. Loading large structures still depends on the visitor's device.

The browser build uses the same five reference datasets and original RDS distributions as the desktop/server app. To avoid including all 110 distributions in the initial `app.json`, reference files are hosted separately and downloaded on first use. Every file is checked against the export's MD5 manifest, and loaded references are cached within the session. The full Plotly library starts downloading when you click **Analyse** instead of delaying the initial form. A small φ/ψ favicon matches the app's teal colour scheme.

For a reproducible deployment, run this from the repository root with the [shinylive R package](https://posit-dev.github.io/r-shinylive/) installed:

```bash
Rscript scripts/export-shinylive.R bikc.be https://bikc.be/RamplotR/reference-data
```

Upload the **contents** of the generated `bikc.be/` directory to the site's document root on one.com, including `RamplotR/reference-data/` and the shared `shinylive/` assets. The export also includes optional, directory-scoped `.htaccess` files for gzip/Brotli (when available) and cautious browser caching; these do not alter the website's root configuration. one.com restricts some Apache features, so check the HTTP response headers after deployment rather than assuming that compression is active.

Ordinary PDB/mmCIF analysis does not require `xml2`. Importing official wwPDB validation XML does require it; compatibility of that optional feature should be checked in the specific exported webR build. The local/server Shiny app and the batch command remain available when a browser package is unsupported.

See the [Shinylive export and one.com deployment guide](docs/shinylive-deployment.md) for hosting checks, caching rules and troubleshooting.

### Analyse many structures

From the repository root, process a structure or a directory containing supported structure files:

```bash
Rscript scripts/ramplotr-batch.R --input structures/ --output results/ --report
```

The offline command produces per-residue CSV, machine-readable JSON and a batch summary; `--report` additionally requests SVG and standalone HTML reports. For a compatible multi-model structure, add `--ensemble-models 20`. Use `--help` for all options, including declared prediction provenance and an optional matching wwPDB validation XML.

See the [batch-analysis instructions](docs/structural-verification.md#offline-batch-mode) for examples and resource limits.

## Scientific interpretation

RamplotR reports two deliberately separate views of backbone geometry. Its **native density regions** evaluate general, glycine, proline and pre-proline residues against the selected bundled RamplotR reference dataset and retain the labels Favoured, Allowed, Generously allowed and Not allowed. The selected plotting background is independent of those residue-specific calculations.

For standard structure validation, RamplotR also evaluates every residue with the current **Rama8000 six-class model** used by cctbx/Phenix: General, Gly, cis-Pro, trans-Pro, pre-Pro and Ile/Val, with Favored/Allowed/Outlier thresholds taken from the same reference score tables. This result is independent of the native RamplotR display background. In the pinned five-structure corpus, **all 3,718 comparable residues receive exactly the same Rama8000 category as the official wwPDB validation reports**, while the calculated φ/ψ angles agree to within 0.055°. See the [Rama8000 implementation and validation](docs/rama8000-validation.md).

Native RamplotR `Not allowed` and Rama8000/wwPDB `Outlier` therefore remain separate concepts in the interface and exports. For deposited experimental structures, an official wwPDB report for the **same structure and model** can additionally be attached as independent evidence.

For predictions, pLDDT and PAE describe model confidence rather than experimental verification. ESMFold normally provides pLDDT but not PAE; experimental thermal B-factors are **never** automatically interpreted as prediction confidence. The native ω, χ1 and Cβ measurements are descriptive, and visualising a density map is not a quantitative map–model fit measurement.

For comparison workflows, RamplotR also assigns finite φ/ψ pairs to one of four coarse backbone states (Alpha-R, Beta, PPII or Alpha-L) by wrapped angular distance to fixed canonical centres, with remote points labelled Other. These state labels make model/ensemble disagreement easier to scan. They are descriptive comparison bins, not secondary-structure assignments, validation categories or evidence of molecular motion.

The default RamplotR teal contour palette provides consistent, recognisable publication figures; changing colours never changes the scientific reference distribution or classification thresholds. Figure and report exports include relevant analysis settings and provenance.

## Documentation

- [Interactive inspection, colours, comparisons and exports](docs/inspection-user-guide.md)
- [AlphaFold, ColabFold and ESMFold confidence analysis](docs/prediction-confidence.md)
- [Geometry, official wwPDB evidence, cryo-EM overlays, ensembles and batch mode](docs/structural-verification.md)
- [Rama8000 standard validation and direct wwPDB comparison](docs/rama8000-validation.md)
- [Prediction ensemble analysis](docs/prediction-ensembles.md)
- [Structure-group conformational comparison](docs/group-conformation-comparison.md)
- [Experimental counterpart discovery for predicted models](docs/experimental-counterparts.md)
- [Experimental Atlas inventory](docs/atlas-inventory.md) · [Conformational Atlas roadmap](docs/conformational-atlas-roadmap.md)
- [Canonical UniProt residue coordinates](docs/canonical-coordinates.md)
- [Conformational Atlas roadmap](docs/conformational-atlas-roadmap.md)
- [Independent wwPDB angle-validation protocol](docs/wwpdb-validation.md) · [Results](docs/validation-results.md)
- [Performance and large-structure benchmarks](docs/benchmark-results.md) · [Scaling results](docs/scaling-results.md)
- [Shinylive browser deployment and one.com caching](docs/shinylive-deployment.md)

Developers can run the pure-R scientific regression suite from the repository root with `Rscript tests/scientific.R`. Additional tests cover structure parsing, validation imports, confidence formats, geometry, ensembles and the batch CLI. GitHub Actions also exercises the application in a real browser and runs the scientific tests on Ubuntu and Windows.

The [`v0.1.0-legacy` tag](https://github.com/BiKC/RamplotR/tree/v0.1.0-legacy) preserves the historical version for reproducing earlier analyses. Record the exact RamplotR revision and reference dataset when publishing results.

## Research use and citation

RamplotR has been used to visualise protein models in the study *[Unveiling Intra-Clonal Diversity of Monkeypox Virus from Brazil's First Outbreak Wave](https://doi.org/10.3390/v18010062)* (Witt et al., *Viruses*, 2026), including supplementary Ramachandran plots of viral polymerase and helicase models. This is an example of its research use, not an independent validation of the software.

When using RamplotR in a publication, cite the software repository and the **exact version or commit** you used, and report the selected reference dataset and analysis mode.

## License

[MIT](LICENSE). The bundled reference datasets and external validation reports retain their respective scientific provenance; consult the linked documentation when reusing those data.
