# RamplotR inspection and publication workflow

The everyday workflow is deliberately short: load a PDB/mmCIF structure, inspect its Ramachandran plot and molecular viewer, select a residue, and review its details. Other workflows are in dedicated **Residue list**, **Sequence**, **Compare** and **Summary** tabs.

## Selection remains visible across views

The inspector beneath the structure input is shared across tabs. Selecting a residue from the Ramachandran plot, DataTable, sequence, or 3D structure updates the inspector, highlights the corresponding point, shows the residue as orange sticks and focuses the NGL camera.

**Show in plot** brings you to the plot and molecular viewer without clearing the selection. **Clear** restores the overview. **Previous issue** and **Next issue** navigate through residues prioritised as outliers, unavailable/terminal backbone angles, and positions within two percentile points of a density cutoff. “Near boundary” is a visual review hint; it is not an additional scientific quality classification.

On a loaded screen with no selected residue, the inspector is reduced to a single review action. After selecting a residue, navigation controls become available. Additional NGL overlays and animation options live in a disclosure below the viewer.

## The residue table

The table uses a single HTML table for headings and rows, without separate DataTables horizontal-scrolling header clones. Columns and angles are readable; displayed angles and density-percentile values are rounded to one decimal place. Exported CSV retains the unrounded values. **Region** and **Review** filters can be combined; selecting a row updates the shared inspector and both visualisations.

## RamplotR colours

The signature default palette, ordered from outside a reference contour to the highest-density region:

| Classification region | Colour |
| --- | --- |
| Not allowed | #FFF8ED |
| Generously allowed | #D4ECE7 |
| Allowed | #7DB9B5 |
| Favoured | #126E74 |

The same palette is used in PNG and SVG exports. This is a recognisable plotting convention **only**: colours do not change reference-density grids or region thresholds. Rampage, PDBsum and custom colours remain available.

## Molecular representations

The NGL style selector offers **Cartoon**, **Ribbon**, **Sticks** (licorice), **Ball & stick**, and **Surface**. Changing styles affects the overall protein representation, not the separate orange residue highlight. Surface calculation can take longer on large structures. **Layers & motion** contains accessible switches for ligands, DNA, RNA, spin and rock. Selecting a residue automatically stops motion so the zoom is stable.

## Structural models

Where Bio3D can retain consistent coordinates across models, a **Structural model** control appears in Reference & validation. It lets you choose the current model while retaining the original structure. Backbone angles are recalculated for the selected coordinates; other reference and presentation settings remain reactive.

If atom records are inconsistent or the parser cannot extract complete multi-model coordinates, the application falls back to the model that was loaded. The analytical model selection and NGL model selector must agree; model selection is not an ensemble quality statistic.

## Pairwise structure comparison

The optional **Compare** tab accepts a second PDB accession or PDB/mmCIF file. Choose one chain from each structure. Residues are paired by a bounded global amino-acid sequence alignment, **not by residue number**. The difference in each angle wraps correctly across ±180°. Gaps remain visible and do not receive invented dihedrals; classification differences are reported only for available classifications.

The paired Ramachandran plot displays the two structures in contrasting colours. Export the aligned table as CSV. The optional 3D viewer superposes the selected chains for qualitative inspection. A large-chain comparison can exceed the alignment size limit; select shorter chains instead.

Identical Ramachandran coordinates do not imply identical Cartesian structure, and an angular difference alone is not evidence of a clinically meaningful change. Distinct models from the same structure are related observations, not independent experiments.

## Publication exports and reproducibility

The Summary tab includes an export disclosure. **Vector SVG** produces an editable figure; **High-resolution PNG** is suitable for manuscripts; **HTML report** includes the plot, summary counts, selected residue data and analysis settings. The residue CSV exports the current table filters; the comparison CSV exports the aligned comparison.

Record the software version, reference dataset, background, classification mode, structural model and residue filters when using RamplotR in a publication. The original paper-era implementation is preserved in the `v0.1.0-legacy` Git tag, which remains independent of the new development branch.

## Regression coverage

- `Rscript tests/inspection.R`: review ordering, sequence identities, insertion/deletion-aware alignment, angular wraparound and model coordinates.
- `Rscript tests/reports.R`: export real vector SVG, PNG and self-contained HTML.
- `Rscript tests/model-integration.R`: an actual multi-model 1D3Z NMR PDB.
- `node tests/ui.test.cjs`: JS message-handler and residue-selection contracts.
- `node tests/ui-browser.cjs`: live Shiny browser workflow, plot-to-NGL/table/sequence selection, default palette, responsive sizing and pairwise self-comparison.

Scientific tests run on Ubuntu and Windows; full-browser checks run on Ubuntu. For changes that affect validation criteria or reference grids, the scientific regression tests must remain unchanged or include explicitly reviewed new reference fixtures.
