# RamplotR inspection and publication workflow

The everyday workflow is deliberately short: load a PDB/mmCIF structure, inspect its Ramachandran plot and molecular viewer, select a residue, and review its details. A compact multi-chain sequence navigator is integrated below the 2D/3D views. Additional workflows use **Residue list**, **Compare** and **Summary** tabs. After loading a structure, the introductory form becomes a slim toolbar; use **Hide settings** for more plotting space on a laptop.

## Selection remains visible across views

The inspector beneath the structure input is shared across tabs. Selecting a residue from the Ramachandran plot, DataTable, sequence, or 3D structure updates the inspector, highlights the corresponding point, shows the residue as orange sticks and focuses the NGL camera.

**Show in plot** brings you to the plot and molecular viewer without clearing the selection. **Clear** restores the overview. **Previous issue** and **Next issue** navigate through residues prioritised as Rama8000 outliers, native RamplotR `Not allowed` regions, unavailable/terminal backbone angles, and positions within two percentile points of a density cutoff. “Near boundary” is a visual review hint; it is not an additional scientific quality classification.

On a loaded screen with no selected residue, the inspector is reduced to a single review action. After selecting a residue, navigation controls become available. Representation buttons and the ligand, DNA, RNA, spin and rock switches remain visible beneath the viewer, so users can discover them without opening a settings menu.

## Why inspect this residue?

Selecting a residue opens a compact evidence explanation in the shared
inspector. RamplotR keeps the underlying signals separate instead of combining
them into an opaque quality score.

Depending on the data available for that residue, the explanation can surface:

- Rama8000 Allowed or Outlier status;
- the native RamplotR `Not allowed` region or proximity to a native contour;
- low or very low pLDDT;
- the combination of very high pLDDT with a Rama8000 outlier;
- cis or twisted peptide geometry;
- an attached official wwPDB Ramachandran or rotamer outlier;
- official local clashes or covalent-geometry outliers;
- nearby non-water hetero residues such as ligands, cofactors or ions, reported
  with their nearest heavy-atom distance when they fall within 6 Å.

Local hetero context is intentionally descriptive. RamplotR reports spatial
proximity but does **not** infer a binding interaction from distance alone.
Water/solvent records and hydrogen atoms are excluded from this context. For
multi-model coordinate files, hetero context is currently reported only for
model 1 so ligand coordinates are never silently reused for another model.

Each item identifies its source and explains why it may be worth examining.
If none of the available evidence is unusual, the panel says so explicitly
rather than inventing a warning.

After a **prediction ensemble** has been analysed, the inspector also adds
model-to-model evidence for the selected residue. It can flag strong circular
phi/psi spread, Rama8000 category disagreement, variable pLDDT, or incomplete
residue coverage across models. A useful discordant case is high mean pLDDT
combined with substantial backbone spread: the predictions are individually
confident but do not converge on one local backbone conformation. The ensemble
table and inspector also report a coarse backbone-state mode and its agreement,
so a residue can be recognized immediately when models move between broad
Alpha-R, Beta, PPII, Alpha-L or Other regions. These are comparison bins, not
DSSP secondary-structure assignments. Ensemble
spread remains prediction uncertainty/heterogeneity and is never presented as
experimental molecular motion.

This is an **inspection aid**, not a residue-quality score. A Rama8000
classification, prediction confidence, local geometry and an official wwPDB
annotation remain scientifically distinct observations.

## Compact multi-chain sequence map

The overview directly under the Ramachandran plot shows every selected
protein chain at the same time. Its small, position-based bars retain the
original residue order; for long chains, bins preserve isolated outliers
and missing-angle positions. Expand the navigator to reveal independently
scrollable, one-letter residue strips for all selected chains. Clicking a
letter updates the same inspector, Ramachandran point and NGL focus.

The map always retains true **PDB residue numbers**, even when amino-acid
or pre-proline filters hide most plotted points. The expanded strip labels
every tenth PDB position (plus its first and last residue) and shows the
currently selected number next to the chain heading. Use **Go to residue**
beside any chain, or press Enter in its number field, to find positions
directly, such as residue 104. A hidden residue is still located, with a
message explaining that its selection is blocked by current plot filters.

For every structure, the expanded navigator keeps the native RamplotR density region as the letter background and shows the independent **Rama8000 standard-validation category** as a small corner marker (Favored, Allowed or Outlier). For predicted models, each residue button additionally shows its numerical **pLDDT** below the amino-acid letter. A separate coloured underline and legend distinguish high, confident, low and very low pLDDT from both geometry classifications. Positions with unavailable confidence show a dash, never an invented zero.

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

The optional **Compare** tab accepts a second PDB accession or PDB/mmCIF file. The comparison-source setup collapses automatically after a successful load, leaving the active structure name and provenance visible in its header; reopen it only when replacing that structure. The **Swap primary ↔ comparison** action remains beside the chain-role controls so role reversal stays available during analysis. Choose one chain from each structure. For uploaded comparison models, declare whether the coordinates are experimental/unknown, AlphaFold 2/ColabFold, AlphaFold 3, ESMFold or another prediction with pLDDT in the B-factor field. AlphaFold 3 comparison models require their matching full-confidence JSON; experimental B-factors are never interpreted as pLDDT. When confidence is available, the aligned table and selected-pair inspector show pLDDT for each side and signed ΔpLDDT, with a direct filter for |ΔpLDDT| ≥20.

The comparison summary groups evidence into **Alignment**, **Backbone**, **Validation** and, when available, **Prediction confidence**, rather than mixing all metrics into one strip. Residues are paired by a bounded global amino-acid sequence alignment, **not by residue number**. The comparison summary reports sequence identity plus coverage of both selected chains, so local conformational differences can be interpreted in the context of alignment quality. A caution appears when identity falls below 50% or either chain has less than 70% aligned coverage. These are interpretation guards, not statistical significance thresholds. The difference in each angle wraps correctly across ±180°. Gaps remain visible and do not receive invented dihedrals; classification and confidence differences are reported only where the corresponding evidence exists.

The paired Ramachandran plot and the **3D superposition are side by side**
on wide screens. The aligned-angle plot uses the **same selected Ramachandran
density background, contour thresholds and palette as the individual plot**,
so changing reference datasets or display palettes stays consistent across
both views. The 3D viewer initially fits both selected chains rather than the
complete uploaded structures, even if hidden chains are very large.
Selecting a plotted point or aligned table row highlights both corresponding
residues in the superposition and focuses the camera on their local
positions. Clicking a residue in either 3D structure finds its aligned
partner. If one structure contains an insertion/deletion at the selected
position, only the available residue is highlighted and the missing
partner is shown as an alignment gap.

Above the paired 2D/3D workspace, the **Conformational change explorer** represents each aligned residue as one compact cell. Its colour ranks the combined wrapped backbone displacement, defined as `sqrt(Δφ² + Δψ²)` after wrapping each angle across ±180°. The bands (<15°, 15–30°, 30–60° and ≥60°) are navigation aids, not statistical significance thresholds or Cartesian distances. Clicking a cell or one of the five largest-shift shortcuts selects that aligned pair everywhere, including both 3D structures.

Use **Find aligned pair** to jump by the true residue number in either
selected chain (for example 104), then inspect the primary/comparison
amino acids, φ/ψ angles, wrapped Δφ/Δψ, combined backbone displacement,
coarse backbone-state transition and Rama8000 categories directly below the
views. The comparison table can be filtered directly to residues that change
broad backbone state.

The same selected-pair card also reports nearby non-water hetero residues
for **both** structures when model-appropriate coordinates are available.
This is useful for apo/holo, cofactor-bound or ion-associated comparisons:
for example, a local backbone shift can be inspected alongside an ATP or
metal ion present on only one side. Distances are minimum heavy-atom
distances within 6 Å and are explicitly structural proximity, not evidence
of biochemical binding.

Use **Swap primary ↔ comparison** when the second structure should become the
coral primary reference for the comparison. This reverses the A/B chain
controls, angle traces, superposition roles and paired-residue navigation
without reloading either structure. The structure loaded in the main RamplotR
workspace is intentionally not replaced, so selections can still synchronize
with its sequence navigator and shared residue inspector.

**Fit both chains** resets the camera without clearing the current
selection. The shared primary-structure inspector and sequence navigator
follow the selected primary residue when it is visible under current plot
filters. Export the aligned table as CSV for reproducibility. A large-chain
comparison can exceed the alignment size limit; select shorter chains instead.

Identical Ramachandran coordinates do not imply identical Cartesian structure, and an angular difference alone is not evidence of a clinically meaningful change. Distinct models from the same structure are related observations, not independent experiments.

## Prediction ensembles

When the loaded structure is a prediction, the **Summary** tab exposes a
prediction-ensemble workflow even if the current coordinate file contains only
one model. Upload additional AF2/ColabFold, ESMFold or compatible
pLDDT-in-B-factor models; the currently loaded compatible model can be included
as one ensemble member. For AlphaFold 3, upload at least two official sample
`*_model.cif` files together with their matching `*_confidences.json`
files and optional `*_summary_confidences.json` files. AF3 pairing uses the
seed/sample filename stem rather than upload order.

RamplotR reports circular phi/psi spread, residue coverage, coarse backbone-state
agreement, Rama8000 agreement and pLDDT spread across the uploaded models. A compact **Prediction variability
map** ranks residues by the larger circular SD of phi or psi. Clicking a map
cell or ensemble-table row selects the corresponding residue in the existing
inspector and linked 2D/3D views.

The variability map is a navigation tool. Model-to-model disagreement is
described as prediction uncertainty or heterogeneity, not molecular dynamics.
Duplicate coordinate files are rejected, and the model export records labels,
declared source and coordinate MD5 hashes.

For AlphaFold 3, atom-level pLDDT is mapped from each sample's matching full
confidence sidecar. pTM, ipTM, ranking score, disordered fraction and clash flag
are retained separately in the model-level table and exports; they do not
change the residue variability map. Missing or ambiguous AF3 file pairs are
rejected instead of being guessed.

## Publication exports and reproducibility

The Summary tab includes an export disclosure. **Vector SVG** produces an editable figure; **High-resolution PNG** is suitable for manuscripts; **HTML report** includes the plot, summary counts, selected residue data and analysis settings. The residue CSV exports the current table filters; the comparison CSV exports the aligned comparison.

Record the software version, reference dataset, background, classification mode, structural model and residue filters when using RamplotR in a publication. The original paper-era implementation is preserved in the `v0.1.0-legacy` Git tag, which remains independent of the new development branch.

## Regression coverage

- `Rscript tests/inspection.R`: review ordering, sequence identities, insertion/deletion-aware alignment, angular wraparound, conformational-displacement ranking and model coordinates.
- `Rscript tests/reports.R`: export real vector SVG, PNG and self-contained HTML.
- `Rscript tests/model-integration.R`: an actual multi-model 1D3Z NMR PDB.
- `node tests/ui.test.cjs`: JS message-handler and residue-selection contracts.
- `node tests/compare-ui.test.cjs`: aligned 3D selection, chain-aware framing, invisible viewers and fit-both reset.
- `node tests/ui-browser.cjs`: live Shiny browser workflow, plot-to-NGL/table/sequence selection, default palette, responsive sizing and pairwise self-comparison.

Scientific tests run on Ubuntu and Windows; full-browser checks run on Ubuntu. For changes that affect validation criteria or reference grids, the scientific regression tests must remain unchanged or include explicitly reviewed new reference fixtures.
