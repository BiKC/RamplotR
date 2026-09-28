# RamplotR inspection and publication upgrade

Branch: `upgrade/inspection-visualization-reporting` (integration PR #12).
Historical reference: `v0.1.0-legacy` is unchanged; all scientific reference grids are preserved.

## Completed work

- [x] Repair DataTables header/body alignment after switching tabs, format on-screen angles and percentiles consistently, add region/review filtering and retain unrounded filtered CSV export.
- [x] Add a shared cross-tab residue inspector, linked plot/table/NGL selection, orange sticks and automatic camera focus, plus outlier and missing-angle review navigation.
- [x] Integrate a compact, collapsible all-chain sequence overview directly below the Ramachandran plot. Overview bins retain isolated outliers and missing angles; the expanded navigator shows each chain in its own horizontal strip with clickable residue letters. Preserve whole-chain sequence positions even when plot filters hide residues.
- [x] Add a distinctive RamplotR default plotting palette while keeping Rampage, PDBSum and custom options.
- [x] Keep NGL's cartoon, ribbon, sticks, ball-and-stick and surface representation buttons and modern ligand, DNA, RNA, spin and rock switches visible underneath the molecule; keep selected-residue highlighting independent from whole-structure representation.
- [x] Add a compact loaded-state layout for standard laptop screens and a reversible Hide settings focus mode. Reframe NGL automatically when a new structure replaces an earlier zoomed structure.
- [x] Export high-resolution PNG, editable SVG, unrounded CSV and self-contained HTML analysis reports including analysis provenance.
- [x] Support appropriate consistent multi-model coordinate sets with a model selector and regression checks, plus sequence-aligned two-structure comparison with gap-aware angle and region comparisons and optional 3D superposition.
- [x] Add isolated regression tests and live-browser tests for residue selection, 3D style switches, responsive layouts, four-chain 1BBB navigation, NGL framing on structure changes and structure/report exports.
- [x] Verify the integration commit `38bebba9be8808b8f20c3cfc6f2162d5a3438b9a` with successful scientific regression tests on Ubuntu and Windows, structure validation and the full live-browser workflow.

## Verification

- [Scientific regression tests](https://github.com/BiKC/RamplotR/actions/runs/36470803767)
- [Structure validation and benchmark](https://github.com/BiKC/RamplotR/actions/runs/36470803613)
- [Live Shiny browser test and visual-preview artifact](https://github.com/BiKC/RamplotR/actions/runs/36470803548)

These runs confirm the implemented regression and browser scenarios. Independent agreement of every RamplotR region classification against MolProbity is separate future validation, not a claim from these tests. Hosted deployments must be updated separately from merging the GitHub repository.

## Scientific boundaries

Compare residues by sequence alignment, not numerical index; missing torsions remain missing. The RamplotR palette is a figure convention, not a new reference distribution or classification method.
