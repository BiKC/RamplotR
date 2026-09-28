# RamplotR inspection and publication upgrade

Branch: upgrade/inspection-visualization-reporting. Preserve v0.1.0-legacy and all scientific reference grids.

## Work packages
- [ ] Fix DT header/body alignment after a hidden tab becomes visible, improve columns and filtered CSV exports.
- [x] Shared inspector, review navigation and a collapsible all-chain sequence navigator beneath the main plot/3D viewer. The overview uses bounded position-based mini-maps, preserving isolated outliers and missing angles when compressed. Expanded chains are independently scrollable and share residue selection with every other view.
- [ ] Publication-quality RamplotR default palette, preserving Rampage, PDBSum and custom settings.
- [ ] Accessible modern NGL viewer switches and cartoon, sticks, ball+stick and surface representations.
- [ ] PNG/SVG and metadata-rich report export.
- [ ] Multi-model support, explicitly retaining model provenance and validating coordinates.
- [ ] Two-structure comparison by sequence alignment, angular differences and region classification changes.
- [ ] Unit and browser tests, scientific regression checks, documentation and integration to main.

## Scientific boundaries
Do not compare residue indices without sequence alignment or silently compare missing torsions. Classifications and reference grids remain untouched. The signature palette is a plotting convention, not a different scientific reference.
