# RamplotR inspection and publication upgrade

Branch: upgrade/inspection-visualization-reporting. Preserve v0.1.0-legacy and all scientific reference grids.

## Work packages
- [ ] Fix DT header/body alignment after a hidden tab becomes visible, improve columns and filtered CSV exports.
- [ ] Shared residue inspector across tabs, outlier navigation, missing-angle inspection, sequence navigation and Show in plot.
- [ ] Publication-quality RamplotR default palette, preserving Rampage, PDBSum and custom settings.
- [ ] Accessible modern NGL viewer switches and cartoon, sticks, ball+stick and surface representations.
- [ ] PNG/SVG and metadata-rich report export.
- [ ] Multi-model support, explicitly retaining model provenance and validating coordinates.
- [ ] Two-structure comparison by sequence alignment, angular differences and region classification changes.
- [ ] Unit and browser tests, scientific regression checks, documentation and integration to main.

## Scientific boundaries
Do not compare residue indices without sequence alignment or silently compare missing torsions. Classifications and reference grids remain untouched. The signature palette is a plotting convention, not a different scientific reference.
