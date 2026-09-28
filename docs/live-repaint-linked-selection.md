# Live controls and linked residue inspection

Structure parsing and backbone torsions are calculated when **Analyze structure**
is pressed. Adjusting the reference dataset, background, residue and chain
filters, validation mode, chain colors, or contour palette then updates the
existing plot and associated statistics automatically. No subsequent Analyze
click is required.

Residue selection is synchronized by chain, residue number and PDB insertion
code. Clicking a Ramachandran point highlights the same residue with an orange
ring in the 2D plot and a ball-and-stick overlay in the NGL viewer. Clicking a
residue table row or a residue in the NGL viewer selects the same residue in
the other views. **Clear selection** removes both highlights.

Selections are ephemeral within each Shiny session. They do not change the
underlying structure or scientific classification. Unclassified residues can
still appear in the table; the 2D plot only shows residues with valid phi/psi
angles.

## Verification

The cross-platform regression workflow checks the existing scientific
calculations and the plot-message handler. A headless browser loads 1CRN
through the file-upload path, changes palettes without clicking Analyze,
selects and clears a plot residue, verifies 3D click messages and residue
table rows, and checks desktop/mobile plot dimensions. This is separate from
independent scientific validation against Bio3D.
