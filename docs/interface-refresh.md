# Interface refresh

The interface refresh changes the layout and plotting presentation without
changing backbone extraction, reference grids or scientific classifications.
The historical application remains available as v0.1.0-legacy.

## What changed

- An explicit structure-source selector shows either the PDB accession field
  or the PDB/mmCIF upload control, with one analysis button.
- Reference dataset, scientific classification, filters and color controls
  are grouped into a sidebar. Detailed palette controls are collapsible.
- The Ramachandran plot and NGL structure viewer share the main analysis tab. On mobile, results appear above advanced settings, with a direct Settings shortcut.
  The residue table and summary have dedicated tabs.
- Plotly uses a responsive container, equal axis scaling, clearer hover labels,
  smaller outlined chain markers, a non-overlapping horizontal legend for up to
  eight chains and image export through the existing mode bar.
- The 3D viewer now uses a soft light background, stronger chain colors and
  controls located directly beneath the molecular structure, rather than
  occupying a separate sidebar panel.
- The summary uses larger totals and an aligned table while preserving the
  exact existing counts. Undefined-angle counts include broken or incomplete
  backbones as well as genuine termini.
- The interface adapts to narrower screens and has visible keyboard focus.

No new R package is required. Styles are local to
shinyRam/www/styles.css and the plotting code is in
shinyRam/www/custom.js.

## Checks

The scientific regression workflow parses app.R, checks the original Shiny
input/output IDs, verifies the responsive-layout contract and runs the Node
plotting test alongside the scientific R tests on Windows and Ubuntu.

## Manual visual smoke test

From the repository root, install dependencies listed in README.md and run:

    shiny::runApp("shinyRam")

In a browser, check both a wide desktop window and a mobile-width window:

1. Load 1CRN by accession, then switch to a local PDB/mmCIF upload.
2. Change background/reference choices, filter chains and amino acids and
   confirm the plot, residue list and summary remain consistent.
3. Confirm the angular axes stay equally scaled when resizing the browser.
4. Switch between Plot, Residue list and Summary tabs. Test the 3D controls
   beneath the viewer. On mobile, confirm the results appear before the
   settings and that the shortcut scrolls to the settings.
5. Test a structure with many chains and a structure with missing angles,
   then download a PNG with the Plotly toolbar.

The old README screenshot is historical; capture a new screenshot from a
running app before replacing it. The shinyapps.io deployment is separate
from merging changes into the GitHub repository.
