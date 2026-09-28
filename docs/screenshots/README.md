# Example screenshots

The screenshots are cropped from real Shiny browser tests, not mock-up renders.
Source: https://github.com/BiKC/RamplotR/actions/runs/36488764177

- Overview and linked inspection: uploaded experimental 1CRN structure.
- Multi-chain sequence: 1BBB multimer.
- Prediction PAE: synthetic AlphaFold 2/ColabFold-style confidence test data (illustrative only).
- Ensemble: 1D3Z NMR multi-model structure.

To regenerate, extract a recent successful interface browser-test artifact
to benchmarks/output/ui-preview, then run
python scripts/curate-readme-screenshots.py.
Review the crops after layout changes.
