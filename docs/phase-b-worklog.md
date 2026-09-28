# Phase B — AlphaFold and ESMFold confidence integration

Base: main `99fd7f38c17e88dce28b043560ecd5a147d647ea`. Scientific Phase A and historical tag stay unchanged.

- [ ] Native prediction provenance and per-residue pLDDT extraction from AlphaFold and ESMFold B-factors; experimental B-factors must never be relabelled confidence.
- [ ] Optional AlphaFold DB accession retrieval and JSON confidence sidecar uploads (AlphaFold 2 PAE and AlphaFold 3 full confidence/summary JSON).
- [ ] Validate AF3 atom-vs-token indexing and multi-chain residue mapping before displaying PAE; do not silently assign unmatched tokens.
- [ ] Compact linked confidence track beneath every protein chain, interactive PAE heatmap, contextual confidence-vs-geometry review in inspector. Keep these panels absent when prediction data are unavailable.
- [ ] Record prediction source, confidence file/provenance and limitations in documentation and publication export.
- [ ] Synthetic tests for AF2/AF3/ESMFold formats, missing confidence, chain alignment, malformed/huge PAE and independent existing R tests. Browser smoke test, CI review, PR/merge when green.

ESMFold PDB stores pLDDT in B-factor fields. AlphaFold3's full confidence JSON has per-atom `atom_plddts` and per-token `pae`; the token indices are **not** atom indices. AF2/DB PAE JSON may use `predicted_aligned_error` and `residue1`/`residue2`. Disallow unsupported formats rather than fabricating confidence values.
