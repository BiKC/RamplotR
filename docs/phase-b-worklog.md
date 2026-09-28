# Phase B — AlphaFold and ESMFold confidence integration

Base: main `99fd7f38c17e88dce28b043560ecd5a147d647ea`. Scientific Phase A and historical tag stay unchanged.

- [x] Native prediction provenance and per-residue pLDDT extraction from AlphaFold and ESMFold B-factors; experimental B-factors must never be relabelled confidence.
- [x] Optional AlphaFold DB accession retrieval and JSON confidence sidecar uploads (AlphaFold 2 PAE and AlphaFold 3 full confidence/summary JSON).
- [x] Validate AF3 atom-vs-token indexing and multi-chain residue mapping before displaying PAE; do not silently assign unmatched tokens.
- [x] Compact linked confidence track beneath every protein chain, interactive PAE heatmap, contextual confidence-vs-geometry review in inspector. Keep these panels absent when prediction data are unavailable.
- [x] Record prediction source, confidence file/provenance and limitations in documentation and publication export.
- [x] Synthetic tests for AF2/AF3/ESMFold formats, absent confidence, insertion-code ambiguity, malformed/oversized PAE and safe token mapping.
- [x] All-chain and multi-model scientific regression tests pass on Ubuntu and Windows; independent wwPDB and structure-validation workflows remain green.
- [x] Real Chromium/Shiny smoke tests pass for ESMFold, AF2 PAE and AF3 full/summary confidence upload, in addition to desktop/mobile interactions.
- [x] Downsampled PAE plots preserve both sides of protein-chain boundaries, and missing phi/psi are never labelled as acceptable geometry.
- [x] Merged PR #14 into main as `dad437d132c58452caa46950d644172e19b2eaa6`. Post-merge scientific regression tests passed on Ubuntu and Windows and all five independent wwPDB checks passed. The final PR-head browser and structure-validation workflows also passed.

ESMFold PDB stores pLDDT in B-factor fields. AlphaFold3's full confidence JSON has per-atom `atom_plddts` and per-token `pae`; the token indices are **not** atom indices. AF2/DB PAE JSON may use `predicted_aligned_error` and `residue1`/`residue2`. Disallow unsupported formats rather than fabricating confidence values.
