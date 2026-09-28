# AlphaFold and ESMFold prediction confidence

RamplotR adds model-confidence evidence to its existing structure-inspection
workflow. Confidence complements stereochemical analysis and experimental
validation; it cannot replace either.

## Loading a prediction

**AlphaFold DB**: choose the AlphaFold DB input option and enter a UniProt
accession. RamplotR requests the official prediction API (using the current
fields), downloads its published PDB when available, and retrieves the
matching PAE JSON if the API provides a supported URL. It does not send user
uploads to AlphaFold DB or any external folding service.

**AlphaFold 2 / ColabFold**: choose Upload file, select the coordinates,
then declare AlphaFold 2 / ColabFold. Per-residue pLDDT comes from PDB
B-factors or parsed mmCIF equivalents. Upload the model's own PAE JSON to
inspect relative-placement uncertainty. The modern AlphaFold DB JSON shape
has a 2D "predicted_aligned_error" array.

**AlphaFold 3**: upload the corresponding mmCIF model and select AlphaFold 3.
The optional full confidence JSON supplies per-atom "atom_plddts" and
per-token "pae". The optional summary confidence JSON supplies pTM/ipTM.
RamplotR checks per-atom confidence against the coordinate model's atom
count and chain order. Because AlphaFold 3's full JSON does not independently
identify every atom by name and residue, **always upload the matching mmCIF and
confidence JSON from the same seed and sample**. The program cannot establish
the identity of a separately reordered or substituted model using the JSON
alone. It displays PAE only when token chain/residue IDs uniquely map to
protein residues. Ligand and nucleic-acid tokens
are not silently assigned to protein residue numbers.

**ESMFold**: upload the PDB returned by ESMFold and explicitly select
ESMFold. Its B-factor field stores per-residue pLDDT. Standard ESMFold
does not supply an AlphaFold PAE matrix; RamplotR does not fabricate one.

Choose Experimental or unknown for a structure with no verified prediction
provenance. Experimental thermal B-factors must never be reinterpreted as
prediction confidence.

## Linked inspection

A small optional pLDDT strip appears under each chain's existing region
track. On predicted structures, a collapsed Confidence panel reports mean
pLDDT and offers a directional PAE heatmap when suitable data exist.
Click a PAE axis position to highlight that residue in the shared inspector,
Ramachandran plot, sequence navigator and 3D structure. Selected high-pLDDT
Ramachandran outliers are explicitly flagged for inspection; confidence does
not automatically make unusual geometry correct.

Extremely large PAE matrices are capped at 1,600 tokens on import. Heatmaps
over 400 mapped protein tokens are downsampled for display with that sampling
clearly labelled. This does not change the underlying PAE matrix. Models
with ambiguous insertion codes or missing AF3 token mapping have no PAE
heatmap rather than a potentially incorrect one.

Confidence is attached to a specific prediction model, not automatically
transferred to unrelated models in an NMR ensemble. For additional models
choose the model and supply compatible confidence data separately.

## Interpretation and provenance

pLDDT is a prediction of local correctness, while PAE measures expected
relative-position error in angstroms and is directional (row aligned, column
target). Local stereochemistry can be plausible despite low prediction
confidence, and vice versa.

The HTML report includes prediction source, model ID, provided confidence
file and PAE availability. For rigorous publication, archive the original
prediction coordinates, all confidence JSON files and tool versions.

Sources:
- ESMFold: https://github.com/facebookresearch/esm
- AFDB PAE JSON: https://alphafold.ebi.ac.uk/faq
- AF3 output format: https://github.com/google-deepmind/alphafold3/blob/main/docs/output.md
- AFDB API field transition (2026): https://www.ebi.ac.uk/pdbe/news/breaking-changes-afdb-predictions-api

For a local uploaded structure, the **Prediction settings** section opens while choosing its prediction type and optional matching confidence JSON. It collapses automatically once analysis finishes to leave more room for the Ramachandran plot. Reopen it any time to change the declared source or confidence sidecars and run the analysis again.
