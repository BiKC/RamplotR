# Experimental counterpart discovery

RamplotR can search for experimentally determined PDB structures related to a
loaded predicted protein model and open a selected hit directly in the linked
Compare workflow.

## Workflow

1. Load an AlphaFold DB, AlphaFold/ColabFold, ESMFold or other declared
   prediction.
2. Expand **Experimental counterparts** below the prediction-confidence panel.
3. Choose the prediction chain and a minimum sequence-identity threshold.
4. Select **Find experimental structures**.
5. Review the returned PDB polymer entities, experimental method and resolution.
6. Open the entry in PDBe for archive context, or choose **Compare in RamplotR**
   to load the experimental structure as the second model.

When a result is opened in RamplotR, the query chain is selected as the primary
comparison chain and the matching PDB chain returned by the archive metadata is
preferred for the experimental structure. The existing Conformational Change
Explorer then provides wrapped delta-phi/delta-psi values, combined local
backbone displacement, Rama8000 categories and linked 3D residue inspection.

## Search implementation

The browser sends the selected protein sequence to the public RCSB PDB Search
API sequence service:

- target: `pdb_protein_sequence`
- return type: `polymer_entity`
- content type: experimental structures only
- default minimum sequence identity: 90%
- maximum displayed hits: 12

Polymer-entity and entry metadata are then obtained from the RCSB Data API.
The network requests are made client-side in `experimental-search.js`, so the
same workflow can be used by the normal Shiny application and the static
Shinylive/webR deployment without adding a server-side HTTP dependency.

Result cards link to the corresponding current PDBe entry page for archive,
annotation and validation context.

## Interpretation

A sequence-similar PDB structure is a **candidate experimental counterpart**,
not proof that it represents the same biochemical state.

Differences can arise from:

- bound ligands or cofactors;
- mutations, constructs or truncations;
- oligomeric state;
- crystallisation or cryo-EM conditions;
- alternate functional states;
- experimental uncertainty;
- genuinely different local conformations.

For that reason RamplotR does not assign a single "agreement score" between a
prediction and an experimental hit. The search only identifies candidates.
Interpretation happens in the existing residue-level comparison workflow.

Lowering the identity threshold can be useful for homologous structures, but
the biological relevance of the comparison becomes increasingly dependent on
domain architecture and functional context.

## Privacy and browser execution

The selected amino-acid sequence is sent to the public RCSB PDB sequence search
service when the user explicitly runs the search. Coordinate files uploaded to
RamplotR are not automatically uploaded to a prediction service by this
feature. Loading a selected PDB hit uses RamplotR's normal public-accession
retrieval path.

## Testing

The browser integration test mocks the RCSB sequence and metadata responses so
continuous integration is deterministic and does not depend on temporary API
availability. It verifies:

- the request targets protein sequence search;
- the chosen identity threshold is transmitted;
- experimental-only search content is requested;
- candidate metadata and PDBe links render;
- **Compare in RamplotR** loads the candidate and selects the matching chain;
- the resulting aligned residues appear in the Conformational Change Explorer.
