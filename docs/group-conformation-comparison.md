# Structure-group conformational comparison

RamplotR can compare two **sets** of related protein structures at residue level.
Typical examples include apo versus ligand-bound structures, wild type versus
mutants, or experimental versus predicted model sets.

The analysis is deliberately conformation-centred. It compares circular
backbone phi/psi distributions after sequence alignment instead of reducing
each structure pair to a single Cartesian RMSD.

## Workflow

1. Load a representative structure in RamplotR.
2. Open **Compare** and expand **Compare groups of structures**.
3. Select the protein chain that defines the reference residue coordinate
   system.
4. Define Group A and Group B labels.
5. Optionally include the currently loaded structure in Group A.
6. Upload one or more PDB/mmCIF files for each group.
7. Run **Analyse groups**.

Each uploaded file contributes model 1. RamplotR chooses the best matching
protein chain independently for every file using sequence identity and
reference-chain coverage. The default acceptance thresholds are 70% identity
and 70% reference coverage; both can be changed before analysis.

## Residue-level calculations

Every accepted candidate chain is globally sequence-aligned to the selected
reference chain. Candidate phi/psi values are then represented on the reference
residue coordinate system. Insertions without a reference residue are not
invented as comparable positions.

For each residue and each group RamplotR calculates:

- number of models contributing phi and psi;
- circular mean phi and psi using only members with **both angles observed**;
- circular standard deviation of phi and psi on that same paired subset;
- separate marginal phi/psi counts and the number of complete pairs;
- modal coarse backbone state and its within-group consistency;
- modal Rama8000 category and its within-group consistency.

The between-group effect is reported as wrapped differences between the two
circular means:

```
delta_phi = wrap(phi_mean_B - phi_mean_A)
delta_psi = wrap(psi_mean_B - psi_mean_A)
backbone_shift = sqrt(delta_phi^2 + delta_psi^2)
```

Wrapping is performed independently across the -180/180-degree boundary. When phi and psi come from different structures but no member has both, the group centroid and between-group displacement are **unavailable**. Partial observations remain recorded for coverage audits. This is also how the local conformational fingerprint overlays calculate their mean crosses.

Click a residue in the between-group track or results table to open its
**local conformational fingerprint**. The overlay contains one measured φ/ψ
point per group member with a complete angle pair. The member table preserves
missing angles rather than treating them as zero; where available, the actual
residue identity in the candidate structure is retained independently from
the reference residue used for sequence alignment.

The inspector reports complete-pair support and coarse backbone-state consensus
separately for both groups. Crosses mark circular group means. The overlay
does not estimate biological populations, conformational transition
probabilities or statistical significance. For Atlas-derived groups the
source is the cached exact UniProt-mapped first-model backbone evidence,
without an additional upload.

### Linking experimental contacts to selected residues

When the selected groups originate in verified Atlas, the local fingerprint
also shows a **Deposited component proximity** table for every selected
experimental member. It joins the current exact UniProt residue to all
first-model non-water nonpolymer HETATM components that have heavy atoms within
4.5 Å of the mapped protein residue. All mapped residues within the cutoff
are retained, even if a different residue is the component's closest atom
contact. Distances belong to the selected residue, not the entire component.

The table distinguishes observed nearby components, completed contact
extraction without a reported local component, incomplete nearest-only legacy
evidence, and unavailable contact evidence. Backbone φ/ψ completeness is shown
separately. It never treats an absent deposited component as a ligand-free
state, makes no causal or statistical claim and does not equate multiple
polymer entities with independent experiments. Download a CSV of the selected
residue's per-member observations for reproducibility. Uploaded structure
group comparisons do not display Atlas-only metadata.

The combined backbone shift is a **navigation effect size**, not a statistical
significance score and not a Cartesian distance.

## Evidence profile and support markers

RamplotR keeps the between-group effect and the evidence supporting it
separate. Every residue reports:

- the combined wrapped backbone shift;
- the maximum within-group circular SD across phi and psi;
- the fraction of structures in each group contributing a **complete phi/psi pair in the same model**;
- the minimum within-group Rama8000 modal-category consistency;
- an evidence profile.

The evidence profile is deliberately descriptive:

- **Small shift**: <15 degrees;
- **Moderate shift**: 15--30 degrees;
- **Large but variable**: >=30 degrees without low within-group dispersion;
- **Low-dispersion shift**: >=30 degrees with maximum within-group SD <=15 degrees;
- **Unreplicated structural difference**: fewer than two members in either group have a complete angle pair;\n- **Sparse coverage**: fewer than 75% of structures in either group contribute\n  a complete angle pair, when both groups have at least two paired observations;
- **Unavailable**: no comparable finite group means.

A residue receives the stronger **high-support shift** marker only when at least two members in each group provide complete phi/psi pairs and it is
a low-dispersion shift, both groups have at least 75% residue coverage, and
the Rama8000 modal category is at least 75% consistent within each group when
that information is available.

These thresholds are transparent navigation criteria, not a hypothesis test,
p-value or claim of biological significance. The exact group means, angular
differences, coverage and within-group dispersion remain visible in the table.

The table also reports whether the modal coarse backbone state differs between
the groups. These Alpha-R/Beta/PPII/Alpha-L/Other labels are broad phi/psi
comparison bins, not secondary-structure assignments. A separate marker
indicates when the modal Rama8000 category differs between the two groups.

## Interpretation

Group-level conformational differences can reflect genuine structural states,
but can also arise from:

- ligands or cofactors;
- construct boundaries and engineered mutations;
- crystal packing;
- cryo-EM classification;
- differences in experimental resolution;
- prediction uncertainty;
- different domain arrangements or oligomeric states.

Sequence alignment alone therefore does not prove that two groups are
biologically exchangeable. The exported member table records the automatically
selected chain, sequence identity and coverage for every input structure so
these assumptions can be audited.

For prediction ensembles, use RamplotR's dedicated prediction-ensemble
workflow when the main question is seed/model uncertainty. Group comparison is
more appropriate when the researcher has already defined biologically
meaningful sets such as apo/holo or WT/mutant.

## Exports

The group comparison can export:

- one residue-level CSV containing circular means, SDs, wrapped differences,
  combined displacement, backbone-state changes and Rama8000 mode changes;
- one long-format fingerprint CSV containing per-model phi/psi angles,
  residue identity, complete-pair flags and coarse backbone states;\n- one member CSV recording each structure label, selected chain, identity,
  reference coverage, candidate coverage and aligned residue count.

## Testing

Unit tests cover:

- automatic best-chain selection;
- rejection of unrelated chains at explicit thresholds;
- circular means across +179/-179 degrees;
- wrapped between-group differences;
- low-dispersion and high-support shift detection;
- sparse-coverage downgrading;
- modal backbone-state changes;
- Rama8000 modal-category changes.

The browser test additionally compares 1CRN against an identical uploaded
1CRN group and requires every residue to remain in the smallest shift band.
