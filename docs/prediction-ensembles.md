# Prediction ensemble analysis

RamplotR can analyse multiple independently generated prediction models as a
single **prediction ensemble**. The purpose is to expose where prediction
seeds/models agree or disagree in local backbone geometry and confidence.

This is deliberately not described as molecular dynamics. Variation between
predicted models can reflect model uncertainty, alternative plausible
solutions, stochastic sampling, different seeds or pipeline settings. It is
not experimental evidence that a protein moves between those conformations.

## Supported inputs

RamplotR accepts one coordinate file per model for:

- AlphaFold 2 / ColabFold;
- AlphaFold 3 sample sets;
- ESMFold;
- other prediction models that store pLDDT in the B-factor field.

For AF2/ColabFold, ESMFold and compatible B-factor predictions, the currently
loaded prediction can optionally be included as one ensemble member.

For **AlphaFold 3**, upload the sample coordinate files together with the
matching full-confidence JSON files. RamplotR pairs files by the official
seed/sample filename stem, not by upload order:

- `*_model.cif`;
- `*_confidences.json`;
- optional `*_summary_confidences.json`.

Every AF3 model must have exactly one matching full-confidence JSON. Missing,
duplicate or ambiguous pairs are rejected. The currently loaded AF3 model is
not silently added because its original sidecar path is not assumed to remain
available.

Each uploaded file must contain exactly one structural model. Up to 30 models
are analysed per run.

## Residue matching

Models are matched by:

- chain;
- residue number;
- insertion code;
- residue identity.

Missing residues remain missing and never receive fabricated coordinates,
angles or confidence scores. The Summary view reports how many residues are
present in every analysed model.

Duplicate coordinate files are rejected using MD5 hashes so repeated copies of
one model cannot artificially inflate ensemble agreement.

## Per-residue statistics

For every matched residue RamplotR reports:

- circular mean phi and psi;
- circular standard deviation of phi and psi;
- number of models contributing each angle;
- native RamplotR density-region agreement;
- coarse backbone-state mode and agreement (Alpha-R, Beta, PPII, Alpha-L or Other);
- Rama8000 Favored/Allowed/Outlier mode and agreement;
- mean, standard deviation, minimum and maximum pLDDT.

Circular statistics are essential because +179 degrees and -179 degrees are
neighbours rather than opposite conformations.

The prediction variability map uses the larger of phi-SD and psi-SD only as a
navigation aid:

| Maximum circular SD | Display |
| --- | --- |
| <5 degrees | stable |
| 5-15 degrees | moderate |
| 15-30 degrees | variable |
| >=30 degrees | high |

These bands are interface thresholds, not statistical significance cutoffs.

The variability map separately marks residues whose coarse backbone state
differs between models and residues whose Rama8000 category differs. The
backbone-state label is assigned from wrapped phi/psi proximity to fixed
canonical centres and is used only as a comparison bin. It is not a DSSP
secondary-structure assignment or a validation result.

## Model-level provenance

The model summary records:

- model/file label;
- residue count;
- number of residues with finite phi/psi;
- Rama8000 outlier count;
- mean and minimum pLDDT;
- declared prediction source;
- coordinate-file MD5 hash when available;
- for AF3 samples, pTM, ipTM, ranking score, disordered fraction and clash flag;
- for AF3 samples, hashes of the matched full-confidence and optional summary
  confidence files.

AF3 model-level ranking metrics are retained as provenance and comparison
context. They do **not** alter the residue-level φ/ψ variability, pLDDT spread
or Rama8000 agreement.

The two downloadable CSV files therefore preserve both residue-level ensemble
results and model-level provenance.

The **HTML ensemble report** packages the same information into a standalone
research record: overview metrics, model labels and coordinate hashes, the 25
largest local backbone-variability positions, all residue-level ensemble
statistics and the analysis settings used for the run. The report repeats the
interpretation warning that prediction disagreement is not experimental
evidence of molecular dynamics.

## Interpretation

A useful pattern is a residue with:

- large phi/psi spread;
- disagreement between broad backbone states;
- changing Rama8000 category;
- and/or large pLDDT variability across seeds.

That combination identifies a position worth inspecting in the 2D plot and 3D
structure. It does **not** by itself identify a biologically relevant
conformational switch.

Conversely, high pLDDT in every model does not imply that all models agree on
the same backbone geometry. RamplotR keeps confidence and geometric agreement
as separate measurements.

## Current limitations

- models are matched by residue identity rather than sequence realignment;
- AF3 sample pairing currently requires the standard filename suffixes so
  sample identity can be established without guessing;
- ensemble members are not superposed or clustered in 3D yet;
- the current map summarizes local backbone variability, not global domain
  motion;
- prediction ensembles are not a substitute for experimental ensemble or
  dynamics data.

Pairwise structural differences can be explored separately in the Compare tab,
including the Conformational Change Explorer.
