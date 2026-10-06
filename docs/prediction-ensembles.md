# Prediction ensemble analysis

RamplotR can analyse multiple independently generated prediction models as a
single **prediction ensemble**. The purpose is to expose where prediction
seeds/models agree or disagree in local backbone geometry and confidence.

This is deliberately not described as molecular dynamics. Variation between
predicted models can reflect model uncertainty, alternative plausible
solutions, stochastic sampling, different seeds or pipeline settings. It is
not experimental evidence that a protein moves between those conformations.

## Supported inputs

The first implementation accepts one coordinate file per model for:

- AlphaFold 2 / ColabFold;
- ESMFold;
- other prediction models that store pLDDT in the B-factor field.

The currently loaded compatible prediction can optionally be included as one
ensemble member.

AlphaFold 3 is intentionally excluded from automatic ensemble confidence
analysis for now. AF3 confidence is atom/token based and needs its matching
confidence JSON sidecar; RamplotR does not silently reinterpret an AF3 model as
AF2.

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

A separate outline marks residues whose Rama8000 category differs between
models.

## Model-level provenance

The model summary records:

- model/file label;
- residue count;
- number of residues with finite phi/psi;
- Rama8000 outlier count;
- mean and minimum pLDDT;
- declared prediction source;
- coordinate-file MD5 hash when available.

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
- AF3 multi-model confidence sidecars are not yet supported;
- ensemble members are not superposed or clustered in 3D yet;
- the current map summarizes local backbone variability, not global domain
  motion;
- prediction ensembles are not a substitute for experimental ensemble or
  dynamics data.

Pairwise structural differences can be explored separately in the Compare tab,
including the Conformational Change Explorer.
