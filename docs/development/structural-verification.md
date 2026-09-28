# Structural verification and batch processing: implementation record

Base: `dad437d132c58452caa46950d644172e19b2eaa6` (the independent-validation and prediction-confidence work already integrated into `main`). Preserve scientific reference grids and `v0.1.0-legacy`.

## Work packages

- [x] **Extended geometry**: compute interpretable omega cis/trans/twisted peptide flags and independent side-chain chi1 and Cβ geometry. Use explicit methodology and missing-value semantics; for authoritative rotamer/clash results, import wwPDB validation results rather than claim an ad hoc calculation is MolProbity-equivalent.
- [x] **Experimental evidence**: allow optional local official wwPDB validation XML and show its residue-specific rotamer/clash/outlier annotations; link to source and show differences against the native RamplotR region. Optional NGL density-map view for local cryo-EM maps with user-controlled contours, not an unvalidated map/model fit score.
- [x] **Ensembles**: compare coherent NMR/prediction models residue by residue using circular angle statistics, coverage and class-change summaries. Identify missing/misaligned atom records, not index-only matching.
- [x] **Batch interface**: a scriptable offline R command, safe output directories, per-residue CSV/JSON and standalone HTML/SVG exports; consistent options and provenance. Structure and reference files remain local.
- [x] Tests on synthetic fixtures and genuine multimeric/NMR samples, CI and docs. Keep advanced controls optional; no extra permanent tabs without need.

## Methodological limits

Independent wwPDB classifications and MolProbity rotamers/clashscores have different algorithms/reference populations and must not be relabelled as new native RamplotR reference classifications. Experimental B factors are not pLDDT. Map visuals are qualitative unless a validated map-fit engine is integrated. Aggregate ensemble measurements require matching residue IDs and explicitly report missing data. AlphaFold/ESMFold predictions remain supported, with their provenance retained.

## Integration

Implemented and verified in [PR #15](https://github.com/BiKC/RamplotR/pull/15). This is an archived implementation checklist; current usage and scientific limits are described in the [structural-verification guide](../structural-verification.md).
