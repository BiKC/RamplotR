# Firefox/Zen startup recovery (fifth bundle)

A local Firefox 159 test reproduced an automatic navigation loop in the
vendor startup loader: serviceWorker.ready resolved before controller became
available, triggering reloads repeatedly and aborting module loads. The worker
was activated but the page never finished gaining control. This was an actual
reproduction in Firefox; the user's Zen profile/URL was not available.

An app-local bootstrap now waits for controllerchange and checks the expected
worker URL before importing Shinylive. It never reloads automatically. A
worker-control timeout/cancellation cleans up listeners; a failed startup
shows an explicit message and a manual reload button. The shared runtime and
scientific code are untouched by this fix. The bootstrap filename is hashed
so a deployment does not reuse an older cached entry script.

Local Firefox testing with the fixed entry completed 28OE initialization,
selected C before the first plot, observed the initial pending indicator
clear, switched to D, and successfully reloaded. Unit checks cover delayed
control, unrelated worker control, timeout and cancellation. The generated
entry replacement is also checked. Existing numerical regressions are
unaffected. The test browser was Firefox, not Zen itself.

# Startup chain selection correction (fourth bundle)

The previous interaction tests waited until loading settled. A new startup
regression test selects chain C as soon as the picker is initialized, before
the first plot or model is ready.

With the previous bundle, NGL representations became available 4.217 s after
the picker event, but the selected-chain filter was not applied until 5.994 s:
an additional 1.777 s showing the old chain selection. With this fix, the
representations became available and were correctly filtered in the same
sampled animation frame (2.754 s after the event). These are local startup
runs with variable network and first-use costs, not controlled total-load
benchmarks. The stable finding is that the extra readiness/acknowledgement
wait is removed; initial geometry and scientific calculations still take time.

Chain identities are now published before classification/viewer readiness.
The picker carries its structure identity and all choices, so an early click
can establish client state before the first server handshake. That handshake
preserves the choice. Newly created NGL representations immediately take the
latest visibility, and NGL lookups are safe before widget initialization.

Regression tests cover early selection, a delayed first handshake and late
representation creation. The exported app was tested on 28OE by selecting C
before its first plot, with both the previous and fixed code. Prior settled
interaction and scientific regression results below remain relevant.

# 28OE performance: third optimized bundle

## Latest measured interactions

These are single local Chromium runs with software WebGL at 1440 × 1000.
Reference files were served locally, with original deployment bytes. Actual
RCSB 28OE coordinates were used. Automated interaction durations include
browser scheduling and Playwright actionability waits; they are indicative,
not a controlled benchmark or a promise for other hardware.

| Interaction | Previous bundle | Third bundle |
| --- | ---: | ---: |
| Cached DNA on while a reference-dataset calculation is running | 4.127 s | 0.029 s |
| Residue click, median of three clicks | 1.245 s | 0.711 s |
| Residue-click range | 0.755–1.360 s | 0.559–0.873 s |
| First residue-list opening | 0.767 s | 1.328 s |

The concurrency test first enabled and disabled DNA, changed the real
reference dataset from original to alphafold, waited for Shiny to become
busy, and enabled DNA again. Previously, DNA appeared after R finished.
With the new client path it appeared while R was still busy. This measures
actual input-event-to-visibility time and excludes Playwright click waits.
The previous run used lightweight R timing instrumentation.

Additional third-bundle measurements:

- First DNA enable: 0.412 s; repeated enable: 0.114 s, no rebuild.
- Plot residue selection: 0.173 s with zero Plotly.react calls.
- Background changes: 0.369–0.863 s; validation-mode changes: 0.425–1.285 s.
- Ribbon/cartoon changes: 0.807–1.870 s, with actual radio input bindings.
- The profiled first validation pipeline took 1.811 s; its components
  included native reference classification 1.449 s, standard validation
  0.153 s, and geometry 0.176 s. These include first-use costs.
- The uninstrumented final run needed 51.428 s for cold WASM startup and
  13.944 s after submitting the structure. Startup/network costs remain.

First residue-list initialization did **not** improve consistently in these
runs. Deferred table DOM creation reduces scaling work, but first-use JS,
layout and initialization still cost about a second. RNA and ligand toggles
were exercised, but 28OE gave zero selected atoms for those layers; their
short timings do not establish performance for large nonempty RNA/ligands.

## Changes in the third bundle

DNA/RNA/ligand controls now update NGL locally and cache representations.
Hiding and showing a layer retains its GPU buffers. They respond even when
R is busy. Layer completion tracks NGL tasks, and rendering-style controls
also keep their progress indicator until representation work finishes.

The residue inspector no longer updates hidden Plotly and NGL views. Returning
to the plot applies the current selection. The main plot selection marker is
a positioned DOM overlay, avoiding a full Plotly restyle per residue click.
Redundant server/client highlight updates are removed. Large structures use
instant camera moves rather than generating animated full-scene frames.
Sequence selection writes only changed accessibility attributes.

Nearby-hetero inspection rejects distant atoms with a conservative bounding
box, then computes distances in vectors rather than copying a data frame for
every hetero atom. Independent all-pairs tests verify distances and nearest
atom identities, including ties and distance boundaries.

Native density validation, six-class standard validation, geometry and
canonical annotations have separate reactive caches. A plotting-background
change in residue-aware mode no longer reclassifies the structure; changing
validation mode/reference does not redo standard validation or geometry.
Tables/exports retain canonical annotations and selected-chain filtering.

The exporter stages exact binary Rama8000 matrices instead of parsing six
CSV grids in WASM. All matrix values and metadata are checked against CSV.
Normal R still supports the original CSV files. The thin app.json decreased
from roughly 3.3 MB to 2.4 MB; reference-data bytes are unchanged.

## Retained earlier optimizations

Chain selection reuses existing plot traces, NGL representations and sequence
sections. Full sequence markup is generated only on expansion and then
cached. Its vectorized renderer produced an equivalent DOM for all 1,428
residues; native R generation/serialization measured 4.19 s versus 0.25 s,
and the earlier WASM opening test measured 1.111 s. Large structures start
without rocking/spinning and use low representation quality. Feedback starts
at the picker event and waits for the latest selection to finish; stale
responses cannot reverse a newer choice.

## Chain C compared with aa3eba2

The deposited PDB and mmCIF files give identical current chain C residue
numbers and names. Author and label identifiers agree for every protein
chain in 28OE. The chain-wise extraction used by aa3eba2 and current code
both return the same 222 protein residues for chain C.

All 217 mutually finite phi values and all 217 mutually finite psi values
agree within 0.01 degrees. Legacy finite angles now marked missing occur at
phi residues 225, 243, 282, 302 and psi residues 136, 240, 276, 297. The current
extractor requires a peptide bond and does not bridge coordinate gaps.
No chain renaming was found.

## Validation

Base R UI, inspection, local-neighbor all-pairs comparison, geometry,
peptide, scientific classification and six-class Rama8000 regressions;
exact CSV/binary matrix checks; actual Shiny cache invalidation and score
preservation tests; vectorized sequence controls/escaping tests; and Node
feedback, plot-marker, layer-cache, style-completion and rapid-input tests.

The final exported app was exercised in Chromium through the residue list,
three residue clicks, a plot selection, DNA/RNA/ligand on/off/on, background
and validation changes, and ribbon/cartoon changes. The separate concurrency
test used the third bundle with the final style/task feedback code.
