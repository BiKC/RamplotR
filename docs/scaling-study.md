# Scaling experiment for RamplotR

Run the GitHub Actions workflow **Large-structure scaling benchmark** on the
scaling pull request or use its manual dispatch option. It fetches the
experimentally determined SARS-CoV-2 spike structure **6VXX** once and saves
the parsed Bio3D object to an RDS file. The workflow uses that same input for
all measurements.

Each scaling factor runs in a **separate R process on the same GitHub runner**:

| Multiplier | Input | Purpose |
| --- | --- | --- |
| 1 | Original 6VXX complex | Real-structure baseline |
| 3 | Three copies of 6VXX with unique chain names | Intermediate synthetic scale |
| 10 | Ten copies of 6VXX with unique chain names | Larger synthetic scale |

The replicated structures are computational stress tests, **not independent
experimental structures**. They preserve chain geometry, peptide continuity
and residue distributions while multiplying the number of atoms and residues.

Within each process, three iterations measure backbone extraction, density
classification and their total. The first iteration builds cached reference
profiles; the next two reuse them. Each multiplier's process is wrapped with
GNU `/usr/bin/time` to record peak resident memory in KiB and process wall
time. Structure fetching and parsing are deliberately excluded from timed
backbone and classification stages. The workflow stores the source-file SHA256,
CSV measurements, session metadata, and memory logs in one Actions artifact.

The script verifies that each replicated chain has identical torsion angles,
region labels and density percentiles to its corresponding first copy. Any
mismatch fails the workflow.

## Local repetition

Install R and Bio3D. From the repository root, prepare a source structure:

```sh
mkdir -p benchmarks/output/scaling
Rscript -e 'source("shinyRam/R/io.R"); saveRDS(ram_load_structure(pdb_id="6VXX"), "benchmarks/output/scaling/source-6vxx.rds")'
```

Then run each size in a separate process:

```sh
for factor in 1 3 10; do
  /usr/bin/time -f 'peak_rss_kb=%M\nwall_seconds=%e' \
    -o "benchmarks/output/scaling/scaling-${factor}.metrics.txt" \
    Rscript benchmarks/scale.R \
      benchmarks/output/scaling/source-6vxx.rds \
      "$factor" 3 original \
      "benchmarks/output/scaling/scaling-${factor}.csv"
done
```

## Interpretation and limits

The goal is to describe scaling within one implementation, not to claim a
speedup over the historical RamplotR version. The scientific classification
algorithm changed between the tagged historical version and current `main`,
so a later head-to-head comparison must separately report equivalent work
and changes in classification semantics.

GitHub-hosted runners have variable loads. Repeat any result intended for the
paper on one documented machine, use more iterations, and report distributions
rather than a single minimum. The peak RSS measurement is **process-wide**;
it includes the R runtime, structure data and cached reference distributions,
and cannot be attributed to one individual stage.
