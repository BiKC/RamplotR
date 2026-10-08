# Reproducible real-crystal benchmark of the Atlas scientific methods.
# Inputs are built by scripts/atlas-adk-fetch.cjs using the SAME SIFTS/mmCIF
# extractor as the browser, not PDB author-number interpolation.
# Usage: Rscript scripts/atlas-adk-benchmark.R benchmarks/output/atlas-adk

source(file.path("shinyRam","R","canonical.R"))
source(file.path("shinyRam","R","backbone.R"))
source(file.path("shinyRam","R","atlas-geometry.R"))
source(file.path("shinyRam","R","atlas-switch.R"))
if(!requireNamespace("jsonlite",quietly=TRUE))
  stop("Install jsonlite before running the live archive benchmark.")
output <- if(length(commandArgs(TRUE))) {
  commandArgs(TRUE)[[1L]]
} else {
  "benchmarks/output/atlas-adk"
}
accession <- "P69441"
load_entry <- function(id) {
  filename <- file.path(output,paste0(id,"_1.json"))
  if(!file.exists(filename)) stop("Missing SIFTS benchmark data: ",filename)
  data <- jsonlite::fromJSON(filename,simplifyVector=TRUE)
  if(!identical(as.character(data$accession),accession) ||
     !identical(as.character(data$pdb_id),id) ||
     !identical(as.character(data$state),"mapped"))
    stop("Wrong or incomplete benchmark entity: ",id)
  if(!is.data.frame(data$mapping) ||
     !is.data.frame(data$ca_points) ||
     !is.data.frame(data$backbone_atoms))
    stop("Malformed benchmark coordinate data for ",id)
  data
}
records <- list("4AKE_1"=load_entry("4AKE"),"1AKE_1"=load_entry("1AKE"))
manifest <- jsonlite::fromJSON(file.path(output,"manifest.json"))
if(!all(c("4AKE_1","1AKE_1") %in% manifest$entries$id))
  stop("Missing archive source-provenance hashes.")
if(any(!grepl("^[a-f0-9]{64}$",manifest$entries$sha256)))
  stop("Benchmark file provenance hashes missing.")

chain_records <- function(record) {
  map <- record$mapping
  table <- sort(table(map$struct_asym_id[!is.na(map$observed) &
                                        map$observed]),decreasing=TRUE)
  selected <- names(table)[as.integer(table)>=100L]
  if(length(selected)<2L)
    stop("Expected at least two sufficiently complete crystal chains.")
  selected[1:2]
}
chain_ids <- lapply(records,chain_records)
print(chain_ids)

exact_ca <- function(record,asym) {
  m <- record$mapping
  m <- m[m$struct_asym_id==asym & !is.na(m$observed) &
         m$observed,c("label_seq_id","uniprot_resi"),drop=FALSE]
  m <- unique(m)
  # Require a one-to-one label<->canonical coordinate relation.
  duplicate <- duplicated(m$label_seq_id) |
    duplicated(m$label_seq_id,fromLast=TRUE) |
    duplicated(m$uniprot_resi) |
    duplicated(m$uniprot_resi,fromLast=TRUE)
  m <- m[!duplicate,,drop=FALSE]
  ca <- record$ca_points
  ca <- ca[ca$struct_asym_id==asym,
    c("label_seq_id","x","y","z"),drop=FALSE]
  duplicate_ca <- duplicated(ca$label_seq_id) |
    duplicated(ca$label_seq_id,fromLast=TRUE)
  ca <- ca[!duplicate_ca,,drop=FALSE]
  result <- merge(m,ca,by="label_seq_id")
  result <- result[order(result$uniprot_resi),,drop=FALSE]
  if(anyDuplicated(result$uniprot_resi) || nrow(result)<100L)
    stop("Insufficient verified C-alpha canonical coverage.")
  result
}

pair <- function(a_name,a_chain,b_name,b_chain,kind) {
  A <- records[[a_name]]
  B <- records[[b_name]]
  caA <- exact_ca(A,a_chain); caB <- exact_ca(B,b_chain)
  common <- sort(intersect(caA$uniprot_resi,caB$uniprot_resi))
  if(length(common)<100L)
    stop("Insufficient pairwise common exact UniProt positions.")
  if(length(common)/nrow(caA)<.6 || length(common)/nrow(caB)<.6)
    stop("Incompatible experimental constructs or missing coordinates.")
  xa <- as.matrix(caA[match(common,caA$uniprot_resi),
                      c("x","y","z"),drop=FALSE])
  xb <- as.matrix(caB[match(common,caB$uniprot_resi),
                      c("x","y","z"),drop=FALSE])
  da <- as.numeric(stats::dist(xa))
  db <- as.numeric(stats::dist(xb))
  drmsd <- sqrt(mean((da-db)^2))
  # Contact-map change per aligned residue describes domain rearrangement,
  # not backbone torsion or per-residue uncertainty.
  dma <- as.matrix(stats::dist(xa))
  dmb <- as.matrix(stats::dist(xb))
  contact_shift <- sqrt(rowMeans((dma-dmb)^2))
  torsA <- ram_atlas_entity_torsions(A,a_chain)
  torsB <- ram_atlas_entity_torsions(B,b_chain)
  difference <- ram_atlas_torsion_delta(torsA,torsB,a_name,a_name,30)
  difference$entity_b <- b_name
  regions <- ram_atlas_candidate_regions(difference)
  label <- paste0(a_name,"-",a_chain,"__",b_name,"-",b_chain)
  utils::write.csv(difference,file.path(output,paste0(label,"_torsions.csv")),
                   row.names=FALSE,na="")
  utils::write.csv(data.frame(uniprot_resi=common,
    ca_distance_map_change=contact_shift),
    file.path(output,paste0(label,"_distance_map_track.csv")),
    row.names=FALSE)
  list(summary=data.frame(
    pair=label,kind=kind,common_positions=length(common),
    coverage_a=round(length(common)/nrow(caA),3),
    coverage_b=round(length(common)/nrow(caB),3),
    ca_distance_map_rmsd_angstrom=round(drmsd,4),
    comparable_phi_psi=sum(difference$comparable),
    above_30deg=sum(difference$candidate),
    candidate_regions=nrow(regions),
    stringsAsFactors=FALSE),
    positions=common,shifts=contact_shift,regions=regions)
}
open_chains <- chain_ids[["4AKE_1"]]
closed_chains <- chain_ids[["1AKE_1"]]
checks <- list(
  pair("4AKE_1",open_chains[[1L]],"1AKE_1",closed_chains[[1L]],
       "known open-versus-closed"),
  pair("4AKE_1",open_chains[[1L]],"4AKE_1",open_chains[[2L]],
       "same crystal, open-state chain copies"),
  pair("1AKE_1",closed_chains[[1L]],"1AKE_1",closed_chains[[2L]],
       "same crystal, closed-state chain copies"),
  pair("4AKE_1",open_chains[[1L]],"4AKE_1",open_chains[[1L]],
       "identity numerical negative control")
)
summary <- do.call(rbind,lapply(checks,function(x)x$summary))
stopifnot(all(is.finite(summary$ca_distance_map_rmsd_angstrom)),
          all(summary$common_positions>=100),
          summary$ca_distance_map_rmsd_angstrom[[4L]]<1e-8,
          summary$above_30deg[[4L]]==0L)
utils::write.csv(summary,file.path(output,"comparison-summary.csv"),
                 row.names=FALSE)
# Literature domain divisions approximate E. coli AK, NOT curated residue
# labels used for validation. Values are exploratory descriptive means only.
target <- checks[[1L]]
domain <- ifelse(target$positions %in% 30:67,"NMP",
           ifelse(target$positions %in% 118:167,"LID","CORE"))
domain_report <- aggregate(target$shifts,list(domain=domain),mean)
names(domain_report)[2L] <- "mean_ca_distance_map_change_angstrom"
utils::write.csv(domain_report,file.path(output,"domain-descriptive.csv"),
                 row.names=FALSE)
summary_text <- c(
  "# Adenylate kinase: experimental Atlas case study",
  "",
  "4AKE (open, apo) versus 1AKE (closed, inhibitor-bound).",
  "Accession P69441. First coordinate model only.",
  "SIFTS polymer label-asym/seq mapping, never author-residue offsets.",
  "The two chains within each crystal are packing/conformation controls,",
  "**not independent biological replicates**.",
  "",
  "No parameter thresholds were trained, selected, or calibrated here.",
  "A 30-degree threshold is a pre-existing navigation rule, not a classifier.",
  "",
  "## Measured comparisons",
  "",
  paste(capture.output(print(summary,row.names=FALSE)),collapse="\n"),
  "",
  "## Per-residue distance-map change, descriptive domains",
  "",
  paste(capture.output(print(domain_report,row.names=FALSE)),collapse="\n"),
  "",
  "## Provenance",
  "",
  paste(manifest$entries$id,manifest$entries$url,
        manifest$entries$sha256,sep=" SHA256=",collapse="\n"),
  "",
  "Interpretation: global domain motion may have few local torsion changes;",
  "those two measurements must remain distinct. Experimental errors,",
  "crystal contacts, ligands, resolution and sequence/construct differences",
  "are possible confounders. One open/closed pair cannot calibrate a",
  "general conformation-state classifier.",
  ""
)
writeLines(summary_text,file.path(output,"benchmark-report.md"))
print(summary,row.names=FALSE)
print(domain_report,row.names=FALSE)
cat("Live experimental Atlas benchmark completed with exact mapping.\n")
