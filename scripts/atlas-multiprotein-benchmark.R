# Multi-protein experimental Atlas case study. Not a state classifier.
# Run from repository root after the node exact-SIFTS fetcher.
source(file.path("shinyRam","R","canonical.R"))
source(file.path("shinyRam","R","backbone.R"))
source(file.path("shinyRam","R","atlas-geometry.R"))
source(file.path("shinyRam","R","atlas-switch.R"))
if(!requireNamespace("jsonlite",quietly=TRUE))
  stop("This optional benchmark requires jsonlite.",call.=FALSE)
argv <- commandArgs(TRUE)
dest <- if(length(argv)) argv[[1L]] else "benchmarks/output/atlas-multiprotein"
conf <- jsonlite::fromJSON("benchmarks/atlas-multiprotein-cases.json",
                           simplifyVector=FALSE)
manifest <- jsonlite::fromJSON(file.path(dest,"manifest.json"),
                                simplifyVector=FALSE)
if(!identical(conf$schema_version,1L) || length(conf$proteins)<3L)
  stop("Missing curated multi-protein benchmark configuration.")
if(length(manifest$entries)<7L) stop("Incomplete live benchmark archive.")
cache <- list()
for(protein in conf$proteins) for(entry in protein$entries) {
  id <- paste0(entry$pdb,"_",entry$entity)
  payload <- jsonlite::fromJSON(file.path(dest,paste0(id,".json")),
                                  simplifyVector=TRUE)
  if(!identical(payload$accession,protein$uniprot) ||
     !identical(payload$pdb_id,entry$pdb) ||
     !identical(payload$state,"mapped"))
    stop("Invalid/mismatched live experiment ",id)
  fields <- c("mapping","ca_points","backbone_atoms","sequence_audit")
  if(!all(vapply(payload[fields],is.data.frame,logical(1L))))
    stop("Missing exact residue/atom tables for ",id)
  audit <- Filter(function(x) identical(x$id,id),manifest$entries)
  if(length(audit)!=1L ||
     !grepl("^[a-f0-9]{64}$",audit[[1L]]$sha256))
    stop("Incomplete source provenance for ",id)
  cache[[id]] <- payload
}
get_chain <- function(id,rank=1L) {
  record <- cache[[id]]
  map <- record$mapping
  count <- sort(table(map$struct_asym_id[!is.na(map$observed) &
                                         map$observed]),decreasing=TRUE)
  names <- names(count)[as.integer(count)>=100L]
  if(length(names)<rank)
    stop(sprintf("Entry %s has no complete chain ranked %d.",id,rank))
  names[[rank]]
}
canonical_ca <- function(id,chain) {
  record <- cache[[id]]
  m <- record$mapping
  m <- m[m$struct_asym_id==chain & !is.na(m$observed) &
           m$observed,c("label_seq_id","uniprot_resi"),drop=FALSE]
  m <- unique(m)
  bad <- duplicated(m$label_seq_id) |
    duplicated(m$label_seq_id,fromLast=TRUE) |
    duplicated(m$uniprot_resi) |
    duplicated(m$uniprot_resi,fromLast=TRUE)
  m <- m[!bad,,drop=FALSE]
  ca <- record$ca_points
  ca <- ca[ca$struct_asym_id==chain,
    c("label_seq_id","x","y","z"),drop=FALSE]
  dup <- duplicated(ca$label_seq_id) |
    duplicated(ca$label_seq_id,fromLast=TRUE)
  ca <- ca[!dup,,drop=FALSE]
  out <- merge(m,ca,by="label_seq_id")
  out <- out[order(out$uniprot_resi),,drop=FALSE]
  if(nrow(out)<100L || anyDuplicated(out$uniprot_resi))
    stop("No complete, unambiguous shared UniProt C-alpha mapping for ",id)
  out
}
canonical_residue_names <- function(id,chain) {
  data <- cache[[id]]$sequence_audit
  data <- data[data$struct_asym_id==chain &
               !is.na(data$uniprot_resi),,drop=FALSE]
  data <- data[!duplicated(data$uniprot_resi) &
               !duplicated(data$uniprot_resi,fromLast=TRUE),,drop=FALSE]
  data
}
pair <- function(case_id,a,b,chain_a,chain_b,kind) {
  aa <- cache[[a]];bb <- cache[[b]]
  if(!identical(aa$accession,bb$accession))
    stop("Cross-UniProt comparison blocked: ",a," and ",b)
  first <- canonical_ca(a,chain_a); second <- canonical_ca(b,chain_b)
  overlap <- sort(intersect(first$uniprot_resi,second$uniprot_resi))
  fraction_a <- length(overlap)/nrow(first)
  fraction_b <- length(overlap)/nrow(second)
  if(length(overlap)<100L || min(fraction_a,fraction_b)<0.60)
    stop(sprintf("Inadequate exact canonical coverage in %s/%s.",a,b))
  pos_a <- as.matrix(first[match(overlap,first$uniprot_resi),
                           c("x","y","z"),drop=FALSE])
  pos_b <- as.matrix(second[match(overlap,second$uniprot_resi),
                            c("x","y","z"),drop=FALSE])
  A <- as.matrix(stats::dist(pos_a))
  B <- as.matrix(stats::dist(pos_b))
  drmsd <- sqrt(mean((A-B)^2)) # uses same residue pair set across compared models
  local_global <- sqrt(rowMeans((A-B)^2))
  torsions_a <- ram_atlas_entity_torsions(aa,chain_a)
  torsions_b <- ram_atlas_entity_torsions(bb,chain_b)
  changes <- ram_atlas_torsion_delta(torsions_a,torsions_b,
      paste(a,chain_a,sep=":"),paste(b,chain_b,sep=":"),30)
  # Only shared complete pairs enter local denominator.
  comparable <- changes$comparable
  change_count <- sum(changes$candidate,na.rm=TRUE)
  if(sum(comparable)<100L)
    stop("Too few complete paired torsions in ",a,"/",b)
  seg <- ram_atlas_candidate_regions(changes)
  names_a <- canonical_residue_names(a,chain_a)
  names_b <- canonical_residue_names(b,chain_b)
  matched_names <- merge(names_a,names_b,by="uniprot_resi",
                          suffixes=c("_a","_b"))
  confirmed <- !is.na(matched_names$residue_name_a) &
               !is.na(matched_names$residue_name_b) &
               nzchar(matched_names$residue_name_a) &
               nzchar(matched_names$residue_name_b)
  mismatches <- sum(matched_names$residue_name_a[confirmed]!=
                    matched_names$residue_name_b[confirmed])
  known_names <- sum(confirmed)
  if(known_names<100L)
    stop("Insufficient cross-entry residue chemistry to audit constructs.")
  fraction_mismatch <- mismatches/known_names
  # A mismatch can reflect a true construct mutation or modified residue.
  # Report all, but do not silently label grossly different constructs same.
  if(fraction_mismatch>0.05)
    stop("Excess residue identity discordance in ",a," / ",b)
  label <- paste(case_id,a,chain_a,b,chain_b,sep="__")
  utils::write.csv(changes,file.path(dest,paste0(label,"_torsions.csv")),
                   row.names=FALSE,na="")
  utils::write.csv(data.frame(uniprot_resi=overlap,
                   distance_map_shift_A=local_global),
                   file.path(dest,paste0(label,"_distance_map.csv")),
                   row.names=FALSE)
  data.frame(protein=case_id,comparison=kind,
    left=paste(a,chain_a,sep=":"),right=paste(b,chain_b,sep=":"),
    positions=length(overlap),coverage_left=round(fraction_a,3),
    coverage_right=round(fraction_b,3),
    drmsd_A=round(drmsd,4),comparable_torsions=sum(comparable),
    angle_changes_30deg=change_count,
    change_fraction=round(change_count/sum(comparable),4),
    candidate_segments=nrow(seg),
    residue_names_checked=known_names,
    residue_name_mismatches=mismatches,
    name_mismatch_fraction=round(fraction_mismatch,4),
    stringsAsFactors=FALSE)
}
all <- list()
for(protein in conf$proteins) {
  entries <- protein$entries
  valid <- vapply(entries,function(x) paste0(x$pdb,"_",x$entity),
                  character(1L))
  names(valid) <- vapply(entries,function(x)x$pdb,character(1L))
  for(cmp in protein$comparisons) {
    a <- valid[[cmp$left]];b <- valid[[cmp$right]]
    if(is.null(a)||is.null(b))
      stop("Invalid case study comparison entry.")
    rank_a <- if(is.null(cmp$left_chain_rank)) 1L
      else as.integer(cmp$left_chain_rank)
    rank_b <- if(is.null(cmp$right_chain_rank)) 1L
      else as.integer(cmp$right_chain_rank)
    row <- pair(protein$id,a,b,get_chain(a,rank_a),
                get_chain(b,rank_b),cmp$class)
    all[[length(all)+1L]] <- row
  }
  # Numerical negative control for each unique protein.
  a <- valid[[1L]]
  chain <- get_chain(a)
  all[[length(all)+1L]] <- pair(protein$id,a,a,chain,chain,
                                "identity_zero_control")
}
results <- do.call(rbind,all)
zero <- results$comparison=="identity_zero_control"
if(sum(zero)!=length(conf$proteins) ||
   any(results$drmsd_A[zero]!=0) ||
   any(results$angle_changes_30deg[zero]!=0L))
  stop("Identity controls must be exact zeros for every protein.")
if(!all(vapply(conf$proteins,function(protein)
  sum(results$protein==protein$id &
      results$comparison=="documented_open_closed")==1L,logical(1L))))
  stop("Every protein must have one documented contrast.")
utils::write.csv(results,file.path(dest,"multiprotein-summary.csv"),
                 row.names=FALSE)
# Contrast report is descriptive. One comparison per protein is not enough to
# learn or evaluate a state classifier or general cutoff.
summary <- c(
  "# Three-protein experimental conformational comparison",
  "",
  "All coordinates are model-1 PDBe updated-mmCIF observations mapped by exact",
  "SIFTS. Sources are hashed and audited in manifest.json.",
  "Labels refer to known experimental forms, not an inferred classifier.",
  "Controls comprise same-crystal chain copies, a separate-crystal",
  "open-state pair and exact identity checks, as available.",
  "",
  "## Per-comparison measurements",
  "",
  paste(capture.output(print(results,row.names=FALSE)),collapse="\n"),
  "",
  "## Scientific limits",
  "",
  "No statistical claim of sensitivity, specificity or general cutoff.",
  "Within-crystal chains are non-independent; periplasmic binding proteins",
  "share related architectures, so the proteins are not a broad fold sample.",
  "Residue-name mismatches and missing coverage must be inspected,",
  "and ligand association here does not establish causality.",
  "Angle changes over 30 degrees are exploratory, not an error or state label.",
  ""
)
writeLines(summary,file.path(dest,"multiprotein-report.md"))
print(results,row.names=FALSE)
cat("Curated multi-protein Atlas live benchmark passed scientific checks.\n")
