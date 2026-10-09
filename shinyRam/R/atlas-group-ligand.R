# Join a selected canonical residue's experimental backbone fingerprint to
# deposited, measured non-water contact evidence. No biological-state inference.
ram_atlas_group_ligand_at <- function(result,resi) {
  if(!is.list(result) || !identical(result$source,"atlas") ||
     !is.data.frame(result$members) ||
     !is.data.frame(result$fingerprint) ||
     !is.list(result$ligand_context) ||
     !is.data.frame(result$ligand_context$entries) ||
     !is.data.frame(result$ligand_context$residue_contacts))
    stop("A verified Atlas group with contact evidence is required.",
      call.=FALSE)
  resi <- suppressWarnings(as.integer(resi))
  if(length(resi)!=1L || is.na(resi) || resi<1L)
    stop("Choose a positive exact UniProt position.",call.=FALSE)
  members <- result$members
  if(!all(c("model","group") %in% names(members)) ||
     anyNA(members$model) || anyDuplicated(as.character(members$model)))
    stop("Ambiguous Atlas group membership.",call.=FALSE)
  entries <- result$ligand_context$entries
  if(anyDuplicated(entries$entity) ||
     any(!members$model %in% entries$entity))
    stop("Ligand evidence does not match group membership.",call.=FALSE)
  sites <- result$ligand_context$residue_contacts
  required <- c("entity","uniprot_resi","comp_id","min_distance_A","scope")
  if(!all(required %in% names(sites)))
    stop("Invalid residue-specific deposited contact evidence.",call.=FALSE)
  out <- lapply(seq_len(nrow(members)),function(i) {
    id <- as.character(members$model[[i]])
    entry <- entries[entries$entity==id,,drop=FALSE]
    if(nrow(entry)!=1L) stop("Unmatched Atlas ligand entry.",call.=FALSE)
    obs <- result$fingerprint[result$fingerprint$member==id &
      result$fingerprint$resi==resi &
      result$fingerprint$chain=="UniProt",,drop=FALSE]
    if(nrow(obs)>1L) stop("Duplicate experimental torsion observation.",call.=FALSE)
    hits <- sites[sites$entity==id & !is.na(sites$uniprot_resi) &
      sites$uniprot_resi==resi,,drop=FALSE]
    all <- sites[sites$entity==id,,drop=FALSE]
    complete <- if(nrow(obs)) isTRUE(obs$paired[[1L]]) else FALSE
    status <- if(!identical(entry$status[[1L]],"measured"))
        "Evidence unavailable"
      else if(nrow(hits)) "Deposited proximity observed"
      else if(any(all$scope=="nearest-only"))
        "Incomplete nearest-only evidence"
      else "No deposited proximity reported"
    if(identical(entry$status[[1L]],"measured") && nrow(hits) &&
       any(!is.finite(hits$min_distance_A) |
           hits$min_distance_A<0 | hits$min_distance_A>4.5))
      stop("Invalid measured contact distance.",call.=FALSE)
    data.frame(group=as.character(members$group[[i]]),member=id,
      uniprot_resi=resi,phi=if(nrow(obs)) obs$phi[[1L]] else NA_real_,
      psi=if(nrow(obs)) obs$psi[[1L]] else NA_real_,
      complete_backbone_pair=complete,
      evidence=status,
      nearby_components=if(identical(entry$status[[1L]],"measured") &&
                           status!="Incomplete nearest-only evidence")
        nrow(hits) else NA_integer_,
      component_codes=if(nrow(hits))
        paste(sort(unique(as.character(hits$comp_id))),collapse=", ")
        else "",
      minimum_distance_A=if(nrow(hits))
        min(hits$min_distance_A) else NA_real_,
      scope=if(nrow(hits)) paste(sort(unique(hits$scope)),collapse=", ")
        else if(any(all$scope=="nearest-only")) "nearest-only"
        else if(identical(entry$status[[1L]],"measured")) "all-mapped-contacts"
        else "unavailable",
      stringsAsFactors=FALSE)
  })
  rows <- do.call(rbind,out)
  rownames(rows) <- NULL
  rows
}
