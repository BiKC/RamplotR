# Group-level conformational comparison.
#
# Structures are sequence-aligned to a single reference chain. Per-group
# circular phi/psi summaries reuse the same statistics as ensemble analysis.
# Between-group shifts are descriptive effect sizes for navigation; no
# significance test is implied.

ram_alignment_quality <- function(reference, candidate) {
  pairing <- ram_align_residues(reference, candidate)
  both <- !is.na(pairing$index_a) & !is.na(pairing$index_b)
  aligned <- sum(both)
  if (!aligned) return(list(
    identity=0, reference_coverage=0, candidate_coverage=0,
    aligned=0L, score=0
  ))
  ra <- toupper(as.character(reference$resn[pairing$index_a[both]]))
  rb <- toupper(as.character(candidate$resn[pairing$index_b[both]]))
  identity <- mean(ra == rb)
  ref_cov <- aligned / max(1L,nrow(reference))
  cand_cov <- aligned / max(1L,nrow(candidate))
  list(
    identity=identity,
    reference_coverage=ref_cov,
    candidate_coverage=cand_cov,
    aligned=as.integer(aligned),
    score=identity * sqrt(ref_cov*cand_cov)
  )
}

ram_best_chain_match <- function(reference, candidate,
                                 min_identity=0.30,
                                 min_reference_coverage=0.50) {
  if (!nrow(reference) || !nrow(candidate))
    stop("Reference and candidate structures need protein residues.")
  chain_values <- unique(as.character(candidate$chain))
  if (!length(chain_values)) stop("Candidate structure has no protein chain.")
  scored <- lapply(chain_values,function(chain) {
    table <- candidate[as.character(candidate$chain)==chain,,drop=FALSE]
    quality <- ram_alignment_quality(reference,table)
    data.frame(
      chain=chain,
      identity=quality$identity,
      reference_coverage=quality$reference_coverage,
      candidate_coverage=quality$candidate_coverage,
      aligned=quality$aligned,
      score=quality$score,
      stringsAsFactors=FALSE
    )
  })
  scores <- do.call(rbind,scored)
  order_index <- order(-scores$score,-scores$identity,
                       -scores$reference_coverage,scores$chain)
  scores <- scores[order_index,,drop=FALSE]
  best <- scores[1L,,drop=FALSE]
  if (!is.finite(best$identity) || best$identity < min_identity ||
      !is.finite(best$reference_coverage) ||
      best$reference_coverage < min_reference_coverage) {
    stop(sprintf(
      "No candidate chain met the minimum alignment criteria (identity %.0f%%, reference coverage %.0f%%). Best chain %s: %.1f%% identity, %.1f%% reference coverage.",
      100*min_identity,100*min_reference_coverage,best$chain,
      100*best$identity,100*best$reference_coverage
    ))
  }
  list(chain=best$chain, metrics=best, all=scores)
}

ram_map_chain_to_reference <- function(reference, candidate,
                                       source_label,
                                       source_chain=unique(candidate$chain)[[1L]]) {
  pairing <- ram_align_residues(reference,candidate)
  both <- !is.na(pairing$index_a) & !is.na(pairing$index_b)
  pairing <- pairing[both,,drop=FALSE]
  if (!nrow(pairing)) return(reference[0,,drop=FALSE])
  ref_idx <- pairing$index_a
  cand_idx <- pairing$index_b
  required <- c("chain","resi","insertion_code","resn","phi","psi","region")
  if (!all(required %in% names(reference)) ||
      !all(required %in% names(candidate)))
    stop("Mapped structures require residue identifiers, phi/psi and regions.")
  out <- reference[ref_idx,c("chain","resi","insertion_code","resn"),
                   drop=FALSE]
  out$phi <- candidate$phi[cand_idx]
  out$psi <- candidate$psi[cand_idx]
  out$region <- candidate$region[cand_idx]
  optional <- c("rama8000_region","rama8000_group","rama8000_score","plddt")
  for(field in optional)
    if(field %in% names(candidate)) out[[field]] <- candidate[[field]][cand_idx]
  out$source_label <- as.character(source_label)
  out$source_chain <- as.character(source_chain)
  out
}

ram_prepare_structure_group <- function(reference, structures, labels=NULL,
                                        min_identity=0.30,
                                        min_reference_coverage=0.50) {
  if (!is.list(structures) || !length(structures))
    stop("Each comparison group needs at least one structure.")
  if (is.null(labels)) labels <- paste("Structure",seq_along(structures))
  labels <- as.character(labels)
  if (length(labels)!=length(structures) || any(!nzchar(labels)))
    stop("Every group structure needs a non-empty label.")
  mapped <- vector("list",length(structures))
  model_info <- vector("list",length(structures))
  for(i in seq_along(structures)) {
    table <- structures[[i]]
    best <- ram_best_chain_match(
      reference,table,min_identity,min_reference_coverage)
    chain <- table[as.character(table$chain)==best$chain,,drop=FALSE]
    mapped[[i]] <- ram_map_chain_to_reference(
      reference,chain,labels[[i]],best$chain)
    model_info[[i]] <- data.frame(
      model=labels[[i]],
      chain=best$chain,
      identity=best$metrics$identity,
      reference_coverage=best$metrics$reference_coverage,
      candidate_coverage=best$metrics$candidate_coverage,
      aligned=best$metrics$aligned,
      stringsAsFactors=FALSE
    )
  }
  list(models=mapped, model_summary=do.call(rbind,model_info))
}

ram_group_conformation_compare <- function(reference, group_a, group_b,
                                           label_a="Group A",
                                           label_b="Group B") {
  if (!is.list(group_a) || !length(group_a) ||
      !is.list(group_b) || !length(group_b))
    stop("Both groups need at least one mapped structure.")
  a <- ram_ensemble_summary(group_a)
  b <- ram_ensemble_summary(group_b)
  key <- function(data)
    paste(data$chain,data$resi,data$insertion_code,toupper(data$resn),sep="\r")
  ka <- key(a); kb <- key(b)
  all_keys <- unique(c(ka,kb))
  ia <- match(all_keys,ka); ib <- match(all_keys,kb)
  template <- rbind(
    a[!duplicated(ka),c("chain","resi","insertion_code","resn"),drop=FALSE],
    b[!duplicated(kb),c("chain","resi","insertion_code","resn"),drop=FALSE]
  )
  kt <- key(template)
  ids <- template[match(all_keys,kt),,drop=FALSE]
  pick <- function(data,index,field,missing=NA_real_) {
    out <- rep(missing,length(index))
    good <- !is.na(index)
    if(any(good) && field %in% names(data)) out[good] <- data[[field]][index[good]]
    out
  }
  out <- ids
  numeric_fields <- c("phi_models","psi_models","phi_mean","phi_sd",
                      "psi_mean","psi_sd","rama8000_models",
                      "rama8000_consistency","basin_models","basin_consistency")
  char_fields <- c("rama8000_mode","basin_mode")
  for(field in numeric_fields) {
    out[[paste0("a_",field)]] <- pick(a,ia,field)
    out[[paste0("b_",field)]] <- pick(b,ib,field)
  }
  for(field in char_fields) {
    out[[paste0("a_",field)]] <- pick(a,ia,field,NA_character_)
    out[[paste0("b_",field)]] <- pick(b,ib,field,NA_character_)
  }
  out$delta_phi <- ram_angular_difference(out$a_phi_mean,out$b_phi_mean)
  out$delta_psi <- ram_angular_difference(out$a_psi_mean,out$b_psi_mean)
  out$angular_displacement <- ram_backbone_angular_displacement(
    out$delta_phi,out$delta_psi)
  out$shift_band <- ram_backbone_shift_band(out$angular_displacement)
  max_finite <- function(...) {
    values <- cbind(...)
    apply(values,1L,function(row) {
      row <- row[is.finite(row)]
      if(!length(row)) NA_real_ else max(row)
    })
  }
  out$max_within_group_sd <- max_finite(
    out$a_phi_sd,out$a_psi_sd,out$b_phi_sd,out$b_psi_sd)
  n_a <- length(group_a)
  n_b <- length(group_b)
  out$a_coverage <- pmin(out$a_phi_models,out$a_psi_models) / max(1L,n_a)
  out$b_coverage <- pmin(out$b_phi_models,out$b_psi_models) / max(1L,n_b)
  out$min_group_coverage <- pmin(out$a_coverage,out$b_coverage)
  out$min_rama8000_consistency <- pmin(
    out$a_rama8000_consistency,out$b_rama8000_consistency,na.rm=TRUE)
  out$min_rama8000_consistency[
    !is.finite(out$a_rama8000_consistency) |
    !is.finite(out$b_rama8000_consistency)] <- NA_real_

  enough <- out$a_phi_models>=2 & out$a_psi_models>=2 &
            out$b_phi_models>=2 & out$b_psi_models>=2
  out$consistent_shift <- enough &
    is.finite(out$angular_displacement) &
    out$angular_displacement>=30 &
    is.finite(out$max_within_group_sd) &
    out$max_within_group_sd<=15
  out$high_support_shift <- out$consistent_shift &
    is.finite(out$min_group_coverage) & out$min_group_coverage>=0.75 &
    (is.na(out$min_rama8000_consistency) |
      out$min_rama8000_consistency>=0.75)
  out$rama8000_mode_changed <- !is.na(out$a_rama8000_mode) &
    !is.na(out$b_rama8000_mode) &
    out$a_rama8000_mode != out$b_rama8000_mode
  out$basin_mode_changed <- !is.na(out$a_basin_mode) &
    !is.na(out$b_basin_mode) &
    out$a_basin_mode != out$b_basin_mode

  out$evidence_profile <- ifelse(
    !is.finite(out$angular_displacement),"Unavailable",
    ifelse(out$min_group_coverage<0.75,"Sparse coverage",
      ifelse(out$angular_displacement>=30 &
               is.finite(out$max_within_group_sd) &
               out$max_within_group_sd<=15,
             "Low-dispersion shift",
        ifelse(out$angular_displacement>=30,"Large but variable",
          ifelse(out$angular_displacement>=15,"Moderate shift","Small shift")))))
  out$group_a <- label_a
  out$group_b <- label_b
  out[order(-as.integer(out$basin_mode_changed),
            -as.integer(out$high_support_shift),
            -as.integer(out$consistent_shift),
            -replace(out$angular_displacement,
                     !is.finite(out$angular_displacement),-Inf),
            out$chain,out$resi,out$insertion_code),,drop=FALSE]
}
