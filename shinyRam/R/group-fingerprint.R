# Per-member residue fingerprints for researcher-defined structure groups.
# This is descriptive backbone evidence, not a test of biological states.
# A model contributes a complete observation only if phi AND psi are finite.

ram_group_fingerprint <- function(group_a,group_b,
                                  label_a="Group A",label_b="Group B") {
  build <- function(models,label) {
    if(!is.list(models) || !length(models))
      stop("A fingerprint requires at least one model per group.",call.=FALSE)
    ids <- names(models)
    if(is.null(ids)) ids <- rep("",length(models))
    ids[is.na(ids) | !nzchar(ids)] <-
      paste0("Structure ",which(is.na(ids) | !nzchar(ids)))
    ids <- make.unique(ids,sep=" #")
    lapply(seq_along(models),function(i) {
      model <- models[[i]]
      required <- c("chain","resi","insertion_code","resn","phi","psi")
      if(!is.data.frame(model) || !all(required %in% names(model)))
        stop("Fingerprint members require mapped residue identities and phi/psi.",call.=FALSE)
      key <- paste(model$chain,model$resi,model$insertion_code,sep="\r")
      if(anyDuplicated(key))
        stop("Fingerprint member contains ambiguous residue identifiers.",call.=FALSE)
      phi <- as.numeric(model$phi); psi <- as.numeric(model$psi)
      paired <- is.finite(phi) & is.finite(psi)
      basin <- ram_backbone_basin(phi,psi)
      # Uploaded groups use reference residue identity for alignment. Keep the
      # candidate's actual residue identity separately where it is available.
      identity <- if("source_resn" %in% names(model))
        as.character(model$source_resn) else as.character(model$resn)
      data.frame(group=label,member=ids[[i]],
        chain=as.character(model$chain),resi=as.integer(model$resi),
        insertion_code=as.character(model$insertion_code),
        amino_acid=identity,phi=phi,psi=psi,
        paired=paired,backbone_state=basin,
        stringsAsFactors=FALSE)
    })
  }
  if(!is.character(label_a) || length(label_a)!=1L ||
     is.na(label_a) || !nzchar(trimws(label_a)) ||
     !is.character(label_b) || length(label_b)!=1L ||
     is.na(label_b) || !nzchar(trimws(label_b)) ||
     identical(label_a,label_b))
    stop("Fingerprint groups require distinct nonempty labels.",call.=FALSE)
  records <- do.call(rbind,c(build(group_a,label_a),build(group_b,label_b)))
  rownames(records) <- NULL
  records[order(records$chain,records$resi,records$insertion_code,
                records$group,records$member),,drop=FALSE]
}

ram_group_fingerprint_at <- function(records,chain,resi,insertion_code="") {
  if(!is.data.frame(records) ||
     !all(c("chain","resi","insertion_code","paired","phi","psi") %in%
          names(records)) ||
     length(chain)!=1L || length(resi)!=1L ||
     length(insertion_code)!=1L)
    stop("Invalid fingerprint residue selection.",call.=FALSE)
  matched <- !is.na(records$chain) & records$chain==as.character(chain) &
    !is.na(records$resi) & records$resi==as.integer(resi) &
    !is.na(records$insertion_code) &
    records$insertion_code==as.character(insertion_code)
  records[matched,,drop=FALSE]
}

ram_group_fingerprint_summary <- function(records,group_names,group_sizes) {
  stopifnot(length(group_names)==2L,length(group_sizes)==2L)
  do.call(rbind,lapply(seq_along(group_names),function(i) {
    group <- group_names[[i]]
    sub <- records[records$group==group,,drop=FALSE]
    paired <- sub[sub$paired %in% TRUE,,drop=FALSE]
    states <- paired$backbone_state[!is.na(paired$backbone_state)]
    mode <- if(length(states)) {
      votes <- sort(table(states),decreasing=TRUE)
      names(votes)[[1L]]
    } else NA_character_
    agreement <- if(length(states)) max(table(states))/length(states)
      else NA_real_
    data.frame(group=group,complete_pairs=nrow(paired),
      members=as.integer(group_sizes[[i]]),
      modal_state=mode,consensus=agreement,
      phi_mean=ram_ensemble_circular(paired$phi)[["mean"]],
      psi_mean=ram_ensemble_circular(paired$psi)[["mean"]],
      phi_sd=ram_ensemble_circular(paired$phi)[["sd"]],
      psi_sd=ram_ensemble_circular(paired$psi)[["sd"]],
      stringsAsFactors=FALSE)
  }))
}

ram_group_fingerprint_plot <- function(records,group_names,contacts=NULL) {
  if(!is.data.frame(records) || length(group_names)!=2L)
    stop("Invalid fingerprint plot data.",call.=FALSE)
  colours <- c("#CE6A4D","#317E9A")
  graphics::plot(NA_real_,NA_real_,xlim=c(-180,180),ylim=c(-180,180),
    xlab=expression(phi~"(degrees)"),ylab=expression(psi~"(degrees)"),
    main="Observed backbone angles",asp=1)
  graphics::grid(nx=6,ny=6,col="#E6EDEF",lty="dotted")
  legend <- character(2L)
  for(i in seq_along(group_names)) {
    data <- records[!is.na(records$group) &
      records$group==group_names[[i]] & records$paired %in% TRUE,,
      drop=FALSE]
    legend[[i]] <- sprintf("%s (%d pairs)",group_names[[i]],nrow(data))
    if(!nrow(data)) next
    graphics::points(data$phi,data$psi,pch=if(i==1L) 16 else 17,
      col=grDevices::adjustcolor(colours[[i]],alpha.f=0.75),cex=1.25)
    # Outline models with a measured deposited nonwater component close
    # to this exact canonical residue; keep group identity in point shape.
    if(is.data.frame(contacts) &&
       all(c("member","evidence") %in% names(contacts))) {
      present <- data$member %in% contacts$member[
        contacts$evidence=="Deposited proximity observed"]
      if(any(present)) graphics::points(data$phi[present],data$psi[present],
        pch=1,col=colours[[i]],cex=1.95,lwd=1.8)
    }
    center <- c(ram_ensemble_circular(data$phi)[["mean"]],
                ram_ensemble_circular(data$psi)[["mean"]])
    if(all(is.finite(center)))
      graphics::points(center[[1L]],center[[2L]],pch=4,
        col=colours[[i]],cex=1.7,lwd=2)
  }
  graphics::legend("bottomleft",legend=legend,pch=c(16,17),
    col=colours,bty="n",cex=0.8)
  if(is.data.frame(contacts) &&
     any(contacts$evidence=="Deposited proximity observed",na.rm=TRUE))
    graphics::legend("bottomright",legend="ring = nearby deposited component",
      pch=1,col="#3B5267",bty="n",cex=0.7)
  invisible(records)
}
