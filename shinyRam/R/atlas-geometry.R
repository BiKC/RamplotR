# Exploratory experimental geometry groups in exact UniProt coordinates.
# Distances compare internal C-alpha pair distances: rigid-body invariant, but
# sensitive to domain movement. The group cutoff is a user-selected heuristic,
# never a validated functional-state or protein population assignment.
#
# Input: on-demand PDBe updated-mmCIF SIFTS verifications, each including
# label_asym_id + label_seq_id C-alpha coordinates for model 1.

ram_atlas_ca_records <- function(records) {
  if(is.null(records)) return(data.frame(
    struct_asym_id=character(),label_seq_id=integer(),
    x=numeric(),y=numeric(),z=numeric()))
  if(is.data.frame(records)) return(records)
  if(!is.list(records)) stop("Invalid coordinate records.",call.=FALSE)
  if(!length(records)) return(ram_atlas_ca_records(NULL))
  if(length(records)>5000L)
    stop("Too many mapped C-alpha coordinates.",call.=FALSE)
  do.call(rbind,lapply(records,function(row) {
    if(!is.list(row)) stop("Malformed C-alpha coordinate.",call.=FALSE)
    data.frame(struct_asym_id=as.character(row$struct_asym_id),
      label_seq_id=suppressWarnings(as.integer(row$label_seq_id)),
      x=suppressWarnings(as.numeric(row$x)),
      y=suppressWarnings(as.numeric(row$y)),
      z=suppressWarnings(as.numeric(row$z)))
  }))
}

ram_atlas_geometry_entities <- function(verified,accession) {
  accession <- ram_uniprot_accession(accession)
  entities <- list()
  if(!length(verified)) return(entities)
  for(key in names(verified)) {
    item <- verified[[key]]
    if(!is.list(item) || !identical(item$state,"mapped")) next
    if(!grepl("^[A-Z0-9]{4}_[1-9][0-9]*$",key)) next
    map <- item$mapping
    if(!is.data.frame(map) || !nrow(map) ||
       !all(c("uniprot_accession","uniprot_resi","label_seq_id",
              "struct_asym_id","observed","chain","resi",
              "insertion_code") %in% names(map))) next
    if(any(is.na(map$uniprot_accession)) ||
       any(map$uniprot_accession != accession))
      stop("Structure mapping contains another UniProt accession.",call.=FALSE)
    ca <- ram_atlas_ca_records(item$ca_points)
    required <- c("struct_asym_id","label_seq_id","x","y","z")
    if(!all(required %in% names(ca)) || !nrow(ca)) next
    ca$struct_asym_id <- as.character(ca$struct_asym_id)
    ca$label_seq_id <- suppressWarnings(as.integer(ca$label_seq_id))
    for(col in c("x","y","z")) ca[[col]] <- suppressWarnings(as.numeric(ca[[col]]))
    good <- is.finite(ca$label_seq_id) & ca$label_seq_id>0L &
      is.finite(ca$x) & is.finite(ca$y) & is.finite(ca$z)
    ca <- ca[good,,drop=FALSE]
    if(!nrow(ca)) next
    # Ambiguous coordinates are excluded, not assigned to the first record.
    ca_key <- paste(ca$struct_asym_id,ca$label_seq_id,sep="\r")
    ca <- ca[!duplicated(ca_key) &
      !duplicated(ca_key,fromLast=TRUE),,drop=FALSE]

    map <- map[!is.na(map$observed) & map$observed &
      is.finite(map$label_seq_id) & is.finite(map$uniprot_resi) &
      map$uniprot_resi>0L,,drop=FALSE]
    if(!nrow(map)) next
    map_key <- paste(map$struct_asym_id,map$label_seq_id,sep="\r")
    target <- as.character(map$uniprot_resi)
    # Remove all records with conflicting label->UniProt OR target->label
    # assignments, including insertion alternatives.
    conflict_label <- duplicated(map_key) | duplicated(map_key,fromLast=TRUE)
    # Perfectly duplicated records should not be considered conflicting.
    map <- unique(map)
    map_key <- paste(map$struct_asym_id,map$label_seq_id,sep="\r")
    target <- paste(map$struct_asym_id,map$uniprot_resi,sep="\r")
    bad <- duplicated(map_key) | duplicated(map_key,fromLast=TRUE) |
      duplicated(target) | duplicated(target,fromLast=TRUE)
    map <- map[!bad,,drop=FALSE]
    if(!nrow(map)) next
    idx <- match(paste(map$struct_asym_id,map$label_seq_id,sep="\r"),
                 paste(ca$struct_asym_id,ca$label_seq_id,sep="\r"))
    map <- map[!is.na(idx),,drop=FALSE]
    idx <- idx[!is.na(idx)]
    if(!nrow(map)) next
    combined <- data.frame(struct_asym_id=as.character(map$struct_asym_id),
      chain=as.character(map$chain),uniprot_resi=as.integer(map$uniprot_resi),
      x=ca$x[idx],y=ca$y[idx],z=ca$z[idx])
    groups <- split(combined,combined$struct_asym_id)
    # A polymer entity can appear in multiple copies; use one complete chain,
    # deterministically, to avoid treating symmetry copies as independent states.
    lengths <- vapply(groups,nrow,integer(1L))
    best <- names(groups)[order(-lengths,names(groups))][[1L]]
    selected <- groups[[best]]
    selected <- selected[order(selected$uniprot_resi),,drop=FALSE]
    entities[[key]] <- list(id=key,struct_asym_id=best,
      chain=selected$chain[[1L]],coordinates=selected)
  }
  entities
}

# A conservative *exploratory* automatic cut of the SAME canonical-core
# dRMSD matrix already used by the manual dendrogram. Average silhouette
# measures geometric cohesion/separation, not biological state validity.
ram_atlas_auto_clusters <- function(distance_matrix,hc,
  min_silhouette=0.50,min_separation=0.35,max_groups=4L) {
  D <- as.matrix(distance_matrix)
  n <- nrow(D)
  if(n<2L || n!=ncol(D) || any(!is.finite(D)) ||
     any(D < -1e-9) || any(abs(D-t(D))>1e-8) ||
     any(abs(diag(D))>1e-8) || !inherits(hc,"hclust"))
    stop("Invalid canonical-core distance matrix for automatic grouping.",
      call.=FALSE)
  groups <- rep(1L,n)
  names(groups) <- rownames(D)
  candidate <- data.frame(k=integer(),silhouette=numeric(),
    separation_A=numeric(),singleton_groups=integer(),
    stringsAsFactors=FALSE)
  reason <- "Fewer than three structures; no automatic split recommended."
  if(n>=3L) {
    for(k in seq.int(2L,min(max_groups,n-1L))) {
      labels <- stats::cutree(hc,k=k)
      memberships <- split(seq_len(n),labels)
      singleton <- sum(lengths(memberships)==1L)
      scores <- numeric(n)
      inside <- numeric()
      outside <- numeric()
      for(i in seq_len(n)) {
        own <- which(labels==labels[[i]] & seq_len(n)!=i)
        other <- which(labels!=labels[[i]])
        if(length(own)) {
          a <- mean(D[i,own])
          nearest <- min(vapply(split(other,labels[other]),
            function(index) mean(D[i,index]),numeric(1L)))
          scores[[i]] <- if(max(a,nearest)>0)
            (nearest-a)/max(a,nearest) else 0
        } else scores[[i]] <- 0 # Singleton silhouette is not 1.
      }
      for(j in seq_len(n-1L)) for(i in seq.int(j+1L,n)) {
        if(labels[[i]]==labels[[j]])
          inside <- c(inside,D[i,j])
        else outside <- c(outside,D[i,j])
      }
      separation <- if(length(inside) && length(outside))
        stats::median(outside)-stats::median(inside)
      else 0
      candidate <- rbind(candidate,data.frame(k=as.integer(k),
        silhouette=round(mean(scores),4),
        separation_A=round(separation,4),
        singleton_groups=as.integer(singleton)))
    }
    accepted <- which(candidate$silhouette>=min_silhouette &
      candidate$separation_A>=min_separation)
    if(length(accepted)) {
      # A small complexity penalty avoids gratuitous subdivision when
      # silhouette scores are close. Exact ties favor fewer groups.
      rank <- candidate$silhouette[accepted] -
        0.02*(candidate$k[accepted]-2L)
      chosen <- accepted[[which.max(rank)]]
      groups <- stats::cutree(hc,k=candidate$k[[chosen]])
      reason <- sprintf(paste0(
        "Suggested %d geometric groups (mean silhouette %.2f; ",
        "median between-minus-within distance %.2f Å)."),
        candidate$k[[chosen]],candidate$silhouette[[chosen]],
        candidate$separation_A[[chosen]])
      if(candidate$singleton_groups[[chosen]]>0)
        reason <- paste0(reason,
          " Some groups are singletons; within-group consistency cannot be assessed.")
    } else reason <- paste0(
      "No sufficiently separated grouping at the exploratory safeguards ",
      sprintf("(mean silhouette ≥%.2f and median distance gap ≥%.2f Å). ",
        min_silhouette,min_separation),
      "Treat the cohort as one geometric group or choose a manual cutoff.")
  }
  k <- length(unique(groups))
  heights <- hc$height
  cutoff <- if(k==1L) max(heights)+max(0.01,0.02*max(heights)) else {
    bottom <- heights[[n-k]]
    top <- heights[[n-k+1L]]
    if(top>bottom) (bottom+top)/2 else top
  }
  list(assignment=as.integer(groups),k=k,cutoff=as.numeric(cutoff),
    reason=reason,candidates=candidate,mode="automatic",
    min_silhouette=min_silhouette,min_separation_A=min_separation)
}

ram_atlas_geometry_groups <- function(verified,accession,selected,
  min_core=30L,min_fraction=0.6,cutoff=1.5,max_core=300L,
  cluster_mode=c("manual","automatic")) {
  cluster_mode <- match.arg(cluster_mode)
  ids <- unique(as.character(selected))
  if(length(ids)<2L || length(ids)>12L ||
     anyNA(ids) || any(!nzchar(ids)))
    stop("Select between 2 and 12 verified experimental entities.",
         call.=FALSE)
  cutoff <- as.numeric(cutoff)
  if(length(cutoff)!=1L || !is.finite(cutoff) ||
     cutoff<=0 || cutoff>10)
    stop("Invalid geometric distance cutoff.",call.=FALSE)
  min_core <- as.integer(min_core)
  max_core <- as.integer(max_core)
  available <- ram_atlas_geometry_entities(verified,accession)
  missing <- setdiff(ids,names(available))
  if(length(missing))
    stop(paste("Missing exact C-alpha coordinates for:",
               paste(missing,collapse=", ")),call.=FALSE)
  data <- available[ids]
  counts <- vapply(data,function(x)nrow(x$coordinates),integer(1L))
  core <- Reduce(intersect,lapply(data,function(x)x$coordinates$uniprot_resi))
  core <- sort(unique(as.integer(core)))
  if(length(core)<min_core)
    stop(sprintf("Only %d common observed UniProt C-alpha positions; at least %d required.",
                 length(core),min_core),call.=FALSE)
  fractions <- length(core)/counts
  if(any(fractions<min_fraction))
    stop(sprintf(paste0("Insufficient common-core coverage (minimum %.0f%% ",
      "of each selected chain); choose more comparable constructs."),
      min_fraction*100),call.=FALSE)
  sampled <- core
  if(length(core)>max_core)
    sampled <- core[unique(as.integer(round(seq(1,length(core),
                                                 length.out=max_core))))]
  # One fixed position set is used for every pair, so the dRMSD distances
  # are comparable and clustering cannot be driven by differing coverage.
  matrices <- lapply(data,function(x) {
    points <- x$coordinates[match(sampled,x$coordinates$uniprot_resi),,
      drop=FALSE]
    as.matrix(points[,c("x","y","z"),drop=FALSE])
  })
  dists <- lapply(matrices,function(coords) as.numeric(
    stats::dist(coords)))
  n <- length(ids)
  D <- matrix(0,nrow=n,ncol=n,dimnames=list(ids,ids))
  for(i in seq_len(n-1L))
    for(j in seq.int(i+1L,n))
      D[i,j] <- D[j,i] <- sqrt(mean((dists[[i]]-dists[[j]])^2))
  hc <- stats::hclust(stats::as.dist(D),method="average")
  automatic <- if(identical(cluster_mode,"automatic"))
    ram_atlas_auto_clusters(D,hc) else NULL
  if(!is.null(automatic)) {
    labels <- stats::setNames(automatic$assignment,ids)
    cutoff <- automatic$cutoff
  } else labels <- stats::cutree(hc,h=cutoff)
  groups <- split(ids,labels)
  representatives <- vapply(groups,function(members) {
    idx <- match(members,ids)
    scores <- vapply(idx,function(i)mean(D[i,idx]),numeric(1L))
    members[order(scores,members)][[1L]]
  },character(1L))
  list(accession=ram_uniprot_accession(accession),
    selected=ids,counts=counts,common_positions=core,
    sampled_positions=sampled,common_fraction=fractions,
    distance_matrix=D,hclust=hc,
    assignment=data.frame(entity=ids,
      chain=vapply(data,function(x)x$chain,character(1L)),
      mapped_ca=unname(counts),
      common_coverage=round(unname(fractions),3),
      geometric_group=as.integer(labels[ids]),
      stringsAsFactors=FALSE),
    representatives=representatives,cutoff=cutoff,
    cluster_mode=cluster_mode,auto_cluster=automatic,
    distance_method="C-alpha intrachain distance-map RMSD (Å), common canonical UniProt positions; average-linkage hierarchical clustering")
}

# Graphics device independent; 2-entry cohorts are a single distance and do
# not need (and should not rely upon) dendrogram conversion.
ram_atlas_geometry_plot <- function(result) {
  if(!is.list(result) || is.null(result$distance_matrix) ||
     is.null(result$selected) || is.null(result$cutoff))
    stop("Invalid Atlas geometry result.",call.=FALSE)
  if(length(result$selected)==2L) {
    distance <- result$distance_matrix[1L,2L]
    graphics::barplot(distance,names.arg=paste(result$selected,
      collapse=" vs "),ylab="C-alpha distance-map RMSD (Å)",
      main="Pairwise experimental geometry",
      ylim=c(0,max(distance,result$cutoff)*1.2))
  } else {
    graphics::plot(result$hclust,main="Experimental geometry similarity",
      xlab="",sub="",ylab="C-alpha distance-map RMSD (Å)")
  }
  graphics::abline(h=result$cutoff,lty=2,col="gray50")
  invisible(result)
}
