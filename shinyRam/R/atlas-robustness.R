# Sensitivity of exploratory Atlas clusters to canonical-core position choice.
# Deterministic contiguous-block deletion of the SAME common, sampled
# UniProt-position matrix. This is not a biological bootstrap and does not
# establish experimental independence or state probabilities.

ram_atlas_cluster_robustness <- function(matrices, ids, positions,
  baseline_labels, cluster_mode=c("automatic","manual"), cutoff=1.5,
  blocks=10L) {
  cluster_mode <- match.arg(cluster_mode)
  ids <- as.character(ids)
  n <- length(ids)
  if(n<2L || n>12L || anyNA(ids) || any(!nzchar(ids)) ||
     anyDuplicated(ids) || !is.list(matrices) || length(matrices)!=n)
    stop("Invalid verified experimental matrices for sensitivity analysis.",
      call.=FALSE)
  positions <- as.integer(positions)
  count <- length(positions)
  if(count<30L || anyNA(positions) || anyDuplicated(positions) ||
     is.unsorted(positions,strictly=TRUE))
    stop("Sensitivity analysis requires a sorted exact canonical core.",
      call.=FALSE)
  clean <- lapply(matrices,function(x) {
    x <- as.matrix(x)
    if(!is.numeric(x) || nrow(x)!=count || ncol(x)!=3L ||
       any(!is.finite(x)))
      stop("Sensitivity matrices need complete, aligned C-alpha XYZ.",
           call.=FALSE)
    x
  })
  baseline <- as.character(baseline_labels)
  if(length(baseline)!=n || anyNA(baseline) ||
     any(!nzchar(baseline)))
    stop("Sensitivity analysis requires one baseline group per entity.",
      call.=FALSE)
  if(!is.null(names(baseline_labels))) {
    if(!setequal(names(baseline_labels),ids))
      stop("Baseline groups must name the same verified entities.",
        call.=FALSE)
    baseline <- as.character(baseline_labels[ids])
  }
  names(baseline) <- ids
  cutoff <- as.numeric(cutoff)
  if(length(cutoff)!=1L || !is.finite(cutoff) ||
     cutoff<=0 || cutoff>10)
    stop("Invalid manual grouping cutoff.",call.=FALSE)
  blocks <- as.integer(blocks)
  if(length(blocks)!=1L || is.na(blocks) ||
     blocks<2L || blocks>20L)
    stop("Use between 2 and 20 position-deletion blocks.",call.=FALSE)

  group_ids <- unique(baseline)
  group_members <- lapply(group_ids,function(g) ids[baseline==g])
  names(group_members) <- group_ids
  pairs <- utils::combn(seq_len(n),2L)
  same_group <- baseline[pairs[1L,]]==baseline[pairs[2L,]]
  same_entry <- substr(ids[pairs[1L,]],1L,4L)==
    substr(ids[pairs[2L,]],1L,4L)
  shared_entry_total <- sum(same_entry)

  # With two entries an automatic split is intentionally forbidden, and
  # leave-position-out agreement would be trivially reassuring.
  if(n<3L) return(list(status="insufficient_structures",
    reason="Two experimental entities cannot establish cluster robustness.",
    iterations=0L,replicates=data.frame(),groups=data.frame(),
    pairs=data.frame(),fraction_unchanged=NA_real_,
    shared_pdb_pairs=shared_entry_total,
    method="Contiguous-block jackknife on exact observed UniProt positions"))

  k <- min(blocks,count %/% 3L)
  # cut() gives contiguous blocks of sorted, sampled UniProt positions,
  # including when the observed canonical sequence contains gaps.
  deleted_block <- cut(seq_len(count),breaks=k,labels=FALSE)
  runs <- vector("list",k)
  replay <- vector("list",k)
  for(b in seq_len(k)) {
    drop <- which(deleted_block==b)
    keep <- which(deleted_block!=b)
    pair_dist <- lapply(clean,function(coords)
      as.numeric(stats::dist(coords[keep,,drop=FALSE])))
    D <- matrix(0,nrow=n,ncol=n,dimnames=list(ids,ids))
    for(i in seq_len(n-1L)) for(j in seq.int(i+1L,n)) {
      value <- sqrt(mean((pair_dist[[i]]-pair_dist[[j]])^2))
      D[i,j] <- D[j,i] <- value
    }
    hc <- stats::hclust(stats::as.dist(D),method="average")
    labels <- if(identical(cluster_mode,"automatic"))
      ram_atlas_auto_clusters(D,hc)$assignment
      else stats::cutree(hc,h=cutoff)
    labels <- as.character(labels)
    together <- labels[pairs[1L,]]==labels[pairs[2L,]]
    replay[[b]] <- list(labels=labels,pairs=together)
    runs[[b]] <- data.frame(
      deleted_block=b,removed_start=positions[min(drop)],
      removed_end=positions[max(drop)],removed_positions=length(drop),
      remaining_positions=length(keep),
      resulting_groups=length(unique(labels)),
      baseline_partition_reproduced=all(together==same_group),
      pair_agreement=mean(together==same_group),
      stringsAsFactors=FALSE)
  }
  runs <- do.call(rbind,runs)
  group_table <- do.call(rbind,lapply(seq_along(group_members),function(i) {
    members <- group_members[[i]]
    unchanged <- vapply(replay,function(run) {
      # Compare exact member sets, never numeric cluster label values.
      any(vapply(unique(run$labels),function(label)
        setequal(ids[run$labels==label],members),logical(1L)))
    },logical(1L))
    index <- which(ids %in% members)
    group_pairs <- which(pairs[1L,] %in% index &
                         pairs[2L,] %in% index)
    data.frame(geometry_group=names(group_members)[[i]],
      size=length(members),members=paste(members,collapse=", "),
      exact_group_recovered=sum(unchanged),
      iterations=k,exact_group_fraction=mean(unchanged),
      within_group_same_pdb_pairs=sum(same_entry[group_pairs]),
      stringsAsFactors=FALSE)
  }))
  pair_table <- data.frame(entity_a=ids[pairs[1L,]],
    entity_b=ids[pairs[2L,]],
    baseline_same_group=as.logical(same_group),
    leave_block_out_same_group=vapply(seq_len(ncol(pairs)),function(i)
      sum(vapply(replay,function(run)run$pairs[[i]],logical(1L))),
      integer(1L)),
    iterations=k,same_group_fraction=vapply(seq_len(ncol(pairs)),function(i)
      mean(vapply(replay,function(run)run$pairs[[i]],logical(1L))),
      numeric(1L)),
    same_pdb_entry=as.logical(same_entry),stringsAsFactors=FALSE)
  list(status="ok",reason=sprintf(paste0(
    "%d/%d contiguous-block omissions reproduced the full ",
    "geometry partition."),sum(runs$baseline_partition_reproduced),k),
    iterations=k,replicates=runs,groups=group_table,pairs=pair_table,
    fraction_unchanged=mean(runs$baseline_partition_reproduced),
    shared_pdb_pairs=shared_entry_total,
    method=paste(
      "Deterministic leave-one-contiguous-block-out of sorted sampled",
      "exact UniProt C-alpha coordinates; recomputed distance-map RMSD,",
      "average linkage and original automatic or manual grouping settings.",
      "Partition comparison ignores arbitrary cluster label numbers.",
      "This is sensitivity to sequence-position selection, not statistical",
      "confidence, a biological state probability, or independent replication."))
}
