# Pure residue inspection and pairwise comparison helpers.
# Nothing in this file changes torsion extraction or scientific classification.

ram_review_queue <- function(data, boundary_margin = 2) {
  if (!nrow(data)) return(data)
  missing <- is.na(data$phi) | is.na(data$psi) | is.na(data$region)
  outlier <- !missing & data$region == "Not allowed"
  # The density percentile is cumulative mass above the residue's density.
  # Closeness to any original contour probability is a review hint only.
  near <- !missing & !is.na(data$density) &
    vapply(data$density, function(value) {
      any(abs(value - c(85, 98, 99.95)) <= boundary_margin)
    }, logical(1))
  data$review_status <- ifelse(missing, "Missing angles",
    ifelse(outlier, "Outlier", ifelse(near, "Near boundary", "Other")))
  priority <- match(data$review_status,
                    c("Outlier", "Missing angles", "Near boundary", "Other"))
  data[order(priority, data$chain, data$resi, data$insertion_code), ,
       drop = FALSE]
}

ram_amino_acid_letters <- c(
  ALA="A", ARG="R", ASN="N", ASP="D", CYS="C", GLU="E", GLN="Q",
  GLY="G", HIS="H", ILE="I", LEU="L", LYS="K", MET="M", PHE="F",
  PRO="P", SER="S", THR="T", TRP="W", TYR="Y", VAL="V"
)

ram_sequence_data <- function(data) {
  if (!nrow(data)) return(data.frame(
    chain=character(), resi=integer(), insertion_code=character(),
    resn=character(), letter=character(), region=character(), plddt=numeric(),
    stringsAsFactors=FALSE))
  key <- paste(data$chain, data$resi, data$insertion_code, sep="\r")
  data <- data[!duplicated(key), , drop=FALSE]
  out <- data.frame(
    chain=as.character(data$chain),
    resi=as.integer(data$resi),
    insertion_code=as.character(data$insertion_code),
    resn=as.character(data$resn),
    letter=unname(ram_amino_acid_letters[data$resn]),
    region=as.character(data$region),
    plddt=if ("plddt" %in% names(data)) as.numeric(data$plddt) else
      rep(NA_real_, nrow(data)),
    stringsAsFactors=FALSE
  )
  out$letter[is.na(out$letter)] <- "X"
  out
}

ram_angular_difference <- function(a, b) {
  delta <- ((as.numeric(b) - as.numeric(a) + 180) %% 360) - 180
  delta[is.na(a) | is.na(b)] <- NA_real_
  delta
}

# Needleman-Wunsch global alignment of one chain from each structure. Rows
# containing gaps are retained for display; differences are NA without both
# measured angles. A modest cell limit prevents unbounded Shiny allocations.
ram_align_residues <- function(a, b, max_cells = 4e6) {
  n <- nrow(a); m <- nrow(b)
  if (n * m > max_cells)
    stop("Comparison exceeds the sequence-alignment size limit; select shorter chains.")
  if (!n || !m) return(data.frame(
    index_a = if (n) seq_len(n) else rep(NA_integer_, m),
    index_b = if (m) seq_len(m) else rep(NA_integer_, n)
  ))
  symbols_a <- unname(ram_amino_acid_letters[a$resn])
  symbols_b <- unname(ram_amino_acid_letters[b$resn])
  symbols_a[is.na(symbols_a)] <- "X"
  symbols_b[is.na(symbols_b)] <- "X"
  # Dynamic programming retains only an integer direction for traceback.
  score <- matrix(0, n + 1L, m + 1L)
  direction <- matrix(0L, n + 1L, m + 1L)
  score[, 1L] <- -2 * (0:n)
  score[1L, ] <- -2 * (0:m)
  direction[-1L, 1L] <- 2L
  direction[1L, -1L] <- 3L
  for (i in seq_len(n)) for (j in seq_len(m)) {
    choice <- c(score[i, j] + if (symbols_a[i] == symbols_b[j]) 2 else -1,
                score[i, j + 1L] - 2, score[i + 1L, j] - 2)
    best <- which.max(choice)
    score[i + 1L, j + 1L] <- choice[[best]]
    direction[i + 1L, j + 1L] <- best
  }
  i <- n; j <- m
  ia <- integer(n + m); ib <- integer(n + m); k <- 0L
  while (i > 0L || j > 0L) {
    k <- k + 1L
    d <- direction[i + 1L, j + 1L]
    if (d == 1L) {
      ia[k] <- i; ib[k] <- j; i <- i - 1L; j <- j - 1L
    } else if (d == 2L) {
      ia[k] <- i; ib[k] <- NA_integer_; i <- i - 1L
    } else {
      ia[k] <- NA_integer_; ib[k] <- j; j <- j - 1L
    }
  }
  data.frame(index_a=rev(ia[seq_len(k)]),
             index_b=rev(ib[seq_len(k)]))
}

ram_compare_torsions <- function(a, b) {
  pairing <- ram_align_residues(a, b)
  value <- function(data, indices, key, missing) {
    out <- rep(missing, nrow(pairing))
    valid <- !is.na(indices)
    if (any(valid)) out[valid] <- data[[key]][indices[valid]]
    out
  }
  result <- data.frame(
    chain_a=value(a, pairing$index_a, "chain", NA_character_),
    residue_a=value(a, pairing$index_a, "resi", NA_integer_),
    insertion_a=value(a, pairing$index_a, "insertion_code", NA_character_),
    amino_a=value(a, pairing$index_a, "resn", NA_character_),
    chain_b=value(b, pairing$index_b, "chain", NA_character_),
    residue_b=value(b, pairing$index_b, "resi", NA_integer_),
    insertion_b=value(b, pairing$index_b, "insertion_code", NA_character_),
    amino_b=value(b, pairing$index_b, "resn", NA_character_),
    phi_a=value(a, pairing$index_a, "phi", NA_real_),
    phi_b=value(b, pairing$index_b, "phi", NA_real_),
    psi_a=value(a, pairing$index_a, "psi", NA_real_),
    psi_b=value(b, pairing$index_b, "psi", NA_real_),
    region_a=value(a, pairing$index_a, "region", NA_character_),
    region_b=value(b, pairing$index_b, "region", NA_character_),
    stringsAsFactors=FALSE
  )
  result$delta_phi <- ram_angular_difference(result$phi_a, result$phi_b)
  result$delta_psi <- ram_angular_difference(result$psi_a, result$psi_b)
  result$class_changed <- !is.na(result$region_a) &
    !is.na(result$region_b) & result$region_a != result$region_b
  result$alignment <- ifelse(is.na(pairing$index_a), "Insertion",
                      ifelse(is.na(pairing$index_b), "Deletion",
                      ifelse(result$amino_a == result$amino_b,
                             "Match", "Substitution")))
  result
}


# Keep all selected chains available in one navigator, in first-appearance
# order. This is intentionally independent of any "active chain" input.
ram_sequence_groups <- function(data) {
  residues <- ram_sequence_data(data)
  if (!nrow(residues)) return(list())
  chain <- ifelse(is.na(residues$chain) | !nzchar(residues$chain),
                  "Unassigned", residues$chain)
  split(residues, factor(chain, levels = unique(chain)), drop = TRUE)
}

ram_sequence_status <- function(region) {
  out <- rep("missing", length(region))
  out[!is.na(region) & region == "Favoured"] <- "favoured"
  out[!is.na(region) & region == "Allowed"] <- "allowed"
  out[!is.na(region) & region == "Generously allowed"] <- "generously-allowed"
  out[!is.na(region) & region == "Not allowed"] <- "outlier"
  out
}

# Small overview strips represent *positions*, not aggregate percentages.
# For long chains each strip cell represents a consecutive residue bin.
# A bin takes the highest-priority review state so isolated outliers remain
# visible instead of vanishing when hundreds of residues are compressed.
ram_sequence_overview_bins <- function(region, max_bins = 180L) {
  if (!length(region)) return(character(0))
  stopifnot(is.numeric(max_bins), length(max_bins) == 1L,
            is.finite(max_bins), max_bins >= 1L)
  status <- ram_sequence_status(region)
  bins <- min(length(status), as.integer(max_bins))
  bucket <- pmin(bins, floor((seq_along(status) - 1) * bins /
                              length(status)) + 1L)
  priority <- c("outlier", "missing", "generously-allowed",
                "allowed", "favoured")
  vapply(seq_len(bins), function(i) {
    candidates <- status[bucket == i]
    priority[match(TRUE, priority %in% candidates)]
  }, character(1))
}


# Per-chain pLDDT strip uses the lowest known confidence in each consecutive
# sequence bin; low-confidence linkers remain visible on long chains.
ram_plddt_overview_bins <- function(plddt, max_bins = 180L) {
  if (!length(plddt)) return(numeric())
  bins <- min(length(plddt), as.integer(max_bins))
  bucket <- pmin(bins, floor((seq_along(plddt) - 1) * bins /
                              length(plddt)) + 1L)
  vapply(seq_len(bins), function(i) {
    values <- plddt[bucket == i]
    if (!any(is.finite(values))) NA_real_ else min(values[is.finite(values)])
  }, numeric(1))
}
