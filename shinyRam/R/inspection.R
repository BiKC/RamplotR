# Pure residue inspection and pairwise comparison helpers.
# Nothing in this file changes torsion extraction or scientific classification.

ram_review_queue <- function(data, boundary_margin = 2) {
  if (!nrow(data)) return(data)
  missing <- is.na(data$phi) | is.na(data$psi) | is.na(data$region)
  not_allowed <- !missing & data$region == "Not allowed"
  standard_outlier <- if ("rama8000_region" %in% names(data))
    !missing & !is.na(data$rama8000_region) & data$rama8000_region == "Outlier"
    else rep(FALSE, nrow(data))
  # The density percentile is cumulative mass above the residue's density.
  # Closeness to any RamplotR contour probability is a review hint only.
  near <- !missing & !is.na(data$density) &
    vapply(data$density, function(value) {
      any(abs(value - c(85, 98, 99.95)) <= boundary_margin)
    }, logical(1))
  data$review_status <- ifelse(missing, "Missing angles",
    ifelse(standard_outlier, "Rama8000 outlier",
      ifelse(not_allowed, "Not allowed",
        ifelse(near, "Near boundary", "Other"))))
  priority <- match(data$review_status,
                    c("Rama8000 outlier", "Not allowed", "Missing angles",
                      "Near boundary", "Other"))
  data[order(priority, data$chain, data$resi, data$insertion_code), ,
       drop = FALSE]
}

ram_amino_acid_letters <- c(
  ALA="A", ARG="R", ASN="N", ASP="D", CYS="C", GLU="E", GLN="Q",
  GLY="G", HIS="H", ILE="I", LEU="L", LYS="K", MET="M", PHE="F",
  PRO="P", SER="S", THR="T", TRP="W", TYR="Y", VAL="V"
)

ram_chain_query_sequence <- function(data, chain) {
  residues <- ram_sequence_data(data)
  if (!nrow(residues)) return(list(sequence="", residues=residues,
                                   known_fraction=NA_real_))
  target <- as.character(chain)
  rows <- residues$chain == target
  residues <- residues[rows, , drop=FALSE]
  if (!nrow(residues)) return(list(sequence="", residues=residues,
                                   known_fraction=NA_real_))
  letters <- toupper(as.character(residues$letter))
  letters[!grepl("^[ACDEFGHIKLMNPQRSTVWY]$", letters)] <- "X"
  known <- letters != "X"
  list(
    sequence=paste0(letters,collapse=""),
    residues=residues,
    known_fraction=mean(known)
  )
}

ram_sequence_data <- function(data) {
  if (!nrow(data)) return(data.frame(
    chain=character(), resi=integer(), insertion_code=character(),
    resn=character(), letter=character(), region=character(),
    rama8000_region=character(), rama8000_group=character(),
    rama8000_score=numeric(), plddt=numeric(),
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
    rama8000_region=if ("rama8000_region" %in% names(data))
      as.character(data$rama8000_region) else rep(NA_character_,nrow(data)),
    rama8000_group=if ("rama8000_group" %in% names(data))
      as.character(data$rama8000_group) else rep(NA_character_,nrow(data)),
    rama8000_score=if ("rama8000_score" %in% names(data))
      as.numeric(data$rama8000_score) else rep(NA_real_,nrow(data)),
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

# A simple local displacement in phi/psi space after each angular component
# has been wrapped independently. This is a navigation/ranking measure, not a
# statistical significance score or a Cartesian structural distance.
ram_backbone_angular_displacement <- function(delta_phi, delta_psi) {
  phi <- as.numeric(delta_phi)
  psi <- as.numeric(delta_psi)
  out <- sqrt(phi^2 + psi^2)
  out[!is.finite(phi) | !is.finite(psi)] <- NA_real_
  out
}

ram_backbone_shift_band <- function(displacement) {
  value <- as.numeric(displacement)
  out <- rep("Unavailable", length(value))
  finite <- is.finite(value)
  out[finite & value < 15] <- "Small"
  out[finite & value >= 15 & value < 30] <- "Moderate"
  out[finite & value >= 30 & value < 60] <- "Large"
  out[finite & value >= 60] <- "Very large"
  out
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
  result$angular_displacement <- ram_backbone_angular_displacement(
    result$delta_phi, result$delta_psi)
  result$shift_band <- ram_backbone_shift_band(result$angular_displacement)
  result$class_changed <- !is.na(result$region_a) &
    !is.na(result$region_b) & result$region_a != result$region_b
  if (all(c("rama8000_region","rama8000_group","rama8000_score") %in% names(a)) &&
      all(c("rama8000_region","rama8000_group","rama8000_score") %in% names(b))) {
    result$rama8000_region_a <- value(
      a, pairing$index_a, "rama8000_region", NA_character_)
    result$rama8000_region_b <- value(
      b, pairing$index_b, "rama8000_region", NA_character_)
    result$rama8000_group_a <- value(
      a, pairing$index_a, "rama8000_group", NA_character_)
    result$rama8000_group_b <- value(
      b, pairing$index_b, "rama8000_group", NA_character_)
    result$rama8000_score_a <- value(
      a, pairing$index_a, "rama8000_score", NA_real_)
    result$rama8000_score_b <- value(
      b, pairing$index_b, "rama8000_score", NA_real_)
    result$rama8000_changed <- !is.na(result$rama8000_region_a) &
      !is.na(result$rama8000_region_b) &
      result$rama8000_region_a != result$rama8000_region_b
  }
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

# Confidence is an independent per-residue signal: never use its colours to
# replace the Ramachandran classification background.
ram_plddt_color <- function(score) {
  vapply(as.numeric(score), function(value) {
    if (!is.finite(value)) "#cbd7db" else if (value < 50) "#d75e56"
    else if (value < 70) "#d6ac52" else if (value < 90) "#7bbcb1"
    else "#126e74"
  }, character(1))
}

# Permanent position labels mark every tenth *PDB residue number*, not every
# tenth item in the sequence. Keep insertion codes on labelled residues.
ram_sequence_position_labels <- function(resi, insertion_code = rep("", length(resi))) {
  stopifnot(length(resi) == length(insertion_code))
  if (!length(resi)) return(character())
  index <- seq_along(resi)
  label <- index == 1L | index == length(resi) |
    (!is.na(resi) & resi %% 10L == 0L)
  out <- rep("", length(resi))
  out[label] <- paste0(resi[label], ifelse(is.na(insertion_code[label]), "",
                                                insertion_code[label]))
  out
}

# Map an NGL/sequence residue to its pair using chain, PDB numbering and
# insertion code. The displayed amino-acid order may contain alignment gaps.
ram_comparison_find <- function(data, side, chain, resi, insertion_code = "") {
  stopifnot(side %in% c("a", "b"))
  if (!nrow(data) || length(chain) != 1L || length(resi) != 1L ||
      length(insertion_code) != 1L || is.na(chain) || is.na(resi) ||
      is.na(insertion_code)) return(NA_integer_)
  number <- suppressWarnings(as.integer(resi))
  if (is.na(number)) return(NA_integer_)
  ix <- which(!is.na(data[[paste0("residue_", side)]]) &
    data[[paste0("chain_", side)]] == as.character(chain) &
    data[[paste0("residue_", side)]] == number &
    data[[paste0("insertion_", side)]] == as.character(insertion_code))
  if (!length(ix)) NA_integer_ else as.integer(ix[[1L]])
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
