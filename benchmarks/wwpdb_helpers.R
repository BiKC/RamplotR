# Independent wwPDB validation report ingestion and residue-wise comparison.
# These data are *not* derived from RamplotR's bundled reference distributions.
# Dependency: xml2 (CI and command-line only; Shiny does not need it).

ram_wwpdb_report_url <- function(accession) {
  id <- tolower(trimws(accession))
  if (length(id) != 1L || is.na(id) ||
      !grepl("^[a-z0-9]{4}$", id)) {
    stop("Expected one four-character PDB accession.", call. = FALSE)
  }
  sprintf("https://files.rcsb.org/pub/pdb/validation_reports/%s/%s/%s_validation.xml.gz",
          substr(id, 2L, 3L), id, id)
}

ram_read_wwpdb_report <- function(path) {
  if (!file.exists(path) || file.info(path)$size <= 0)
    stop("The wwPDB validation XML is missing or empty.", call. = FALSE)
  raw <- readBin(path, what = "raw", n = file.info(path)$size)
  if (grepl("\\.gz$", path, ignore.case = TRUE))
    raw <- memDecompress(raw, type = "gzip")
  doc <- xml2::read_xml(raw)
  nodes <- xml2::xml_find_all(doc, ".//*[local-name()='ModelledSubgroup']")
  if (!length(nodes)) {
    stop("No ModelledSubgroup entries found in the wwPDB validation XML.",
         call. = FALSE)
  }
  attr <- function(key) {
    value <- xml2::xml_attr(nodes, key)
    value[is.na(value)] <- ""
    trimws(value)
  }
  result <- data.frame(
    model = attr("model"), chain = attr("chain"),
    resi = suppressWarnings(as.integer(attr("resnum"))),
    insertion_code = attr("icode"), resn = toupper(attr("resname")),
    altcode = attr("altcode"),
    phi = suppressWarnings(as.numeric(attr("phi"))),
    psi = suppressWarnings(as.numeric(attr("psi"))),
    wwpdb_region = tolower(attr("rama")),
    stringsAsFactors = FALSE
  )
  result$wwpdb_region[result$wwpdb_region == "favoured"] <- "favored"
  result$wwpdb_region[!result$wwpdb_region %in%
    c("favored", "allowed", "outlier")] <- NA_character_
  if (!any(!is.na(result$wwpdb_region)))
    stop("No residue-level Ramachandran labels found in the wwPDB report.",
         call. = FALSE)
  result
}

# The existing original four-region RamplotR contours use 85, 98 and
# 99.95% cutoffs. For an *illustrative* 3-way crosswalk, Favoured + Allowed
# map to favored; Generously allowed maps to allowed; Not allowed to outlier.
# This harmonizes labels for contingency tables but DOES NOT make the
# underlying density grids or residue-specific criteria equivalent.
ram_wwpdb_label_crosswalk <- function(region) {
  out <- rep(NA_character_, length(region))
  out[!is.na(region) & region %in% c("Favoured", "Allowed")] <- "favored"
  out[!is.na(region) & region == "Generously allowed"] <- "allowed"
  out[!is.na(region) & region == "Not allowed"] <- "outlier"
  out
}

ram_circular_angle_difference <- function(a, b) {
  abs(((a - b + 180) %% 360) - 180)
}

ram_validation_key <- function(chain, resi, insertion_code, resn) {
  paste(as.character(chain), as.character(resi),
        as.character(insertion_code), toupper(as.character(resn)),
        sep = "\r")
}

ram_compare_wwpdb <- function(ours, external, model = 1L) {
  own_columns <- c("chain", "resi", "insertion_code", "resn", "phi",
                   "psi", "region", "reference_group")
  ref_columns <- c("model", "chain", "resi", "insertion_code", "resn",
                   "altcode", "phi", "psi", "wwpdb_region")
  if (!all(own_columns %in% names(ours)) ||
      !all(ref_columns %in% names(external)))
    stop("Validation comparison requires residue identifiers, angles and labels.")
  ref <- external[external$model == as.character(model) &
                   !is.na(external$resi), , drop = FALSE]
  if (!nrow(ref)) stop("Selected model is absent from the wwPDB report.")
  # Prefer non-alternate wwPDB rows; allow alternate A only if no unlabelled
  # row exists. Ambiguous remaining labels cannot be compared reliably.
  ref <- ref[order(!ref$altcode %in% c("", ".", "?"),
                   ref$altcode != "A"), , drop = FALSE]
  ref_key <- ram_validation_key(ref$chain, ref$resi,
                                 ref$insertion_code, ref$resn)
  ref <- ref[!duplicated(ref_key), , drop = FALSE]
  ref_key <- unique(ref_key)

  our_key <- ram_validation_key(ours$chain, ours$resi,
                                ours$insertion_code, ours$resn)
  if (anyDuplicated(our_key))
    stop("RamplotR has ambiguous chain/residue/insertion identifiers.")
  match_index <- match(our_key, ref_key)
  observed <- ref[match_index, , drop = FALSE]
  if (nrow(observed) != nrow(ours)) stop("Validation join size mismatch.")

  out <- data.frame(
    model = model,
    chain = ours$chain, resi = ours$resi,
    insertion_code = ours$insertion_code, resn = ours$resn,
    reference_group = ours$reference_group,
    ram_phi = ours$phi, ram_psi = ours$psi,
    wwpdb_phi = observed$phi, wwpdb_psi = observed$psi,
    phi_abs_degrees = ram_circular_angle_difference(ours$phi, observed$phi),
    psi_abs_degrees = ram_circular_angle_difference(ours$psi, observed$psi),
    ram_region_four = ours$region,
    ram_region_three = ram_wwpdb_label_crosswalk(ours$region),
    wwpdb_region = observed$wwpdb_region,
    identifier_matched = !is.na(match_index),
    stringsAsFactors = FALSE
  )
  out$labels_comparable <- !is.na(out$wwpdb_region) &
    !is.na(out$ram_region_three) &
    is.finite(out$ram_phi) & is.finite(out$ram_psi)
  out$labels_agree <- rep(NA, nrow(out))
  ok <- which(out$labels_comparable)
  out$labels_agree[ok] <-
    out$ram_region_three[ok] == out$wwpdb_region[ok]
  out
}

ram_wwpdb_summary <- function(comparison, accession,
                               min_coverage = 0.90,
                               max_angle_deviation = 1.5) {
  if (!nrow(comparison)) stop("No residues for wwPDB comparison.")
  finite_ours <- is.finite(comparison$ram_phi) &
    is.finite(comparison$ram_psi)
  fully_matched <- finite_ours & is.finite(comparison$wwpdb_phi) &
    is.finite(comparison$wwpdb_psi)
  comparable <- comparison$labels_comparable & fully_matched
  coverage <- sum(fully_matched) / max(1L, sum(finite_ours))
  agree <- sum(comparison$labels_agree[comparable], na.rm = TRUE)
  report <- data.frame(
    accession = accession,
    ram_residues = nrow(comparison),
    ram_finite_phi_psi = sum(finite_ours),
    wwpdb_identifier_matches = sum(comparison$identifier_matched),
    matched_phi_psi = sum(fully_matched),
    matched_angle_coverage = coverage,
    max_phi_difference = if (any(fully_matched))
      max(comparison$phi_abs_degrees[fully_matched]) else NA_real_,
    max_psi_difference = if (any(fully_matched))
      max(comparison$psi_abs_degrees[fully_matched]) else NA_real_,
    comparable_labels = sum(comparable),
    matching_labels = agree,
    differing_labels = sum(comparable) - agree,
    label_agreement = agree / max(1L, sum(comparable)),
    stringsAsFactors = FALSE
  )
  if (coverage < min_coverage)
    stop(sprintf("Only %.1f%% of finite angles matched wwPDB in %s.",
                 100 * coverage, accession), call. = FALSE)
  if (!any(fully_matched) ||
      any(comparison$phi_abs_degrees[fully_matched] > max_angle_deviation) ||
      any(comparison$psi_abs_degrees[fully_matched] > max_angle_deviation)) {
    stop(sprintf("Independent wwPDB angular comparison failed in %s.",
                 accession), call. = FALSE)
  }
  if (!sum(comparable))
    stop("No comparable independent residue classifications.", call. = FALSE)
  report
}

ram_wwpdb_contingency <- function(comparison) {
  groups <- unique(comparison$reference_group)
  data <- comparison[comparison$labels_comparable, , drop = FALSE]
  # Every category remains visible even if absent in this structure.
  levels <- c("favored", "allowed", "outlier")
  if (!nrow(data)) stop("No comparable labels for the contingency matrix.")
  do.call(rbind, lapply(groups, function(group) {
    group_data <- data[data$reference_group == group, , drop = FALSE]
    cross <- table(
      ramplotr = factor(group_data$ram_region_three, levels = levels),
      wwpdb = factor(group_data$wwpdb_region, levels = levels)
    )
    result <- as.data.frame(cross, responseName = "residues")
    names(result) <- c("ramplotr", "wwpdb", "residues")
    result$reference_group <- group
    result[, c("reference_group", "ramplotr", "wwpdb", "residues")]
  }))
}
