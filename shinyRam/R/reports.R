# Export scientific figures and a self-contained, editable HTML analysis report.
# All figures use the selected reference *background*, the same contour
# thresholds and the displayed coordinates; numerical data is not rounded
# until writing a human-readable table.

ram_draw_figure <- function(data, reference, colors, chain_colors,
                            title = "RamplotR | Ramachandran analysis",
                            use_raster = FALSE) {
  stopifnot(length(colors) == 4L, is.matrix(reference$z))
  limits <- ram_density_thresholds(reference)
  z <- reference$z
  code <- matrix(1L, nrow(z), ncol(z))
  code[z > limits[[3L]]] <- 2L
  code[z > limits[[2L]]] <- 3L
  code[z > limits[[1L]]] <- 4L
  graphics::par(mar = c(4.7, 4.7, 2.8, 2.0), bg = "white")
  graphics::image(
    x = reference$x, y = reference$y, z = t(code),
    col = as.character(colors), breaks = seq(0.5, 4.5, by = 1),
    xlim = c(-180, 180), ylim = c(-180, 180), asp = 1,
    xlab = expression(phi~"(degrees)"), ylab = expression(psi~"(degrees)"),
    main = title, axes = FALSE, useRaster = use_raster
  )
  graphics::axis(1, at = c(-180, -90, 0, 90, 180))
  graphics::axis(2, at = c(-180, -90, 0, 90, 180), las = 1)
  graphics::box(col = "#9DB7B8")
  graphics::abline(h = 0, v = 0, col = "#B5CBCA", lty = 2)
  valid <- is.finite(data$phi) & is.finite(data$psi)
  present <- unique(as.character(data$chain[valid]))
  for (i in seq_along(present)) {
    rows <- which(valid & data$chain == present[[i]])
    chain_index <- match(present[[i]], names(chain_colors))
    color <- if (!is.na(chain_index)) chain_colors[[chain_index]] else
      c("#D57359", "#307D99", "#80649A", "#BF8E40")[[((i - 1L) %% 4L) + 1L]]
    graphics::points(data$phi[rows], data$psi[rows],
                     pch = 21, cex = 0.65, bg = color,
                     col = "#FFFFFF", lwd = 0.35)
  }
  if (length(present) && length(present) <= 8L) {
    graphics::legend("topright", legend = paste("Chain", present),
                     pch = 21, pt.bg = vapply(present, function(chain) {
                       if (chain %in% names(chain_colors))
                         chain_colors[[chain]] else "#4D8F96"
                     }, character(1)), bty = "n", cex = 0.7)
  }
}

ram_save_figure <- function(path, data, reference, colors, chain_colors,
                            format = c("svg", "png"), title = "RamplotR") {
  format <- match.arg(format)
  if (format == "svg") {
    grDevices::svg(path, width = 7.1, height = 7.1, pointsize = 10)
  } else {
    grDevices::png(path, width = 2200, height = 2200, res = 300,
                   type = if (capabilities("cairo")) "cairo" else "Xlib")
  }
  on.exit(grDevices::dev.off(), add = TRUE)
  ram_draw_figure(data, reference, colors, chain_colors, title,
                  use_raster = identical(format, "png"))
  invisible(path)
}

ram_report_metadata <- function(name, reference_set, background,
                                validation_mode, model, reference_file) {
  list(
    structure = name,
    reference_set = reference_set,
    displayed_background = background,
    classification = validation_mode,
    model = as.character(model),
    reference_md5 = if (file.exists(reference_file))
      unname(tools::md5sum(reference_file)) else "unavailable",
    generated_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
    R_version = R.version.string,
    bio3d_version = as.character(utils::packageVersion("bio3d"))
  )
}

ram_save_html_report <- function(path, data, metadata, svg_path,
                                 max_report_rows=10000L, ensemble=NULL) {
  h <- function(label,value) htmltools::tags$tr(
    htmltools::tags$th(label),
    htmltools::tags$td(paste(as.character(value),collapse=", ")))
  table_for <- function(frame) {
    rows <- lapply(seq_len(nrow(frame)),function(i)
      htmltools::tags$tr(lapply(frame[i,,drop=FALSE],function(value)
        htmltools::tags$td(if(is.na(value[[1L]])) "" else
                             as.character(value[[1L]])))))
    htmltools::tags$div(style="overflow-x:auto",
      htmltools::tags$table(
        htmltools::tags$thead(htmltools::tags$tr(
          lapply(names(frame),htmltools::tags$th))),
        htmltools::tags$tbody(rows)))
  }
  tally <- vapply(c("Favoured","Allowed","Generously allowed","Not allowed"),
    function(region) sum(data$region==region,na.rm=TRUE),integer(1))
  standard_tally <- if ("rama8000_region" %in% names(data))
    vapply(c("Favored","Allowed","Outlier"),
      function(region) sum(data$rama8000_region==region,na.rm=TRUE),integer(1))
    else NULL
  counts <- lapply(names(tally),function(label) h(label,tally[[label]]))
  counts <- c(counts,list(
    h("Missing angles",sum(is.na(data$phi)|is.na(data$psi))),
    h("All selected residues",nrow(data))))
  meta <- lapply(names(metadata),function(name)
    h(gsub("_"," ",name),metadata[[name]]))
  shown <- utils::head(data[,intersect(c(
    "chain","resi","insertion_code","resn","phi","psi","region",
    "density","rama8000_region","rama8000_group","rama8000_score",
    "omega","omega_status","chi1","cb_ca_distance",
    "cb_signed_volume","plddt","confidence_category"),
    names(data)),drop=FALSE],max_report_rows)
  for(field in intersect(c("phi","psi","density","rama8000_score",
                       "omega","chi1","cb_ca_distance",
                       "cb_signed_volume","plddt"),
                         names(shown)))
    shown[[field]] <- round(shown[[field]],2L)
  official <- NULL
  if("wwpdb_matched" %in% names(data)) {
    fields <- intersect(c("chain","resi","insertion_code","resn",
      "wwpdb_rama","wwpdb_rotamer","wwpdb_clashes",
      "wwpdb_bond_outliers","wwpdb_angle_outliers","wwpdb_rscc","wwpdb_rsrz"),
      names(data))
    official <- utils::head(data[data$wwpdb_matched %in% TRUE,
                                fields,drop=FALSE],max_report_rows)
  }
  if(!is.null(ensemble)) {
    ensemble <- utils::head(ensemble,max_report_rows)
    for(field in intersect(c("phi_mean","phi_sd","psi_mean","psi_sd",
                             "class_consistency"),names(ensemble)))
      ensemble[[field]] <- round(ensemble[[field]],2L)
  }
  lines <- readLines(svg_path,warn=FALSE)
  chart <- htmltools::HTML(paste(lines[!grepl("^<\\?xml",lines)],
                                collapse="\n"))
  doc <- htmltools::tags$html(
    htmltools::tags$head(
      htmltools::tags$meta(charset="UTF-8"),
      htmltools::tags$title("RamplotR structure analysis"),
      htmltools::tags$style(htmltools::HTML(
        "body{font:15px system-ui,sans-serif;margin:3em auto;max-width:1100px;color:#18323b}h1,h2{color:#146a70}table{border-collapse:collapse;width:100%}th,td{padding:6px 12px;border-bottom:1px solid #dbe6e6;text-align:left}th{background:#eef6f4}svg{width:min(100%,760px);height:auto}p.note{color:#516b74;line-height:1.5}section{margin:2em 0}"))
    ),
    htmltools::tags$body(
      htmltools::tags$h1("RamplotR structure analysis"),
      htmltools::tags$p(class="note",
        "Portable report: figure, counts, geometry and provenance. Native and independent results are reported separately."),
      htmltools::tags$section(htmltools::tags$h2("Ramachandran plot"),chart),
      htmltools::tags$section(htmltools::tags$h2("RamplotR density regions"),
        htmltools::tags$p(class="note",
          "Native RamplotR density regions. Not allowed is not a MolProbity/wwPDB outlier label."),
        table_for(data.frame(
          Region=c(names(tally),"Missing angles","All selected residues"),
          Residues=c(as.integer(tally),
            sum(is.na(data$phi)|is.na(data$psi)),nrow(data))))),
      if(!is.null(standard_tally))
        htmltools::tags$section(
          htmltools::tags$h2("Rama8000 standard validation"),
          htmltools::tags$p(class="note",
            "Independent six-class Rama8000 evaluation using the current cctbx/Phenix reference tables and Favored/Allowed/Outlier thresholds."),
          table_for(data.frame(
            Region=names(standard_tally),
            Residues=as.integer(standard_tally)
          ))),
      if("omega_status" %in% names(data))
        htmltools::tags$section(
          htmltools::tags$h2("Additional descriptive geometry"),
          htmltools::tags$p(class="note",
            "Omega (peptide dihedral) and chi1 (first side-chain dihedral) are calculated locally. No rotamer or clash outlier status is inferred from these measurements."),
          htmltools::tags$p(sprintf(
            "%d cis and %d twisted peptide bonds; %d measurable chi1 angles.",
            sum(data$omega_status=="Cis",na.rm=TRUE),
            sum(data$omega_status=="Twisted",na.rm=TRUE),
            sum(is.finite(data$chi1))))
        ),
      if(!is.null(official))
        htmltools::tags$section(
          htmltools::tags$h2("Independent wwPDB validation"),
          htmltools::tags$p(class="note",
            "Official residue-level validation values from the attached report, not computed by RamplotR. Missing records are not scored; underlying contour classifications may differ."),
          if(sum(data$wwpdb_matched %in% TRUE)>max_report_rows)
            htmltools::tags$p(class="note",
              "Independent table truncated; use the detailed CSV for every matched residue."),
          table_for(official)),
      if(!is.null(ensemble))
        htmltools::tags$section(
          htmltools::tags$h2("Structural ensemble"),
          htmltools::tags$p(class="note",
            "Residue-matched circular dihedral means and model consistency. Missing and one-model observations do not imply zero variation."),
          table_for(ensemble)),
      htmltools::tags$section(htmltools::tags$h2("Analysis provenance"),
        table_for(data.frame(Setting=gsub("_"," ",names(metadata)),
          Value=vapply(metadata,function(x)paste(as.character(x),
                                                 collapse=", "),character(1))))),
      htmltools::tags$section(htmltools::tags$h2("Residue details"),
        if(nrow(data)>max_report_rows)
          htmltools::tags$p(class="note",
            sprintf("Showing the first %d of %d residues; the CSV contains all rows.",
              max_report_rows,nrow(data))),
        table_for(shown))
    )
  )
  writeLines(as.character(doc),path,useBytes=TRUE)
  invisible(path)
}


# Standalone report for independently generated prediction ensembles. This is
# intentionally separate from the single-structure report: model disagreement
# is prediction uncertainty/heterogeneity, not evidence of molecular dynamics.
ram_save_prediction_ensemble_report <- function(path, result, metadata=list(),
                                                max_report_rows=10000L) {
  if(!is.list(result) || is.null(result$summary) ||
     !is.data.frame(result$summary) || is.null(result$model_summary) ||
     !is.data.frame(result$model_summary))
    stop("Prediction ensemble report requires residue and model summaries.")
  if(!is.numeric(max_report_rows) || length(max_report_rows)!=1L ||
     !is.finite(max_report_rows) || max_report_rows<1L)
    stop("max_report_rows must be a positive finite number.")
  max_report_rows <- as.integer(max_report_rows)

  table_for <- function(frame) {
    if(!nrow(frame))
      return(htmltools::tags$p(class="note","No rows available."))
    rows <- lapply(seq_len(nrow(frame)),function(i)
      htmltools::tags$tr(lapply(frame[i,,drop=FALSE],function(value)
        htmltools::tags$td(if(length(value)==0L || is.na(value[[1L]])) ""
          else as.character(value[[1L]])))))
    htmltools::tags$div(class="table-wrap",
      htmltools::tags$table(
        htmltools::tags$thead(htmltools::tags$tr(
          lapply(names(frame),htmltools::tags$th))),
        htmltools::tags$tbody(rows)))
  }
  metric <- function(label,value)
    htmltools::tags$div(class="metric",
      htmltools::tags$span(label),htmltools::tags$strong(value))

  residues <- result$summary
  spread <- pmax(residues$phi_sd,residues$psi_sd,na.rm=TRUE)
  spread[!is.finite(residues$phi_sd) & !is.finite(residues$psi_sd)] <- NA_real_
  rama_changes <- if("rama8000_changes" %in% names(residues))
    sum(residues$rama8000_changes,na.rm=TRUE) else NA_integer_
  plddt_variable <- if("plddt_sd" %in% names(residues))
    sum(is.finite(residues$plddt_sd) & residues$plddt_sd>=10) else NA_integer_
  angular_variable <- sum(is.finite(spread) & spread>=20)
  common <- if(!is.null(result$common_residues)) result$common_residues else
    sum(residues$models_present==result$analyzed_models,na.rm=TRUE)

  ranked <- residues
  ranked$max_angular_sd <- spread
  ranked <- ranked[order(-replace(ranked$max_angular_sd,
                                  !is.finite(ranked$max_angular_sd),-Inf)),
                   ,drop=FALSE]
  ranked <- utils::head(ranked,25L)
  rank_fields <- intersect(c("chain","resi","insertion_code","resn",
    "models_present","phi_sd","psi_sd","max_angular_sd",
    "rama8000_mode","rama8000_consistency",
    "plddt_mean","plddt_sd","plddt_min","plddt_max"),names(ranked))
  ranked <- ranked[,rank_fields,drop=FALSE]

  all_fields <- intersect(c("chain","resi","insertion_code","resn",
    "models_present","phi_models","psi_models","phi_mean","phi_sd",
    "psi_mean","psi_sd","rama8000_models","rama8000_mode",
    "rama8000_consistency","rama8000_changes","plddt_models","plddt_mean",
    "plddt_sd","plddt_min","plddt_max"),names(residues))
  shown <- utils::head(residues[,all_fields,drop=FALSE],max_report_rows)
  numeric_fields <- c("phi_mean","phi_sd","psi_mean","psi_sd",
    "rama8000_consistency","plddt_mean","plddt_sd","plddt_min","plddt_max",
    "max_angular_sd")
  for(field in intersect(numeric_fields,union(names(shown),names(ranked)))) {
    if(field %in% names(shown)) shown[[field]] <- round(shown[[field]],2L)
    if(field %in% names(ranked)) ranked[[field]] <- round(ranked[[field]],2L)
  }

  models <- result$model_summary
  if(!is.null(result$provenance) && is.data.frame(result$provenance)) {
    provenance <- result$provenance
    if(nrow(provenance)!=nrow(models))
      stop("Prediction ensemble provenance no longer matches model order.")
    for(field in setdiff(names(provenance),"model"))
      models[[field]] <- provenance[[field]]
  }
  for(field in intersect(c("plddt_mean","plddt_min"),names(models)))
    models[[field]] <- round(models[[field]],2L)
  for(field in intersect(c("ptm","iptm","ranking_score",
                           "fraction_disordered"),names(models)))
    models[[field]] <- round(models[[field]],4L)

  provenance <- c(metadata,list(
    prediction_source=if(is.null(result$source)) "unknown" else result$source,
    analyzed_models=result$analyzed_models,
    available_models=result$available_models,
    generated_utc=format(Sys.time(),"%Y-%m-%dT%H:%M:%SZ",tz="UTC"),
    R_version=R.version.string
  ))
  provenance_frame <- data.frame(
    Setting=gsub("_"," ",names(provenance)),
    Value=vapply(provenance,function(x)
      paste(as.character(x),collapse=", "),character(1)),
    stringsAsFactors=FALSE
  )

  doc <- htmltools::tags$html(
    htmltools::tags$head(
      htmltools::tags$meta(charset="UTF-8"),
      htmltools::tags$title("RamplotR Prediction ensemble"),
      htmltools::tags$style(htmltools::HTML(
        "body{font:15px system-ui,sans-serif;margin:3em auto;max-width:1180px;color:#18323b}h1,h2{color:#146a70}.lead,.note{color:#516b74;line-height:1.55}.warning{padding:12px 14px;border-left:4px solid #c36b4f;background:#fff5ee}.metrics{display:flex;gap:10px;flex-wrap:wrap}.metric{border:1px solid #dbe7e5;border-radius:8px;padding:8px 11px;min-width:145px;background:#f8fbfa}.metric span{display:block;color:#657c82;font-size:12px}.metric strong{font-size:20px;color:#164f57}.table-wrap{overflow-x:auto}table{border-collapse:collapse;width:100%}th,td{padding:6px 10px;border-bottom:1px solid #dbe6e6;text-align:left;white-space:nowrap}th{background:#eef6f4}section{margin:2em 0}"))
    ),
    htmltools::tags$body(
      htmltools::tags$h1("RamplotR Prediction ensemble"),
      htmltools::tags$p(class="lead",
        "Residue-level agreement across independently generated prediction models/seeds."),
      htmltools::tags$p(class="warning",
        htmltools::tags$strong("Interpretation: "),
        "model-to-model spread represents prediction uncertainty or heterogeneity. It is not experimental evidence of molecular motion or conformational dynamics."),
      htmltools::tags$section(
        htmltools::tags$h2("Ensemble overview"),
        htmltools::tags$div(class="metrics",
          metric("Models analysed",result$analyzed_models),
          metric("Common residues",common),
          metric("Rama8000 disagreements",
            if(is.na(rama_changes)) "n/a" else rama_changes),
          metric("Angular SD >=20°",angular_variable),
          metric("pLDDT SD >=10",
            if(is.na(plddt_variable)) "n/a" else plddt_variable)
        )
      ),
      htmltools::tags$section(
        htmltools::tags$h2("Model provenance"),
        htmltools::tags$p(class="note",
          "Coordinate hashes identify the exact uploaded models when available."),
        table_for(models)
      ),
      htmltools::tags$section(
        htmltools::tags$h2("Largest local backbone variability"),
        htmltools::tags$p(class="note",
          "Ranked by max(circular SD of phi, circular SD of psi). This is a navigation ranking, not a significance statistic."),
        table_for(ranked)
      ),
      htmltools::tags$section(
        htmltools::tags$h2("All residue-level ensemble statistics"),
        if(nrow(residues)>max_report_rows)
          htmltools::tags$p(class="note",
            sprintf("Showing the first %d of %d rows.",max_report_rows,nrow(residues))),
        table_for(shown)
      ),
      htmltools::tags$section(
        htmltools::tags$h2("Analysis provenance"),
        table_for(provenance_frame)
      )
    )
  )
  writeLines(as.character(doc),path,useBytes=TRUE)
  invisible(path)
}
