#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#


# shiny related packages
library(shiny)
library(shinyWidgets)
library(colourpicker)

# used for pdb files
library(bio3d)
library(NGLVieweR)

# Used for processing data

source(file.path("R", "ramachandran.R"), local = TRUE)
source(file.path("R", "backbone.R"), local = TRUE)
source(file.path("R", "io.R"), local = TRUE)

color_set <- c(
  "#137C79", "#8662A8", "#C47B36", "#4E86B3",
  "#648E5E", "#B95873", "#9C7545", "#62798D"
)

rampage<-c("#F1EEF6","#BDC9E1","#74A9CF","#0570B0")
pdbsum<-c("#FEFFB2","#F3F300","#C37800","#F30100")

# load bg image
#fig <- readRDS("static/bgfig.RDS")

# densMatGeneral <- readRDS("shinyRam/static/alphafold/densMatGeneral")
# densMatGeneral$z<-t(densMatGeneral$z)
# saveRDS(densMatGeneral,"shinyRam/static/alphafold/densMatGeneral")
# 
# densMatGly <- readRDS("shinyRam/static/alphafold/densMatGly")
# densMatGly$z<-t(densMatGly$z)
# saveRDS(densMatGly,"shinyRam/static/alphafold/densMatGly")
# 
# densMatPrepro <- readRDS("shinyRam/static/alphafold/densMatprePro")
# densMatPrepro$z<-t(densMatPrepro$z)
# saveRDS(densMatPrepro,"shinyRam/static/alphafold/densMatPrepro")
# 
# densMatPro <- readRDS("shinyRam/static/alphafold/densMatPro")
# densMatPro$z<-t(densMatPro$z)
# saveRDS(densMatPro,"shinyRam/static/alphafold/densMatPro")

allAA <- c(
  "ALA",
  "ARG",
  "ASN",
  "ASP",
  "CYS",
  "GLU",
  "GLN",
  "GLY",
  "HIS",
  "ILE",
  "LEU",
  "LYS",
  "MET",
  "PHE",
  "PRO",
  "SER",
  "THR",
  "TRP",
  "TYR",
  "VAL"
)

# Presentation stays separate from scientific analysis. Input and output IDs
# are preserved for compatibility with plotting and the molecular viewer.
ui <- fluidPage(
  tags$head(
    tags$meta(name = "viewport", content = "width=device-width, initial-scale=1"),
    tags$title("RamplotR | Ramachandran analysis"),
    tags$link(rel = "stylesheet", type = "text/css", href = "styles.css"),
    tags$script(src = "https://cdn.plot.ly/plotly-2.14.0.min.js")
  ),
  tags$div(
    class = "ram-app",
    tags$header(
      class = "ram-header",
      tags$div(
        class = "ram-inner",
        tags$div(
          class = "ram-brand",
          tags$span(class = "ram-brand-mark", "φψ", "aria-hidden" = "true"),
          tags$div(
            tags$div(class = "ram-brand-name", "RamplotR"),
            tags$div(class = "ram-brand-subtitle", "Protein structure analysis")
          )
        ),
        tags$div(class = "ram-header-badge", "Interactive structural validation")
      )
    ),
    tags$div(
      class = "ram-inner",
      tags$div(
        class = "ram-intro",
        tags$p(class = "ram-overline", "Ramachandran analysis"),
        tags$h1("Explore backbone geometry."),
        tags$p("Inspect residue conformations, compare reference distributions and examine the protein in 3D.")
      ),
      tags$section(
        class = "ram-source", "aria-label" = "Structure input",
        tags$div(
          class = "ram-section-heading",
          tags$div(
            tags$h2("Load a structure"),
            tags$p("Enter a PDB accession or use a local PDB/mmCIF file.")
          )
        ),
        tags$div(
          class = "ram-source-controls",
          tags$div(
            class = "ram-source-choice",
            radioButtons(
              "inputSource", "Structure source",
              choices = c("PDB accession" = "pdb", "Uploaded file" = "upload"),
              selected = "pdb", inline = TRUE
            )
          ),
          tags$div(
            id = "ram-pdb-wrap", class = "ram-source-picker",
            textInput("PDB", "PDB accession", value = "1BBB",
                      placeholder = "e.g. 1CRN")
          ),
          tags$div(
            id = "ram-upload-wrap", class = "ram-source-picker is-hidden",
            fileInput(
              "structfile", "Structure file",
              accept = c(".pdb", ".ent", ".cif", ".mmcif", ".mcif")
            )
          ),
          tags$div(
            class = "ram-submit",
            actionButton("submit", "Analyze structure", class = "btn-primary")
          )
        ),
        tags$p(
          class = "ram-tip",
          "Classification uses the selected reference dataset; the plotted background can be changed independently."
        )
      ),
      tags$div(
        class = "ram-workspace",
        tags$aside(
          id = "ram-settings", class = "ram-sidebar", "aria-label" = "Analysis settings",
          tags$section(
            class = "ram-panel",
            tags$div(
              class = "ram-section-heading",
              tags$div(
                tags$h3("Reference & validation"),
                tags$p("Choose how residues and plot regions are interpreted.")
              )
            ),
            selectInput(
              "validationMode", "Residue classification",
              c("Residue-aware (recommended)" = "residue",
                "Selected background (legacy)" = "legacy")
            ),
            selectInput(
              "bgtype", "Reference dataset",
              c("original", "alphafold", "alphafold_filtered",
                "astral2.08", "custom_high_resolution")
            ),
            selectInput(
              "background", "Plot background",
              c("General", "Glycine", "Preproline", "Proline")
            ),
            tags$p(class = "ram-field-hint",
                   "Residue-aware mode evaluates glycine, proline and pre-proline against their own reference distributions.")
          ),
          tags$section(
            class = "ram-panel",
            tags$div(
              class = "ram-section-heading",
              tags$div(
                tags$h3("Residue selection"),
                tags$p("Focus on individual amino acids and protein chains.")
              )
            ),
            pickerInput(
              "AA", "Amino acids", choices = allAA, multiple = TRUE,
              options = list("actions-box" = TRUE),
              selected = allAA
            ),
            uiOutput("chains"),
            tags$p(class = "ram-field-hint",
                   "Selections apply to the plot, residue list and summary.")
          ),
          tags$section(
            class = "ram-panel",
            tags$div(
              class = "ram-section-heading",
              tags$div(
                tags$h3("Appearance"),
                tags$p("Keep the default palette or customize the contours.")
              )
            ),
            selectInput(
              "colorscheme", "Contour palette",
              choices = c("Rampage", "PDBSum", "custom"),
              selected = "Rampage"
            ),
            tags$details(
              class = "ram-details",
              tags$summary("Customize region colors"),
              tags$p(class = "ram-field-hint",
                     "Select Custom as your contour palette to use these colours."),
              tags$div(
                class = "ram-details-body ram-appearance-controls",
                colourpicker::colourInput("bg1", "Not allowed", value = "#F1EEF6"),
                colourpicker::colourInput("bg2", "Generously allowed", value = "#BDC9E1"),
                colourpicker::colourInput("bg3", "Allowed", value = "#74A9CF"),
                colourpicker::colourInput("bg4", "Favoured", value = "#0570B0")
              )
            ),
            tags$div(class = "ram-panel-divider"),
            tags$div(class = "ram-chain-controls", uiOutput("chainColors"))
          )
        ),
        tags$main(
          class = "ram-main",
          tabsetPanel(
            id = "analysisTabs",
            tabPanel(
              title = "Ramachandran plot", value = "plot",
              tags$div(
                class = "ram-result-head",
                tags$div(
                  tags$h2("Conformation overview"),
                  tags$p("Use the plot to inspect phi (φ) and psi (ψ) backbone angles.")
                ),
                tags$div(
                  class = "ram-result-actions",
                  tags$a(class = "ram-mobile-settings-link",
                         href = "#ram-settings", "Analysis settings"),
                  tags$span(id = "ram-current-structure",
                            class = "ram-status", "No structure loaded")
                )
              ),
              tags$div(
                class = "ram-selected-residue",
                uiOutput("selectedResidueInfo"),
                actionButton("clearResidue", "Clear selection", class = "btn-default btn-sm")
              ),
              tags$div(
                class = "ram-charts",
                tags$section(
                  class = "ram-chart-card", "aria-label" = "Ramachandran plot",
                  tags$h3(class = "ram-chart-label", "Residue distribution"),
                  tags$p(class = "ram-chart-help",
                         "Click a point to highlight its residue in the table and 3D structure."),
                  tags$div(
                    id = "plot-empty", class = "ram-plot-empty",
                    tags$span(class = "ram-empty-mark",
                              "aria-hidden" = "true", "φψ"),
                    tags$strong("Your plot will appear here"),
                    tags$p("Enter an accession and select Analyze structure.")
                  ),
                  tags$div(id = "plotly", class = "ram-plot",
                           role = "img", "aria-label" = "Interactive Ramachandran plot")
                ),
                tags$section(
                  class = "ram-chart-card", "aria-label" = "3D molecular viewer",
                  tags$h3(class = "ram-chart-label", "Molecular structure"),
                  tags$p(class = "ram-chart-help",
                         "Click a residue to highlight it in the plot and table. Drag to rotate."),
                  tags$div(class = "ram-ngl",
                           NGLVieweR::NGLVieweROutput("NGL")),
                  tags$div(
                    class = "ram-viewer-options",
                    tags$div(class = "ram-viewer-options-title", "Display"),
                    tags$div(
                      class = "ram-toggles",
                      checkboxInput("ligands", "Ligands"),
                      checkboxInput("dna", "DNA"),
                      checkboxInput("rna", "RNA"),
                      checkboxInput("spinning", "Spin"),
                      checkboxInput("rocking", "Rock", value = TRUE)
                    )
                  )
                )
              )
            ),
            tabPanel(
              title = "Residue list", value = "residues",
              tags$div(
                class = "ram-subtab-content",
                tags$div(
                  class = "ram-result-head",
                  tags$div(
                    tags$h2("Residue details"),
                    tags$p("Filter results by region and search individual residues.")
                  )
                ),
                selectInput(
                  "regionselect", "Region",
                  c("All", "Not allowed", "Generously allowed",
                    "Allowed", "Favoured"),
                  selected = "All"
                ),
                DT::DTOutput("regions")
              )
            ),
            tabPanel(
              title = "Summary", value = "summary",
              tags$div(
                class = "ram-subtab-content",
                tags$div(
                  class = "ram-result-head",
                  tags$div(
                    tags$h2("Classification summary"),
                    tags$p("Counts and percentages for the selected residues and chains.")
                  )
                ),
                htmlOutput("summary")
              )
            )
          )
        )
      )
    ),
    tags$footer(
      class = "ram-foot",
      tags$div(
        class = "ram-inner",
        "RamplotR · Interactive Ramachandran analysis · ",
        tags$a(href = "https://github.com/BiKC/RamplotR",
               "Source code", target = "_blank", rel = "noopener noreferrer")
      )
    )
  ),
  tags$script(src = "custom.js")
)
# Structure parsing is deliberately triggered by the Analyse button. Every
# downstream result is a reactive expression, so adjusting settings never
# refetches the structure or recomputes backbone torsions.
server <- function(input, output, session) {
  loaded <- reactiveVal(NULL)
  selected_residue <- reactiveVal(NULL)
  viewer_ready <- reactiveVal(FALSE)

  default_color <- function(index) color_set[((index - 1L) %% length(color_set)) + 1L]
  residue_key <- function(chain, resi, insertion_code = "") {
    paste(as.character(chain), as.integer(resi),
          if (is.na(insertion_code)) "" else as.character(insertion_code),
          sep = "\r")
  }
  selection_string <- function(row) {
    chain <- as.character(row$chain[[1L]])
    insertion <- as.character(row$insertion_code[[1L]])
    # NGL supports insertion codes with ^ and chain IDs with :. These fields
    # come from parsed structure atoms, never from an arbitrary user string.
    if (!grepl("^[A-Za-z0-9_-]*$", chain) ||
        !grepl("^[A-Za-z0-9]*$", insertion)) return("none")
    paste0(as.integer(row$resi[[1L]]),
           if (nzchar(insertion)) paste0("^", insertion) else "",
           if (nzchar(chain)) paste0(":", chain) else "")
  }

  observeEvent(input$bgtype, {
    req(input$bgtype)
    files <- list.files(file.path("static", input$bgtype))
    if (!length(files)) return()
    choices <- list(
      "Commonly used" = files[!files %in% allAA],
      "Per amino acid" = files[files %in% allAA]
    )
    current <- isolate(input$background)
    chosen <- if (!is.null(current) && current %in% files) current else
      if ("General" %in% files) "General" else files[[1L]]
    updateSelectInput(session, "background", choices = choices, selected = chosen)
  })

  update_color_inputs <- function(colors) {
    updateColourInput(session, "bg1", value = colors[[1L]])
    updateColourInput(session, "bg2", value = colors[[2L]])
    updateColourInput(session, "bg3", value = colors[[3L]])
    updateColourInput(session, "bg4", value = colors[[4L]])
  }
  # A preset is applied atomically in the plotted data; picker updates are
  # presentation-only. Custom colours apply when Custom is selected, avoiding
  # circular observers that used to undo a preset during asynchronous updates.
  active_palette <- reactive({
    if (identical(input$colorscheme, "Rampage")) return(unname(rampage))
    if (identical(input$colorscheme, "PDBSum")) return(unname(pdbsum))
    values <- c(input$bg1, input$bg2, input$bg3, input$bg4)
    if (length(values) != 4L || anyNA(values)) return(unname(rampage))
    unname(values)
  })
  observeEvent(input$colorscheme, {
    colors <- switch(input$colorscheme, Rampage = rampage,
                     PDBSum = pdbsum, NULL)
    if (!is.null(colors)) update_color_inputs(colors)
  })
  observeEvent(input$background, {
    # Preserve the existing convenience behaviour for per-amino-acid plots.
    if (!is.null(input$background) && input$background %in% allAA)
      updatePickerInput(session, "AA", selected = input$background)
  }, ignoreInit = TRUE)

  observeEvent(input$submit, {
    type <- if (identical(input$inputSource, "upload")) "file" else "code"
    if (type == "file" && (is.null(input$structfile) ||
                           is.null(input$structfile$datapath))) {
      showNotification("Choose a PDB or mmCIF file first.", type = "error")
      return()
    }
    source_id <- if (type == "file") input$structfile$datapath else
      toupper(trimws(input$PDB))
    key <- paste(type, source_id, sep = ":")
    previous <- isolate(loaded())
    if (!is.null(previous) && identical(previous$key, key)) {
      # The user may click Analyse again, but no expensive reloading is needed.
      return()
    }
    withProgress(message = "Analysing structure", value = 0, {
      incProgress(0.25, detail = "Loading coordinates")
      pdb <- tryCatch(
        ram_load_structure(
          path = if (type == "file") source_id else NULL,
          original_name = if (type == "file") input$structfile$name else NULL,
          pdb_id = if (type == "code") source_id else NULL
        ),
        error = function(e) {
          showNotification(conditionMessage(e), type = "error", duration = 12)
          NULL
        }
      )
      if (is.null(pdb)) return()
      incProgress(0.5, detail = "Calculating backbone geometry")
      torsions <- tryCatch(ram_extract_torsions(pdb),
        error = function(e) {
          showNotification(conditionMessage(e), type = "error", duration = 12)
          NULL
        }
      )
      if (is.null(torsions)) return()
      chains <- unique(torsions$chain)
      name <- if (type == "file")
        tools::file_path_sans_ext(basename(input$structfile$name)) else source_id

      # Invalidate selections before changing the 3D stage, even if a prior
      # structure used the same chain and residue numbering.
      selected_residue(NULL)
      viewer_ready(FALSE)
      viewer_format <- if (type == "file")
        ram_detect_format(input$structfile$name) else NULL
      widget <- NGLVieweR(data = source_id, format = viewer_format) %>%
        NGLVieweR::stageParameters(backgroundColor = "#f7fafb") %>%
        setRock()
      for (k in seq_along(chains)) {
        chain <- chains[[k]]
        widget <- addRepresentation(
          widget, "cartoon",
          param = list(
            name = paste0("ram-chain-", chain),
            sele = paste0(":", chain, " and protein"),
            color = default_color(k)
          )
        )
      }
      # Keep the highlight representation alive, initially selecting nothing.
      # updateSelection changes it without rebuilding the molecular viewer.
      widget <- addRepresentation(widget, "ball+stick",
        param = list(name = "ram-highlight", sele = "none",
                     color = "#ff9e2c", scale = 1.5))
      output$NGL <- NGLVieweR::renderNGLVieweR(widget)

      output$chainColors <- renderUI({
        widgets <- lapply(seq_along(chains), function(k) {
          colourpicker::colourInput(
            paste0("chain", chains[[k]]),
            label = paste("Chain", chains[[k]]), value = default_color(k)
          )
        })
        shinyWidgets::dropdown(widgets, label = "Chain color settings")
      })
      output$chains <- renderUI({
        shinyWidgets::pickerInput(
          "chainselection", "Chain selection", choices = chains,
          multiple = TRUE, options = list("actions-box" = TRUE),
          selected = chains
        )
      })
      loaded(list(key = key, name = name, torsions = torsions, chains = chains))
      incProgress(0.25, detail = "Preparing interactive views")
    })
  }, ignoreInit = TRUE)

  plot_reference <- reactive({
    req(loaded(), input$bgtype, input$background)
    choice <- input$background
    refname <- if (identical(choice, "preProline")) "preProline" else
      if (choice %in% allAA) choice else "General"
    ram_read_reference(file.path("static", input$bgtype, refname))
  })
  classified <- reactive({
    data <- req(loaded())
    req(input$validationMode, input$bgtype)
    ram_classify_torsions(
      data$torsions,
      reference_dir = file.path("static", input$bgtype),
      selected_reference = plot_reference(),
      mode = input$validationMode,
      threshold_fn = ram_density_thresholds
    )
  })
  displayed <- reactive({
    data <- classified()
    # Unlike structure loading, this cheap filtering reacts immediately.
    aa <- input$AA
    data <- data[data$resn %in% aa, , drop = FALSE]
    if (!is.null(input$chainselection))
      data <- data[data$chain %in% input$chainselection, , drop = FALSE]
    if (identical(input$background, "preProline"))
      data <- data[!is.na(data$bonded_to_next) & data$bonded_to_next &
                     !is.na(data$next_resn) & data$next_resn == "PRO",
                   , drop = FALSE]
    data
  })
  table_rows <- reactive({
    data <- displayed()
    region <- input$regionselect
    if (!is.null(region) && !identical(region, "All"))
      data <- data[!is.na(data$region) & data$region == region, , drop = FALSE]
    data
  })

  output$regions <- DT::renderDT({
    data <- table_rows()
    # Keep a previous plot or 3D pick highlighted when the table tab is
    # opened for the first time, before its DT proxy has been initialized.
    selected <- selected_residue()
    marked <- if (is.null(selected)) integer(0) else which(
      data$chain == selected$chain & data$resi == selected$resi &
      data$insertion_code == selected$insertion_code
    )
    columns <- c("chain", "resi", "insertion_code", "resn",
                 "phi", "psi", "region", "density")
    DT::datatable(
      data[, columns, drop = FALSE],
      rownames = FALSE,
      selection = list(mode = "single",
                       selected = if (length(marked)) marked[[1L]] else integer(0)),
      options = list(
        pageLength = 12, scrollX = TRUE, autoWidth = TRUE,
        dom = "ftip", order = list(list(0, "asc"), list(1, "asc")),
        language = list(emptyTable = "No residues match the current filters.")
      ),
      class = "compact stripe hover"
    )
  }, server = FALSE)

  output$summary <- renderUI({
    data <- displayed()
    eligible <- data[!is.na(data$region) & !data$resn %in% c("GLY", "PRO"),
                     , drop = FALSE]
    count <- nrow(eligible)
    count_region <- function(region) sum(eligible$region == region)
    percentage <- function(n) if (!count) "n/a" else sprintf("%.2f%%", 100 * n/count)
    metric <- function(label, value, hint) {
      tags$div(class = "ram-summary-metric",
        tags$span(class = "ram-summary-metric-label", label),
        tags$strong(format(value, big.mark = ",")), tags$small(hint))
    }
    region_row <- function(label, n) {
      tags$tr(tags$td(label),
        tags$td(class = "ram-numeric", format(n, big.mark = ",")),
        tags$td(class = "ram-numeric", percentage(n)))
    }
    outlier <- count_region("Not allowed")
    tags$div(class = "ram-summary",
      tags$div(class = "ram-summary-metrics",
        metric("Selected residues", nrow(data), "Across selected chains"),
        metric("Classified, excluding Gly/Pro", count,
               "Residues with defined backbone angles"),
        metric("Outliers", outlier, "Outside the selected reference regions")
      ),
      tags$h3("Region breakdown"),
      tags$p(class = "ram-summary-note",
        "Percentages use classified residues other than glycine and proline as the denominator."),
      tags$div(class = "ram-summary-table-wrap",
        tags$table(class = "ram-summary-table",
          tags$thead(tags$tr(tags$th("Region"), tags$th("Residues"), tags$th("Share"))),
          tags$tbody(
            region_row("Favoured", count_region("Favoured")),
            region_row("Allowed", count_region("Allowed")),
            region_row("Generously allowed", count_region("Generously allowed")),
            region_row("Not allowed", outlier),
            region_row("Total classified", count)
          ))),
      tags$div(class = "ram-summary-footnotes",
        tags$div(
          tags$strong(format(sum(is.na(data$region) &
            !data$resn %in% c("GLY", "PRO")), big.mark = ",")),
          tags$span("Missing or terminal angles (excluding Gly/Pro)")
        ),
        tags$div(tags$strong(format(sum(data$resn == "GLY"), big.mark = ",")),
                 tags$span("Glycine residues")),
        tags$div(tags$strong(format(sum(data$resn == "PRO"), big.mark = ",")),
                 tags$span("Proline residues"))
      )
    )
  })

  # Dynamic chain colour controls only exist after loading a structure.
  # Tracking them here makes both the plot and molecular viewer repaint.
  current_chain_colors <- reactive({
    data <- req(loaded())
    vapply(seq_along(data$chains), function(k) {
      value <- input[[paste0("chain", data$chains[[k]])]]
      if (is.null(value) || length(value) != 1L ||
          is.na(value) || !grepl("^#[[:xdigit:]]{6}$", value))
        default_color(k) else value
    }, character(1))
  })
  observe({
    data <- req(loaded())
    req(viewer_ready())
    colors <- current_chain_colors()
    for (k in seq_along(data$chains))
      NGLVieweR_proxy("NGL") %>% updateColor(
        name = paste0("ram-chain-", data$chains[[k]]), color = colors[[k]]
      )
  })

  plot_state <- reactive({
    data <- req(loaded())
    list(
      df = displayed(),
      matrix = plot_reference(),
      name = data$name,
      backgroundColors = active_palette(),
      chainColors = {
        colors <- current_chain_colors()
        colors[match(unique(displayed()$chain), data$chains)]
      },
      limits = ram_density_thresholds(plot_reference())
    )
  }) %>% shiny::debounce(100)
  observeEvent(plot_state(), {
    session$sendCustomMessage("process", plot_state())
  }, ignoreInit = FALSE)

  resolve_selection <- function(pick) {
    data <- isolate(displayed())
    if (!is.list(pick) || is.null(pick$chain) || is.null(pick$resi) ||
        length(pick$chain) != 1L || length(pick$resi) != 1L)
      return(NULL)
    number <- suppressWarnings(as.integer(pick$resi))
    if (is.na(number)) return(NULL)
    chain <- as.character(pick$chain)
    insertion <- if (is.null(pick$insertion_code)) "" else
      as.character(pick$insertion_code)
    if (length(insertion) != 1L || is.na(insertion)) return(NULL)
    matches <- which(data$chain == chain & data$resi == number &
                       data$insertion_code == insertion)
    if (!length(matches)) return(NULL)
    data[matches[[1L]], , drop = FALSE]
  }
  select_from <- function(pick) {
    row <- resolve_selection(pick)
    if (is.null(row)) return()
    selected_residue(list(
      chain = as.character(row$chain[[1L]]),
      resi = as.integer(row$resi[[1L]]),
      insertion_code = as.character(row$insertion_code[[1L]])
    ))
  }
  observeEvent(input$ramPlotPick, select_from(input$ramPlotPick))
  observeEvent(input$ramNglPick, select_from(input$ramNglPick))
  observeEvent(input$regions_row_last_clicked, {
    rows <- table_rows()
    index <- suppressWarnings(as.integer(input$regions_row_last_clicked))
    if (length(index) != 1L || is.na(index) || index < 1L ||
        index > nrow(rows)) return()
    row <- rows[index, , drop = FALSE]
    select_from(list(chain = row$chain[[1L]], resi = row$resi[[1L]],
                     insertion_code = row$insertion_code[[1L]]))
  })
  observeEvent(input$clearResidue, selected_residue(NULL))

  # Selection is a single source of truth for the plot, table and 3D viewer.
  selected_row <- reactive({
    selected <- selected_residue()
    if (is.null(selected)) return(NULL)
    data <- displayed()
    ix <- which(data$chain == selected$chain & data$resi == selected$resi &
                  data$insertion_code == selected$insertion_code)
    if (!length(ix)) return(NULL)
    data[ix[[1L]], , drop = FALSE]
  })
  observe({
    selection <- selected_residue()
    if (!is.null(selection) && is.null(selected_row())) selected_residue(NULL)
  })
  output$selectedResidueInfo <- renderUI({
    row <- selected_row()
    if (is.null(row)) return(tags$span("Click a residue in any view to inspect it."))
    tags$strong(
      sprintf("%s %d%s · %s · %s",
        if (nzchar(row$chain[[1L]])) paste("Chain", row$chain[[1L]]) else "Chain",
        as.integer(row$resi[[1L]]), row$insertion_code[[1L]],
        row$resn[[1L]],
        if (is.na(row$region[[1L]])) "Angles unavailable" else row$region[[1L]])
    )
  })
  observe({
    row <- selected_row()
    session$sendCustomMessage("ram-selection", if (is.null(row)) NULL else list(
      chain = as.character(row$chain[[1L]]),
      resi = as.integer(row$resi[[1L]]),
      insertion_code = as.character(row$insertion_code[[1L]]),
      resn = as.character(row$resn[[1L]]),
      region = as.character(row$region[[1L]]),
      phi = as.numeric(row$phi[[1L]]),
      psi = as.numeric(row$psi[[1L]])
    ))
    if (!is.null(isolate(loaded())) && isTRUE(viewer_ready())) {
      sele <- if (is.null(row)) "none" else selection_string(row)
      NGLVieweR_proxy("NGL") %>% updateSelection(
        name = "ram-highlight", sele = sele)
    }
    rows <- table_rows()
    index <- if (is.null(row)) integer(0) else which(
      rows$chain == row$chain[[1L]] & rows$resi == row$resi[[1L]] &
      rows$insertion_code == row$insertion_code[[1L]])
    DT::selectRows(DT::dataTableProxy("regions", session = session),
                   if (length(index)) index[[1L]] else integer(0))
  })

  observeEvent(input$NGL_rendering, {
    if (!identical(input$NGL_rendering, FALSE) ||
        is.null(isolate(loaded()))) return()
    viewer_ready(TRUE)
    session$sendCustomMessage("ram-bind-ngl", list())
  })

  # Existing molecular toggles operate independently of plot redraws.
  observeEvent(input$ligands, {
    if (is.null(loaded()) || !isTRUE(viewer_ready())) return()
    if (isTRUE(input$ligands))
      NGLVieweR_proxy("NGL") %>% addSelection(
        "ball+stick", param = list(name = "ligand", sele = "ligand"))
    else NGLVieweR_proxy("NGL") %>% removeSelection("ligand")
  })
  observeEvent(input$dna, {
    if (is.null(loaded()) || !isTRUE(viewer_ready())) return()
    if (isTRUE(input$dna))
      NGLVieweR_proxy("NGL") %>% addSelection(
        "cartoon", param = list(name = "dna", sele = "dna"))
    else NGLVieweR_proxy("NGL") %>% removeSelection("dna")
  })
  observeEvent(input$rna, {
    if (is.null(loaded()) || !isTRUE(viewer_ready())) return()
    if (isTRUE(input$rna))
      NGLVieweR_proxy("NGL") %>% addSelection(
        "cartoon", param = list(name = "rna", sele = "rna"))
    else NGLVieweR_proxy("NGL") %>% removeSelection("rna")
  })
  observeEvent(input$rocking, {
    if (is.null(loaded()) || !isTRUE(viewer_ready())) return()
    NGLVieweR_proxy("NGL") %>% updateRock(rock = isTRUE(input$rocking))
    if (isTRUE(input$rocking) && isTRUE(input$spinning))
      updateCheckboxInput(session, "spinning", value = FALSE)
  })
  observeEvent(input$spinning, {
    if (is.null(loaded()) || !isTRUE(viewer_ready())) return()
    NGLVieweR_proxy("NGL") %>% updateSpin(spin = isTRUE(input$spinning))
    if (isTRUE(input$spinning) && isTRUE(input$rocking))
      updateCheckboxInput(session, "rocking", value = FALSE)
  })
}

shinyApp(ui = ui, server = server)
