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
                class = "ram-charts",
                tags$section(
                  class = "ram-chart-card", "aria-label" = "Ramachandran plot",
                  tags$h3(class = "ram-chart-label", "Residue distribution"),
                  tags$p(class = "ram-chart-help", "Click a residue to highlight it in the 3D structure and residue list."),
                  tags$p(class = "ram-chart-help",
                         "Hover over a point to identify its chain and residue."),
                  tags$div(
                    id = "plot-empty", class = "ram-plot-empty",
                    tags$span(class = "ram-empty-mark",
                              "aria-hidden" = "true", "φψ"),
                    tags$strong("Your plot will appear here"),
                    tags$p("Enter an accession and select Analyze structure.")
                  ),
                  tags$div(id = "plotly", class = "ram-plot",
                           role = "img", "aria-label" = "Interactive Ramachandran plot"),
                  tags$div(class = "ram-selection-bar",
                    tags$span(id = "ram-selected-residue", "Click a point to select a residue"),
                    actionButton("clearSelection", "Clear selection", class = "btn-default btn-sm")
                  )
                ),
                tags$section(
                  class = "ram-chart-card", "aria-label" = "3D molecular viewer",
                  tags$h3(class = "ram-chart-label", "Molecular structure"),
                  tags$p(class = "ram-chart-help",
                         "Drag to rotate, scroll to zoom."),
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
                tags$p(class = "ram-field-hint", "Click a row to highlight that residue in the plot and molecular viewer."),
                dataTableOutput("regions")
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
# Define server logic required to draw a histogram
server <- function(input, output, session) {
  # Structure parsing and backbone torsions are intentionally independent of
  # visual controls. Only the Analyze button loads a new structure.
  structure <- reactiveVal(NULL)
  torsions <- reactiveVal(NULL)
  structure_name <- reactiveVal(NULL)
  selected_residue <- reactiveVal(NULL)
  viewer_highlight_added <- reactiveVal(FALSE)
  available_chains <- reactiveVal(character())

  output$chains <- renderUI({
    chains <- available_chains()
    if (!length(chains)) return(NULL)
    pickerInput("chainselection", "Chain selection", choices = chains,
      multiple = TRUE, options = list("actions-box" = TRUE), selected = chains)
  })
  output$chainColors <- renderUI({
    chains <- available_chains()
    if (!length(chains)) return(NULL)
    dropdown(lapply(seq_along(chains), function(k) {
      colourpicker::colourInput(
        paste0("chain", chains[[k]]),
        label = paste("Chain", chains[[k]]),
        value = color_set[((k - 1L) %% length(color_set)) + 1L])
    }), label = "Chain color settings")
  })

  observeEvent(input$submit, {
    source_type <- if (identical(input$inputSource, "upload")) "upload" else "pdb"
    if (source_type == "upload" && is.null(input$structfile$datapath)) {
      showNotification("Choose a PDB or mmCIF file first.", type = "error")
      return()
    }
    accession <- toupper(trimws(input$PDB))
    if (source_type == "pdb" && !nzchar(accession)) {
      showNotification("Enter a PDB accession first.", type = "error")
      return()
    }
    withProgress(message = "Loading structure", value = 0, {
      loaded <- tryCatch(ram_load_structure(
        path = if (source_type == "upload") input$structfile$datapath else NULL,
        original_name = if (source_type == "upload") input$structfile$name else NULL,
        pdb_id = if (source_type == "pdb") accession else NULL
      ), error = function(e) {
        showNotification(conditionMessage(e), type = "error", duration = 12)
        NULL
      })
      if (is.null(loaded)) return()
      incProgress(.5, detail = "Calculating backbone torsions")
      angles <- tryCatch(ram_extract_torsions(loaded), error = function(e) {
        showNotification(conditionMessage(e), type = "error", duration = 12)
        NULL
      })
      if (is.null(angles)) return()
      label <- if (source_type == "pdb") accession else
        tools::file_path_sans_ext(input$structfile$name)
      viewer_source <- if (source_type == "pdb") accession else
        input$structfile$datapath
      viewer_format <- if (source_type == "upload")
        ram_detect_format(input$structfile$name) else NULL
      chains <- unique(as.character(angles$chain))
      molecule <- NGLVieweR(data = viewer_source, format = viewer_format) %>%
        NGLVieweR::stageParameters(backgroundColor = "#f7fafb") %>%
        setRock() %>%
        NGLVieweR::selectionParameters(3, "residue")
      for (k in seq_along(chains)) {
        molecule <- addRepresentation(molecule, "cartoon",
          param = list(
            name = paste0("chain-", k),
            sele = paste0(":", chains[k], " and protein"),
            color = color_set[((k - 1L) %% length(color_set)) + 1L]))
      }
      # Clear selection and signal the browser before switching structures.
      selected_residue(NULL)
      viewer_highlight_added(FALSE)
      session$sendCustomMessage("ramplotr-select", NULL)
      structure(loaded)
      torsions(angles)
      available_chains(chains)
      structure_name(label)
      output$NGL <- NGLVieweR::renderNGLVieweR(molecule)
      incProgress(.5)
    })
  })

  # Only the selected background controls the plotted density. All region
  # classifications use per-residue reference grids unless legacy is selected.
  observeEvent(input$bgtype, {
    req(input$bgtype)
    files <- list.files(file.path("static", input$bgtype))
    choices <- list("Commonly used" = files[!files %in% allAA],
                    "Per amino acid" = files[files %in% allAA])
    choice <- if (!is.null(input$background) && input$background %in% files)
      input$background else if ("General" %in% files) "General" else files[1]
    updateSelectInput(session, "background", choices = choices, selected = choice)
  })
  observeEvent(input$background, {
    if (!is.null(input$background) && input$background %in% allAA) {
      updatePickerInput(session, "AA", selected = input$background)
    }
  })
  updateColorInputs <- function(colors) {
    for (k in seq_len(4L))
      colourpicker::updateColourInput(session, paste0("bg", k), value = colors[k])
  }
  observeEvent(input$colorscheme, {
    if (identical(input$colorscheme, "Rampage")) updateColorInputs(rampage)
    if (identical(input$colorscheme, "PDBSum")) updateColorInputs(pdbsum)
  })
  observe({
    colors <- c(input$bg1, input$bg2, input$bg3, input$bg4)
    if (length(colors) != 4L || anyNA(colors)) return()
    name <- if (identical(colors, unname(rampage))) "Rampage" else
      if (identical(colors, unname(pdbsum))) "PDBSum" else "custom"
    if (!identical(input$colorscheme, name))
      updateSelectInput(session, "colorscheme", selected = name)
  })

  selected_matrix <- reactive({
    req(input$bgtype, input$background)
    background <- input$background
    # Some reference sets name this group preProline.
    if (identical(background, "preProline") ||
        identical(background, "Preproline")) background <- "preProline"
    if (!background %in% allAA && background != "preProline")
      background <- "General"
    ram_read_reference(file.path("static", input$bgtype, background))
  })

  classified <- reactive({
    req(torsions(), input$bgtype, input$validationMode)
    ram_classify_torsions(torsions(), file.path("static", input$bgtype),
      selected_matrix(), mode = input$validationMode,
      threshold_fn = ram_density_thresholds)
  })

  visible_rows <- reactive({
    rows <- classified()
    # An explicit empty selection should display no residues; NULL is used
    # only before a dynamic chain-picker has been mounted.
    if (!is.null(input$AA)) rows <- rows[rows$resn %in% input$AA, , drop = FALSE]
    if (!is.null(input$chainselection))
      rows <- rows[rows$chain %in% input$chainselection, , drop = FALSE]
    if (identical(input$background, "preProline"))
      rows <- rows[!is.na(rows$bonded_to_next) & rows$bonded_to_next &
        !is.na(rows$next_resn) & rows$next_resn == "PRO", , drop = FALSE]
    rows
  })

  # Every visual change invalidates this observer; parsing/torsions remain
  # cached. Plotly.react handles repainting without resetting the viewport.
  observe({
    req(structure(), selected_matrix())
    rows <- visible_rows()
    groups <- unique(as.character(rows$chain))
    colors <- vapply(groups, function(chain) {
      widget <- input[[paste0("chain", chain)]]
      if (!is.null(widget) && nzchar(widget)) widget else
        color_set[((match(chain, available_chains()) - 1L) %% length(color_set)) + 1L]
    }, character(1))
    palette <- c(input$bg1, input$bg2, input$bg3, input$bg4)
    req(length(palette) == 4L, all(nzchar(palette)))
    matrix <- selected_matrix()
    session$sendCustomMessage("process", list(
      df = rows, matrix = matrix, name = structure_name(),
      backgroundColors = palette, chainColors = unname(colors),
      limits = ram_density_thresholds(matrix),
      selected = selected_residue()))
  })

  region_rows <- reactive({
    rows <- visible_rows()
    if (!identical(input$regionselect, "All") && !is.null(input$regionselect))
      rows <- rows[!is.na(rows$region) &
        rows$region == input$regionselect, , drop = FALSE]
    rows
  })
  output$regions <- renderDataTable({
    rows <- region_rows()
    data.frame(
      Chain = rows$chain, Residue = rows$resi,
      Insertion = rows$insertion_code, Amino_acid = rows$resn,
      Phi = round(rows$phi, 2), Psi = round(rows$psi, 2),
      Region = rows$region, Density_percentile = round(rows$density, 2),
      check.names = FALSE)
  }, selection = "single", rownames = FALSE, options = list(pageLength = 15,
     scrollX = TRUE))

  # All summary counts use the same live selection as the scatter plot.
  output$summary <- renderUI({
    rows <- visible_rows()
    exclude <- rows$resn %in% c("GLY", "PRO")
    eligible <- rows[!exclude & !is.na(rows$region), , drop = FALSE]
    counts <- vapply(c("Favoured", "Allowed", "Generously allowed", "Not allowed"),
      function(x) sum(eligible$region == x), integer(1))
    percentage <- function(n) if (nrow(eligible)) sprintf("%.2f%%", 100 * n /
      nrow(eligible)) else "n/a"
    region_row <- function(label, number) tags$tr(
      tags$td(label), tags$td(class = "ram-numeric", number),
      tags$td(class = "ram-numeric", percentage(number)))
    metric <- function(label, value, hint) tags$div(class = "ram-summary-metric",
      tags$span(class = "ram-summary-metric-label", label),
      tags$strong(as.character(value)), tags$small(hint))
    tags$div(class = "ram-summary",
      tags$div(class = "ram-summary-metrics",
        metric("Selected residues", nrow(rows), "Across selected chains"),
        metric("Classified, excluding Gly/Pro", nrow(eligible),
          "Residues with defined backbone angles"),
        metric("Outliers", counts[["Not allowed"]],
          "Outside the selected reference regions")),
      tags$h3("Region breakdown"),
      tags$p(class = "ram-summary-note",
        "Percentages use classified residues other than glycine and proline as the denominator."),
      tags$div(class = "ram-summary-table-wrap",
        tags$table(class = "ram-summary-table",
          tags$thead(tags$tr(tags$th("Region"), tags$th("Residues"), tags$th("Share"))),
          tags$tbody(lapply(names(counts), function(label)
            region_row(label, counts[[label]])),
            region_row("Total classified", nrow(eligible))))),
      tags$div(class = "ram-summary-footnotes",
        tags$div(tags$strong(sum(!exclude & is.na(rows$region))),
          tags$span("Missing or terminal angles (excluding Gly/Pro)")),
        tags$div(tags$strong(sum(rows$resn == "GLY")),
          tags$span("Glycine residues")),
        tags$div(tags$strong(sum(rows$resn == "PRO")),
          tags$span("Proline residues"))))
  })

  # A single selected residue drives plot emphasis, 3D highlighting and
  # table selection, whichever representation the user clicked first.
  valid_selection <- function(x) {
    if (is.null(x) || is.null(x$chain) || is.null(x$resi)) return(NULL)
    row <- torsions()
    if (is.null(row)) return(NULL)
    index <- which(as.character(row$chain) == as.character(x$chain) &
       row$resi == suppressWarnings(as.integer(x$resi)) &
       as.character(row$insertion_code) ==
         if (is.null(x$insertion_code)) "" else as.character(x$insertion_code))
    if (!length(index)) return(NULL)
    row[index[1L], , drop = FALSE]
  }
  select_row <- function(row) {
    if (is.null(row) || !nrow(row)) return()
    selected_residue(list(
      chain = as.character(row$chain[1L]),
      resi = as.integer(row$resi[1L]),
      insertion_code = as.character(row$insertion_code[1L]),
      resn = as.character(row$resn[1L])))
  }
  observeEvent(input$ramplotr_point, {
    select_row(valid_selection(input$ramplotr_point))
  })
  observeEvent(input$regions_rows_selected, {
    index <- input$regions_rows_selected
    rows <- region_rows()
    if (length(index) != 1L || index < 1L || index > nrow(rows)) return()
    select_row(rows[index, , drop = FALSE])
  })
  observeEvent(input$NGL_selection, {
    # NGLVieweR emits values such as [ALA]12:A.CA; ignore atom suffixes.
    selection <- input$NGL_selection
    if (!is.character(selection) || !length(selection)) return()
    token <- regmatches(selection[1L],
      regexec("\\[[^]]+\\](-?[0-9]+)(?:\\^([A-Za-z0-9]))?:([^\\.\\s]+)",
              selection[1L], perl = TRUE))[[1L]]
    if (length(token) < 4L) return()
    row <- valid_selection(list(resi = as.integer(token[2L]),
      chain = token[4L], insertion_code = if (length(token) > 2L) token[3L] else ""))
    select_row(row)
  })

  observeEvent(input$clearSelection, {
    selected_residue(NULL)
  })
  observeEvent(selected_residue(), {
    selected <- selected_residue()
    session$sendCustomMessage("ramplotr-select", selected)
    if (is.null(selected) || is.null(torsions())) {
      if (viewer_highlight_added()) {
        NGLVieweR_proxy("NGL") %>% removeSelection("ramplotr-highlight")
        viewer_highlight_added(FALSE)
      }
      return()
    }
    chain <- selected$chain
    if (!grepl("^[A-Za-z0-9]$", chain)) return()
    insertion <- selected$insertion_code
    sele <- paste0(selected$resi,
      if (nzchar(insertion)) paste0("^", insertion) else "",
      ":", chain)
    if (viewer_highlight_added()) {
      NGLVieweR_proxy("NGL") %>% updateSelection("ramplotr-highlight", sele)
    } else {
      NGLVieweR_proxy("NGL") %>% addSelection("ball+stick", param = list(
        name = "ramplotr-highlight", sele = sele, colorValue = "#ff5630",
        colorScheme = "uniform", radiusScale = 1.3))
      viewer_highlight_added(TRUE)
    }
  }, ignoreNULL = FALSE)

  observeEvent(input$ligands, {
    if (isTRUE(input$ligands))
      NGLVieweR_proxy("NGL") %>%
        addSelection("ball+stick", param = list(name = "ligand", sele = "ligand"))
    else NGLVieweR_proxy("NGL") %>% removeSelection("ligand")
  })
  observeEvent(input$dna, {
    if (isTRUE(input$dna))
      NGLVieweR_proxy("NGL") %>%
        addSelection("cartoon", param = list(name = "dna", sele = "dna"))
    else NGLVieweR_proxy("NGL") %>% removeSelection("dna")
  })
  observeEvent(input$rna, {
    if (isTRUE(input$rna))
      NGLVieweR_proxy("NGL") %>%
        addSelection("cartoon", param = list(name = "rna", sele = "rna"))
    else NGLVieweR_proxy("NGL") %>% removeSelection("rna")
  })
  observeEvent(input$rocking, {
    if (isTRUE(input$rocking)) {
      NGLVieweR_proxy("NGL") %>% updateRock()
      if (isTRUE(input$spinning)) updateCheckboxInput(session, "spinning", value = FALSE)
    } else NGLVieweR_proxy("NGL") %>% updateRock(rock = FALSE)
  })
  observeEvent(input$spinning, {
    if (isTRUE(input$spinning)) {
      NGLVieweR_proxy("NGL") %>% updateSpin()
      if (isTRUE(input$rocking)) updateCheckboxInput(session, "rocking", value = FALSE)
    } else NGLVieweR_proxy("NGL") %>% updateSpin(spin = FALSE)
  })
}

shinyApp(ui = ui, server = server)
