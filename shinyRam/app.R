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

# Confidence JSON files can be substantially larger than Shiny's 5 MB
# default upload limit. The parser separately rejects JSON above 32 MB.
options(shiny.maxRequestSize = 40 * 1024^2)

# Used for processing data

source(file.path("R", "reference-loader.R"), local = TRUE)
source(file.path("R", "ramachandran.R"), local = TRUE)
source(file.path("R", "backbone.R"), local = TRUE)
source(file.path("R", "io.R"), local = TRUE)
source(file.path("R", "inspection.R"), local = TRUE)
source(file.path("R", "reports.R"), local = TRUE)
source(file.path("R", "predictions.R"), local = TRUE)
source(file.path("R", "geometry.R"), local = TRUE)
source(file.path("R", "experimental.R"), local = TRUE)
source(file.path("R", "ensemble.R"), local = TRUE)

# Chain colours and contour colours are designed together for a recognisable
# RamplotR publication identity. Region meaning is encoded by ordered contrast,
# not by hue alone; legacy presets remain available.
color_set <- c(
  "#CE6A4D", "#317E9A", "#8065A3", "#B68A3E",
  "#498777", "#B65B7A", "#5070A0", "#826B4B"
)
ramplotr_palette <- c("#FFF8ED", "#D4ECE7", "#7DB9B5", "#126E74")

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
    tags$link(rel = "icon", type = "image/svg+xml", href = "favicon.svg")
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
            tags$h2(
              tags$span(class = "ram-heading-full", "Load a structure"),
              tags$span(class = "ram-heading-compact", "Structure")
            ),
            tags$p("Enter a PDB accession or use a local PDB/mmCIF file.")
          )
        ),
        tags$div(
          class = "ram-source-controls",
          tags$div(
            class = "ram-source-choice",
            radioButtons(
              "inputSource", "Structure source",
              choices = c("PDB ID" = "pdb", "Upload file" = "upload",
                           "AlphaFold DB" = "afdb"),
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
            id = "ram-afdb-wrap", class = "ram-source-picker is-hidden",
            textInput("afdbAccession", "AlphaFold DB UniProt accession",
                      placeholder = "e.g. P69905")
          ),
          tags$div(
            class = "ram-submit",
            actionButton("submit", "Analyze structure", class = "btn-primary")
          )
        ),
        tags$details(id = "ram-prediction-upload", class = "ram-prediction-upload is-hidden",
          tags$summary(class = "ram-prediction-summary", "Prediction settings (optional)"),
          tags$div(class = "ram-prediction-upload-fields",
            selectInput("predictionSource", "Uploaded structure type",
              choices = c("Experimental or unknown (no confidence)" = "experimental",
                          "AlphaFold 2 / ColabFold" = "alphafold2",
                          "AlphaFold 3" = "alphafold3",
                          "ESMFold" = "esmfold",
                          "Other predicted model (B-factor pLDDT)" = "other_prediction"),
              selected = "experimental", selectize = FALSE),
            tags$div(id = "ram-confidence-sidecars",
              class = "ram-prediction-sidecars is-hidden",
              fileInput("predictionJson", "PAE / full confidence JSON",
                accept = c(".json")),
              fileInput("predictionSummaryJson", "Summary JSON (AF3, optional)",
                accept = c(".json"))
            )
          ),
          tags$p(class = "ram-field-hint",
            "For declared predictions, pLDDT is read from B-factors or verified AF3 atom confidence. AlphaFold supports optional PAE; ESMFold provides local pLDDT but not native PAE.")
        )
      ),
      tags$section(
        class = "ram-global-inspector is-empty", "aria-label" = "Residue inspection",
        tags$div(class = "ram-inspector-copy", "aria-live" = "polite",
                 uiOutput("selectedResidueInfo")),
        tags$div(class = "ram-inspector-actions",
          actionButton("prevReview", "Previous issue", class = "btn-default btn-sm"),
          actionButton("nextReview", "Review issues", class = "btn-default btn-sm"),
          actionButton("showInPlot", "Show in plot", class = "btn-primary btn-sm"),
          actionButton("clearResidue", "Clear", class = "btn-default btn-sm")
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
            uiOutput("modelControl"),
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
                   "Residue-aware mode uses separate RamplotR reference distributions."),
            tags$details(class = "ram-details",
              tags$summary("How do these regions differ from MolProbity?"),
              tags$p(class = "ram-field-hint",
                "RamplotR reference groups and contours are not identical to official wwPDB/MolProbity criteria. Avoid interpreting these labels as official structural-validation results."),
              tags$a(href = "https://github.com/BiKC/RamplotR/blob/main/docs/phase-a-results.md",
                     target = "_blank", rel = "noopener noreferrer",
                     "Independent validation results")
            )
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
                tags$p("The RamplotR palette is designed for consistent, identifiable publication figures.")
              )
            ),
            selectInput(
              "colorscheme", "Contour palette",
              choices = c("RamplotR" = "RamplotR", "Rampage" = "Rampage",
                          "PDBsum" = "PDBSum", "Custom colours" = "custom"),
              selected = "RamplotR"
            ),
            tags$details(
              class = "ram-details",
              tags$summary("Customize region colors"),
              tags$p(class = "ram-field-hint",
                     "Select Custom as your contour palette to use these colours."),
              tags$div(
                class = "ram-details-body ram-appearance-controls",
                colourpicker::colourInput("bg1", "Not allowed", value = "#FFF8ED"),
                colourpicker::colourInput("bg2", "Generously allowed", value = "#D4ECE7"),
                colourpicker::colourInput("bg3", "Allowed", value = "#7DB9B5"),
                colourpicker::colourInput("bg4", "Favoured", value = "#126E74")
              )
            ),
            tags$div(class = "ram-panel-divider"),
            tags$div(class = "ram-chain-controls", uiOutput("chainColors"))
          )
        ),
        tags$main(
          class = "ram-main",
          tags$div(class = "ram-workspace-toolbar",
            tags$button(id = "ram-toggle-settings", type = "button",
                        class = "ram-settings-toggle",
                        "Hide settings", "aria-controls" = "ram-settings",
                        "aria-expanded" = "true",
                        title = "Expand the plots by hiding the settings sidebar")
          ),
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
                           role = "img", "aria-label" = "Interactive Ramachandran plot"),
                  # Available in the main analysis view for every chain. The
                  # collapsed position map is compact; expand it for residue
                  # letters and synchronized 2D / 3D selection.
                  tags$details(
                    id = "ram-sequence-panel", class = "ram-sequence-panel",
                    tags$summary(
                      tags$div(class = "ram-sequence-summary-title",
                        tags$strong("Sequence navigator"),
                        tags$span(class = "ram-sequence-summary-hint",
                          "Every selected chain · expand to inspect residues")
                      ),
                      uiOutput("sequenceOverview")
                    ),
                    tags$div(class = "ram-sequence-detail",
                      tags$p(class = "ram-sequence-instruction",
                        "Permanent labels show actual PDB residue numbers every ten positions. Enter a number beside a chain to jump directly to it. The colour behind each letter shows Ramachandran classification; the separate coloured underline and number indicate pLDDT, when available."),
                      tags$div(class = "ram-sequence-legend",
                        tags$span(class="ram-swatch ram-sw-favoured", "Favoured"),
                        tags$span(class="ram-swatch ram-sw-allowed", "Allowed"),
                        tags$span(class="ram-swatch ram-sw-generously-allowed", "Generously allowed"),
                        tags$span(class="ram-swatch ram-sw-outlier", "Outlier"),
                        tags$span(class="ram-swatch ram-sw-missing", "Missing angles")
                      ),
                      tags$div(class="ram-sequence-confidence-key",
                        tags$strong("Model confidence · pLDDT"),
                        tags$span(class="ram-confidence-key-high", "≥90"),
                        tags$span(class="ram-confidence-key-good", "70–89"),
                        tags$span(class="ram-confidence-key-low", "50–69"),
                        tags$span(class="ram-confidence-key-poor", "<50"),
                        tags$span("The number below each amino acid is its pLDDT score.")
                      ),
                      uiOutput("sequenceView")
                    )
                  )
                ),
                tags$section(
                  class = "ram-chart-card", "aria-label" = "3D molecular viewer",
                  tags$div(class = "ram-viewer-header",
                    tags$div(class = "ram-viewer-heading",
                      tags$h3(class = "ram-chart-label", "Molecular structure"),
                      tags$p(class = "ram-chart-help",
                        "Click a residue to inspect it. Drag to rotate, scroll to zoom.")
                    )
                  ),
                  tags$div(class = "ram-ngl",
                    NGLVieweR::NGLVieweROutput("NGL")),
                  tags$div(class = "ram-viewer-options", "aria-label" = "3D view controls",
                    tags$div(class = "ram-viewer-control-heading",
                      tags$span("Representation"),
                      tags$span(class = "ram-viewer-control-note",
                        "Selected residues stay highlighted in orange")
                    ),
                    tags$div(class = "ram-representation",
                      radioButtons("nglRepresentation", label = NULL,
                        choices = c("Cartoon" = "cartoon", "Ribbon" = "ribbon",
                          "Sticks" = "licorice", "Ball & stick" = "ball+stick",
                          "Surface" = "surface"),
                        selected = "cartoon", inline = TRUE)
                    ),
                    tags$div(class = "ram-viewer-layers",
                      tags$div(class = "ram-viewer-control-heading ram-layer-heading",
                        tags$span("Layers & motion")
                      ),
                      tags$div(class = "ram-toggles",
                        checkboxInput("ligands", "Ligands"),
                        checkboxInput("dna", "DNA"),
                        checkboxInput("rna", "RNA"),
                        checkboxInput("spinning", "Spin"),
                        checkboxInput("rocking", "Rock", value = TRUE)
                      ),
                      tags$p(class = "ram-viewer-hint",
                        "Surface rendering can take longer for large structures.")
                    ),
                    tags$details(class="ram-density-panel",
                      tags$summary("Local cryo-EM density map"),
                      tags$p(class="ram-field-hint",
                        "Overlay a local CCP4/MRC map. This is a visual aid, not an experimental map-fit score."),
                      # Insert the file picker after Shiny's initial input
                      # binding so the map stays client-side in the browser.
                      tags$div(id="ram-density-file-slot",
                               class="ram-density-picker"),
                      tags$div(class="ram-density-level",
                        tags$label("Map threshold (σ)", `for`="ram-density-level"),
                        tags$input(type="range",id="ram-density-level",
                                   min="0.5",max="5",step="0.25",value="2"),
                        tags$span(id="ram-density-value","2.0σ")
                      ),
                      tags$div(class="ram-density-buttons",
                        tags$button(type="button",id="ram-density-load",
                                    class="btn btn-primary btn-sm","Show map"),
                        tags$button(type="button",id="ram-density-clear",
                                    class="btn btn-default btn-sm","Remove")
                      ),
                      tags$p(id="ram-density-status",role="status",
                             "No map loaded.")
                    )
                  )
                )
              ),
              uiOutput("predictionPanel"),
              uiOutput("geometryPanel")
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
                tags$div(class = "ram-table-toolbar",
                  selectInput(
                    "regionselect", "Region",
                    c("All", "Not allowed", "Generously allowed",
                      "Allowed", "Favoured"), selected = "All"
                  ),
                  selectInput("reviewFilter", "Review",
                    c("All residues" = "All", "Outliers" = "Outlier",
                      "Missing angles" = "Missing angles",
                      "Near a contour boundary" = "Near boundary"),
                    selected = "All"
                  ),
                  downloadButton("downloadResidues", "Export filtered CSV")
                ),
                tags$p(class = "ram-table-hint",
                  "Click any residue to inspect it. Your selection remains available above every tab."),
                tags$div(class = "ram-residue-table", DT::DTOutput("regions"))
              )
            ),
            tabPanel(
              title = "Compare", value = "compare",
              tags$div(class = "ram-subtab-content",
                tags$div(class = "ram-result-head",
                  tags$div(tags$h2("Compare protein conformations"),
                    tags$p("Align one chain from each structure by amino-acid sequence. Compare angles and classifications, not residue numbers alone."))
                ),
                tags$div(class = "ram-compare-source",
                  radioButtons("compareInputSource", "Comparison input",
                    choices = c("PDB accession" = "pdb", "Uploaded file" = "upload"),
                    selected = "pdb", inline = TRUE),
                  tags$div(id = "ram-compare-pdb", textInput(
                    "comparePDB", "Second structure", value = "1CRN",
                    placeholder = "e.g. 1CRN")),
                  tags$div(id = "ram-compare-upload", class = "is-hidden",
                    fileInput("compareFile", "Second PDB/mmCIF file",
                      accept = c(".pdb", ".ent", ".cif", ".mmcif", ".mcif"))),
                  actionButton("compareSubmit", "Load comparison", class="btn-primary")
                ),
                uiOutput("compareChainControls"),
                tags$div(class = "ram-compare-status", uiOutput("compareSummary")),
                tags$div(class = "ram-compare-toolbar",
                  selectInput("compareJumpSide", "Locate in", c(
                    "Primary chain" = "a", "Comparison chain" = "b")),
                  numericInput("compareJumpResidue", "Residue number",
                    value = NA, min = 1, step = 1, width = "135px"),
                  actionButton("compareJump", "Find aligned pair",
                    class = "btn-primary"),
                  tags$p(class = "ram-compare-toolbar-hint",
                    "Uses actual residue numbers, including alignment gaps.")
                ),
                tags$div(class = "ram-compare-workspace",
                  tags$section(class = "ram-compare-card",
                    tags$div(class = "ram-compare-card-head",
                      tags$h3("Aligned backbone angles"),
                      tags$p("Select either colour to inspect that aligned residue pair.")
                    ),
                    tags$div(id = "comparePlot", class = "ram-compare-plot")
                  ),
                  tags$section(class = "ram-compare-card ram-compare-viewer",
                    tags$div(class = "ram-compare-card-head",
                      tags$h3("3D superposition"),
                      tags$p("Primary chain in coral, comparison chain in blue. Click either structure to inspect aligned residues.")
                    ),
                    tags$div(class = "ram-compare-viewer-controls",
                      checkboxInput("showComparison3D", "Show 3D", value = TRUE),
                      actionButton("compareResetView", "Fit both chains",
                        class = "btn-default btn-sm")
                    ),
                    conditionalPanel(condition = "input.showComparison3D",
                      tags$div(class = "ram-compare-ngl",
                        NGLVieweR::NGLVieweROutput("NGLCompare",
                          height = "410px")))
                  )
                ),
                uiOutput("compareSelectionInfo"),
                tags$div(class = "ram-table-toolbar",
                  selectInput("compareFilter", "Show comparison",
                    choices = c("All aligned residues" = "All",
                      "Changed classification" = "changed",
                      "Angle difference ≥ 30°" = "large",
                      "Insertions / deletions" = "gaps"), selected = "All"),
                  downloadButton("downloadComparison", "Export comparison CSV")
                ),
                tags$p(class = "ram-table-hint",
                  "Select a row to highlight its corresponding residues in both 3D structures and the angle plot."),
                tags$div(class = "ram-residue-table", DT::DTOutput("comparison"))
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
                htmlOutput("summary"),
                uiOutput("ensemblePanel"),
                tags$details(class = "ram-details ram-export-panel",
                  tags$summary("Export figures and a reproducible report"),
                  tags$p(class = "ram-field-hint",
                    "The SVG/PNG figure uses your current palette, chains and selected reference; the HTML report includes all selected residues and analysis provenance."),
                  tags$div(class = "ram-export-actions",
                    downloadButton("downloadSVG", "Vector SVG"),
                    downloadButton("downloadPNG", "High-resolution PNG"),
                    downloadButton("downloadReport", "HTML report")
                  )
                )
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
  tags$script(src = "plotly-loader.js"),
  tags$script(src = "custom.js"),
  tags$script(src = "compare.js"),
  tags$script(src = "prediction.js"),
  tags$script(src = "density.js")
)
# Structure parsing is deliberately triggered by the Analyse button. Every
# downstream result is a reactive expression, so adjusting settings never
# refetches the structure or recomputes backbone torsions.
server <- function(input, output, session) {
  loaded <- reactiveVal(NULL)
  comparison_loaded <- reactiveVal(NULL)
  external_validation <- reactiveVal(NULL)
  ensemble_results <- reactiveVal(NULL)
  selected_residue <- reactiveVal(NULL)
  selected_comparison <- reactiveVal(NULL)
  viewer_ready <- reactiveVal(FALSE)
  current_model <- reactive({
    value <- input$modelChoice
    if (is.null(value) || !nzchar(value)) return(1L)
    suppressWarnings(as.integer(value))
  })
  output$modelControl <- renderUI({
    data <- req(loaded())
    if (data$nmodels <= 1L) return(NULL)
    selectInput("modelChoice", "Structural model",
      choices = stats::setNames(as.character(seq_len(data$nmodels)),
        paste("Model", seq_len(data$nmodels))),
      selected = "1")
  })

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
    files <- ram_reference_choices(file.path("static", input$bgtype))
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
    if (identical(input$colorscheme, "RamplotR")) return(unname(ramplotr_palette))
    if (identical(input$colorscheme, "Rampage")) return(unname(rampage))
    if (identical(input$colorscheme, "PDBSum")) return(unname(pdbsum))
    values <- c(input$bg1, input$bg2, input$bg3, input$bg4)
    if (length(values) != 4L || anyNA(values)) return(unname(ramplotr_palette))
    unname(values)
  })
  observeEvent(input$colorscheme, {
    colors <- switch(input$colorscheme, RamplotR = ramplotr_palette,
                     Rampage = rampage, PDBSum = pdbsum, NULL)
    if (!is.null(colors)) update_color_inputs(colors)
  })
  observeEvent(input$background, {
    # Preserve the existing convenience behaviour for per-amino-acid plots.
    if (!is.null(input$background) && input$background %in% allAA)
      updatePickerInput(session, "AA", selected = input$background)
  }, ignoreInit = TRUE)

  prediction_downloads <- character()
  session$onSessionEnded(function() unlink(prediction_downloads))
  observeEvent(input$submit, {
    source_type <- input$inputSource
    if (!source_type %in% c("pdb", "upload", "afdb")) return()
    is_upload <- identical(source_type, "upload")
    is_afdb <- identical(source_type, "afdb")
    if (is_upload && (is.null(input$structfile) ||
                      is.null(input$structfile$datapath))) {
      showNotification("Choose a PDB or mmCIF file first.", type = "error")
      return()
    }
    source_label <- if (is_upload) input$structfile$datapath else if (is_afdb)
      toupper(trimws(input$afdbAccession)) else toupper(trimws(input$PDB))
    declared_source <- if (is_afdb) "alphafold_db" else if (is_upload)
      input$predictionSource else "experimental"
    if (is.null(declared_source) || !nzchar(declared_source))
      declared_source <- "experimental"
    sidecar <- if (is_upload && !is.null(input$predictionJson))
      input$predictionJson$datapath else ""
    summary_file <- if (is_upload && !is.null(input$predictionSummaryJson))
      input$predictionSummaryJson$datapath else ""
    key <- paste(source_type, source_label, declared_source,
                 sidecar, summary_file, sep = ":")
    previous <- isolate(loaded())
    if (!is.null(previous) && identical(previous$key, key)) return()
    withProgress(message = "Analysing structure", value = 0, {
      incProgress(0.15, detail = "Loading coordinates")
      afdb_files <- NULL
      if (is_afdb) {
        afdb_files <- tryCatch({
          if (!requireNamespace("jsonlite", quietly = TRUE))
            stop("Install jsonlite to retrieve AlphaFold DB structures.")
          ram_download_afdb(ram_afdb_entry(source_label))
        }, error = function(e) {
          showNotification(conditionMessage(e), type = "error", duration = 12)
          NULL
        })
        if (is.null(afdb_files)) return()
        prediction_downloads <<- c(prediction_downloads,
                                     afdb_files$structure, afdb_files$pae)
      }
      source_id <- if (is_afdb) afdb_files$structure else source_label
      original_name <- if (is_afdb) afdb_files$original_name else if (is_upload)
        input$structfile$name else NULL
      pdb <- tryCatch(
        ram_load_structure(
          path = if (is_upload || is_afdb) source_id else NULL,
          original_name = original_name,
          pdb_id = if (identical(source_type, "pdb")) source_id else NULL
        ),
        error = function(e) {
          showNotification(conditionMessage(e), type = "error", duration = 12)
          NULL
        }
      )
      if (is.null(pdb)) return()
      incProgress(0.5, detail = "Calculating backbone geometry")
      torsions <- tryCatch(ram_extract_torsions(ram_model_at(pdb, 1L)),
        error = function(e) {
          showNotification(conditionMessage(e), type = "error", duration = 12)
          NULL
        }
      )
      if (is.null(torsions)) return()
      chains <- unique(torsions$chain)
      name <- if (is_afdb) paste0("AF-", source_label) else if (is_upload)
        tools::file_path_sans_ext(basename(input$structfile$name)) else source_id
      prediction <- NULL
      if (!identical(declared_source, "experimental")) {
        incProgress(0.15, detail = "Reading prediction confidence")
        confidence_file <- if (is_afdb) afdb_files$pae else sidecar
        prediction <- tryCatch(
          ram_prepare_prediction(
            ram_model_at(pdb, 1L), torsions, declared_source,
            sidecar = confidence_file, summary_file = summary_file,
            notes = if (is_afdb) afdb_files$notes else character(),
            model_id = name
          ), error = function(e) {
            showNotification(paste("Prediction confidence:",
              conditionMessage(e)), type = "warning", duration = 15)
            NULL
          }
        )
        if (!is.null(prediction) && is_afdb)
          prediction$confidence_file <- afdb_files$pae_source
        if (!is.null(prediction) && is_upload && !is.null(input$predictionJson))
          prediction$confidence_file <- input$predictionJson$name
        if (!is.null(prediction) && !is.null(confidence_file) &&
            nzchar(confidence_file) && file.exists(confidence_file)) {
          prediction$confidence_md5 <- unname(tools::md5sum(confidence_file))
        }
        if (!is.null(prediction) && length(prediction$notes))
          showNotification(paste(prediction$notes, collapse = " "),
                           type = "warning", duration = 14)
      }

      # Invalidate selections before changing the 3D stage, even if a prior
      # structure used the same chain and residue numbering.
      selected_residue(NULL)
      viewer_ready(FALSE)
      viewer_format <- if (is_upload || is_afdb)
        ram_detect_format(original_name) else NULL
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
      session$sendCustomMessage("ram-clear-density",list())
      loaded(list(key = key, name = name, torsions = torsions, chains = chains,
                  pdb = pdb, nmodels = ram_model_count(pdb), source_id = source_id,
                  viewer_format = viewer_format, prediction = prediction,
                  declared_source = declared_source, input_source = source_type))
      incProgress(0.25, detail = "Preparing interactive views")
    })
  }, ignoreInit = TRUE)

  # Comparison loading is a separate, deliberate action, so changing plot
  # settings does not repeatedly refetch the secondary structure.
  observeEvent(input$compareSubmit, {
    is_upload <- identical(input$compareInputSource, "upload")
    if (is_upload && (is.null(input$compareFile) ||
                      is.null(input$compareFile$datapath))) {
      showNotification("Choose a second PDB or mmCIF file.", type="error")
      return()
    }
    secondary <- if (is_upload) input$compareFile$datapath else
      toupper(trimws(input$comparePDB))
    source_name <- if (is_upload) input$compareFile$name else secondary
    data <- tryCatch(
      ram_load_structure(
        path = if (is_upload) secondary else NULL,
        original_name = if (is_upload) source_name else NULL,
        pdb_id = if (is_upload) NULL else secondary
      ), error = function(e) {
        showNotification(conditionMessage(e), type="error", duration=12)
        NULL
      }
    )
    if (is.null(data)) return()
    torsions <- tryCatch(ram_extract_torsions(ram_model_at(data, 1L)),
      error = function(e) {
        showNotification(conditionMessage(e), type="error",duration=12)
        NULL
      })
    if (is.null(torsions)) return()
    comparison_loaded(list(pdb=data, torsions=torsions,
      name=tools::file_path_sans_ext(basename(source_name)),
      source_id=secondary,
      viewer_format=if (is_upload) ram_detect_format(source_name) else NULL,
      nmodels=ram_model_count(data)))
  }, ignoreInit=TRUE)
  output$compareChainControls <- renderUI({
    first <- req(loaded())
    second <- req(comparison_loaded())
    tags$div(class="ram-compare-chains",
      selectInput("compareChainA", paste("Chain in", first$name),
        choices=first$chains, selected=first$chains[[1L]]),
      selectInput("compareChainB", paste("Chain in", second$name),
        choices=unique(second$torsions$chain),
        selected=unique(second$torsions$chain)[[1L]]),
      if (second$nmodels > 1L)
        selectInput("compareModel", "Second structure model",
          choices=as.character(seq_len(second$nmodels)), selected="1")
    )
  })

  plot_reference <- reactive({
    req(loaded(), input$bgtype, input$background)
    choice <- input$background
    refname <- if (identical(choice, "preProline")) "preProline" else
      if (choice %in% allAA) choice else "General"
    ram_read_reference(file.path("static", input$bgtype, refname))
  })
  model_torsions <- reactive({
    data <- req(loaded())
    model <- current_model()
    if (length(model) != 1L || is.na(model) || model < 1L ||
        model > data$nmodels) return(data$torsions)
    if (model == 1L) return(data$torsions)
    ram_extract_torsions(ram_model_at(data$pdb, model))
  })
  # Extra geometry depends on the selected model, not on plot palettes.
  model_geometry <- reactive({
    structure <- req(loaded())
    ram_extra_geometry(ram_model_at(structure$pdb,current_model()),
                       model_torsions())
  })
  observeEvent(loaded(), external_validation(NULL), ignoreInit=TRUE)
  observeEvent(input$attachValidation, {
    structure <- req(loaded())
    if (!identical(structure$declared_source, "experimental") ||
        identical(structure$input_source, "afdb")) {
      showNotification(
        "Official wwPDB reports apply to their deposited experimental structure, not a predicted model.",
        type = "error", duration = 14)
      return()
    }
    if (!isTRUE(input$confirmValidationSource)) {
      showNotification(
        "Confirm that the official wwPDB report belongs to this exact deposited structure and model.",
        type = "warning", duration = 12)
      return()
    }
    file <- req(input$validationXml)
    record <- tryCatch(ram_external_validation_read(file$datapath),
      error=function(e) {
        showNotification(conditionMessage(e),type="error",duration=12)
        NULL
      })
    if(is.null(record)) return()
    external_validation(list(key=structure$key,records=record,
      name=file$name,md5=unname(tools::md5sum(file$datapath))))
    updateCheckboxInput(session, "confirmValidationSource", value = FALSE)
    showNotification("Official wwPDB annotations attached. Check model and residue coverage.",
                     type="message")
  },ignoreInit=TRUE)
  observeEvent(input$clearValidation, {
    external_validation(NULL)
  },ignoreInit=TRUE)

  classified <- reactive({
    structure <- req(loaded(), input$validationMode, input$bgtype)
    result <- ram_classify_torsions(
      model_torsions(),
      reference_dir = file.path("static", input$bgtype),
      selected_reference = plot_reference(),
      mode = input$validationMode,
      threshold_fn = ram_density_thresholds
    )
    if (current_model() == 1L)
      result <- ram_apply_prediction(result, structure$prediction)
    result <- ram_join_geometry(result,model_geometry())
    official <- external_validation()
    if(!is.null(official) && identical(official$key,structure$key))
      result <- ram_external_validation_join(result,official$records,
                                               model=current_model())
    result
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
  review_queue <- reactive({
    queue <- ram_review_queue(displayed())
    queue[queue$review_status != "Other", , drop = FALSE]
  })
  table_rows <- reactive({
    data <- displayed()
    region <- input$regionselect
    if (!is.null(region) && !identical(region, "All"))
      data <- data[!is.na(data$region) & data$region == region, , drop = FALSE]
    review <- input$reviewFilter
    if (!is.null(review) && !identical(review, "All")) {
      data <- ram_review_queue(data)
      data <- data[data$review_status == review, , drop = FALSE]
    }
    data
  })

  output$regions <- DT::renderDT({
    data <- table_rows()
    # The proxy handles selection without rebuilding DT while a user clicks.
    selected <- isolate(selected_residue())
    marked <- if (is.null(selected)) integer(0) else which(
      data$chain == selected$chain & data$resi == selected$resi &
      data$insertion_code == selected$insertion_code
    )
    columns <- c("chain", "resi", "insertion_code", "resn",
                 "phi", "psi", "region", "density")
    if ("plddt" %in% names(data))
      columns <- c(columns, "plddt", "confidence_category")
    shown <- data[, columns, drop = FALSE]
    if ("plddt" %in% names(shown)) shown$plddt <- round(shown$plddt, 1L)
    shown$phi <- round(shown$phi, 1L)
    shown$psi <- round(shown$psi, 1L)
    shown$density <- round(shown$density, 1L)
    widget <- DT::datatable(
      shown, rownames = FALSE,
      colnames = c("Chain", "Residue", "Ins.", "AA", "Phi (°)", "Psi (°)",
                   "Region", "Percentile",
                   if ("plddt" %in% names(shown)) c("pLDDT", "Confidence")),
      selection = list(mode = "single",
                       selected = if (length(marked)) marked[[1L]] else integer(0)),
      options = list(
        pageLength = 15, scrollX = FALSE, autoWidth = FALSE,
        dom = "ftip", order = list(list(0, "asc"), list(1, "asc")),
        language = list(emptyTable = "No residues match the selected filters.")
      ),
      class = "compact stripe hover"
    )
    DT::formatStyle(widget, "region",
      backgroundColor = DT::styleEqual(
        c("Favoured", "Allowed", "Generously allowed", "Not allowed"),
        c("#D4ECE7", "#E8F3F1", "#FFF5E1", "#FCE5DD")
      ),
      fontWeight = "600"
    )
  }, server = FALSE)
  output$downloadResidues <- downloadHandler(
    filename = function() {
      name <- if (is.null(loaded())) "RamplotR" else loaded()$name
      paste0(gsub("[^A-Za-z0-9_-]", "_", name), "_filtered_residues.csv")
    },
    content = function(file) {
      utils::write.csv(isolate(table_rows()), file, row.names = FALSE, na = "")
    }
  )
  # Publication figures are drawn independently of the browser's plot size;
  # exported data always matches the selected chains, reference and palette.
  safe_filename <- function(extension) {
    current <- isolate(loaded())
    title <- if (is.null(current)) "RamplotR" else current$name
    paste0(gsub("[^A-Za-z0-9_-]", "_", title), "_ramplotr.", extension)
  }
  export_plot <- function(path, format) {
    data <- req(displayed())
    structure <- req(loaded())
    ram_save_figure(
      path, data, plot_reference(), active_palette(),
      stats::setNames(current_chain_colors(), structure$chains),
      format = format, title = paste0(structure$name, " · RamplotR")
    )
  }
  output$downloadSVG <- downloadHandler(
    filename = function() safe_filename("svg"),
    content = function(file) export_plot(file, "svg")
  )
  output$downloadPNG <- downloadHandler(
    filename = function() safe_filename("png"),
    content = function(file) export_plot(file, "png")
  )
  output$downloadReport <- downloadHandler(
    filename = function() safe_filename("html"),
    content = function(file) {
      data <- req(displayed())
      structure <- req(loaded())
      image <- tempfile(fileext = ".svg")
      on.exit(unlink(image), add = TRUE)
      export_plot(image, "svg")
      reference_file <- file.path("static", input$bgtype,
        if (input$background %in% allAA ||
            identical(input$background, "preProline")) input$background
        else "General")
      provenance <- ram_report_metadata(
        structure$name, input$bgtype, input$background,
        input$validationMode, current_model(), reference_file
      )
      if (!is.null(structure$prediction)) {
        provenance$prediction_source <- structure$prediction$source
        provenance$prediction_model <- structure$prediction$model_id
        provenance$confidence_file <- structure$prediction$confidence_file
        if (!is.null(structure$prediction$confidence_md5))
          provenance$confidence_file_md5 <- structure$prediction$confidence_md5
        if (is.finite(structure$prediction$ptm))
          provenance$prediction_pTM <- structure$prediction$ptm
        if (is.finite(structure$prediction$iptm))
          provenance$prediction_ipTM <- structure$prediction$iptm
        provenance$PAE_available <- !is.null(structure$prediction$pae)
        provenance$confidence_limitations <- paste(
          structure$prediction$notes, collapse = "; ")
      }
      ext <- external_validation()
      if(!is.null(ext) && identical(ext$key,structure$key)) {
        provenance$official_wwPDB_report <- ext$name
        provenance$official_wwPDB_md5 <- ext$md5
        provenance$official_wwPDB_model <- current_model()
      }
      provenance$extended_native_geometry <- "Omega and descriptive chi1; not MolProbity-equivalent"
      ram_save_html_report(file, data, provenance, image)
    }
  )

  sequence_groups <- reactive({
    req(loaded())
    # A sequence must retain its true residue positions even if an amino-acid
    # or pre-proline filter limits the points currently drawn in the plot.
    data <- classified()
    if (!is.null(input$chainselection))
      data <- data[data$chain %in% input$chainselection, , drop=FALSE]
    ram_sequence_groups(data)
  })
  # Keep the condensed position maps reactive even while the full sequence
  # navigator is collapsed; all selected chains remain visible.
  output$sequenceOverview <- renderUI({
    groups <- sequence_groups()
    if (!length(groups)) return(tags$span(class="ram-sequence-empty",
      "No residues match the current filters."))
    tags$div(class="ram-sequence-overview", role="group",
      "aria-label"="Selected protein chains and residue classifications",
      lapply(seq_along(groups), function(k) {
        chain <- groups[[k]]
        title <- names(groups)[[k]]
        status <- ram_sequence_status(chain$region)
        bins <- ram_sequence_overview_bins(chain$region)
        tags$div(class="ram-sequence-overview-chain",
          tags$span(class="ram-sequence-chain-name",
            if (identical(title, "Unassigned")) title else paste("Chain", title)),
          tags$div(class="ram-sequence-mini", role="img",
            "aria-label"=sprintf("%s: %d residues, %d outliers, %d missing angles.",
              title, nrow(chain), sum(status=="outlier"), sum(status=="missing")),
            lapply(bins, function(value) {
              tags$span(class=paste("ram-sequence-mini-cell",
                 paste0("ram-seq-",value)), "aria-hidden"="true")
            })
          ),
          if (any(is.finite(chain$plddt))) tags$div(
            class = "ram-confidence-mini", role = "img",
            "aria-label" = sprintf("%s prediction confidence; teal is high, amber/red is low.", title),
            lapply(ram_plddt_overview_bins(chain$plddt), function(value) {
              color <- if (!is.finite(value)) "#cbd7db" else if (value < 50)
                "#d75e56" else if (value < 70) "#d6ac52" else if (value < 90)
                "#7bbcb1" else "#126e74"
              tags$span(class = "ram-confidence-mini-cell",
                style = paste0("background:", color),
                title = if (is.finite(value)) sprintf("Minimum pLDDT %.1f", value)
                        else "Confidence unavailable")
            })
          ),
          tags$span(class="ram-sequence-chain-count",
            sprintf("%s aa",format(nrow(chain),big.mark=",")))
        )
      })
    )
  })
  outputOptions(output,"sequenceOverview",suspendWhenHidden=FALSE)

  output$sequenceView <- renderUI({
    groups <- sequence_groups()
    if (!length(groups)) return(tags$p("No chains match the current selection."))
    shown <- displayed()
    shown_keys <- paste(shown$chain, shown$resi,
                        shown$insertion_code, sep="\r")
    tags$div(class="ram-sequence-chains", role="group",
      "aria-label"="Residue navigation for all selected protein chains",
      lapply(seq_along(groups), function(k) {
        chain <- groups[[k]]
        chain_name <- names(groups)[[k]]
        statuses <- ram_sequence_status(chain$region)
        selectable <- paste(chain$chain, chain$resi,
                            chain$insertion_code, sep="\r") %in% shown_keys
        tags$section(class="ram-sequence-chain",
          tags$div(class="ram-sequence-chain-heading",
            tags$strong(if (identical(chain_name, "Unassigned"))
              chain_name else paste("Chain",chain_name)),
            tags$span(sprintf("%s residues",
              format(nrow(chain),big.mark=","))),
            tags$div(class="ram-sequence-jump",
              tags$label("Go to", class="sr-only"),
              tags$input(type="number", class="ram-seq-jump-input",
                min="1", step="1", placeholder="Residue #",
                "aria-label"=paste("Jump to residue number in",chain_name)),
              tags$button(type="button", class="ram-seq-jump",
                "data-chain"=if (identical(chain_name,"Unassigned")) "" else chain_name,
                "aria-label"=paste("Go to residue in",chain_name), "Go")
            ),
            tags$span(class="ram-sequence-scroll-hint","Scroll sideways →")
          ),
          tags$div(class="ram-sequence-grid", role="group",
            "aria-label"=paste("Select a residue in",chain_name),
            {
              labels <- ram_sequence_position_labels(chain$resi, chain$insertion_code)
              show_confidence <- any(is.finite(chain$plddt))
              lapply(seq_len(nrow(chain)), function(i) {
                residue <- chain[i,,drop=FALSE]
                score <- residue$plddt[[1L]]
                position <- paste0(residue$resi[[1L]],
                                   residue$insertion_code[[1L]])
                tags$div(class="ram-seq-slot",
                  tags$span(class="ram-seq-position",
                    if (nzchar(labels[[i]])) labels[[i]] else "\u00a0",
                    "aria-hidden"="true"),
                  tags$button(type="button",
                    class=paste("ram-seq-res",
                      paste0("ram-seq-",statuses[[i]]),
                      if (show_confidence) "ram-seq-with-confidence" else ""),
                    style=if (is.finite(score))
                      paste0("--ram-plddt-color:",ram_plddt_color(score)) else NULL,
                    disabled=if (!selectable[[i]]) "disabled" else NULL,
                    "data-chain"=residue$chain[[1L]],
                    "data-resi"=residue$resi[[1L]],
                    "data-insertion"=residue$insertion_code[[1L]],
                    title=paste0(residue$resn[[1L]], " ",
                      residue$chain[[1L]], position, " · ",
                      if (is.na(residue$region[[1L]])) "Missing angles"
                      else residue$region[[1L]],
                      if (is.finite(score)) sprintf(" · pLDDT %.1f",score)
                      else "",
                      if (!selectable[[i]]) " · Hidden by current filters"
                      else ""),
                    "aria-pressed"="false",
                    tags$span(class="ram-seq-aa",residue$letter[[1L]]),
                    if (show_confidence) tags$span(class="ram-seq-plddt",
                      if (is.finite(score)) sprintf("%.0f",score) else "\u2014")
                  )
                )
              })
            }
          )
        )
      })
    )
  })

  output$geometryPanel <- renderUI({
    structure <- req(loaded())
    tags$details(id="ram-geometry-panel",class="ram-confidence-panel",
      tags$summary(
        tags$span(class="ram-confidence-title","Extended structure verification"),
        tags$span(class="ram-confidence-subtitle",
          "Peptide and side-chain diagnostics · independent wwPDB evidence")
      ),
      tags$div(class="ram-confidence-body",
        tags$p(class="ram-confidence-explainer",
          "Omega and chi1 are descriptive measurements. Rotamer, clash, bond-angle and experimental-fit assessments are imported from the official wwPDB report when you attach one."),
        uiOutput("geometrySummary"),
        tags$div(class="ram-phase-c-attach",
          fileInput("validationXml","Attach wwPDB validation XML (.xml or .xml.gz)",
                    accept=c(".xml",".gz")),
          checkboxInput("confirmValidationSource",
            "This official report belongs to the loaded deposited experimental structure.",
            value = FALSE),
          actionButton("attachValidation","Attach report",class="btn-primary btn-sm"),
          actionButton("clearValidation","Clear",class="btn-default btn-sm")
        ),
        uiOutput("officialSummary"),
        tags$p(class="ram-confidence-explainer",
          "Independent validation reports are for deposited experimental structures. They cannot validate an unpublished AlphaFold or ESMFold prediction."),
        downloadButton("downloadGeometry","Export detailed residue CSV")
      )
    )
  })
  output$geometrySummary <- renderUI({
    data <- displayed()
    if(!nrow(data)) return(NULL)
    metric <- function(label,count)
      tags$span(class="ram-confidence-metric",paste(label,format(count,big.mark=",")))
    tags$div(class="ram-confidence-metrics",
      metric("Cis peptide bonds",sum(data$omega_status=="Cis",na.rm=TRUE)),
      metric("Twisted peptide bonds",sum(data$omega_status=="Twisted",na.rm=TRUE)),
      metric("Measured χ1 angles",sum(is.finite(data$chi1))),
      metric("Measured Cβ bonds",sum(is.finite(data$cb_ca_distance))),
      metric("Missing ω",sum(!is.finite(data$omega)))
    )
  })
  output$officialSummary <- renderUI({
    ext <- external_validation()
    structure <- req(loaded())
    if(is.null(ext) || !identical(ext$key,structure$key))
      return(tags$p(class="ram-field-hint","No official validation report attached."))
    values <- ram_external_validation_summary(classified())
    tags$div(class="ram-official-summary",
      tags$strong(paste("Official report:",ext$name)),
      tags$p(sprintf("Matched %s of %s residues in model %s. %s independent Ramachandran outliers, %s rotamer outliers and %s residues with local clashes.",
        values$matched,values$total,current_model(),
        values$official_rama_outliers,values$official_rotamer_outliers,
        values$residues_with_clashes)),
      if(values$matched==0L)
        tags$p(class="ram-confidence-warning",
          "No report residues match this model's chain, numbering, insertion codes and residue types."),
      tags$p(class="ram-field-hint",
        "Independent wwPDB values may disagree with RamplotR's reference-specific region labels. Provenance and report checksum are preserved in exports.")
    )
  })
  output$downloadGeometry <- downloadHandler(
    filename=function() safe_filename("extended.csv"),
    content=function(file) utils::write.csv(displayed(),file,row.names=FALSE,na="")
  )

  output$predictionPanel <- renderUI({
    structure <- req(loaded())
    prediction <- structure$prediction
    if (is.null(prediction)) return(NULL)
    if (current_model() > 1L)
      return(tags$section(class = "ram-confidence-panel",
        tags$p("Confidence applies to model 1 only. Switch back to model 1 to view its predictions.")))
    known <- prediction$residues$plddt
    available <- known[is.finite(known)]
    title <- switch(prediction$source,
      alphafold_db = "AlphaFold DB", alphafold2 = "AlphaFold / ColabFold",
      alphafold3 = "AlphaFold 3", esmfold = "ESMFold",
      other_prediction = "Predicted structure", "Predicted structure")
    label <- if (length(available))
      sprintf("Mean pLDDT %.1f · %d / %d residues",
              mean(available), length(available), length(known))
      else "No usable pLDDT values"
    metric <- function(title, value) {
      if (is.finite(value)) tags$span(class = "ram-confidence-metric",
        paste0(title, " ", sprintf("%.2f", value))) else NULL
    }
    tags$details(id = "ram-confidence-panel",
      class = "ram-confidence-panel",
      tags$summary(
        tags$span(class = "ram-confidence-title",
          paste(title, "confidence")),
        tags$span(class = "ram-confidence-subtitle", label),
        if (!is.null(prediction$pae))
          tags$span(class = "ram-confidence-available", "PAE available")
      ),
      tags$div(class = "ram-confidence-body",
        tags$p(class = "ram-confidence-explainer",
          "pLDDT estimates local prediction confidence. PAE estimates uncertainty in relative residue placement. Neither replaces experimental or stereochemical validation."),
        tags$div(class = "ram-confidence-metrics",
          metric("pTM", prediction$ptm),
          metric("ipTM", prediction$iptm),
          uiOutput("predictionReviewMetrics")),
        if (!is.null(prediction$pae)) tagList(
          tags$div(class = "ram-confidence-map-title",
            tags$strong("Predicted aligned error (PAE)"),
            tags$span("Click an axis residue to inspect it in 2D and 3D.")),
          tags$div(id = "ram-pae-plot", role = "img",
            "aria-label" = "Interactive predicted aligned error heatmap"),
          tags$p(class = "ram-pae-note", id = "ram-pae-note")
        ) else tags$p(class = "ram-confidence-explainer",
          "No matching PAE matrix was provided for this model."),
        if (length(prediction$notes)) tags$p(
          class = "ram-confidence-warning",
          paste(prediction$notes, collapse = " "))
      )
    )
  })
  output$predictionReviewMetrics <- renderUI({
    if (is.null(req(loaded())$prediction) || current_model() != 1L)
      return(NULL)
    data <- classified()
    if (!"plddt" %in% names(data)) return(NULL)
    high_outliers <- sum(is.finite(data$plddt) & data$plddt >= 90 &
      !is.na(data$region) & data$region == "Not allowed")
    lower_inrange <- sum(is.finite(data$plddt) & data$plddt < 70 &
      !is.na(data$region) & data$region != "Not allowed")
    tagList(
      tags$span(class = if (high_outliers > 0L)
        "ram-confidence-metric ram-review-high" else "ram-confidence-metric",
        paste(high_outliers, "high-confidence Ramachandran outliers")),
      tags$span(class = "ram-confidence-metric",
        paste(lower_inrange, "lower-confidence residues with in-range geometry"))
    )
  })
  observe({
    structure <- req(loaded())
    prediction <- structure$prediction
    if (is.null(prediction) || current_model() != 1L) {
      session$sendCustomMessage("ram-confidence", list(clear = TRUE))
    } else {
      plot_data <- tryCatch(
        ram_pae_plot_data(prediction, structure$torsions),
        error = function(e) NULL
      )
      session$sendCustomMessage("ram-confidence",
        if (is.null(plot_data)) list(clear = TRUE)
        else plot_data)
    }
  })

  comparison_torsions <- reactive({
    second <- req(comparison_loaded())
    choice <- if (is.null(input$compareModel)) 1L else
      suppressWarnings(as.integer(input$compareModel))
    if (length(choice) != 1L || is.na(choice) || choice <= 1L ||
        choice > second$nmodels) return(second$torsions)
    ram_extract_torsions(ram_model_at(second$pdb, choice))
  })
  comparison_data <- reactive({
    first <- req(loaded())
    second <- req(comparison_loaded())
    req(input$compareChainA, input$compareChainB, input$bgtype,
        input$validationMode)
    original <- classified()
    original <- original[original$chain == input$compareChainA, , drop=FALSE]
    secondary <- ram_classify_torsions(comparison_torsions(),
      reference_dir = file.path("static", input$bgtype),
      selected_reference=plot_reference(),
      mode=input$validationMode,
      threshold_fn=ram_density_thresholds)
    secondary <- secondary[secondary$chain == input$compareChainB, , drop=FALSE]
    if (!nrow(original) || !nrow(secondary))
      return(data.frame())
    result <- ram_compare_torsions(original, secondary)
    result$row_id <- seq_len(nrow(result))
    result
  })
  # A selection identifies a whole aligned pair, not just a numeric residue.
  # It remains stable when the comparison table is filtered or re-ordered.
  choose_comparison <- function(index) {
    data <- isolate(comparison_data())
    row_index <- suppressWarnings(as.integer(index))
    if (length(row_index) != 1L || is.na(row_index) ||
        row_index < 1L || row_index > nrow(data)) return(invisible(FALSE))
    row <- data[row_index, , drop=FALSE]
    selected_comparison(row$row_id[[1L]])
    # Make the shared inspector and main sequence navigator follow the
    # primary chain, without selecting residues hidden by main plot filters.
    if (!is.na(row$residue_a[[1L]])) {
      visible <- isolate(displayed())
      matches <- which(visible$chain == row$chain_a[[1L]] &
        visible$resi == row$residue_a[[1L]] &
        visible$insertion_code == row$insertion_a[[1L]])
      if (length(matches)) selected_residue(list(
        chain = row$chain_a[[1L]],
        resi = as.integer(row$residue_a[[1L]]),
        insertion_code = row$insertion_a[[1L]]
      ))
    }
    invisible(TRUE)
  }
  observeEvent(list(input$compareChainA, input$compareChainB,
                    input$compareModel, comparison_loaded()), {
    selected_comparison(NULL)
  }, ignoreInit=TRUE)
  observeEvent(input$ramComparePlotPick, {
    choose_comparison(input$ramComparePlotPick)
  }, ignoreInit=TRUE)
  observeEvent(input$comparison_row_last_clicked, {
    rows <- filtered_comparison()
    i <- suppressWarnings(as.integer(input$comparison_row_last_clicked))
    if (length(i) != 1L || is.na(i) || i < 1L || i > nrow(rows)) return()
    choose_comparison(rows$row_id[[i]])
  }, ignoreInit=TRUE)
  observeEvent(input$ramCompareNglPick, {
    item <- input$ramCompareNglPick
    if (!is.list(item) || is.null(item$side) ||
        !item$side %in% c("a", "b")) return()
    data <- isolate(comparison_data())
    index <- ram_comparison_find(data, item$side, item$chain,
                                 item$resi, if (is.null(item$insertion_code))
                                   "" else item$insertion_code)
    if (!is.na(index)) choose_comparison(index)
  }, ignoreInit=TRUE)
  observeEvent(input$compareJump, {
    value <- suppressWarnings(as.integer(input$compareJumpResidue))
    if (length(value) != 1L || is.na(value)) {
      showNotification("Enter a valid residue number.", type="warning")
      return()
    }
    data <- isolate(comparison_data())
    side <- isolate(input$compareJumpSide)
    chain <- if (identical(side,"b")) isolate(input$compareChainB) else
      isolate(input$compareChainA)
    index <- ram_comparison_find(data, if (identical(side,"b")) "b" else "a",
                                  chain, value)
    if (is.na(index)) {
      # A jump to 104 should also find insertion-only 104A when necessary.
      positions <- data[[paste0("residue_", if (identical(side,"b")) "b" else "a")]]
      candidates <- which(!is.na(positions) & positions == value &
        data[[paste0("chain_", if (identical(side,"b")) "b" else "a")]] == chain)
      index <- if (length(candidates)) candidates[[1L]] else NA_integer_
    }
    if (is.na(index)) {
      showNotification("That number is not present in the selected chain.",
                       type="warning")
      return()
    }
    choose_comparison(index)
  })
  # Selecting a residue in the primary sequence/plot also locates its aligned
  # partner when the comparison is available.
  observeEvent(selected_residue(), {
    item <- selected_residue()
    if (is.null(item) || is.null(isolate(input$compareChainA)) ||
        !identical(item$chain, isolate(input$compareChainA)) ||
        is.null(isolate(comparison_loaded()))) return()
    index <- ram_comparison_find(isolate(comparison_data()), "a",
               item$chain, item$resi, item$insertion_code)
    if (!is.na(index) &&
        !identical(isolate(selected_comparison()), index))
      selected_comparison(index)
  }, ignoreNULL=TRUE)
  filtered_comparison <- reactive({
    result <- comparison_data()
    if (!nrow(result)) return(result)
    criterion <- input$compareFilter
    if (identical(criterion, "changed"))
      result <- result[result$class_changed, , drop=FALSE]
    else if (identical(criterion, "large"))
      result <- result[(!is.na(result$delta_phi) & abs(result$delta_phi)>=30) |
                       (!is.na(result$delta_psi) & abs(result$delta_psi)>=30),
                       , drop=FALSE]
    else if (identical(criterion, "gaps"))
      result <- result[result$alignment %in% c("Insertion","Deletion"),
                       , drop=FALSE]
    result
  })
  output$compareSummary <- renderUI({
    result <- req(comparison_data())
    if (!nrow(result)) return(tags$p("Select two nonempty protein chains."))
    aligned <- result$alignment %in% c("Match", "Substitution")
    tags$div(class="ram-compare-metrics",
      tags$span(tags$strong(sum(aligned)), " aligned residues"),
      tags$span(tags$strong(sum(result$class_changed)), " region changes"),
      tags$span(tags$strong(sum(!aligned)), " insertions / deletions"),
      tags$span("Angular differences account for the -180° / +180° boundary.")
    )
  })
  output$comparison <- DT::renderDT({
    result <- filtered_comparison()
    fields <- c("chain_a", "residue_a", "insertion_a", "amino_a",
      "chain_b", "residue_b", "insertion_b", "amino_b",
      "delta_phi", "delta_psi", "class_changed", "alignment")
    if (!all(fields %in% names(result)))
      return(DT::datatable(data.frame()))
    shown <- result[, fields, drop=FALSE]
    shown$pos_a <- ifelse(is.na(shown$residue_a), "—",
      paste0(shown$residue_a, shown$insertion_a))
    shown$pos_b <- ifelse(is.na(shown$residue_b), "—",
      paste0(shown$residue_b, shown$insertion_b))
    shown$delta_phi <- round(shown$delta_phi, 1)
    shown$delta_psi <- round(shown$delta_psi, 1)
    shown$class_changed <- ifelse(shown$class_changed, "Yes", "No")
    shown <- shown[, c("chain_a", "pos_a", "amino_a",
      "chain_b", "pos_b", "amino_b", "delta_phi", "delta_psi",
      "class_changed", "alignment"), drop=FALSE]
    DT::datatable(shown, rownames=FALSE,
      colnames=c("Chain A", "Pos A", "AA A", "Chain B", "Pos B", "AA B",
                 "Δφ (°)", "Δψ (°)", "Region changed", "Alignment"),
      selection="single",
      options=list(pageLength=15,scrollX=FALSE,autoWidth=FALSE,dom="ftip"),
      class="compact stripe hover")
  }, server=FALSE)
  observe({
    rows <- req(filtered_comparison())
    index <- selected_comparison()
    pos <- if (is.null(index)) integer() else match(index, rows$row_id)
    DT::selectRows(DT::dataTableProxy("comparison", session = session),
      if (length(pos) && !is.na(pos)) pos else integer())
  })
  output$downloadComparison <- downloadHandler(
    filename=function() "RamplotR_structure_comparison.csv",
    content=function(file) utils::write.csv(
      isolate(filtered_comparison()), file,row.names=FALSE,na="")
  )
  observeEvent(comparison_data(), {
    result <- comparison_data()
    if (!nrow(result)) return()
    session$sendCustomMessage("ram-comparison", list(
      nameA=req(loaded())$name,
      nameB=req(comparison_loaded())$name,
      phiA=result$phi_a, psiA=result$psi_a,
      phiB=result$phi_b, psiB=result$psi_b,
      rowIds=result$row_id,
      chainA=result$chain_a, posA=result$residue_a,
      insA=result$insertion_a, aminoA=result$amino_a,
      chainB=result$chain_b, posB=result$residue_b,
      insB=result$insertion_b, aminoB=result$amino_b,
      deltaPhi=result$delta_phi, deltaPsi=result$delta_psi,
      alignment=result$alignment
    ))
    session$sendCustomMessage("ram-compare-config", list(
      chainA=input$compareChainA, chainB=input$compareChainB,
      modelA=current_model(), modelB=if (is.null(input$compareModel)) 1L
        else as.integer(input$compareModel),
      multipleA=req(loaded())$nmodels > 1L,
      multipleB=req(comparison_loaded())$nmodels > 1L
    ))
  })
  output$compareSelectionInfo <- renderUI({
    data <- req(comparison_data())
    id <- selected_comparison()
    if (!length(id) || is.null(id) || !nrow(data))
      return(tags$div(class="ram-compare-selection ram-compare-selection-empty",
        tags$strong("Inspect an aligned pair"),
        tags$p("Click a point, table row or residue in either 3D structure. Use the position finder for residues such as 104.")))
    match_index <- match(id, data$row_id)
    if (is.na(match_index)) return(NULL)
    row <- data[match_index,,drop=FALSE]
    label <- function(side) {
      number <- row[[paste0("residue_",side)]][[1L]]
      if (is.na(number)) return("Alignment gap")
      paste0(row[[paste0("amino_",side)]][[1L]], " ",
        row[[paste0("chain_",side)]][[1L]], ":",
        number, row[[paste0("insertion_",side)]][[1L]])
    }
    angle <- function(x) if (is.finite(x)) sprintf("%.1f°",x) else "N/A"
    tags$div(class="ram-compare-selection",
      tags$div(class="ram-compare-selection-pair",
        tags$span(class="ram-compare-primary",
          tags$small("Primary"), tags$strong(label("a")),
          tags$span(paste("φ",angle(row$phi_a[[1L]]),
                          "· ψ",angle(row$psi_a[[1L]])))),
        tags$span(class="ram-compare-pair-arrow", "↔", "aria-hidden"="true"),
        tags$span(class="ram-compare-secondary",
          tags$small("Comparison"), tags$strong(label("b")),
          tags$span(paste("φ",angle(row$phi_b[[1L]]),
                          "· ψ",angle(row$psi_b[[1L]]))))
      ),
      tags$div(class="ram-compare-selection-deltas",
        tags$span(paste("Δφ",angle(row$delta_phi[[1L]]))),
        tags$span(paste("Δψ",angle(row$delta_psi[[1L]]))),
        tags$span(row$alignment[[1L]]),
        if (isTRUE(row$class_changed[[1L]])) tags$span(
          class="ram-compare-change", "Classification changed")
      )
    )
  })
  outputOptions(output, "compareSelectionInfo", suspendWhenHidden=FALSE)
  observe({
    id <- selected_comparison()
    data <- comparison_data()
    if (!length(id) || is.null(id) || !nrow(data)) {
      session$sendCustomMessage("ram-comparison-selected", list(clear=TRUE))
      session$sendCustomMessage("ram-compare-pair", list(clear=TRUE))
      return()
    }
    pos <- match(id,data$row_id)
    if (is.na(pos)) return()
    row <- data[pos,,drop=FALSE]
    pair <- function(side,model,multiple) {
      number <- row[[paste0("residue_",side)]][[1L]]
      if (is.na(number)) return(NULL)
      list(chain=as.character(row[[paste0("chain_",side)]][[1L]]),
        resi=as.integer(number),
        insertion_code=as.character(row[[paste0("insertion_",side)]][[1L]]),
        modelIndex=model, multipleModels=multiple)
    }
    second <- req(comparison_loaded())
    session$sendCustomMessage("ram-comparison-selected",
      list(rowId=id))
    session$sendCustomMessage("ram-compare-pair", list(
      a=pair("a",current_model(),req(loaded())$nmodels > 1L),
      b=pair("b",if (is.null(input$compareModel)) 1L else
                      as.integer(input$compareModel),second$nmodels > 1L)
    ))
  })
  output$NGLCompare <- NGLVieweR::renderNGLVieweR({
    req(input$showComparison3D,input$compareChainA,input$compareChainB)
    req(comparison_data())
    first <- req(loaded()); second <- req(comparison_loaded())
    model_a <- if (first$nmodels > 1L)
      paste0(" and /",current_model()-1L) else ""
    model_b <- if (second$nmodels > 1L)
      paste0(" and /",if (is.null(input$compareModel)) 0L else
                     as.integer(input$compareModel)-1L) else ""
    sel_a <- paste0(":", input$compareChainA, model_a, " and protein")
    sel_b <- paste0(":", input$compareChainB, model_b, " and protein")
    widget <- NGLVieweR(data=first$source_id,format=first$viewer_format) %>%
      NGLVieweR::stageParameters(backgroundColor="#f7fafb") %>%
      addRepresentation("cartoon",param=list(
        sele=sel_a,color="#CE6A4D",name="ram-compare-chain-a")) %>%
      addRepresentation("ball+stick",param=list(
        sele="none",color="#ffc04a",scale=1.5,name="ram-compare-highlight-a"))
    widget <- NGLVieweR::addStructure(widget, data=second$source_id,
                                     format=second$viewer_format) %>%
      addRepresentation("cartoon",param=list(
        sele=sel_b,color="#317E9A",name="ram-compare-chain-b")) %>%
      addRepresentation("ball+stick",param=list(
        sele="none",color="#83e6f5",scale=1.5,name="ram-compare-highlight-b"))
    NGLVieweR::setSuperpose(widget, reference=1,
      sele_reference=sel_a, sele_target=sel_b)
  })
  observeEvent(input$NGLCompare_PDB, {
    if (is.null(isolate(comparison_loaded()))) return()
    session$sendCustomMessage("ram-compare-ready", list())
  }, ignoreInit=TRUE)
  observeEvent(input$NGLCompare_rendering, {
    if (!identical(input$NGLCompare_rendering, FALSE) ||
        is.null(isolate(comparison_loaded()))) return()
    # Some NGL versions do not emit a changed PDB input when users switch
    # chains of the same structure. The JS readiness guard checks both
    # structure objects before reframing.
    session$sendCustomMessage("ram-compare-ready", list())
  }, ignoreInit=TRUE)

  output$ensemblePanel <- renderUI({
    structure <- req(loaded())
    if(structure$nmodels<=1L) return(NULL)
    tags$details(id="ram-ensemble-panel",class="ram-confidence-panel",
      tags$summary(
        tags$span(class="ram-confidence-title","Ensemble analysis"),
        tags$span(class="ram-confidence-subtitle",
          paste(structure$nmodels,"structural models · circular φ/ψ variation and region consistency"))
      ),
      tags$div(class="ram-confidence-body",
        tags$p(class="ram-confidence-explainer",
          "Model variation is matched by chain, residue and insertion code. Circular statistics correctly handle the -180°/180° boundary; models with missing coordinates contribute only observed angles."),
        tags$div(class="ram-ensemble-actions",
          actionButton("calculateEnsemble","Analyse ensemble",
                       class="btn-primary btn-sm"),
          downloadButton("downloadEnsemble","Export ensemble CSV")
        ),
        uiOutput("ensembleResultSummary"),
        tags$div(class="ram-residue-table",DT::DTOutput("ensembleRows"))
      )
    )
  })
  ensemble_matches <- reactive({
    value <- ensemble_results()
    if(is.null(value)) return(NULL)
    if(!identical(value$key,req(loaded())$key) ||
       !identical(value$mode,input$validationMode) ||
       !identical(value$reference,input$bgtype) ||
       !identical(value$background,input$background)) return(NULL)
    value$result
  })
  observeEvent(input$calculateEnsemble, {
    structure <- req(loaded())
    if(structure$nmodels<=1L) return()
    withProgress(message="Analysing compatible ensemble models",value=0.2,{
      result <- tryCatch(
        ram_ensemble_analyze(structure$pdb,max_models=min(30L,structure$nmodels),
          classifier=function(torsions)
            ram_classify_torsions(torsions,
              reference_dir=file.path("static",input$bgtype),
              selected_reference=plot_reference(),mode=input$validationMode,
              threshold_fn=ram_density_thresholds)),
        error=function(e) {
          showNotification(conditionMessage(e),type="error",duration=12)
          NULL
        })
      if(!is.null(result))
        ensemble_results(list(key=structure$key,mode=input$validationMode,
          reference=input$bgtype,background=input$background,result=result))
      incProgress(0.8)
    })
  },ignoreInit=TRUE)
  output$ensembleResultSummary <- renderUI({
    result <- ensemble_matches()
    if(is.null(result)) return(tags$p(class="ram-field-hint",
      "Run the ensemble analysis. Results are recalculated on request after changing the reference or classification settings."))
    data <- result$summary
    tags$div(class="ram-confidence-metrics",
      tags$span(class="ram-confidence-metric",
        sprintf("%s of %s models analysed",result$analyzed_models,
                result$available_models)),
      tags$span(class="ram-confidence-metric",
        sprintf("%s residues with classification changes",
                sum(data$changes_class,na.rm=TRUE))),
      tags$span(class="ram-confidence-metric",
        sprintf("%s residues with ≥20° angular spread",
                sum(pmax(data$phi_sd,data$psi_sd,na.rm=TRUE)>=20,
                    na.rm=TRUE))),
      if(result$limited)
        tags$span(class="ram-confidence-warning",
          "Only the first 30 models are included; export records this limit.")
    )
  })
  output$ensembleRows <- DT::renderDT({
    result <- ensemble_matches()
    req(result)
    data <- result$summary
    if(!nrow(data)) return(DT::datatable(data,rownames=FALSE))
    fields <- c("chain","resi","insertion_code","resn",
      "phi_mean","phi_sd","psi_mean","psi_sd","models_present",
      "class_consistency","changes_class")
    shown <- data[,fields,drop=FALSE]
    for(field in c("phi_mean","phi_sd","psi_mean","psi_sd"))
      shown[[field]] <- round(shown[[field]],1L)
    shown$class_consistency <- round(shown$class_consistency*100,1L)
    shown$changes_class <- ifelse(shown$changes_class,"Changed","Stable")
    DT::datatable(shown,rownames=FALSE,selection="single",
      colnames=c("Chain","Residue","Ins.","AA","φ mean","φ SD",
        "ψ mean","ψ SD","Models","Class agreement (%)","Class"),
      options=list(pageLength=12,autoWidth=FALSE,scrollX=TRUE,
                   dom="ftip",order=list(list(10,"asc"))),
      class="compact stripe hover")
  },server=FALSE)
  observeEvent(input$ensembleRows_rows_selected, {
    data <- req(ensemble_matches())$summary
    ix <- input$ensembleRows_rows_selected[[1L]]
    if(!is.finite(ix) || ix<1L || ix>nrow(data)) return()
    row <- data[ix,,drop=FALSE]
    selected_residue(list(chain=as.character(row$chain[[1L]]),
      resi=as.integer(row$resi[[1L]]),
      insertion_code=as.character(row$insertion_code[[1L]])))
  })
  output$downloadEnsemble <- downloadHandler(
    filename=function() safe_filename("ensemble.csv"),
    content=function(file) utils::write.csv(req(ensemble_matches())$summary,
                                               file,row.names=FALSE,na="")
  )

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
  observeEvent(input$ramSeqPick, select_from(input$ramSeqPick))
  observeEvent(input$ramPaePick, select_from(input$ramPaePick))
  observeEvent(input$showInPlot, {
    updateTabsetPanel(session, "analysisTabs", selected = "plot")
  })
  advance_review <- function(direction) {
    queue <- isolate(review_queue())
    if (!nrow(queue)) {
      showNotification("No residues need review in the current selection.", type="message")
      return()
    }
    current <- isolate(selected_residue())
    index <- if (is.null(current)) integer(0) else which(
      queue$chain == current$chain & queue$resi == current$resi &
      queue$insertion_code == current$insertion_code
    )
    next_index <- if (!length(index)) 1L else
      ((index[[1L]] - 1L + direction + nrow(queue)) %% nrow(queue)) + 1L
    row <- queue[next_index, , drop = FALSE]
    selected_residue(list(chain = row$chain[[1L]], resi = row$resi[[1L]],
                         insertion_code = row$insertion_code[[1L]]))
  }
  observeEvent(input$nextReview, advance_review(1L))
  observeEvent(input$prevReview, advance_review(-1L))
  observeEvent(input$clearResidue, selected_residue(NULL))

  # Pause viewer motion when the user starts inspecting a particular
  # residue; keep the 3D control checkboxes in sync with that camera state.
  observeEvent(selected_residue(), {
    if (isTRUE(isolate(input$rocking)))
      updateCheckboxInput(session, "rocking", value = FALSE)
    if (isTRUE(isolate(input$spinning)))
      updateCheckboxInput(session, "spinning", value = FALSE)
  }, ignoreNULL = TRUE)

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
    if (is.null(row)) return(tags$span(class = "ram-inspector-empty",
      "Select a residue in the plot, table, sequence or 3D model. Review controls navigate to the next issue."))
    angle <- function(value) if (is.finite(value)) sprintf("%.1f°", value) else "Unavailable"
    tags$div(class = "ram-inspector-data",
      tags$div(tags$strong(sprintf("%s %d%s · %s",
        if (nzchar(row$chain[[1L]])) paste("Chain", row$chain[[1L]]) else "Chain",
        as.integer(row$resi[[1L]]), row$insertion_code[[1L]],
        row$resn[[1L]])),
        tags$span(class = "ram-inspector-classification",
          if (is.na(row$region[[1L]])) "Missing angles" else row$region[[1L]])),
      tags$div(class = "ram-inspector-angles",
        tags$span(paste("φ", angle(row$phi[[1L]]))),
        tags$span(paste("ψ", angle(row$psi[[1L]]))),
        tags$span(if (is.finite(row$density[[1L]]))
          sprintf("Density percentile %.1f", row$density[[1L]]) else ""),
        if ("omega" %in% names(row) && is.finite(row$omega[[1L]]))
          tags$span(sprintf("ω %.1f° · %s",row$omega[[1L]],
                            row$omega_status[[1L]])),
        if ("chi1" %in% names(row) && is.finite(row$chi1[[1L]]))
          tags$span(sprintf("χ1 %.1f°",row$chi1[[1L]])),
        if ("wwpdb_rotamer" %in% names(row) &&
            !is.na(row$wwpdb_rotamer[[1L]]))
          tags$span(class="ram-inspector-plddt",
                    paste("wwPDB rotamer",row$wwpdb_rotamer[[1L]])),
        if ("wwpdb_clashes" %in% names(row) &&
            is.finite(row$wwpdb_clashes[[1L]]) &&
            row$wwpdb_clashes[[1L]]>0)
          tags$span(class="ram-inspector-warning",
                    paste("Official wwPDB clashes",row$wwpdb_clashes[[1L]])),
        if ("plddt" %in% names(row) && is.finite(row$plddt[[1L]]))
          tags$span(class = "ram-inspector-plddt",
            sprintf("pLDDT %.1f · %s", row$plddt[[1L]],
                    row$confidence_category[[1L]])),
        if ("plddt" %in% names(row) && is.finite(row$plddt[[1L]]) &&
            row$plddt[[1L]] >= 90 &&
            identical(as.character(row$region[[1L]]), "Not allowed"))
          tags$span(class = "ram-inspector-warning",
            "High model confidence, unusual backbone geometry; inspect locally.")
      )
    )
  })
  # Keep the shared selection label current while the plot tab is hidden
  # (for example when a residue is selected from the DataTable).
  outputOptions(output, "selectedResidueInfo", suspendWhenHidden = FALSE)
  observe({
    row <- selected_row()
    session$sendCustomMessage("ram-confidence-selected",
      if (is.null(row)) list(clear = TRUE) else list(
        chain = as.character(row$chain[[1L]]),
        resi = as.integer(row$resi[[1L]]),
        insertion_code = as.character(row$insertion_code[[1L]])))
    session$sendCustomMessage("ram-selection", if (is.null(row)) list(clear = TRUE) else list(
      chain = as.character(row$chain[[1L]]),
      resi = as.integer(row$resi[[1L]]),
      insertion_code = as.character(row$insertion_code[[1L]]),
      resn = as.character(row$resn[[1L]]),
      region = as.character(row$region[[1L]]),
      phi = as.numeric(row$phi[[1L]]),
      psi = as.numeric(row$psi[[1L]]),
      modelIndex = current_model(),
      multipleModels = isolate(loaded())$nmodels > 1L
    ))
    if (!is.null(isolate(loaded())) && isTRUE(viewer_ready())) {
      sele <- if (is.null(row)) "none" else selection_string(row)
      if (!is.null(row) && isolate(loaded())$nmodels > 1L)
        sele <- paste0(sele, " and /", current_model() - 1L)
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
    first_ready <- !isTRUE(isolate(viewer_ready()))
    viewer_ready(TRUE)
    # The NGL widget reuses its stage when a second PDB is loaded. Reframe
    # the new structure after it has finished loading, not at submit time.
    session$sendCustomMessage("ram-bind-ngl",
      list(resetView = first_ready))
  })

  # Swap only named chain representations. The named orange highlight
  # representation survives changes to the whole-structure rendering mode.
  observe({
    data <- req(loaded())
    req(viewer_ready(), input$nglRepresentation)
    style <- input$nglRepresentation
    model <- current_model()
    if (is.na(model)) model <- 1L
    suffix <- if (data$nmodels > 1L)
      paste0(" and /", model - 1L) else ""
    if (!style %in% c("cartoon", "ribbon", "licorice", "ball+stick", "surface"))
      return()
    colors <- isolate(current_chain_colors())
    for (k in seq_along(data$chains)) {
      chain <- data$chains[[k]]
      name <- paste0("ram-chain-", chain)
      proxy <- NGLVieweR_proxy("NGL")
      proxy %>% removeSelection(name)
      proxy %>% addSelection(style, param = list(
        name = name, sele = paste0(":", chain, " and protein", suffix),
        color = colors[[k]],
        opacity = if (identical(style, "surface")) 0.8 else 1
      ))
    }
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
