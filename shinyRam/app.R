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

# Confidence JSON and multi-file prediction ensembles can exceed Shiny's
# 5 MB default upload limit. JSON parsing still enforces its own 32 MB limit.
options(shiny.maxRequestSize = 120 * 1024^2)

# Used for processing data

source(file.path("R", "reference-loader.R"), local = TRUE)
source(file.path("R", "ramachandran.R"), local = TRUE)
source(file.path("R", "rama8000.R"), local = TRUE)
source(file.path("R", "backbone.R"), local = TRUE)
source(file.path("R", "conformation.R"), local = TRUE)
source(file.path("R", "canonical.R"), local = TRUE)
source(file.path("R", "atlas-sifts.R"), local = TRUE)
source(file.path("R", "atlas-cohort.R"), local = TRUE)
source(file.path("R", "atlas-geometry.R"), local = TRUE)
source(file.path("R", "atlas-robustness.R"), local = TRUE)
source(file.path("R", "atlas-construct.R"), local = TRUE)
source(file.path("R", "atlas-experimental-context.R"), local = TRUE)
source(file.path("R", "atlas-ligand-evidence.R"), local = TRUE)
source(file.path("R", "atlas-group-ligand.R"), local = TRUE)
source(file.path("R", "atlas-switch.R"), local = TRUE)
source(file.path("R", "atlas.R"), local = TRUE)
source(file.path("R", "io.R"), local = TRUE)
source(file.path("R", "inspection.R"), local = TRUE)
source(file.path("R", "reports.R"), local = TRUE)
source(file.path("R", "predictions.R"), local = TRUE)
source(file.path("R", "geometry.R"), local = TRUE)
source(file.path("R", "experimental.R"), local = TRUE)
source(file.path("R", "ensemble.R"), local = TRUE)
source(file.path("R", "group-comparison.R"), local = TRUE)
source(file.path("R", "group-fingerprint.R"), local = TRUE)
source(file.path("R", "atlas-group-handoff.R"), local = TRUE)
source(file.path("R", "guide.R"), local = TRUE)

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
        tags$div(class="ram-source-guide",
          actionLink("openGuide","New here? Explore the workflow guide")),
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
            tags$div(class="ram-standard-validation-note",
              tags$strong("Standard validation · Rama8000"),
              tags$p(class="ram-field-hint",
                "Always calculated with the current six-class cctbx/Phenix reference model. Favored, Allowed and Outlier are shown separately from RamplotR density regions.")
            ),
            selectInput(
              "validationMode", "RamplotR density classification",
              c("Residue-aware" = "residue",
                "Selected plotting background" = "legacy")
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
                        "Permanent labels show actual PDB residue numbers every ten positions. Enter a number beside a chain to jump directly to it. The letter background shows the native RamplotR density region, the small corner dot shows Rama8000 standard validation, and the separate underline/number shows pLDDT when available."),
                      tags$div(class = "ram-sequence-legend",
                        tags$span(class="ram-swatch ram-sw-favoured", "Favoured"),
                        tags$span(class="ram-swatch ram-sw-allowed", "Allowed"),
                        tags$span(class="ram-swatch ram-sw-generously-allowed", "Generously allowed"),
                        tags$span(class="ram-swatch ram-sw-outlier", "Not allowed"),
                        tags$span(class="ram-swatch ram-sw-missing", "Missing angles")
                      ),
                      tags$div(class="ram-sequence-standard-key",
                        tags$strong("Standard validation · Rama8000"),
                        tags$span(class="ram-standard-key-favored", "Favored"),
                        tags$span(class="ram-standard-key-allowed", "Allowed"),
                        tags$span(class="ram-standard-key-outlier", "Outlier"),
                        tags$span("Corner dot on each residue.")
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
              uiOutput("experimentalCounterpartPanel"),
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
                    c("All residues" = "All",
                      "Rama8000 outliers" = "Rama8000 outlier",
                      "RamplotR: Not allowed" = "Not allowed",
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
                tags$details(
                  id="ram-compare-source-panel",
                  class="ram-details ram-compare-source-panel",
                  open=NA,
                  tags$summary(
                    tags$span(class="ram-compare-source-title",
                      "Comparison structure"),
                    uiOutput("compareSourceSummary",inline=TRUE)
                  ),
                  tags$div(class = "ram-compare-source",
                    radioButtons("compareInputSource", "Comparison input",
                      choices = c("PDB accession" = "pdb",
                                  "Uploaded file" = "upload"),
                      selected = "pdb", inline = TRUE),
                    tags$div(id = "ram-compare-pdb", textInput(
                      "comparePDB", "Second structure", value = "1CRN",
                      placeholder = "e.g. 1CRN")),
                    tags$div(id = "ram-compare-upload", class = "is-hidden",
                      fileInput("compareFile", "Second PDB/mmCIF file",
                        accept = c(".pdb", ".ent", ".cif", ".mmcif", ".mcif")),
                      selectInput("comparePredictionSource",
                        "Comparison structure type",
                        choices=c(
                          "Experimental / unknown"="experimental",
                          "AlphaFold 2 / ColabFold"="alphafold2",
                          "AlphaFold 3"="alphafold3",
                          "ESMFold"="esmfold",
                          "Other prediction with pLDDT in B-factor"="other_prediction"
                        ),selected="experimental",selectize=FALSE),
                      conditionalPanel(
                        condition="input.comparePredictionSource == 'alphafold3'",
                        fileInput("comparePredictionJson",
                          "Matching AF3 full confidences JSON",
                          accept=c(".json")),
                        fileInput("comparePredictionSummaryJson",
                          "AF3 summary confidences JSON (optional)",
                          accept=c(".json")),
                        tags$p(class="ram-field-hint",
                          "Use confidence files from the same AF3 seed/sample as the uploaded coordinates.")
                      ),
                      tags$p(class="ram-field-hint",
                        "Prediction confidence is interpreted only after an explicit prediction source is selected. Experimental B-factors are never treated as pLDDT.")
                    ),
                    actionButton("compareSubmit", "Load comparison",
                      class="btn-primary")
                  )
                ),
                uiOutput("compareChainControls"),
                uiOutput("atlasCanonicalNotice"),
                tags$div(class = "ram-compare-status", uiOutput("compareSummary")),
                uiOutput("compareChangeTrack"),
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
                      "Changed RamplotR region" = "changed",
                      "Changed Rama8000 category" = "standard_changed",
                      "Changed backbone state" = "basin_changed",
                      "Rama8000 outlier in either structure" = "standard_outlier",
                      "Combined backbone shift ≥ 30°" = "shift_large",
                      "Either angle difference ≥ 30°" = "large",
                      "|ΔpLDDT| ≥ 20" = "confidence_large",
                      "Insertions / deletions" = "gaps"), selected = "All"),
                  downloadButton("downloadComparison", "Export comparison CSV")
                ),
                tags$p(class = "ram-table-hint",
                  "Select a row to highlight its corresponding residues in both 3D structures and the angle plot."),
                tags$div(class = "ram-residue-table", DT::DTOutput("comparison")),
                tags$details(id="ram-group-comparison-panel",
                  class="ram-details ram-group-comparison-panel",
                  tags$summary("Compare groups of structures"),
                  tags$div(class="ram-group-intro",
                    tags$p("Compare experimental states, mutants or prediction sets using uploaded structures or verified Atlas entries. You choose and label both groups; Atlas clusters are exploratory, not state assignments."),
                    actionLink("groupGuide","How do group comparisons work?")
                  ),
                  radioButtons("groupInputMode","Structure source",
                    choices=c("Upload structures"="uploads",
                      "Use verified Atlas selection"="atlas"),
                    selected="uploads",inline=TRUE),
                  uiOutput("groupAtlasSelection"),
                  tags$div(class="ram-group-labels",
                    textInput("groupALabel","Group A label",value="Group A"),
                    textInput("groupBLabel","Group B label",value="Group B")),
                  conditionalPanel(condition="input.groupInputMode !== 'atlas'",
                  tags$div(class="ram-group-step-heading",
                    tags$span("1"), tags$strong("Choose a reference and chain-match criteria")),
                  tags$div(class="ram-group-compare-controls",
                    tags$div(class="ram-group-field",
                      selectInput("groupReferenceChain","Reference chain",
                        choices=character(),selectize=FALSE)),
                    tags$div(class="ram-group-field",
                      numericInput("groupMinIdentity","Minimum chain identity (%)",
                        value=70,min=20,max=100,step=5)),
                    tags$div(class="ram-group-field",
                      numericInput("groupMinCoverage","Minimum reference coverage (%)",
                        value=70,min=20,max=100,step=5))
                  ),
                  tags$p(class="ram-field-hint",
                    "Identity and coverage refer to matching uploaded chains to the selected reference, not to a threshold for biological state changes."),
                  tags$div(class="ram-group-step-heading",
                    tags$span("2"), tags$strong("Upload structures for both conditions")),
                  tags$div(class="ram-group-upload-grid",
                    tags$section(class="ram-group-upload-card",
                      checkboxInput("groupIncludeLoadedA",
                        "Include loaded structure in Group A",value=TRUE),
                      fileInput("groupAFiles","Additional Group A structures",
                        multiple=TRUE,
                        accept=c(".pdb",".ent",".cif",".mmcif",".mcif")),
                      tags$p(class="ram-field-hint",
                        "The structure loaded at the top of RamplotR can be your first Group A member.")
                    ),
                    tags$section(class="ram-group-upload-card",
                      fileInput("groupBFiles","Group B structures",
                        multiple=TRUE,
                        accept=c(".pdb",".ent",".cif",".mmcif",".mcif")),
                      tags$p(class="ram-field-hint",
                        "Upload at least one structure for your comparison condition.")
                    )
                  ),
                  tags$p(class="ram-field-hint",
                    "For meaningful within-group variation, add several independent structures per condition where possible. Members use their first model.")),
                  tags$div(class="ram-group-step-heading",
                    tags$span("3"), tags$strong("Analyse and inspect residue-level differences")),
                  tags$div(class="ram-group-run-row",
                    actionButton("runGroupComparison","Analyse groups",
                      class="btn-primary btn-sm"),
                    tags$span("Results and CSV exports appear after a successful analysis.")
                  ),
                  uiOutput("groupComparisonSummary"),
                  uiOutput("groupComparisonTrack"),
                  uiOutput("groupComparisonSelectionInfo"),
                  uiOutput("groupFingerprintPanel"),
                  tags$div(class="ram-residue-table",
                    DT::DTOutput("groupComparisonRows")),
                  uiOutput("groupComparisonExports")
                )
              )
            ),
            tabPanel(
              title = "Atlas", value = "atlas",
              tags$div(class="ram-subtab-content",
                tags$div(class="ram-result-head",
                  tags$div(tags$h2("Conformational Atlas"),
                    tags$p("Find experimental PDB counterparts, verify exact UniProt residue mapping, and explore geometric groups and candidate backbone changes."))),
                tags$section(class="ram-panel",
                  actionLink("atlasGuide","Read the experimental Atlas walkthrough"),
                  tags$p(class="ram-field-hint",
                    "Searches public RCSB PDB by UniProt cross-reference, restricted to experimental entries. Results are candidate polymer entities, not distinct conformational states."),
                  tags$details(class="ram-details ram-atlas-network-check",
                    tags$summary("Check browser access to RCSB and PDBe"),
                    tags$p(class="ram-field-hint",
                      "Run an on-demand connectivity check from this browser. It tests a known experimental protein through RCSB Search, RCSB metadata and PDBe updated mmCIF with exact SIFTS mapping. No loaded structures or analyses are modified."),
                    actionButton("atlasConnectivityRun","Test archive connections",
                      class="btn-default btn-sm"),
                    uiOutput("atlasConnectivityStatus"),
                    tags$p(class="ram-field-hint",
                      "The check runs from your browser, not a GitHub Actions server. A successful check does not guarantee every other archive entry will download."),
                    tags$a("Troubleshooting and source endpoints",
                      href="https://github.com/BiKC/RamplotR/blob/main/docs/atlas-browser-connectivity.md",
                      target="_blank",rel="noopener noreferrer")),
                  tags$div(class="ram-atlas-controls",
                    textInput("atlasAccession","UniProt accession",
                      placeholder="e.g. P00533"),
                    actionButton("atlasDiscover","Find experimental structures",
                      class="btn-primary")),
                  uiOutput("atlasStatus"),
                  uiOutput("atlasResults"),
                  uiOutput("atlasGeometryPanel"),
                  uiOutput("atlasSwitchPanel"),
                  tags$p(class="ram-field-hint",
                    "Results load 50 experimental entities at a time. Canonical residue coverage, construct equivalence, and conformational-state identity must be verified before interpreting structural states.")
                )
              )
            ),
            tabPanel(
              title = "Guide", value = "guide",
              ram_guide_ui()
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
  tags$script(src = "canonical-mapping.js"),
  tags$script(src = "atlas-sifts-exact.js"),
  tags$script(src = "atlas-ligand-contacts.js"),
  tags$script(src = "atlas-discovery.js"),
  tags$script(src = "atlas-connectivity.js"),
  tags$script(src = "experimental-search.js"),
  tags$script(src = "compare.js"),
  tags$script(src = "prediction.js"),
  tags$script(src = "density.js")
)
# Structure parsing is deliberately triggered by the Analyse button. Every
# downstream result is a reactive expression, so adjusting settings never
# refetches the structure or recomputes backbone torsions.
server <- function(input, output, session) {
  # Simple task navigation. These links never modify the loaded analysis.
  guide_targets <- c(
    openGuide="guide",groupGuide="guide",atlasGuide="guide",
    guideGoPlot="plot",guideGoResidues="residues",
    guideGoSummary="summary",guideGoSummary2="summary",
    guideGoCompare="compare",guideGoGroups="compare",guideGoAtlas="atlas"
  )
  for (id in names(guide_targets)) {
    local({
      trigger <- id
      target <- unname(guide_targets[[id]])
      observeEvent(input[[trigger]], {
        updateTabsetPanel(session, "analysisTabs", selected=target)
      }, ignoreInit=TRUE)
    })
  }

  loaded <- reactiveVal(NULL)
  comparison_loaded <- reactiveVal(NULL)
  external_validation <- reactiveVal(NULL)
  ensemble_results <- reactiveVal(NULL)
  prediction_ensemble_results <- reactiveVal(NULL)
  group_comparison_results <- reactiveVal(NULL)
  experimental_search_results <- reactiveVal(NULL)
  experimental_search_status <- reactiveVal(NULL)
  experimental_search_request <- reactiveVal(0L)
  atlas_request <- reactiveVal(0L)
  atlas_payload <- reactiveVal(NULL)
  atlas_status <- reactiveVal(NULL)
  atlas_connectivity_request <- reactiveVal(0L)
  atlas_connectivity_status <- reactiveVal(NULL)
  atlas_exact_request <- reactiveVal(0L)
  atlas_exact_selected <- reactiveVal(NULL)
  atlas_exact_results <- reactiveVal(list())
  atlas_candidate_picks <- reactiveVal(character())
  atlas_verify_queue <- reactiveVal(character())
  atlas_verify_progress <- reactiveVal(NULL)
  atlas_geometry_result <- reactiveVal(NULL)
  atlas_switch_result <- reactiveVal(NULL)
  atlas_selected_position <- reactiveVal(NULL)
  atlas_pair_handoff <- reactiveVal(NULL)
  atlas_alignment_context <- reactiveVal(NULL)
  atlas_group_transfer <- reactiveVal(NULL)
  canonical_segments <- reactiveVal(ram_canonical_empty_segments())
  canonical_mapping <- reactiveVal(ram_canonical_empty_map())
  canonical_status <- reactiveVal(NULL)
  canonical_request <- reactiveVal(0L)
  selected_residue <- reactiveVal(NULL)
  selected_comparison <- reactiveVal(NULL)
  compare_swapped <- reactiveVal(FALSE)
  compare_swap_chains <- reactiveVal(NULL)
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
  # Shared loader: user submissions and explicit Atlas representative-pair
  # inspection both use the same primary PDB parsing and NGL preparation.
  load_primary_structure <- function(source_type,source_label,
    declared_source="experimental",sidecar="",summary_file="") {
    if (length(source_type)!=1L || is.na(source_type) ||
        !source_type %in% c("pdb","upload","afdb")) return(invisible(FALSE))
    is_upload <- identical(source_type,"upload")
    is_afdb <- identical(source_type,"afdb")
    if (is_upload && (is.null(input$structfile) ||
                      is.null(input$structfile$datapath))) {
      showNotification("Choose a PDB or mmCIF file first.",type="error")
      return(invisible(FALSE))
    }
    key <- paste(source_type,source_label,declared_source,sidecar,
                 summary_file,sep=":")
    previous <- isolate(loaded())
    if (!is.null(previous) && identical(previous$key,key))
      return(invisible(TRUE))
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
                  declared_source = declared_source, input_source = source_type,
                  pdb_accession = if (identical(source_type,"pdb"))
                    toupper(source_label) else NULL,
                  uniprot_accession = if (is_afdb)
                    ram_uniprot_accession(source_label) else NULL))
      incProgress(0.25, detail = "Preparing interactive views")
    })
    invisible(TRUE)
  }
  observeEvent(input$submit, {
    atlas_alignment_context(NULL)
    source_type <- input$inputSource
    is_upload <- identical(source_type,"upload")
    is_afdb <- identical(source_type,"afdb")
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
    load_primary_structure(source_type,source_label,declared_source,
      sidecar,summary_file)
  }, ignoreInit=TRUE)

  # Canonical UniProt coordinates are an additive annotation. PDB author
  # numbering remains the local coordinate system used by selection and NGL.
  # Public PDB entries request SIFTS mappings client-side so Shinylive keeps
  # working without a server HTTP dependency. AlphaFold DB models already use
  # the requested UniProt sequence numbering.
  observeEvent(loaded(), {
    structure <- loaded()
    canonical_segments(ram_canonical_empty_segments())
    canonical_mapping(ram_canonical_empty_map())
    request_id <- isolate(canonical_request()) + 1L
    canonical_request(request_id)

    if (is.null(structure)) {
      canonical_status(NULL)
      return()
    }
    if (identical(structure$input_source,"afdb") &&
        !is.null(structure$uniprot_accession)) {
      mapping <- tryCatch(
        ram_afdb_canonical_map(structure$torsions,
                               structure$uniprot_accession),
        error=function(e) {
          canonical_status(list(state="error",
            message=conditionMessage(e)))
          NULL
        }
      )
      if (!is.null(mapping)) {
        canonical_mapping(mapping)
        canonical_status(list(
          state="mapped",source="AlphaFold DB",
          expanded_mapping_rows=nrow(mapping),
          accessions=unique(mapping$uniprot_accession)
        ))
      }
      return()
    }
    if (identical(structure$input_source,"pdb") &&
        !is.null(structure$pdb_accession)) {
      canonical_status(list(state="searching",source="PDBe SIFTS",
        pdb_id=structure$pdb_accession))
      session$sendCustomMessage("ram-canonical-map",list(
        request_id=as.character(request_id),
        pdb_id=structure$pdb_accession
      ))
      return()
    }
    canonical_status(list(state="unavailable",
      message="Canonical mapping is not inferred automatically for uploaded structures."))
  },ignoreInit=TRUE)

  observeEvent(input$ramCanonicalMapping, {
    value <- input$ramCanonicalMapping
    if (!is.list(value) || is.null(value$request_id)) return()
    if (!identical(as.character(value$request_id),
                   as.character(isolate(canonical_request())))) return()
    structure <- isolate(loaded())
    if (is.null(structure) || !identical(structure$input_source,"pdb"))
      return()
    if (!is.null(value$pdb_id) &&
        !identical(toupper(as.character(value$pdb_id)),
                   toupper(as.character(structure$pdb_accession)))) return()

    if (!identical(as.character(value$state),"ok")) {
      canonical_status(list(state="error",source="PDBe SIFTS",
        message=if(is.null(value$message)) "Canonical mapping unavailable."
          else as.character(value$message)))
      return()
    }

    segments <- tryCatch(
      ram_sifts_normalize_segments(value$segments,
                                   pdb_id=structure$pdb_accession),
      error=function(e) {
        canonical_status(list(state="error",source="PDBe SIFTS",
          message=conditionMessage(e)))
        NULL
      }
    )
    if (is.null(segments)) return()
    mapping <- ram_sifts_expand_safe(segments)
    canonical_segments(segments)
    canonical_mapping(mapping)
    safe_n <- sum(segments$safe_linear,na.rm=TRUE)
    canonical_status(list(
      state=if(nrow(segments) && safe_n==nrow(segments) && nrow(mapping))
        "mapped" else "partial",
      source="PDBe SIFTS",
      endpoint=if(is.null(value$endpoint)) "" else as.character(value$endpoint),
      segments=nrow(segments),
      safe_segments=safe_n,
      expanded_mapping_rows=nrow(mapping),
      accessions=unique(segments$uniprot_accession)
    ))
  },ignoreInit=TRUE)

  # Comparison loading is a separate, deliberate action, so changing plot
  # settings does not repeatedly refetch the secondary structure. The helper is
  # shared by the manual Compare tab and prediction-to-experiment discovery.
  load_comparison_structure <- function(path = NULL, original_name = NULL,
                                        pdb_id = NULL,
                                        preferred_chain_a = NULL,
                                        preferred_chain_b = NULL,
                                        declared_source = "experimental",
                                        sidecar = NULL,
                                        summary_file = NULL) {
    is_upload <- !is.null(path)
    source_id <- if (is_upload) path else toupper(trimws(pdb_id))
    source_name <- if (is_upload) original_name else source_id
    data <- tryCatch(
      ram_load_structure(
        path = path,
        original_name = original_name,
        pdb_id = if (is_upload) NULL else source_id
      ), error = function(e) {
        showNotification(conditionMessage(e), type="error", duration=12)
        NULL
      }
    )
    if (is.null(data)) return(FALSE)
    torsions <- tryCatch(ram_extract_torsions(ram_model_at(data, 1L)),
      error = function(e) {
        showNotification(conditionMessage(e), type="error",duration=12)
        NULL
      })
    if (is.null(torsions)) return(FALSE)
    if (!is_upload) declared_source <- "experimental"
    permitted_sources <- c("experimental","alphafold2","alphafold3",
                           "esmfold","other_prediction")
    if (length(declared_source)!=1L || !declared_source %in% permitted_sources)
      declared_source <- "experimental"
    if (identical(declared_source,"alphafold3") &&
        (is.null(sidecar) || !nzchar(sidecar) || !file.exists(sidecar))) {
      showNotification(
        "AlphaFold 3 comparison requires the matching full confidences JSON.",
        type="error",duration=14)
      return(FALSE)
    }
    comparison_name <- tools::file_path_sans_ext(basename(source_name))
    prediction <- NULL
    if (!identical(declared_source,"experimental")) {
      prediction <- tryCatch(
        ram_prepare_prediction(
          ram_model_at(data,1L),torsions,declared_source,
          sidecar=if(!is.null(sidecar) && nzchar(sidecar)) sidecar else NULL,
          summary_file=if(!is.null(summary_file) && nzchar(summary_file))
            summary_file else NULL,
          model_id=comparison_name
        ),
        error=function(e) {
          showNotification(paste("Comparison confidence:",
            conditionMessage(e)),type="error",duration=15)
          NULL
        }
      )
      if (is.null(prediction)) return(FALSE)
      if (length(prediction$notes))
        showNotification(paste(prediction$notes,collapse=" "),
          type="warning",duration=14)
    }
    comparison_loaded(list(
      pdb=data, torsions=torsions,
      name=comparison_name,
      source_id=source_id,
      viewer_format=if (is_upload) ram_detect_format(source_name) else NULL,
      nmodels=ram_model_count(data),
      prediction=prediction,
      declared_source=declared_source,
      preferred_chain_a=preferred_chain_a,
      preferred_chain_b=preferred_chain_b
    ))
    compare_swapped(FALSE)
    compare_swap_chains(NULL)
    session$sendCustomMessage("ram-compare-source-state",list(open=FALSE))
    TRUE
  }

  observeEvent(input$compareSubmit, {
    atlas_alignment_context(NULL)
    is_upload <- identical(input$compareInputSource, "upload")
    if (is_upload && (is.null(input$compareFile) ||
                      is.null(input$compareFile$datapath))) {
      showNotification("Choose a second PDB or mmCIF file.", type="error")
      return()
    }
    if (is_upload) {
      declared_source <- input$comparePredictionSource
      if (is.null(declared_source) || !nzchar(declared_source))
        declared_source <- "experimental"
      sidecar <- if (identical(declared_source,"alphafold3") &&
                        !is.null(input$comparePredictionJson))
        input$comparePredictionJson$datapath else NULL
      summary_file <- if (identical(declared_source,"alphafold3") &&
                             !is.null(input$comparePredictionSummaryJson))
        input$comparePredictionSummaryJson$datapath else NULL
      load_comparison_structure(
        path=input$compareFile$datapath,
        original_name=input$compareFile$name,
        declared_source=declared_source,
        sidecar=sidecar,
        summary_file=summary_file
      )
    } else {
      pdb_id <- toupper(trimws(input$comparePDB))
      if (!nzchar(pdb_id)) {
        showNotification("Enter a PDB accession.", type="error")
        return()
      }
      load_comparison_structure(pdb_id=pdb_id)
    }
  }, ignoreInit=TRUE)

  observeEvent(input$compareSwap, {
    req(comparison_loaded(), loaded())
    old_a <- isolate(input$compareChainA)
    old_b <- isolate(input$compareChainB)
    compare_swap_chains(list(a=old_b,b=old_a))
    compare_swapped(!isTRUE(isolate(compare_swapped())))
    selected_comparison(NULL)
  }, ignoreInit=TRUE)

  output$compareSourceSummary <- renderUI({
    comparison <- comparison_loaded()
    if (is.null(comparison))
      return(tags$span(class="ram-compare-source-summary",
        "Load a second structure to begin"))

    source_label <- switch(
      comparison$declared_source,
      experimental="Experimental / unknown",
      alphafold2="AlphaFold 2 / ColabFold",
      alphafold3="AlphaFold 3",
      esmfold="ESMFold",
      other_prediction="Predicted model",
      "Experimental / unknown"
    )
    origin <- if (is.null(comparison$viewer_format))
      "PDB accession" else "Uploaded file"
    tags$span(class="ram-compare-source-summary is-loaded",
      tags$strong(comparison$name),
      tags$span(paste(origin,source_label,sep=" · "))
    )
  })

  output$compareChainControls <- renderUI({
    main <- req(loaded())
    comparison <- req(comparison_loaded())
    swapped <- isTRUE(compare_swapped())
    first <- if (swapped) comparison else main
    second <- if (swapped) main else comparison
    first_chains <- if (!is.null(first$chains)) first$chains else
      unique(first$torsions$chain)
    second_chains <- if (!is.null(second$chains)) second$chains else
      unique(second$torsions$chain)
    preferred_first <- if (swapped) comparison$preferred_chain_b
      else comparison$preferred_chain_a
    preferred_second <- if (swapped) comparison$preferred_chain_a
      else comparison$preferred_chain_b
    pending <- compare_swap_chains()
    selected_a <- if (!is.null(pending) && !is.null(pending$a) &&
                      pending$a %in% first_chains)
      pending$a else if (!is.null(preferred_first) &&
                         preferred_first %in% first_chains)
      preferred_first else first_chains[[1L]]
    selected_b <- if (!is.null(pending) && !is.null(pending$b) &&
                      pending$b %in% second_chains)
      pending$b else if (!is.null(preferred_second) &&
                         preferred_second %in% second_chains)
      preferred_second else second_chains[[1L]]

    tags$div(class="ram-compare-chains",
      tags$div(class="ram-compare-role",
        tags$span(class="ram-compare-role-label","Primary · coral"),
        selectInput("compareChainA", paste("Chain in", first$name),
          choices=first_chains, selected=selected_a)),
      tags$div(class="ram-compare-role",
        tags$span(class="ram-compare-role-label","Comparison · blue"),
        selectInput("compareChainB", paste("Chain in", second$name),
          choices=second_chains, selected=selected_b)),
      if (comparison$nmodels > 1L)
        selectInput("compareModel", paste("Model in", comparison$name),
          choices=as.character(seq_len(comparison$nmodels)), selected="1"),
      tags$div(class="ram-compare-chain-actions",
        actionButton("compareSwap", "Swap primary ↔ comparison",
          class="btn-default btn-sm",
          title="Swap A/B roles without reloading either structure"),
        tags$span(class="ram-compare-swap-state",
          if (swapped) "Roles swapped" else "Loaded structure is primary")
      )
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
    # Compute the six-class Rama8000 result in parallel. It is independent
    # of the selected RamplotR plotting background and is exposed explicitly
    # rather than relabelling the native four-region result.
    result <- ram_rama8000_classify(result, file.path("static", "rama8000"))
    if (current_model() == 1L)
      result <- ram_apply_prediction(result, structure$prediction)
    result <- ram_join_geometry(result,model_geometry())
    official <- external_validation()
    if(!is.null(official) && identical(official$key,structure$key))
      result <- ram_external_validation_join(result,official$records,
                                               model=current_model())
    result <- ram_canonical_join(result,canonical_mapping())
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
                 "phi", "psi", "region", "density",
                 "rama8000_region", "rama8000_group", "rama8000_score")
    has_canonical <- all(c("uniprot_accession","uniprot_resi") %in% names(data)) &&
      any(!is.na(data$uniprot_accession) & is.finite(data$uniprot_resi))
    if (has_canonical)
      columns <- c(columns,"uniprot_accession","uniprot_resi")
    if ("plddt" %in% names(data))
      columns <- c(columns, "plddt", "confidence_category")
    shown <- data[, columns, drop = FALSE]
    if ("plddt" %in% names(shown)) shown$plddt <- round(shown$plddt, 1L)
    shown$phi <- round(shown$phi, 1L)
    shown$psi <- round(shown$psi, 1L)
    shown$density <- round(shown$density, 1L)
    shown$rama8000_score <- round(100 * shown$rama8000_score, 2L)
    widget <- DT::datatable(
      shown, rownames = FALSE,
      colnames = c("Chain", "Residue", "Ins.", "AA", "Phi (°)", "Psi (°)",
                   "RamplotR region", "Percentile", "Rama8000", "Rama8000 class",
                   "Rama8000 score (%)",
                   if (has_canonical) c("UniProt","UniProt residue"),
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
    widget <- DT::formatStyle(widget, "region",
      backgroundColor = DT::styleEqual(
        c("Favoured", "Allowed", "Generously allowed", "Not allowed"),
        c("#D4ECE7", "#E8F3F1", "#FFF5E1", "#FCE5DD")
      ),
      fontWeight = "600"
    )
    DT::formatStyle(widget, "rama8000_region",
      backgroundColor = DT::styleEqual(
        c("Favored", "Allowed", "Outlier"),
        c("#D4ECE7", "#FFF5E1", "#FCE5DD")
      ),
      fontWeight = "650"
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
      mapping <- canonical_status()
      if(!is.null(mapping)) {
        provenance$canonical_mapping_state <- mapping$state
        provenance$canonical_mapping_source <- if(is.null(mapping$source))
          "none" else mapping$source
        if(!is.null(mapping$accessions) && length(mapping$accessions))
          provenance$canonical_uniprot_accessions <-
            paste(mapping$accessions,collapse=",")
        provenance$canonical_mapped_residues <-
          if("canonical_status" %in% names(data))
            sum(data$canonical_status=="mapped",na.rm=TRUE) else 0L
        if(!is.null(mapping$segments))
          provenance$canonical_mapping_segments <- mapping$segments
        if(!is.null(mapping$safe_segments))
          provenance$canonical_safe_linear_segments <- mapping$safe_segments
        if(!is.null(mapping$endpoint) && nzchar(mapping$endpoint))
          provenance$canonical_mapping_endpoint <- mapping$endpoint
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
            tags$span(class="ram-seq-current", role="status",
              "Select a residue"),
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
                standard <- residue$rama8000_region[[1L]]
                standard_class <- if (is.na(standard)) "ram-seq-standard-missing"
                  else paste0("ram-seq-standard-",tolower(standard))
                position <- paste0(residue$resi[[1L]],
                                   residue$insertion_code[[1L]])
                tags$div(class="ram-seq-slot",
                  tags$span(class="ram-seq-position",
                    if (nzchar(labels[[i]])) labels[[i]] else "\u00a0",
                    "aria-hidden"="true"),
                  tags$button(type="button",
                    class=paste("ram-seq-res",
                      paste0("ram-seq-",statuses[[i]]),
                      standard_class,
                      if (show_confidence) "ram-seq-with-confidence" else ""),
                    style=if (is.finite(score))
                      paste0("--ram-plddt-color:",ram_plddt_color(score)) else NULL,
                    disabled=if (!selectable[[i]]) "disabled" else NULL,
                    "data-chain"=residue$chain[[1L]],
                    "data-resi"=residue$resi[[1L]],
                    "data-insertion"=residue$insertion_code[[1L]],
                    "data-plddt"=if (is.finite(score))
                      sprintf("%.1f",score) else "",
                    "data-rama8000"=if (!is.na(standard)) standard else "",
                    title=paste0(residue$resn[[1L]], " ",
                      residue$chain[[1L]], position, " · ",
                      if (is.na(residue$region[[1L]])) "Missing angles"
                      else residue$region[[1L]],
                      if (!is.na(standard)) paste0(" · Rama8000 ",standard,
                        " (",residue$rama8000_group[[1L]],")") else "",
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

  # On-demand checks are executed from the visitor's browser origin;
  # this never mutates any scientific result or loaded structure.
  observeEvent(input$atlasConnectivityRun, {
    id <- isolate(atlas_connectivity_request())+1L
    atlas_connectivity_request(id)
    atlas_connectivity_status(list(state="running",passed=0L,total=3L,
      checks=list()))
    session$sendCustomMessage("ram-atlas-connectivity",
      list(request_id=as.character(id)))
  },ignoreInit=TRUE)
  observeEvent(input$ramAtlasConnectivity, {
    value <- input$ramAtlasConnectivity
    if(!is.list(value) ||
       !identical(as.character(value$request_id),
                  as.character(isolate(atlas_connectivity_request()))))
      return()
    atlas_connectivity_status(value)
  },ignoreInit=TRUE)
  output$atlasConnectivityStatus <- renderUI({
    status <- atlas_connectivity_status()
    if(is.null(status)) return(NULL)
    checks <- status$checks
    if(!is.list(checks)) checks <- list()
    state <- as.character(status$state)
    if(identical(state,"error"))
      return(tags$p(class="ram-confidence-warning",
        if(!is.null(status$message)) as.character(status$message)
        else "The browser diagnostic could not run."))
    if(!state %in% c("ok","partial_failure","running"))
      return(NULL)
    passed <- suppressWarnings(as.integer(status$passed))
    total <- suppressWarnings(as.integer(status$total))
    if(length(passed)!=1L || is.na(passed)) passed <- 0L
    if(length(total)!=1L || is.na(total)) total <- 3L
    tags$div(class="ram-atlas-connectivity-results",
      tags$p(class=if(identical(state,"partial_failure"))
          "ram-confidence-warning" else "ram-field-hint",
        if(identical(state,"running"))
          sprintf("Checking archive connections: %d/%d passed so far...",
            passed,total)
        else if(identical(state,"ok"))
          sprintf("All %d endpoints passed in this browser.",total)
        else sprintf("%d/%d endpoints passed. Check failures below.",
          passed,total)),
      if(!is.null(status$tested_origin))
        tags$p(class="ram-field-hint",
          paste("Browser origin:",as.character(status$tested_origin))),
      tags$ul(class="ram-atlas-connectivity-list",
        lapply(checks,function(check) {
          if(!is.list(check)) return(NULL)
          good <- identical(check$state,"ok")
          name <- as.character(check$name)
          value <- if(good && !is.null(check$summary))
            as.character(check$summary)
          else if(!is.null(check$detail)) as.character(check$detail)
          else "Result unavailable."
          reason <- if(!good && !is.null(check$reason))
            paste0(" [",as.character(check$reason),"]") else ""
          duration <- if(!is.null(check$elapsed_ms))
            sprintf(" (%d ms)",as.integer(check$elapsed_ms)) else ""
          tags$li(class=if(good) "ram-atlas-probe-ok"
              else "ram-atlas-probe-failed",
            tags$strong(paste0(if(good) "Passed: " else "Failed: ",name)),
            tags$span(paste0(" — ",value,reason,duration)))
        }))
    )
  })

  # Prefill from AlphaFold DB without issuing a network request.
  observeEvent(loaded(), {
    structure <- loaded()
    if(!is.null(structure) && !is.null(structure$uniprot_accession))
      updateTextInput(session,"atlasAccession",
        value=structure$uniprot_accession)
  },ignoreInit=TRUE)

  # Avoid showing stale results after the requested accession changes.
  observeEvent(input$atlasAccession, {
    atlas_payload(NULL)
    atlas_status(NULL)
    atlas_exact_selected(NULL)
    atlas_exact_results(list())
    atlas_candidate_picks(character())
    atlas_verify_queue(character())
    atlas_verify_progress(NULL)
    atlas_exact_request(isolate(atlas_exact_request())+1L)
  },ignoreInit=TRUE)

  # Atlas inventory is independent of the loaded structure. It is an
  # experimental polymer-entity discovery view, not yet a state classification.
  observeEvent(input$atlasDiscover, {
    accession <- tryCatch(ram_uniprot_accession(input$atlasAccession),
      error=function(e) {
        atlas_status(list(state="error",message=conditionMessage(e)))
        NULL
      })
    if(is.null(accession)) return()
    id <- isolate(atlas_request())+1L
    atlas_request(id)
    atlas_payload(NULL)
    atlas_exact_selected(NULL)
    atlas_exact_results(list())
    atlas_candidate_picks(character())
    atlas_verify_queue(character())
    atlas_verify_progress(NULL)
    atlas_exact_request(isolate(atlas_exact_request())+1L)
    atlas_status(list(state="searching",
      message=paste("Searching experimental PDB entities for",accession,"...")))
    session$sendCustomMessage("ram-atlas-discover",list(
      request_id=as.character(id),accession=accession,
      rows=50L,start=0L))
  },ignoreInit=TRUE)

  observeEvent(input$atlasLoadMore, {
    previous <- isolate(atlas_payload())
    status <- isolate(atlas_status())
    if(is.null(previous) || !isTRUE(previous$has_more) ||
       identical(status$state,"searching")) return()
    offset <- as.integer(previous$next_offset)
    id <- isolate(atlas_request())+1L
    atlas_request(id)
    atlas_status(list(state="searching",
      message=sprintf("Loading experimental entities %d–%d...",
        offset+1L,min(previous$total_count,offset+50L))))
    session$sendCustomMessage("ram-atlas-discover",list(
      request_id=as.character(id),accession=previous$accession,
      rows=50L,start=offset))
  },ignoreInit=TRUE)

  observeEvent(input$ramAtlasResults, {
    data <- input$ramAtlasResults
    if(!is.list(data) || is.null(data$request_id) ||
       !identical(as.character(data$request_id),
                  as.character(isolate(atlas_request())))) return()
    if(is.null(data$accession) ||
       !identical(toupper(trimws(as.character(data$accession))),
                  toupper(trimws(as.character(isolate(input$atlasAccession))))))
      return()
    if(identical(as.character(data$state),"error")) {
      # Keep successfully retrieved pages when a later page fails.
      atlas_status(list(state="error",
        message=if(is.null(data$message)) "Atlas search failed; retry the page."
          else as.character(data$message)))
    } else {
      merged <- tryCatch(
        ram_atlas_merge_page(isolate(atlas_payload()),data),
        error=function(e) {
          atlas_status(list(state="error",message=conditionMessage(e)))
          NULL
        }
      )
      if(is.null(merged)) return()
      atlas_payload(merged)
      atlas_status(list(state=if(isTRUE(merged$stalled)) "error" else "done",
        message=if(isTRUE(merged$stalled))
          "RCSB returned an empty page before its reported total; restart the search."
        else sprintf("Loaded %d of %d returned search hits across %d page%s.",
          merged$returned_count,merged$total_count,merged$pages,
          if(merged$pages==1L) "" else "s")))
    }
  },ignoreInit=TRUE)

  # Select experimental entities directly from cards. The choice is
  # deliberately separate from verification: no archive file is downloaded
  # merely because a search card is selected.
  observeEvent(input$ramAtlasCandidatePick, {
    change <- input$ramAtlasCandidatePick
    cohort <- isolate(atlas_payload())
    if(is.null(cohort) || !is.list(change) ||
       !identical(as.character(change$accession),cohort$accession))
      return()
    key <- as.character(change$key)
    available <- vapply(cohort$results,ram_atlas_record_key,character(1L))
    if(length(key)!=1L || is.na(key) || !key %in% available) return()
    selected <- isolate(atlas_candidate_picks())
    if(isTRUE(change$selected)) {
      if(!key %in% selected && length(selected)>=12L) {
        showNotification("Atlas supports up to 12 selected entities per comparison.",
          type="warning",duration=10)
        session$sendCustomMessage("ram-atlas-set-selection",
          list(ids=selected))
        return()
      }
      selected <- unique(c(selected,key))
    } else selected <- setdiff(selected,key)
    atlas_candidate_picks(selected)
  },ignoreInit=TRUE)
  observeEvent(input$atlasClearSelection, {
    atlas_candidate_picks(character())
    atlas_verify_queue(character())
    atlas_verify_progress(NULL)
    session$sendCustomMessage("ram-atlas-set-selection",
      list(ids=character()))
  },ignoreInit=TRUE)
  observeEvent(atlas_candidate_picks(), {
    # A previous grouping must never be presented as if it reflects the
    # newly chosen experimental subset.
    atlas_geometry_result(NULL)
    atlas_switch_result(NULL)
    atlas_group_transfer(NULL)
  },ignoreInit=TRUE)

  start_atlas_verification <- function(key,cohort) {
    available <- vapply(cohort$results,ram_atlas_record_key,character(1L))
    if(length(key)!=1L || is.na(key) || !key %in% available)
      return(FALSE)
    split <- strsplit(key,"_",fixed=TRUE)[[1L]]
    if(length(split)!=2L || !grepl("^[A-Z0-9]{4}$",split[[1L]]) ||
       !grepl("^[1-9][0-9]*$",split[[2L]])) return(FALSE)
    seq <- isolate(atlas_exact_request())+1L
    atlas_exact_request(seq)
    atlas_exact_selected(list(key=key,request_id=as.character(seq),
                              accession=cohort$accession))
    results <- isolate(atlas_exact_results())
    results[[key]] <- list(state="searching")
    atlas_exact_results(results)
    session$sendCustomMessage("ram-atlas-sifts-exact",list(
      request_id=as.character(seq),pdb_id=split[[1L]],
      entity_id=split[[2L]],accession=cohort$accession))
    TRUE
  }

  # Individual verification remains available. The sequential batch path
  # avoids overlapping browser requests overwriting the request-id guard.
  observeEvent(input$ramAtlasVerifyPick, {
    pick <- input$ramAtlasVerifyPick
    cohort <- isolate(atlas_payload())
    if(is.null(cohort) || !is.list(pick)) return()
    if(length(isolate(atlas_verify_queue())) ||
       !is.null(isolate(atlas_verify_progress()))) {
      showNotification("Selected-structure verification is in progress.",
        type="message",duration=6)
      return()
    }
    pdb <- toupper(trimws(as.character(pick$pdb_id)))
    entity <- as.character(pick$entity_id)
    if(length(pdb)!=1L || length(entity)!=1L ||
       is.na(pdb) || is.na(entity)) return()
    key <- paste0(pdb,"_",entity)
    if(!key %in% vapply(cohort$results,ram_atlas_record_key,character(1L)))
      return()
    atlas_candidate_picks(unique(c(isolate(atlas_candidate_picks()),key)))
    start_atlas_verification(key,cohort)
  },ignoreInit=TRUE)

  observeEvent(input$atlasVerifySelection, {
    cohort <- isolate(atlas_payload())
    if(is.null(cohort)) return()
    selected <- isolate(atlas_candidate_picks())
    if(!length(selected)) {
      showNotification("Select structures from the cards first.",
        type="warning",duration=9)
      return()
    }
    if(!is.null(isolate(atlas_verify_progress()))) return()
    available <- vapply(cohort$results,ram_atlas_record_key,character(1L))
    selected <- intersect(selected,available)
    current <- isolate(atlas_exact_results())
    needed <- selected[!vapply(selected,function(key)
      identical(current[[key]]$state,"mapped"),logical(1L))]
    if(!length(needed)) {
      showNotification("All selected structures are already verified.",
        type="message",duration=8)
      return()
    }
    atlas_verify_progress(list(total=length(needed),completed=0L,
                               accession=cohort$accession))
    atlas_verify_queue(tail(needed,-1L))
    start_atlas_verification(needed[[1L]],cohort)
  },ignoreInit=TRUE)

  observeEvent(input$ramAtlasSiftsExact, {
    value <- input$ramAtlasSiftsExact
    chosen <- isolate(atlas_exact_selected())
    cohort <- isolate(atlas_payload())
    if(!is.list(value) || is.null(chosen) || is.null(cohort) ||
       !identical(as.character(value$request_id),chosen$request_id) ||
       !identical(as.character(value$accession),chosen$accession) ||
       !identical(paste0(as.character(value$pdb_id),"_",
                         as.character(value$entity_id)),chosen$key)) return()
    outcomes <- isolate(atlas_exact_results())
    if(!identical(as.character(value$state),"ok")) {
      outcomes[[chosen$key]] <- list(state="error",
        message=if(is.null(value$message))
          "Exact SIFTS mapping could not be retrieved."
          else as.character(value$message))
    } else {
      mapped <- tryCatch(
        ram_atlas_exact_sifts_map(value, value$pdb_id,
          value$entity_id,chosen$accession),
        error=function(e) {
          outcomes[[chosen$key]] <<- list(state="error",
            message=conditionMessage(e))
          NULL
        })
      if(!is.null(mapped)) {
        summary <- ram_atlas_sifts_summary(mapped)
        outcomes[[chosen$key]] <- c(list(state="mapped",
          source=as.character(value$source),endpoint=as.character(value$endpoint),
          matched_sifts_rows=as.integer(value$matched_sifts_rows),
          unlinked_sifts_rows=as.integer(value$unlinked_sifts_rows)),
          summary,list(mapping=mapped,
            backbone_atoms=if(is.null(value$backbone_atoms)) list()
              else value$backbone_atoms,
            ca_points=if(is.null(value$ca_points)) list() else value$ca_points,
            ca_warning=if(is.null(value$ca_warning)) "" else
              as.character(value$ca_warning),
            ligand_contacts=if(is.null(value$ligand_contacts))
              list(status="unavailable",warning="Contact evidence missing.")
              else value$ligand_contacts))
      }
    }
    atlas_exact_results(outcomes)
    progress <- isolate(atlas_verify_progress())
    if(!is.null(progress)) {
      progress$completed <- progress$completed+1L
      remaining <- isolate(atlas_verify_queue())
      if(length(remaining) && identical(cohort$accession,progress$accession)) {
        atlas_verify_queue(tail(remaining,-1L))
        atlas_verify_progress(progress)
        start_atlas_verification(remaining[[1L]],cohort)
      } else {
        atlas_verify_queue(character())
        atlas_verify_progress(NULL)
        showNotification(sprintf(
          "Verified selection: %d completed. Choose at least two mapped structures to cluster.",
          progress$completed),type="message",duration=10)
      }
    }
  },ignoreInit=TRUE)

  output$atlasStatus <- renderUI({
    item <- atlas_status()
    if(is.null(item)) return(NULL)
    tags$p(class=if(identical(item$state,"error"))
      "ram-confidence-warning" else "ram-field-hint",item$message)
  })

  output$atlasSelectionToolbar <- renderUI({
    payload <- atlas_payload()
    if(is.null(payload)) return(NULL)
    chosen <- atlas_candidate_picks()
    progress <- atlas_verify_progress()
    tags$div(class="ram-atlas-selection-toolbar",
      tags$strong(sprintf("%d selected for analysis (maximum 12)",
        length(chosen))),
      tags$span("Tick entries below, then verify the selected structures before clustering."),
      actionButton("atlasVerifySelection",
        sprintf("Verify %d selected structure%s",
          length(chosen),if(length(chosen)==1L) "" else "s"),
        class="btn-primary btn-sm"),
      actionLink("atlasClearSelection","Clear selection"),
      if(!is.null(progress)) tags$p(class="ram-field-hint",
        sprintf("Verifying selection: %d/%d completed.",
          progress$completed,progress$total)))
  })
  output$atlasResults <- renderUI({
    payload <- atlas_payload()
    if(is.null(payload)) return(NULL)
    results <- payload$results
    exact <- atlas_exact_results()
    chosen <- isolate(atlas_candidate_picks())
    verified <- ram_atlas_cohort_summary(exact,payload$accession)
    get <- function(item,key,default="") {
      value <- item[[key]]
      if(is.null(value) || !length(value) || is.na(value[[1L]])) default
      else as.character(value[[1L]])
    }
    pdbs <- unique(vapply(results,get,character(1L),key="pdb_id"))
    total <- suppressWarnings(as.integer(payload$total_count))
    returned <- suppressWarnings(as.integer(payload$returned_count))
    unresolved <- suppressWarnings(as.integer(payload$incomplete_metadata))
    tags$div(class="ram-atlas-inventory",
      uiOutput("atlasSelectionToolbar"),
      tags$p(class="ram-field-hint",
        sprintf("%d enriched polymer entities in %d PDB entries. RCSB reports %s matching entities%s.",
          length(results),length(pdbs),
          if(is.finite(total)) as.character(total) else "an unknown number of",
          if(is.finite(total) && is.finite(returned) && total>returned)
            sprintf("; %d of %d search hits retrieved so far",returned,total)
            else "; all currently reported hits retrieved"),
        if(is.finite(unresolved) && unresolved>0L)
          sprintf(" Metadata unavailable for %d returned entities.",unresolved) else "",
        if(isTRUE(payload$duplicate_count>0L))
          sprintf(" %d duplicate entity IDs collapsed.",payload$duplicate_count) else ""),
      if(verified$verified_entities>0L)
        tags$div(class="ram-atlas-verified-cohort",
          tags$strong(sprintf(
            "%d verified experimental entities · %d distinct observed UniProt positions",
            verified$verified_entities,verified$observed_positions)),
          tags$p(class="ram-field-hint",
            sprintf(paste0(
              "%d exact SIFTS residue rows; %d ambiguous local residue IDs. ",
              "Counts are verified-entity support, not distinct conformational ",
              "states, independent replicates or sequence completeness."),
              verified$exact_rows,verified$ambiguous_local_residues)),
          if(verified$exact_rows>0L)
            tags$div(class="ram-export-actions",
              downloadButton("downloadAtlasCanonical",
                "Export verified residue mapping CSV"),
              downloadButton("downloadAtlasSupport",
                "Export UniProt position support CSV"))),
      if(length(payload$failed_entity_ids))
        tags$details(class="ram-details",
          tags$summary(sprintf("Show %d entity IDs with unavailable metadata",
            length(payload$failed_entity_ids))),
          tags$p(class="ram-field-hint",
            paste(payload$failed_entity_ids,collapse=", "))),
      if(!length(results))
        tags$p(class="ram-counterpart-empty",
          "No enriched entity metadata were returned on the pages checked so far."),
      tags$div(class="ram-counterpart-results",
        lapply(results,function(item) {
          id <- get(item,"pdb_id")
          entity <- get(item,"entity_id")
          verify <- exact[[paste0(id,"_",entity)]]
          chain <- get(item,"chain")
          resolution <- suppressWarnings(as.numeric(get(item,"resolution",NA_character_)))
          reference_coverage <- suppressWarnings(as.numeric(
            get(item,"reference_sequence_coverage",NA_character_)))
          entity_coverage <- suppressWarnings(as.numeric(
            get(item,"entity_sequence_coverage",NA_character_)))
          tags$article(class=paste("ram-counterpart-card",
              if(paste0(id,"_",entity) %in% chosen) "is-selected" else ""),
            tags$div(class="ram-counterpart-card-main",
              tags$label(class="ram-atlas-pick-control",
                tags$input(type="checkbox",class="ram-atlas-pick",
                  "data-key"=paste0(id,"_",entity),
                  "data-accession"=payload$accession,
                  checked=if(paste0(id,"_",entity) %in% chosen) "checked" else NULL),
                tags$span("Select for analysis")),
              tags$div(class="ram-counterpart-id",
                tags$strong(id),
                tags$span(paste("Entity",entity)),
                tags$span(if(nzchar(chain)) paste("Chain",chain) else "Chain unknown")),
              tags$div(class="ram-counterpart-copy",
                tags$strong(get(item,"description","Protein entity")),
                tags$p(get(item,"title")),
                tags$div(class="ram-counterpart-meta",
                  tags$span(get(item,"method","Unknown method")),
                  tags$span(if(is.finite(resolution))
                    sprintf("%.2f Å",resolution) else "Resolution n/a"),
                  tags$span(if(is.finite(reference_coverage))
                    sprintf("%.1f%% UniProt sequence coverage",
                      100*reference_coverage)
                    else "UniProt coverage unreported"),
                  if(is.finite(entity_coverage))
                    tags$span(sprintf("%.1f%% entity sequence aligned",
                      100*entity_coverage)),
                  tags$span(get(item,"release_date"))),
              if(!is.null(verify))
                tags$p(class="ram-field-hint",
                  if(identical(verify$state,"searching"))
                    "Retrieving exact SIFTS residue mapping..."
                  else if(identical(verify$state,"error"))
                    paste("SIFTS unavailable:",verify$message)
                  else sprintf(paste0(
                    "Verified exact SIFTS: %d distinct PDB residues; ",
                    "%d observed; %d conflicting positions; %d unmatched ",
                    "SIFTS rows. Not a complete structure-state assessment."),
                    verify$unique_residues,verify$observed_residues,
                    verify$conflicting_residues,verify$unlinked_sifts_rows)),
                  if(identical(verify$state,"mapped"))
                    tags$span(if(length(verify$ca_points))
                      sprintf(" %d mapped C-alpha coordinates available for geometry comparison.",
                        length(verify$ca_points))
                      else if(nzchar(verify$ca_warning))
                        paste(" Geometry unavailable:",verify$ca_warning)
                      else " No C-alpha coordinates available."))),
            tags$div(class="ram-counterpart-actions",
              tags$a("RCSB entry",
                href=paste0("https://www.rcsb.org/structure/",id),
                target="_blank",rel="noopener noreferrer",
                class="btn btn-default btn-sm"),
              tags$button(type="button",
                class="btn btn-default btn-sm ram-atlas-verify",
                "data-pdb"=id,"data-entity"=entity,
                "Verify SIFTS mapping"),
              tags$button(type="button",class="btn btn-primary btn-sm ram-atlas-compare",
                "data-pdb"=id,"data-chain"=chain,"data-entity"=entity,
                "Compare with loaded structure")
            ))
        })),
      if(isTRUE(payload$has_more))
        actionButton("atlasLoadMore",
          sprintf("Load next %d experimental entities",
            min(50L,as.integer(payload$total_count-payload$next_offset))),
          class="btn-default btn-sm"),
      if(!isTRUE(payload$has_more) && !isTRUE(payload$stalled))
        tags$p(class="ram-field-hint",
          "The currently reported experimental search cohort has been retrieved.")
    )
  })

  # An experimental atlas geometry comparison is deliberately user-triggered.
  # The same canonical core is used for every distance and no group is given
  # a biological state name by the software.
  observeEvent(atlas_exact_results(), {
    atlas_geometry_result(NULL)
  },ignoreInit=TRUE)

  output$atlasGeometryPanel <- renderUI({
    cohort <- req(atlas_payload())
    verified <- atlas_exact_results()
    candidates <- tryCatch(
      ram_atlas_geometry_entities(verified,cohort$accession),
      error=function(e) list())
    eligible <- intersect(atlas_candidate_picks(),
      names(candidates)[vapply(candidates,function(x)
        nrow(x$coordinates)>=30L,logical(1L))])
    if(length(eligible)<2L) {
      if(!length(atlas_candidate_picks())) return(NULL)
      return(tags$p(class="ram-field-hint",
        "Select and verify at least two experimental entities to enable clustering. Each needs 30 observed, unambiguous C-alpha positions."))
    }
    prior <- isolate(input$atlasGeometryEntities)
    prior <- intersect(prior,eligible)
    if(length(prior)<2L) prior <- head(eligible,6L)
    tags$section(class="ram-panel",
      tags$h3("Experimental geometry groups"),
      tags$p(class="ram-field-hint",paste0(
        "Compare common canonical C-alpha distance maps. Groups are ",
        "exploratory geometric similarities, not verified functional states. ",
        "Only observed residues with exact SIFTS coordinates count.")),
      selectInput("atlasGeometryEntities","Verified experimental entities",
        choices=eligible,selected=prior,multiple=TRUE),
      uiOutput("atlasConstructReview"),
      radioButtons("atlasClusterMode","How should structures be grouped?",
        choices=c("Suggest groups automatically"="automatic",
          "Choose a distance cutoff myself"="manual"),
        selected="automatic",inline=TRUE),
      conditionalPanel(condition="input.atlasClusterMode === 'manual'",
        numericInput("atlasGeometryCutoff","Distance-map RMSD group cutoff (Å)",
          value=1.5,min=0.1,max=10,step=0.1)),
      tags$p(class="ram-field-hint",
        "Automatic mode evaluates geometry-only cohesion and separation on the same verified UniProt core. It can recommend one group when no clear split is supported. This does not infer functional states."),
      actionButton("atlasRunGeometry","Cluster experimental structures",
        class="btn-primary btn-sm"),
      uiOutput("atlasGeometrySummary"),
      plotOutput("atlasGeometryPlot",height="260px"),
      tableOutput("atlasGeometryTable"),
      uiOutput("atlasExperimentalContextPanel"),
      uiOutput("atlasLigandContextPanel"),
      uiOutput("atlasAutoQuality"),
      uiOutput("atlasRobustnessPanel"),
      uiOutput("atlasGroupsHandoffPanel"))
  })

  # Researcher-assigned group labels; the optional geometric clusters are
  # proposed starting points, never inferred functional-state assignments.
  output$atlasGroupsHandoffPanel <- renderUI({
    geometry <- atlas_geometry_result()
    if(is.null(geometry) || !is.null(geometry$error) ||
       length(geometry$selected)<2L) return(NULL)
    members <- split(as.character(geometry$assignment$entity),
      geometry$assignment$geometric_group)
    # A one-group automatic suggestion must not select the entire cohort
    # as Group A and leave Group B empty. Offer an explicit editable pair.
    default_a <- if(length(members)>1L) members[[1L]]
      else geometry$selected[[1L]]
    default_b <- if(length(members)>1L) members[[2L]]
      else geometry$selected[[2L]]
    tags$section(class="ram-atlas-group-handoff",
      tags$h4("Compare Atlas structures as groups"),
      tags$p(class="ram-field-hint",
        "Review the suggested clusters, then freely edit Group A and Group B. Geometric grouping is not a functional-state label. The analysis reuses exact SIFTS-mapped backbone atoms without reuploading files."),
      if(length(members)<2L)
        tags$p(class="ram-confidence-warning",
          "Atlas did not identify two well-separated clusters. These default single-entry groups are only a starting point for an explicitly researcher-defined comparison."),
      if(length(members)>2L)
        tags$p(class="ram-field-hint",
          sprintf("%d geometric clusters found. Group Compare accepts two sets: select which clusters or individual structures to compare.",
            length(members))),
      tags$div(class="ram-atlas-group-handoff-grid",
        selectInput("atlasGroupAEntities","Group A experimental entries",
          choices=geometry$selected,selected=default_a,multiple=TRUE),
        selectInput("atlasGroupBEntities","Group B experimental entries",
          choices=geometry$selected,selected=default_b,multiple=TRUE)),
      tags$p(class="ram-field-hint",
        "Each entry contributes one verified first-model polymer chain. Same-crystal copies are not independent biological replicates. Residues with unknown or differing chemistry have no classification label."),
      actionButton("atlasSendGroups","Use these entries in Compare Groups",
        class="btn-primary btn-sm"),
      uiOutput("atlasGroupHandoffStatus")
    )
  })
  output$atlasGroupHandoffStatus <- renderUI({
    current <- atlas_group_transfer()
    if(is.null(current)) return(NULL)
    tags$p(class="ram-field-hint",sprintf(
      "Selected %d Group A and %d Group B entries. Open Compare Groups to analyse.",
      length(current$group_a),length(current$group_b)))
  })
  observeEvent(atlas_geometry_result(), {
    if(!is.null(isolate(atlas_group_transfer()))) {
      atlas_group_transfer(NULL)
      if(identical(isolate(group_comparison_results())$source,"atlas"))
        group_comparison_results(NULL)
    }
  },ignoreInit=TRUE)
  observeEvent(input$atlasSendGroups, {
    geometry <- isolate(atlas_geometry_result())
    picked <- tryCatch(ram_atlas_group_select(geometry,
      isolate(input$atlasGroupAEntities),
      isolate(input$atlasGroupBEntities)),error=function(e)e)
    if(inherits(picked,"error")) {
      showNotification(conditionMessage(picked),type="warning",duration=14)
      return()
    }
    atlas_group_transfer(picked)
    group_comparison_results(NULL)
    updateRadioButtons(session,"groupInputMode",selected="atlas")
    updateTabsetPanel(session,"analysisTabs",selected="compare")
    session$sendCustomMessage("ram-open-group-panel",list())
    showNotification(sprintf(
      "Transferred %d + %d verified entries. Set group labels and choose Analyse groups.",
      length(picked$group_a),length(picked$group_b)),
      type="message",duration=12)
  },ignoreInit=TRUE)
  output$groupAtlasSelection <- renderUI({
    if(!identical(input$groupInputMode,"atlas")) return(NULL)
    picked <- atlas_group_transfer()
    if(is.null(picked))
      return(tags$p(class="ram-confidence-warning",
        "Select structures in Atlas, then choose Use these entries in Compare Groups."))
    tags$div(class="ram-panel",
      tags$strong(sprintf("Verified Atlas cohort: %s",picked$accession)),
      tags$p(paste0("Group A: ",paste(picked$group_a,collapse=", "))),
      tags$p(paste0("Group B: ",paste(picked$group_b,collapse=", "))),
      tags$p(class="ram-field-hint",
        "Exact observed UniProt positions define residue correspondence. Uploaded-file chain identity thresholds are not applied. Rama8000 classification is unavailable in the Atlas-only group analysis."),
      actionLink("groupBackToAtlas","Change group membership in Atlas"))
  })
  observeEvent(input$groupBackToAtlas, {
    updateTabsetPanel(session,"analysisTabs",selected="atlas")
  },ignoreInit=TRUE)

  atlas_construct_audit <- reactive({
    cohort <- req(atlas_payload())
    ids <- input$atlasGeometryEntities
    if(length(ids)<2L) return(NULL)
    ram_atlas_construct_audit(atlas_exact_results(),cohort$accession,ids)
  })
  output$atlasConstructReview <- renderUI({
    audit <- tryCatch(atlas_construct_audit(),
      error=function(e) list(error=conditionMessage(e)))
    if(is.null(audit)) return(NULL)
    if(!is.null(audit$error))
      return(tags$p(class="ram-confidence-warning",audit$error))
    tags$div(class="ram-atlas-construct-review",
      tags$h4("Construct and sequence-chemistry review"),
      tags$p(class="ram-field-hint",sprintf(
        "%d experimental-entity pairs compared using their shared observed UniProt positions.",
        audit$pair_count)),
      if(audit$differs)
        tags$p(class="ram-confidence-warning",
          "Verified polymer residue chemistry differs between some structures. Different residue chemistry can reflect mutations, modifications or construct differences; the geometry groups cannot isolate their effects."),
      if(audit$incomplete)
        tags$p(class="ram-field-hint",
          "Some experimental residues lack comparable chemistry metadata or occur outside the shared observed core. Unknown does not mean equivalent."),
      DT::DTOutput("atlasConstructTable"),
      downloadButton("downloadAtlasConstruct",
        "Export construct-comparability CSV"),
      if(audit$differs)
        checkboxInput("atlasConfirmDifferentChemistry",
          "Include structures with verified residue-chemistry differences in exploratory geometry grouping",
          value=FALSE),
      tags$p(class="ram-field-hint",
        "This is not a biological-state equivalence test. Structure preparation, missing regions, ligands, sequence changes and experimental conditions can affect the geometry."))
  })
  output$atlasConstructTable <- DT::renderDT({
    audit <- atlas_construct_audit()
    req(!is.null(audit))
    data <- audit$pairs[,c("entity_a","entity_b","common_observed",
      "chemistry_checked","chemistry_unknown","chemistry_differences",
      "unmatched_a","unmatched_b","difference_examples"),drop=FALSE]
    names(data) <- c("Structure A","Structure B","Common C-alpha",
      "Chemistry checked","Chemistry unknown","Different residues",
      "A-only","B-only","Differences (UniProt: A/B)")
    DT::datatable(data,rownames=FALSE,selection="none",
      options=list(pageLength=8,scrollX=TRUE,dom="tip"),
      class="compact stripe hover")
  },server=FALSE)
  output$downloadAtlasConstruct <- downloadHandler(
    filename=function() "ramplotr_atlas_construct_review.csv",
    content=function(file) {
      audit <- req(atlas_construct_audit())
      utils::write.csv(audit$pairs,file,row.names=FALSE,na="")
    }
  )
  atlas_chemistry_acknowledged <- reactiveVal(NULL)
  observeEvent(input$atlasGeometryEntities, {
    atlas_chemistry_acknowledged(NULL)
  },ignoreInit=TRUE)
  observeEvent(atlas_exact_results(), {
    atlas_chemistry_acknowledged(NULL)
  },ignoreInit=TRUE)
  observeEvent(input$atlasConfirmDifferentChemistry, {
    ids <- isolate(input$atlasGeometryEntities)
    atlas_chemistry_acknowledged(
      if(isTRUE(input$atlasConfirmDifferentChemistry) && length(ids)>=2L)
        sort(as.character(ids)) else NULL)
  },ignoreInit=TRUE)
  observeEvent(input$atlasRunGeometry, {
    cohort <- isolate(atlas_payload())
    if(is.null(cohort)) return()
    output <- tryCatch({
      verified <- isolate(atlas_exact_results())
      ids <- isolate(input$atlasGeometryEntities)
      audit <- ram_atlas_construct_audit(verified,cohort$accession,ids)
      if(audit$differs &&
         !identical(isolate(atlas_chemistry_acknowledged()),
                    sort(as.character(ids))))
        stop("Review the observed residue-chemistry differences and explicitly acknowledge their confounding effects before exploratory grouping.",
          call.=FALSE)
      geometry <- ram_atlas_geometry_groups(
        verified,cohort$accession,ids,
        min_core=30L,min_fraction=0.6,
        cutoff=isolate(input$atlasGeometryCutoff),max_core=300L,
        cluster_mode=isolate(input$atlasClusterMode))
      geometry$construct_audit <- audit
      geometry$experimental_context <- ram_atlas_experimental_context(
        cohort,geometry,audit)
      geometry$ligand_context <- ram_atlas_observed_ligand_context(
        verified,geometry)
      geometry
    },error=function(e) list(error=conditionMessage(e)))
    atlas_geometry_result(output)
  },ignoreInit=TRUE)

  output$atlasGeometrySummary <- renderUI({
    result <- atlas_geometry_result()
    if(is.null(result)) return(NULL)
    if(!is.null(result$error))
      return(tags$p(class="ram-confidence-warning",result$error))
    automatic <- result$auto_cluster
    count <- length(unique(result$assignment$geometric_group))
    tags$div(class="ram-field-hint",
      tags$p(sprintf(paste0(
        "%d experimental entities, %d common observed UniProt C-alpha ",
        "positions (%d sampled). %d exploratory geometry group(s)."),
        length(result$selected),length(result$common_positions),
        length(result$sampled_positions),count)),
      if(!is.null(automatic))
        tags$p(tags$strong("Automatic suggestion: "),
          automatic$reason)
      else tags$p(sprintf(
        "Manual distance-map RMSD cutoff: %.2f Å.",result$cutoff)),
      tags$p("Average-linkage clustering of C-alpha internal distance differences. All comparisons use the same verified canonical core. Neither automated clusters nor manual groups establish biological states, independent experimental replicates or ligand conditions."),
      tags$p(paste("Group representatives:",
        paste(unname(result$representatives),collapse=", "))))
  })
  output$atlasExperimentalContextPanel <- renderUI({
    result <- atlas_geometry_result()
    if(is.null(result) || !is.null(result$error) ||
       is.null(result$experimental_context)) return(NULL)
    evidence <- result$experimental_context
    tags$details(class="ram-details ram-atlas-context",
      tags$summary("Experimental method, construct and deposition context"),
      tags$p(class="ram-field-hint",sprintf(paste0(
        "%d selected PDB entities; %d missing method(s), %d missing ",
        "resolution(s), %d same-deposition pair(s)."),
        nrow(evidence$entries),evidence$missing_methods,
        evidence$missing_resolution,evidence$shared_pdb_pairs)),
      if(evidence$shared_pdb_pairs>0L)
        tags$p(class="ram-confidence-warning",
          "Some selected entities belong to the same PDB entry. Multiple chains from one deposition are not independent experimental measurements."),
      if(evidence$differing_chemistry_pairs>0L)
        tags$p(class="ram-confidence-warning",
          sprintf("%d pair(s) contain documented monomer-chemistry differences. Modified residues can differ without a gene mutation.",evidence$differing_chemistry_pairs)),
      tags$p(class="ram-field-hint",evidence$ligand_status),
      tags$p(class="ram-field-hint",
        "Experimental method, resolution and release date are from RCSB discovery. Exact observed residue chemistry and coverage are from the verified SIFTS/mmCIF audit. Different methods or resolution do not establish functional states."),
      tags$h5("Selected experimental entries"),
      DT::DTOutput("atlasContextEntries"),
      tags$h5("Pairwise metadata and construct caveats"),
      DT::DTOutput("atlasContextPairs"),
      tags$div(class="ram-ensemble-actions",
        downloadButton("downloadAtlasContextEntries",
          "Export experimental entries CSV"),
        downloadButton("downloadAtlasContextPairs",
          "Export experimental pair evidence CSV"))
    )
  })
  output$atlasContextEntries <- DT::renderDT({
    result <- req(atlas_geometry_result())
    req(is.null(result$error),!is.null(result$experimental_context))
    shown <- result$experimental_context$entries[,
      c("entity","group","method","resolution_A","initial_release_date",
        "observed_canonical_coverage","monomer_known","monomer_observed"),
      drop=FALSE]
    shown$resolution_A <- round(shown$resolution_A,2L)
    shown$observed_canonical_coverage <-
      round(100*shown$observed_canonical_coverage,1L)
    DT::datatable(shown,rownames=FALSE,
      colnames=c("PDB entity","Geometry group","Method","Resolution (Å)",
        "Released","Canonical core (%)","Monomers known",
        "Monomers observed"),
      options=list(pageLength=8,scrollX=TRUE,dom="ftip"),
      class="compact stripe")
  },server=FALSE)
  output$atlasContextPairs <- DT::renderDT({
    result <- req(atlas_geometry_result())
    req(is.null(result$error),!is.null(result$experimental_context))
    evidence <- result$experimental_context$pairs
    shown <- evidence[,c("entity_a","entity_b","same_geometric_group",
      "same_pdb_entry","method_a","method_b","chemistry_differences",
      "chemistry_unknown","unmatched_a","unmatched_b",
      "context_warning"),drop=FALSE]
    shown$same_geometric_group <- ifelse(shown$same_geometric_group,
      "Yes","No")
    shown$same_pdb_entry <- ifelse(shown$same_pdb_entry,"Yes","No")
    DT::datatable(shown,rownames=FALSE,
      colnames=c("PDB A","PDB B","Same group","Same deposition",
        "Method A","Method B","Chemistry differences",
        "Unknown chemistry","Unmatched A","Unmatched B","Caveats"),
      options=list(pageLength=8,scrollX=TRUE,dom="ftip"),
      class="compact stripe")
  },server=FALSE)
  output$downloadAtlasContextEntries <- downloadHandler(
    filename=function() "ramplotr_atlas_experimental_entries.csv",
    content=function(file) utils::write.csv(
      req(atlas_geometry_result())$experimental_context$entries,
      file,row.names=FALSE,na=""))
  output$downloadAtlasContextPairs <- downloadHandler(
    filename=function() "ramplotr_atlas_experimental_pairs.csv",
    content=function(file) utils::write.csv(
      req(atlas_geometry_result())$experimental_context$pairs,
      file,row.names=FALSE,na=""))

  output$atlasLigandContextPanel <- renderUI({
    result <- atlas_geometry_result()
    if(is.null(result) || !is.null(result$error) ||
       is.null(result$ligand_context)) return(NULL)
    context <- result$ligand_context
    tags$details(class="ram-details ram-atlas-ligand-context",
      tags$summary("Observed non-water components near the protein"),
      tags$p(class="ram-field-hint",sprintf(
        "Proximity extraction completed for %d of %d verified structures; %d have one or more nearby non-water components.",
        context$measured,nrow(context$entries),context$with_proximity)),
      if(context$unavailable>0L) tags$p(class="ram-confidence-warning",
        "Some structures have unavailable atom-level contact evidence. These are unknown, not ligand-free."),
      tags$p(class="ram-field-hint",
        "Reports first-model deposited non-water HETATM heavy atoms within 4.5 Å of the verified protein's mapped heavy atoms. Solvent waters and modified polymer residues are excluded. This is proximity, not proven binding, occupancy or an apo/holo classification."),
      tags$h5("Per-structure contact availability"),
      DT::DTOutput("atlasLigandEntries"),
      if(nrow(context$contacts)) tags$h5("Nearby deposited components"),
      if(nrow(context$contacts)) DT::DTOutput("atlasLigandContacts"),
      tags$div(class="ram-ensemble-actions",
        downloadButton("downloadAtlasLigandEntries",
          "Export contact availability CSV"),
        downloadButton("downloadAtlasLigandContacts",
          "Export observed proximity CSV")),
      tags$p(class="ram-field-hint",
        "No detected component within this distance cutoff cannot prove an apo state. Bound molecules may be unmodelled, absent from the deposition, too distant from mapped residues, or represented by polymer components. Review the experimental publication before assigning a ligand state.")
    )
  })
  output$atlasLigandEntries <- DT::renderDT({
    result <- req(atlas_geometry_result())
    req(is.null(result$error),!is.null(result$ligand_context))
    shown <- result$ligand_context$entries
    DT::datatable(shown,rownames=FALSE,
      colnames=c("PDB entity","Geometry group","Extraction",
        "Nearby components","Non-water sites","Radius (Å)","Caveat"),
      options=list(pageLength=8,scrollX=TRUE,dom="ftip"),
      class="compact stripe")
  },server=FALSE)
  output$atlasLigandContacts <- DT::renderDT({
    result <- req(atlas_geometry_result())
    req(is.null(result$error),!is.null(result$ligand_context))
    shown <- result$ligand_context$contacts
    DT::datatable(shown,rownames=FALSE,
      colnames=c("PDB entity","Deposited code","Asym ID","Author position",
        "Closest UniProt residue","Min heavy-atom distance (Å)",
        "Component heavy atoms"),
      options=list(pageLength=8,scrollX=TRUE,dom="ftip"),
      class="compact stripe")
  },server=FALSE)
  output$downloadAtlasLigandEntries <- downloadHandler(
    filename=function() "ramplotr_atlas_hetero_contact_availability.csv",
    content=function(file) utils::write.csv(
      req(atlas_geometry_result())$ligand_context$entries,
      file,row.names=FALSE,na=""))
  output$downloadAtlasLigandContacts <- downloadHandler(
    filename=function() "ramplotr_atlas_observed_hetero_proximity.csv",
    content=function(file) utils::write.csv(
      req(atlas_geometry_result())$ligand_context$contacts,
      file,row.names=FALSE,na=""))

  output$atlasAutoQuality <- renderUI({
    result <- atlas_geometry_result()
    if(is.null(result) || !is.null(result$error) ||
       is.null(result$auto_cluster)) return(NULL)
    quality <- result$auto_cluster$candidates
    if(!nrow(quality)) return(NULL)
    tags$details(class="ram-details",
      tags$summary("Why did Atlas suggest these groups?"),
      tags$p(class="ram-field-hint",
        "Mean silhouette compares each structure with its own and other clusters (range -1 to 1; singletons score 0). The distance gap is median between-cluster minus within-cluster C-alpha dRMSD in Å. Exploratory safeguards require silhouette ≥0.50 and gap ≥0.35 Å; these values are not biologically calibrated thresholds."),
      tags$table(class="table table-condensed",
        tags$thead(tags$tr(lapply(
          c("Groups","Mean silhouette","Distance gap (Å)","Singleton groups"),
          tags$th))),
        tags$tbody(lapply(seq_len(nrow(quality)),function(i)
          tags$tr(lapply(as.list(quality[i,c("k","silhouette",
            "separation_A","singleton_groups")]),tags$td))))),
      tags$p(class="ram-field-hint",
        "For two structures, no automatic split is proposed. Switch to manual mode to inspect a two-structure comparison.")
    )
  })

  output$atlasRobustnessPanel <- renderUI({
    result <- atlas_geometry_result()
    if(is.null(result) || !is.null(result$error) ||
       is.null(result$robustness)) return(NULL)
    stability <- result$robustness
    if(identical(stability$status,"insufficient_structures"))
      return(tags$p(class="ram-field-hint",
        "Cluster robustness: two experimental structures cannot provide a meaningful sensitivity estimate. You can still inspect their measured structural difference."))
    tags$details(class="ram-details ram-atlas-robustness",
      tags$summary("Cluster robustness and shared-entry cautions"),
      tags$p(class="ram-field-hint",
        sprintf("%d of %d leave-block-out analyses reproduced the full geometric grouping.",
          sum(stability$replicates$baseline_partition_reproduced),
          stability$iterations)),
      tags$p(class="ram-field-hint",
        paste("Contiguous blocks of sampled canonical UniProt positions are",
          "omitted one at a time. Atlas recalculates C-alpha distance maps",
          "and reruns the same automatic or manual clustering.",
          "A stable partition shows robustness to this particular position",
          "selection, not statistical confidence or biological states.")),
      if(stability$shared_pdb_pairs>0L)
        tags$p(class="ram-confidence-warning",
          sprintf(paste("%d selected entity pair(s) originate in the same",
            "PDB entry. Polymer chains from a single deposition cannot be",
            "treated as independent experimental replicates."),
            stability$shared_pdb_pairs)),
      tags$h5("Exact group recovery"),
      tableOutput("atlasRobustnessGroups"),
      tags$h5("Pairwise grouping under position deletion"),
      DT::DTOutput("atlasRobustnessPairs"),
      tags$div(class="ram-ensemble-actions",
        downloadButton("downloadAtlasRobustnessGroups",
          "Export group robustness CSV"),
        downloadButton("downloadAtlasRobustnessPairs",
          "Export pair robustness CSV"),
        downloadButton("downloadAtlasRobustnessRuns",
          "Export omitted-block runs CSV")),
      tags$p(class="ram-field-hint",
        "Group recovery requires the exact same member set, regardless of numeric group labels. Pairwise co-assignment is the fraction of block omissions keeping two entries together. It is not a transition probability.")
    )
  })
  output$atlasRobustnessGroups <- renderTable({
    result <- req(atlas_geometry_result())
    req(is.null(result$error),identical(result$robustness$status,"ok"))
    shown <- result$robustness$groups
    shown$exact_group_fraction <- sprintf("%.0f%%",
      100*shown$exact_group_fraction)
    names(shown) <- c("Group","Size","PDB entities","Recovered",
      "Runs","Recovery","Same-PDB pairs")
    shown
  },striped=TRUE,spacing="xs",rownames=FALSE)
  output$atlasRobustnessPairs <- DT::renderDT({
    result <- req(atlas_geometry_result())
    req(is.null(result$error),identical(result$robustness$status,"ok"))
    shown <- result$robustness$pairs
    shown$same_group_fraction <- round(100*shown$same_group_fraction,1L)
    shown$baseline_same_group <- ifelse(shown$baseline_same_group,"Yes","No")
    shown$same_pdb_entry <- ifelse(shown$same_pdb_entry,"Yes","No")
    DT::datatable(shown,rownames=FALSE,
      colnames=c("Entry A","Entry B","Same baseline group",
        "Together after omissions","Runs","Together (%)","Same PDB"),
      options=list(pageLength=8,scrollX=TRUE,dom="ftip"),
      class="compact stripe")
  },server=FALSE)
  output$downloadAtlasRobustnessGroups <- downloadHandler(
    filename=function() "ramplotr_atlas_cluster_group_robustness.csv",
    content=function(file) utils::write.csv(
      req(atlas_geometry_result())$robustness$groups,
      file,row.names=FALSE,na=""))
  output$downloadAtlasRobustnessPairs <- downloadHandler(
    filename=function() "ramplotr_atlas_cluster_pair_robustness.csv",
    content=function(file) utils::write.csv(
      req(atlas_geometry_result())$robustness$pairs,
      file,row.names=FALSE,na=""))
  output$downloadAtlasRobustnessRuns <- downloadHandler(
    filename=function() "ramplotr_atlas_cluster_position_omissions.csv",
    content=function(file) utils::write.csv(
      req(atlas_geometry_result())$robustness$replicates,
      file,row.names=FALSE,na=""))

  output$atlasGeometryPlot <- renderPlot({
    result <- atlas_geometry_result()
    req(!is.null(result),is.null(result$error))
    ram_atlas_geometry_plot(result)
  })

  output$atlasGeometryTable <- renderTable({
    result <- atlas_geometry_result()
    req(!is.null(result),is.null(result$error))
    table <- result$assignment
    names(table) <- c("PDB entity","Chain","Mapped C-alpha",
      "Common core fraction","Geometry group")
    group_counts <- base::table(table[["Geometry group"]])
    table[["Group size"]] <- as.integer(group_counts[
      as.character(table[["Geometry group"]])])
    table
  },striped=TRUE,spacing="xs",rownames=FALSE)

  # Local torsion comparisons are secondary evidence for the previously
  # computed exploratory global geometric groups. They do not assign a
  # functional state or claim ligand-driven transitions.
  observeEvent(atlas_geometry_result(), {
    atlas_switch_result(NULL)
    atlas_selected_position(NULL)
  },ignoreInit=TRUE)

  output$atlasSwitchPanel <- renderUI({
    geometry <- atlas_geometry_result()
    if(is.null(geometry) || !is.null(geometry$error)) return(NULL)
    available <- as.character(geometry$selected)
    if(length(available)<2L) return(NULL)
    defaults <- unname(geometry$representatives)
    if(length(defaults)<2L) defaults <- available[1:2]
    tags$div(class="ram-panel",
      tags$h4("Local backbone-change candidates"),
      tags$p(class="ram-field-hint",
        paste("Compare backbone φ/ψ between any two verified experimental",
        "structures, even when automatic clustering suggests one group.",
        "Exact SIFTS positions, complete N/CA/C atoms and continuous peptide",
        "bonds are required. These are exploratory changes, not validated functional states.")),
      if(length(geometry$representatives)<2L)
        tags$p(class="ram-confidence-warning",
          "No convincing geometric split was suggested. You can inspect two selected structures, but do not treat them as separate structural states."),
      selectInput("atlasSwitchRepresentativeA","Reference experimental structure",
        choices=available,selected=defaults[[1L]]),
      selectInput("atlasSwitchRepresentativeB","Other experimental structure",
        choices=available,selected=defaults[[2L]]),
      numericInput("atlasSwitchThreshold",
        "Combined circular φ/ψ change threshold (degrees)",
        value=30,min=5,max=180,step=5),
      actionButton("atlasRunSwitch","Find local backbone changes",
        class="btn-primary btn-sm"),
      uiOutput("atlasSwitchSummary"),
      plotOutput("atlasSwitchPlot",height="230px"),
      tableOutput("atlasSwitchRegions"),
      tags$p(class="ram-field-hint","Click a residue to inspect the angles and exact experimental PDB identifiers."),
      DT::DTOutput("atlasSwitchResidues"),
      uiOutput("atlasSwitchResidueInspector"),
      uiOutput("atlasSwitchExport"))
  })

  observeEvent(input$atlasRunSwitch, {
    geometry <- isolate(atlas_geometry_result())
    if(is.null(geometry) || !is.null(geometry$error)) return()
    result <- tryCatch(ram_atlas_group_switches(
      isolate(atlas_exact_results()),geometry,
      isolate(input$atlasSwitchThreshold),
      representative_ids=c(isolate(input$atlasSwitchRepresentativeA),
        isolate(input$atlasSwitchRepresentativeB))),
      error=function(e) list(error=conditionMessage(e)))
    atlas_switch_result(result)
    atlas_selected_position(NULL)
  },ignoreInit=TRUE)

  output$atlasSwitchSummary <- renderUI({
    result <- atlas_switch_result()
    if(is.null(result))return(NULL)
    if(!is.null(result$error))
      return(tags$p(class="ram-confidence-warning",result$error))
    tags$div(class="ram-field-hint",
      tags$p(sprintf(
        "%s versus %s: %d of %d canonical positions have both valid φ/ψ pairs; %d candidate residue(s) at ≥%.0f°.",
        result$representatives[[1L]],result$representatives[[2L]],
        result$comparable,result$total_positions,
        sum(result$residues$candidate),result$threshold)),
      tags$p(sprintf(
        "%d contiguous candidate segment(s), including isolated positions. Missing and broken-backbone angles are excluded.",
        nrow(result$regions))))
  })

  output$atlasSwitchPlot <- renderPlot({
    result <- atlas_switch_result()
    req(!is.null(result),is.null(result$error))
    data <- result$residues
    if(!nrow(data) || !any(is.finite(data$angular_shift))) {
      graphics::plot.new()
      graphics::text(0.5,0.5,"No comparable complete backbone torsions",
        cex=0.95)
    } else {
      shift <- data$angular_shift
      graphics::plot(data$uniprot_resi,shift,type="h",lwd=2,
        xlab="UniProt residue",ylab="Wrapped φ/ψ displacement (°)",
        main="Local backbone differences")
      graphics::abline(h=result$threshold,lty=2,col="gray50")
      graphics::points(data$uniprot_resi[data$candidate],
        shift[data$candidate],pch=16)
    }
  })

  output$atlasSwitchRegions <- renderTable({
    result <- atlas_switch_result()
    req(!is.null(result),is.null(result$error))
    if(!nrow(result$regions))return(NULL)
    regions <- result$regions
    names(regions) <- c("Start UniProt","End UniProt","Residues",
      "Mean angular change (°)","Peak change (°)")
    regions
  },striped=TRUE,spacing="xs",rownames=FALSE)

  output$atlasSwitchResidues <- DT::renderDT({
    result <- atlas_switch_result()
    req(!is.null(result),is.null(result$error))
    rows <- result$residues
    if(!nrow(rows)) return(DT::datatable(data.frame()))
    view <- rows[,c("uniprot_resi","chain_a","resi_a","insertion_a",
      "chain_b","resi_b","insertion_b","phi_a","psi_a",
      "phi_b","psi_b","angular_shift","candidate"),drop=FALSE]
    for(col in c("phi_a","psi_a","phi_b","psi_b","angular_shift"))
      view[[col]] <- round(view[[col]],1)
    view$candidate <- ifelse(view$candidate,"Candidate","")
    DT::datatable(view,rownames=FALSE,selection="single",
      colnames=c("UniProt","Chain A","PDB res. A","Ins. A",
        "Chain B","PDB res. B","Ins. B","φ A","ψ A",
        "φ B","ψ B","Shift (°)","Review"),
      options=list(pageLength=12,scrollX=TRUE,dom="ftip",
        order=list(list(11,"desc"))),
      class="compact stripe hover")
  },server=FALSE)
  observeEvent(input$atlasSwitchResidues_rows_selected, {
    result <- isolate(atlas_switch_result())
    if(is.null(result) || !is.null(result$error)) return()
    ix <- suppressWarnings(as.integer(input$atlasSwitchResidues_rows_selected))
    if(length(ix)!=1L || is.na(ix) || ix<1L ||
       ix>nrow(result$residues)) return()
    atlas_selected_position(result$residues$uniprot_resi[[ix]])
  },ignoreInit=TRUE)
  output$atlasSwitchResidueInspector <- renderUI({
    result <- atlas_switch_result()
    pos <- atlas_selected_position()
    if(is.null(result) || !is.null(result$error) || is.null(pos))
      return(NULL)
    rows <- result$residues[result$residues$uniprot_resi==pos,,drop=FALSE]
    if(nrow(rows)!=1L) return(NULL)
    row <- rows[1L,,drop=FALSE]
    format_angle <- function(x) if(is.finite(x)) sprintf("%.1f°",x)
      else "Unavailable"
    local_id <- function(side) {
      chain <- row[[paste0("chain_",side)]][[1L]]
      resi <- row[[paste0("resi_",side)]][[1L]]
      insertion <- row[[paste0("insertion_",side)]][[1L]]
      if(is.na(chain) || is.na(resi)) return("Not mapped/observed")
      paste0("Chain ",chain,", residue ",resi,
        if(!is.na(insertion) && nzchar(insertion)) insertion else "")
    }
    tags$section(class="ram-panel ram-atlas-residue-inspector",
      tags$h4(sprintf("UniProt residue %d",pos)),
      tags$p(class="ram-field-hint",
        if(isTRUE(row$comparable[[1L]]))
          sprintf("Circular φ/ψ displacement: %.1f°%s",
            row$angular_shift[[1L]],
            if(isTRUE(row$candidate[[1L]])) " · above review threshold" else "")
        else "Incomplete torsions: no paired angle displacement can be calculated."),
      tags$div(class="ram-atlas-inspection-pair",
        tags$div(tags$strong(result$representatives[[1L]]),
          tags$p(local_id("a")),
          tags$p(paste("φ",format_angle(row$phi_a[[1L]]),
            "· ψ",format_angle(row$psi_a[[1L]])))),
        tags$div(tags$strong(result$representatives[[2L]]),
          tags$p(local_id("b")),
          tags$p(paste("φ",format_angle(row$phi_b[[1L]]),
            "· ψ",format_angle(row$psi_b[[1L]]))))),
      tags$p(class="ram-field-hint",
        "These are first-model experimental torsions at the same exact SIFTS UniProt position. A large shift is an inspection candidate, not proof of a functional transition."),
      tags$p(class="ram-field-hint",
        "Opens both selected experimental representatives in Compare. This replaces the current primary analysis."),
      actionButton("atlasInspectCompare","Inspect both structures in 2D/3D",
        class="btn-default btn-sm"))
  })
  observeEvent(input$atlasInspectCompare, {
    result <- isolate(atlas_switch_result())
    pos <- isolate(atlas_selected_position())
    if(is.null(result) || !is.null(result$error) || is.null(pos)) return()
    row <- result$residues[result$residues$uniprot_resi==pos,,drop=FALSE]
    if(nrow(row)!=1L || anyNA(row[,c("chain_a","resi_a",
        "insertion_a","chain_b","resi_b","insertion_b"),drop=FALSE])) {
      showNotification("Both experimental residues must have exact author identifiers before paired 3D inspection.",
        type="warning",duration=12)
      return()
    }
    ids <- as.character(result$representatives)
    if(length(ids)!=2L || any(!grepl("^[A-Z0-9]{4}_[1-9][0-9]*$",ids)))
      return()
    pdb_a <- substr(ids[[1L]],1L,4L)
    pdb_b <- substr(ids[[2L]],1L,4L)
    geometry <- isolate(atlas_geometry_result())
    verified <- isolate(atlas_exact_results())
    available <- tryCatch(
      ram_atlas_geometry_entities(verified,geometry$accession),
      error=function(e) list())
    if(any(!ids %in% names(available))) {
      showNotification("Exact verified Atlas representatives are unavailable.",
        type="error",duration=12)
      return()
    }
    if(!identical(as.character(row$chain_a[[1L]]),
                  as.character(available[[ids[[1L]]]]$chain)) ||
       !identical(as.character(row$chain_b[[1L]]),
                  as.character(available[[ids[[2L]]]]$chain))) {
      showNotification("The Atlas residue no longer matches the verified representative chains.",
        type="error",duration=12)
      return()
    }
    context <- list(pdb_a=pdb_a,pdb_b=pdb_b,
      entity_a=ids[[1L]],entity_b=ids[[2L]],
      accession=geometry$accession,
      chain_a=as.character(row$chain_a[[1L]]),
      chain_b=as.character(row$chain_b[[1L]]),
      asym_a=available[[ids[[1L]]]]$struct_asym_id,
      asym_b=available[[ids[[2L]]]]$struct_asym_id,
      map_a=verified[[ids[[1L]]]]$mapping,
      map_b=verified[[ids[[2L]]]]$mapping)
    atlas_pair_handoff(NULL)
    atlas_alignment_context(NULL)
    # An unrelated primary protein cannot serve as the Atlas representative.
    first_ok <- isTRUE(load_primary_structure("pdb",pdb_a))
    if(!first_ok) {
      showNotification(paste("Could not load Atlas representative",ids[[1L]]),
        type="error",duration=12)
      return()
    }
    second_ok <- isTRUE(load_comparison_structure(pdb_id=pdb_b,
      preferred_chain_a=as.character(row$chain_a[[1L]]),
      preferred_chain_b=as.character(row$chain_b[[1L]])))
    if(!second_ok) {
      showNotification(paste("Primary loaded, but comparison representative",
        ids[[2L]],"could not be loaded."),type="error",duration=12)
      return()
    }
    updateRadioButtons(session,"inputSource",selected="pdb")
    updateTextInput(session,"PDB",value=pdb_a)
    atlas_alignment_context(context)
    atlas_pair_handoff(list(
      pdb_a=pdb_a,pdb_b=pdb_b,entity_a=ids[[1L]],entity_b=ids[[2L]],
      position=as.integer(pos),residue=row))
    updateTabsetPanel(session,"analysisTabs",selected="compare")
  },ignoreInit=TRUE)

  output$atlasSwitchExport <- renderUI({
    result <- atlas_switch_result()
    if(is.null(result) || !is.null(result$error) ||
       !nrow(result$residues))return(NULL)
    downloadButton("downloadAtlasSwitch",
      "Export canonical φ/ψ differences CSV")
  })
  output$downloadAtlasSwitch <- downloadHandler(
    filename=function() {
      result <- req(atlas_switch_result())
      paste0("ramplotr_atlas_",result$representatives[[1L]],"_",
        result$representatives[[2L]],"_backbone_changes.csv")
    },
    content=function(file) {
      result <- req(atlas_switch_result())
      if(!is.null(result$error))
        stop("Backbone change comparison failed.")
      utils::write.csv(result$residues,file,row.names=FALSE,na="")
    }
  )

  output$downloadAtlasCanonical <- downloadHandler(
    filename=function() {
      cohort <- req(atlas_payload())
      paste0("ramplotr_atlas_",ram_uniprot_accession(cohort$accession),
             "_exact_residues.csv")
    },
    content=function(file) {
      cohort <- req(atlas_payload())
      rows <- ram_atlas_cohort_table(isolate(atlas_exact_results()),
                                     cohort$accession)
      if(!nrow(rows)) stop("No verified SIFTS residue data to export.")
      utils::write.csv(rows,file,row.names=FALSE,na="")
    }
  )
  output$downloadAtlasSupport <- downloadHandler(
    filename=function() {
      cohort <- req(atlas_payload())
      paste0("ramplotr_atlas_",ram_uniprot_accession(cohort$accession),
             "_observed_position_support.csv")
    },
    content=function(file) {
      cohort <- req(atlas_payload())
      positions <- ram_atlas_cohort_position_support(
        isolate(atlas_exact_results()),cohort$accession)
      utils::write.csv(positions,file,row.names=FALSE,na="")
    }
  )

  observeEvent(input$ramAtlasComparePick, {
    candidate <- input$ramAtlasComparePick
    primary <- isolate(loaded())
    if(is.null(primary)) {
      showNotification("Load a primary structure first, then choose an Atlas counterpart.",
        type="warning",duration=10)
      return()
    }
    if(!is.list(candidate) || is.null(candidate$pdb_id)) return()
    id <- toupper(trimws(as.character(candidate$pdb_id)))
    if(!grepl("^[A-Z0-9]{4}$",id)) return()
    chain_b <- if(is.null(candidate$chain)) "" else
      as.character(candidate$chain)
    chain_a <- if(length(primary$chains)) primary$chains[[1L]] else ""
    if(!is.null(primary$uniprot_accession) &&
       !is.null(isolate(atlas_payload())) &&
       !identical(as.character(primary$uniprot_accession),
                  as.character(isolate(atlas_payload())$accession)))
      showNotification(
        "Loaded AlphaFold model and Atlas query refer to different UniProt accessions. Check biological comparability.",
        type="warning",duration=12)
    segments <- isolate(canonical_segments())
    searched <- isolate(atlas_payload())
    if(!is.null(searched) && nrow(segments)) {
      matches <- segments$chain[segments$uniprot_accession ==
        as.character(searched$accession)]
      matches <- intersect(matches,primary$chains)
      if(length(matches)) chain_a <- matches[[1L]]
    }
    success <- load_comparison_structure(pdb_id=id,
      preferred_chain_a=chain_a,preferred_chain_b=chain_b)
    if(isTRUE(success))
      updateTabsetPanel(session,"analysisTabs",selected="compare")
  },ignoreInit=TRUE)

  output$experimentalCounterpartPanel <- renderUI({
    structure <- req(loaded())
    prediction <- structure$prediction
    if (is.null(prediction) || current_model() != 1L) return(NULL)
    tags$details(id="ram-experimental-counterparts",
      class="ram-confidence-panel ram-counterpart-panel",
      tags$summary(
        tags$span(class="ram-confidence-title",
          "Experimental counterparts"),
        tags$span(class="ram-confidence-subtitle",
          "Find related PDB structures and compare local backbone conformations")
      ),
      tags$div(class="ram-confidence-body",
        tags$p(class="ram-confidence-explainer",
          "Search experimental PDB polymer entities by sequence similarity, then open a candidate directly in RamplotR's linked comparison view. A sequence match is not evidence that two structures represent the same functional or biochemical state."),
        tags$div(class="ram-counterpart-controls",
          selectInput("experimentalSearchChain","Prediction chain",
            choices=structure$chains,selected=structure$chains[[1L]],
            selectize=FALSE),
          selectInput("experimentalIdentity","Minimum sequence identity",
            choices=c("100%"="1","95%"="0.95","90%"="0.90",
                      "70%"="0.70","50%"="0.50"),
            selected="0.90",selectize=FALSE),
          actionButton("findExperimentalStructures",
            "Find experimental structures",class="btn-primary btn-sm")
        ),
        uiOutput("experimentalCounterpartStatus"),
        uiOutput("experimentalCounterpartResults"),
        tags$p(class="ram-field-hint",
          "Search uses the public RCSB PDB sequence service and requests experimental entries only. Entry links open the corresponding PDBe page.")
      )
    )
  })

  output$experimentalCounterpartStatus <- renderUI({
    status <- experimental_search_status()
    if (is.null(status)) return(NULL)
    state <- if (is.null(status$state)) "" else as.character(status$state)
    message <- if (is.null(status$message)) "" else as.character(status$message)
    cls <- if (identical(state,"error")) "ram-confidence-warning"
      else if (identical(state,"searching")) "ram-counterpart-searching"
      else "ram-field-hint"
    tags$p(class=cls,message)
  })

  output$experimentalCounterpartResults <- renderUI({
    payload <- experimental_search_results()
    if (is.null(payload)) return(NULL)
    results <- payload$results
    if (is.null(results) || !length(results))
      return(tags$div(class="ram-counterpart-empty",
        tags$strong("No experimental matches at this threshold."),
        tags$p("Try a lower sequence-identity threshold if a more distant structural homologue would still be informative.")))
    total <- suppressWarnings(as.integer(payload$total_count))
    card <- function(item) {
      get <- function(name, default="") {
        value <- item[[name]]
        if (is.null(value) || !length(value) || is.na(value[[1L]]))
          default else as.character(value[[1L]])
      }
      pdb_id <- toupper(get("pdb_id"))
      entity <- get("entity_id")
      chain <- get("chain")
      title <- get("title")
      description <- get("description","Protein polymer entity")
      method <- get("method","Experimental structure")
      resolution <- suppressWarnings(as.numeric(get("resolution",NA_character_)))
      resolution_label <- if (is.finite(resolution))
        sprintf("%.2f Å",resolution) else "Resolution n/a"
      chains <- item[["chains"]]
      chain_label <- if (!is.null(chains) && length(chains))
        paste(as.character(unlist(chains)),collapse=", ") else
        if (nzchar(chain)) chain else "n/a"
      pdbe <- paste0("https://www.ebi.ac.uk/pdbe/entry/pdb/",
                     tolower(pdb_id))
      tags$article(class="ram-counterpart-card",
        tags$div(class="ram-counterpart-card-main",
          tags$div(class="ram-counterpart-id",
            tags$strong(pdb_id),
            tags$span(paste("Entity",entity)),
            tags$span(paste("Chain",chain_label))
          ),
          tags$div(class="ram-counterpart-copy",
            tags$strong(description),
            if (nzchar(title) && !identical(title,description))
              tags$p(title),
            tags$div(class="ram-counterpart-meta",
              tags$span(method),tags$span(resolution_label))
          )
        ),
        tags$div(class="ram-counterpart-actions",
          tags$a("PDBe entry",href=pdbe,target="_blank",
            rel="noopener noreferrer",class="btn btn-default btn-sm"),
          tags$button(type="button",
            class="btn btn-primary btn-sm ram-experimental-compare",
            "data-pdb"=pdb_id,"data-entity"=entity,
            "data-chain"=chain,"data-title"=title,
            "Compare in RamplotR")
        )
      )
    }
    tags$div(class="ram-counterpart-results",
      tags$div(class="ram-counterpart-result-head",
        tags$strong(sprintf("%d candidate%s shown",
          length(results),if(length(results)==1L) "" else "s")),
        if (is.finite(total) && total > length(results))
          tags$span(sprintf("%d total hits met the sequence threshold",total))
      ),
      lapply(results,card)
    )
  })

  observeEvent(loaded(), {
    experimental_search_results(NULL)
    experimental_search_status(NULL)
    experimental_search_request(experimental_search_request()+1L)
  }, ignoreInit=TRUE)

  observeEvent(input$findExperimentalStructures, {
    structure <- req(loaded())
    if (is.null(structure$prediction) || current_model() != 1L) return()
    chain <- req(input$experimentalSearchChain)
    query <- ram_chain_query_sequence(classified(),chain)
    if (nchar(query$sequence) < 20L) {
      experimental_search_status(list(state="error",
        message="The selected chain is too short for a useful sequence search."))
      return()
    }
    if (!is.finite(query$known_fraction) || query$known_fraction < 0.60) {
      experimental_search_status(list(state="error",
        message="Too much of this chain is unresolved or non-standard for a reliable protein-sequence search."))
      return()
    }
    identity <- suppressWarnings(as.numeric(input$experimentalIdentity))
    if (!is.finite(identity)) identity <- 0.90
    request_id <- experimental_search_request()+1L
    experimental_search_request(request_id)
    experimental_search_results(NULL)
    experimental_search_status(list(state="searching",
      message=sprintf("Searching experimental structures related to chain %s…",chain)))
    session$sendCustomMessage("ram-experimental-search",list(
      request_id=as.character(request_id),
      sequence=query$sequence,
      identity_cutoff=identity,
      rows=12L
    ))
  },ignoreInit=TRUE)

  observeEvent(input$ramExperimentalSearchStatus, {
    value <- input$ramExperimentalSearchStatus
    if (!is.list(value) || is.null(value$request_id)) return()
    if (!identical(as.character(value$request_id),
                   as.character(isolate(experimental_search_request())))) return()
    experimental_search_status(value)
  },ignoreInit=TRUE)

  observeEvent(input$ramExperimentalSearchResults, {
    value <- input$ramExperimentalSearchResults
    if (!is.list(value) || is.null(value$request_id)) return()
    if (!identical(as.character(value$request_id),
                   as.character(isolate(experimental_search_request())))) return()
    experimental_search_results(value)
  },ignoreInit=TRUE)

  observeEvent(input$ramExperimentalComparePick, {
    item <- input$ramExperimentalComparePick
    if (!is.list(item) || is.null(item$pdb_id)) return()
    pdb_id <- toupper(trimws(as.character(item$pdb_id)))
    chain_b <- if (is.null(item$chain)) "" else as.character(item$chain)
    chain_a <- isolate(input$experimentalSearchChain)
    if (!grepl("^[A-Za-z0-9]{4}$",pdb_id)) {
      showNotification("Invalid PDB identifier returned by the search.",
                       type="error")
      return()
    }
    updateRadioButtons(session,"compareInputSource",selected="pdb")
    updateTextInput(session,"comparePDB",value=pdb_id)
    ok <- load_comparison_structure(
      pdb_id=pdb_id,
      preferred_chain_a=chain_a,
      preferred_chain_b=chain_b
    )
    if (isTRUE(ok)) {
      updateTabsetPanel(session,"analysisTabs",selected="compare")
      showNotification(
        paste("Loaded",pdb_id,"for experimental comparison."),
        type="message",duration=5)
    }
  },ignoreInit=TRUE)

  output$predictionReviewMetrics <- renderUI({
    if (is.null(req(loaded())$prediction) || current_model() != 1L)
      return(NULL)
    data <- classified()
    if (!"plddt" %in% names(data)) return(NULL)
    high_standard_outliers <- sum(is.finite(data$plddt) & data$plddt >= 90 &
      !is.na(data$rama8000_region) & data$rama8000_region == "Outlier")
    high_not_allowed <- sum(is.finite(data$plddt) & data$plddt >= 90 &
      !is.na(data$region) & data$region == "Not allowed")
    lower_inrange <- sum(is.finite(data$plddt) & data$plddt < 70 &
      !is.na(data$rama8000_region) & data$rama8000_region != "Outlier")
    tagList(
      tags$span(class = if (high_standard_outliers > 0L)
        "ram-confidence-metric ram-review-high" else "ram-confidence-metric",
        paste(high_standard_outliers, "high-confidence Rama8000 outliers")),
      tags$span(class = "ram-confidence-metric",
        paste(high_not_allowed, "high-confidence RamplotR Not allowed")),
      tags$span(class = "ram-confidence-metric",
        paste(lower_inrange, "lower-confidence residues without a Rama8000 outlier"))
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
    req(loaded(), comparison_loaded(), input$compareChainA, input$compareChainB,
        input$bgtype, input$validationMode)
    main_data <- classified()
    secondary_data <- ram_classify_torsions(comparison_torsions(),
      reference_dir = file.path("static", input$bgtype),
      selected_reference=plot_reference(),
      mode=input$validationMode,
      threshold_fn=ram_density_thresholds)
    secondary_data <- ram_rama8000_classify(
      secondary_data, file.path("static", "rama8000"))
    comparison <- comparison_loaded()
    comparison_model <- if (is.null(input$compareModel)) 1L else
      suppressWarnings(as.integer(input$compareModel))
    if (length(comparison_model)!=1L || is.na(comparison_model) ||
        comparison_model<1L || comparison_model>comparison$nmodels)
      comparison_model <- 1L
    if (!is.null(comparison$prediction) && comparison_model==1L)
      secondary_data <- ram_apply_prediction(
        secondary_data,comparison$prediction)

    if (isTRUE(compare_swapped())) {
      first <- secondary_data[
        secondary_data$chain == input$compareChainA, , drop=FALSE]
      second <- main_data[
        main_data$chain == input$compareChainB, , drop=FALSE]
    } else {
      first <- main_data[
        main_data$chain == input$compareChainA, , drop=FALSE]
      second <- secondary_data[
        secondary_data$chain == input$compareChainB, , drop=FALSE]
    }
    if (!nrow(first) || !nrow(second)) return(data.frame())
    context <- atlas_alignment_context()
    comparison_model <- if(is.null(input$compareModel)) 1L else
      suppressWarnings(as.integer(input$compareModel))
    canonical <- !is.null(context) && !isTRUE(compare_swapped()) &&
      identical(loaded()$pdb_accession,context$pdb_a) &&
      identical(comparison$source_id,context$pdb_b) &&
      identical(input$compareChainA,context$chain_a) &&
      identical(input$compareChainB,context$chain_b) &&
      identical(current_model(),1L) &&
      identical(comparison_model,1L)
    if(canonical) {
      mapped <- ram_atlas_canonical_pairing(first,second,
        context$map_a,context$map_b,context$accession,
        context$chain_a,context$asym_a,
        context$chain_b,context$asym_b)
      result <- if(nrow(mapped$pairing))
        ram_compare_torsions(first,second,pairing=mapped$pairing)
      else data.frame()
      attr(result,"atlas_canonical") <- mapped
    } else {
      result <- ram_compare_torsions(first,second)
    }
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
    # The global inspector belongs to the structure loaded in the main app.
    # If comparison roles are swapped, that structure is side B rather than A.
    main_side <- if (isTRUE(isolate(compare_swapped()))) "b" else "a"
    number <- row[[paste0("residue_",main_side)]][[1L]]
    if (!is.na(number)) {
      chain <- row[[paste0("chain_",main_side)]][[1L]]
      insertion <- row[[paste0("insertion_",main_side)]][[1L]]
      visible <- isolate(displayed())
      matches <- which(visible$chain == chain &
        visible$resi == number & visible$insertion_code == insertion)
      if (length(matches)) selected_residue(list(
        chain=chain, resi=as.integer(number), insertion_code=insertion))
    }
    invisible(TRUE)
  }
  # Once the correct models and chains are aligned, use the existing
  # shared comparison selection to focus the same exact PDB residue pair.
  observe({
    handoff <- atlas_pair_handoff()
    req(handoff)
    primary <- req(loaded())
    secondary <- req(comparison_loaded())
    if(!identical(primary$pdb_accession,handoff$pdb_a) ||
       !identical(secondary$name,handoff$pdb_b) ||
       isTRUE(compare_swapped())) return()
    chain_a <- req(input$compareChainA)
    chain_b <- req(input$compareChainB)
    if(!identical(chain_a,as.character(handoff$residue$chain_a[[1L]])) ||
       !identical(chain_b,as.character(handoff$residue$chain_b[[1L]])))
      return()
    aligned <- req(comparison_data())
    if(!nrow(aligned)) return()
    matched <- ram_atlas_comparison_pair_index(aligned,handoff$residue)
    atlas_pair_handoff(NULL)
    if(is.na(matched)) {
      showNotification(sprintf(paste0("UniProt position %d could not be ",
        "matched to both exact PDB residues in the sequence alignment. ",
        "Inspect the alignment before drawing conclusions."),
        handoff$position),type="warning",duration=18)
      return()
    }
    choose_comparison(matched)
  },priority=-2)

  observeEvent(list(input$compareChainA, input$compareChainB,
                    input$compareModel, comparison_loaded(), compare_swapped()), {
    selected_comparison(NULL)
  }, ignoreInit=TRUE)
  observeEvent(input$ramComparePlotPick, {
    choose_comparison(input$ramComparePlotPick)
  }, ignoreInit=TRUE)
  observeEvent(input$ramCompareTrackPick, {
    data <- isolate(comparison_data())
    row_id <- suppressWarnings(as.integer(input$ramCompareTrackPick))
    index <- match(row_id, data$row_id)
    if (length(index) == 1L && !is.na(index)) choose_comparison(index)
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
    if (is.null(item) || is.null(isolate(comparison_loaded()))) return()
    main_side <- if (isTRUE(isolate(compare_swapped()))) "b" else "a"
    main_chain <- if (identical(main_side,"a"))
      isolate(input$compareChainA) else isolate(input$compareChainB)
    if (is.null(main_chain) || !identical(item$chain, main_chain)) return()
    data <- isolate(comparison_data())
    index <- ram_comparison_find(data, main_side,
               item$chain, item$resi, item$insertion_code)
    if (!is.na(index) &&
        !identical(isolate(selected_comparison()), data$row_id[[index]]))
      selected_comparison(data$row_id[[index]])
  }, ignoreNULL=TRUE)
  filtered_comparison <- reactive({
    result <- comparison_data()
    if (!nrow(result)) return(result)
    criterion <- input$compareFilter
    if (identical(criterion, "changed"))
      result <- result[result$class_changed, , drop=FALSE]
    else if (identical(criterion, "standard_changed") &&
             "rama8000_changed" %in% names(result))
      result <- result[result$rama8000_changed, , drop=FALSE]
    else if (identical(criterion, "basin_changed") &&
             "basin_changed" %in% names(result))
      result <- result[result$basin_changed, , drop=FALSE]
    else if (identical(criterion, "standard_outlier") &&
             all(c("rama8000_region_a","rama8000_region_b") %in% names(result)))
      result <- result[
        result$rama8000_region_a == "Outlier" |
        result$rama8000_region_b == "Outlier", , drop=FALSE]
    else if (identical(criterion, "shift_large"))
      result <- result[is.finite(result$angular_displacement) &
                       result$angular_displacement >= 30, , drop=FALSE]
    else if (identical(criterion, "large"))
      result <- result[(!is.na(result$delta_phi) & abs(result$delta_phi)>=30) |
                       (!is.na(result$delta_psi) & abs(result$delta_psi)>=30),
                       , drop=FALSE]
    else if (identical(criterion, "confidence_large")) {
      if ("delta_plddt" %in% names(result))
        result <- result[is.finite(result$delta_plddt) &
                         abs(result$delta_plddt)>=20,,drop=FALSE]
      else result <- result[0,,drop=FALSE]
    }
    else if (identical(criterion, "gaps"))
      result <- result[result$alignment %in% c("Insertion","Deletion"),
                       , drop=FALSE]
    result
  })
  output$atlasCanonicalNotice <- renderUI({
    ctx <- atlas_alignment_context()
    if(is.null(ctx)) return(NULL)
    result <- comparison_data()
    meta <- attr(result,"atlas_canonical")
    if(!is.null(meta)) {
      tags$div(class="ram-panel",
        tags$strong("Verified UniProt comparison"),
        tags$p(class="ram-field-hint",
          sprintf(paste0("%d common observed canonical positions in %s. ",
            "Verified mapping covers %d/%d primary and %d/%d comparison residues. ",
            "Other residues are excluded, not called insertions or deletions."),
            meta$matched,ctx$accession,meta$mapped_a,meta$total_a,
            meta$mapped_b,meta$total_b)),
        if(meta$matched==0L)
          tags$p(class="ram-confidence-warning",
            "No verified common pairs. RamplotR will not substitute a guessed alignment."),
        actionButton("atlasUseSequenceAlignment","Switch to sequence alignment",
          class="btn-default btn-sm"))
    } else {
      tags$div(class="ram-field-hint",
        tags$p("Atlas matching applies only to the originally selected chains and first structural models. The current view uses sequence alignment."),
        actionButton("atlasUseSequenceAlignment","Discard Atlas mapping",
          class="btn-default btn-sm"))
    }
  })
  observeEvent(input$atlasUseSequenceAlignment, {
    atlas_alignment_context(NULL)
    atlas_pair_handoff(NULL)
    selected_comparison(NULL)
  },ignoreInit=TRUE)

  output$compareSummary <- renderUI({
    result <- req(comparison_data())
    if (!nrow(result)) return(tags$p("Select two nonempty protein chains."))
    aligned <- result$alignment %in% c("Match", "Substitution")
    quality <- ram_comparison_alignment_quality(result)
    pct <- function(value) if(is.finite(value))
      sprintf("%.1f%%",100*value) else "n/a"
    standard_changes <- if ("rama8000_changed" %in% names(result))
      sum(result$rama8000_changed, na.rm=TRUE) else 0L
    basin_changes <- if ("basin_changed" %in% names(result))
      sum(result$basin_changed, na.rm=TRUE) else 0L
    standard_outliers <- if (all(c("rama8000_region_a","rama8000_region_b") %in%
                                 names(result)))
      sum(result$rama8000_region_a == "Outlier" |
          result$rama8000_region_b == "Outlier", na.rm=TRUE) else 0L
    shifts <- result$angular_displacement[is.finite(result$angular_displacement)]
    confidence_pairs <- if ("delta_plddt" %in% names(result))
      sum(is.finite(result$delta_plddt)) else 0L
    confidence_large <- if ("delta_plddt" %in% names(result))
      sum(is.finite(result$delta_plddt) & abs(result$delta_plddt)>=20) else 0L
    weak_alignment <- (is.finite(quality$identity) && quality$identity < 0.50) ||
      (is.finite(quality$coverage_a) && quality$coverage_a < 0.70) ||
      (is.finite(quality$coverage_b) && quality$coverage_b < 0.70)

    metric <- function(value,label)
      tags$span(class="ram-compare-summary-metric",
        tags$strong(value),tags$span(label))
    group <- function(label,...)
      tags$section(class="ram-compare-summary-group",
        tags$h4(label),
        tags$div(class="ram-compare-summary-values",...))

    tags$div(
      class="ram-compare-summary-block",
      tags$div(class="ram-compare-summary-groups",
        group("Alignment",
          metric(quality$aligned,
            if(!is.null(attr(result,"atlas_canonical"))) "exact UniProt pairs"
            else "aligned residues"),
          metric(pct(quality$identity),
            if(!is.null(attr(result,"atlas_canonical")))
              "identity among paired residues" else "sequence identity"),
          if(is.null(attr(result,"atlas_canonical")))
            metric(pct(quality$coverage_a),"primary coverage"),
          if(is.null(attr(result,"atlas_canonical")))
            metric(pct(quality$coverage_b),"comparison coverage"),
          if(is.null(attr(result,"atlas_canonical")))
            metric(sum(!aligned),"insertions / deletions")
        ),
        group("Backbone",
          metric(sum(shifts>=30),"pairs with ≥30° combined shift"),
          metric(basin_changes,"broad backbone-state changes"),
          metric(if(length(shifts)) sprintf("%.1f°",max(shifts)) else "n/a",
                 "largest combined shift")
        ),
        group("Validation",
          metric(sum(result$class_changed),"RamplotR region changes"),
          metric(standard_changes,"Rama8000 category changes"),
          metric(standard_outliers,"pairs with a Rama8000 outlier")
        ),
        if(confidence_pairs>0L)
          group("Prediction confidence",
            metric(confidence_pairs,"pairs with pLDDT on both sides"),
            metric(confidence_large,"pairs with |ΔpLDDT| ≥20")
          )
      ),
      tags$p(class="ram-compare-summary-note",
        "Angular differences wrap across the -180° / +180° boundary."),
      if(weak_alignment)
        tags$p(class="ram-compare-alignment-warning",
          "Alignment identity or coverage is limited. Interpret local conformational shifts cautiously and inspect the aligned sequence context.")
    )
  })

  output$compareChangeTrack <- renderUI({
    data <- comparison_data()
    if (!nrow(data) || !"angular_displacement" %in% names(data)) return(NULL)
    finite <- which(is.finite(data$angular_displacement) &
                    data$alignment %in% c("Match","Substitution"))
    if (!length(finite)) return(NULL)
    selected <- selected_comparison()
    band_class <- function(value)
      paste0("ram-change-",tolower(gsub(" ","-",value,fixed=TRUE)))
    residue_label <- function(side, i) {
      chain <- data[[paste0("chain_",side)]][[i]]
      resi <- data[[paste0("residue_",side)]][[i]]
      ins <- data[[paste0("insertion_",side)]][[i]]
      aa <- data[[paste0("amino_",side)]][[i]]
      paste0(aa," ",chain,resi,ifelse(is.na(ins),"",ins))
    }
    cells <- lapply(seq_along(finite),function(k) {
      i <- finite[[k]]
      shift <- data$angular_displacement[[i]]
      number <- data$residue_a[[i]]
      insertion <- data$insertion_a[[i]]
      show_number <- k==1L || k==length(finite) ||
        (!is.na(number) && number %% 10L == 0L)
      tags$div(class="ram-change-slot",
        tags$span(class="ram-change-position",
          if(show_number) paste0(number,ifelse(is.na(insertion),"",insertion))
          else "\u00a0",
          "aria-hidden"="true"),
        tags$button(
          type="button",
          class=paste("ram-change-cell","ram-change-pick",
            band_class(data$shift_band[[i]]),
            if ("delta_plddt" %in% names(data) &&
                is.finite(data$delta_plddt[[i]]) &&
                abs(data$delta_plddt[[i]])>=20)
              "has-confidence-shift" else "",
            if (!is.null(selected) && identical(data$row_id[[i]],selected))
              "is-selected" else ""),
          "data-row-id"=data$row_id[[i]],
          title=sprintf("%s ↔ %s · Δφ %.1f° · Δψ %.1f° · combined %.1f°",
            residue_label("a",i),residue_label("b",i),
            data$delta_phi[[i]],data$delta_psi[[i]],shift),
          "aria-label"=sprintf(
            "Inspect aligned residue pair with %.1f degree backbone shift",shift)
        )
      )
    })
    ranked <- finite[order(data$angular_displacement[finite],decreasing=TRUE)]
    ranked <- head(ranked,5L)
    tags$section(class="ram-change-explorer",
      tags$div(class="ram-change-head",
        tags$div(
          tags$h3("Conformational change explorer"),
          tags$p("Each cell is one aligned residue. Colour ranks the combined wrapped φ/ψ displacement; it is a navigation measure, not a significance score.")
        ),
        tags$div(class="ram-change-legend",
          tags$span(class="ram-change-small","<15°"),
          tags$span(class="ram-change-moderate","15–30°"),
          tags$span(class="ram-change-large","30–60°"),
          tags$span(class="ram-change-very-large","≥60°"),
          if ("delta_plddt" %in% names(data) &&
              any(is.finite(data$delta_plddt)))
            tags$span(class="ram-change-confidence-key",
              "outline = |ΔpLDDT| ≥20")
        )
      ),
      tags$div(class="ram-change-track",role="group",
        "aria-label"="Aligned residue conformational-change track",cells),
      tags$div(class="ram-change-top",
        tags$strong("Largest local shifts"),
        lapply(ranked,function(i) tags$button(
          type="button",class="ram-change-top-item ram-change-pick",
          "data-row-id"=data$row_id[[i]],
          sprintf("%s ↔ %s · %.1f°",
            residue_label("a",i),residue_label("b",i),
            data$angular_displacement[[i]])
        ))
      )
    )
  })
  outputOptions(output,"compareChangeTrack",suspendWhenHidden=FALSE)

  output$comparison <- DT::renderDT({
    result <- filtered_comparison()
    required <- c("chain_a","residue_a","insertion_a","amino_a",
      "chain_b","residue_b","insertion_b","amino_b",
      "delta_phi","delta_psi","angular_displacement","shift_band",
      "class_changed","rama8000_region_a","rama8000_region_b",
      "rama8000_changed","alignment")
    if (!all(required %in% names(result)))
      return(DT::datatable(data.frame()))
    show_conf_a <- "plddt_a" %in% names(result) &&
      any(is.finite(result$plddt_a))
    show_conf_b <- "plddt_b" %in% names(result) &&
      any(is.finite(result$plddt_b))
    show_delta <- "delta_plddt" %in% names(result) &&
      any(is.finite(result$delta_plddt))
    fields <- c("chain_a","residue_a","insertion_a","amino_a",
      "chain_b","residue_b","insertion_b","amino_b",
      "delta_phi","delta_psi","angular_displacement","shift_band")
    if(all(c("basin_a","basin_b","basin_changed") %in% names(result)))
      fields <- c(fields,"basin_a","basin_b","basin_changed")
    if(show_conf_a) fields <- c(fields,"plddt_a","confidence_a")
    if(show_conf_b) fields <- c(fields,"plddt_b","confidence_b")
    if(show_delta) fields <- c(fields,"delta_plddt")
    fields <- c(fields,"class_changed","rama8000_region_a","rama8000_region_b",
      "rama8000_changed","alignment")
    if("uniprot_resi" %in% names(result)) fields <- c("uniprot_resi",fields)
    shown <- result[,fields,drop=FALSE]
    shown$pos_a <- ifelse(is.na(shown$residue_a),"—",
      paste0(shown$residue_a,shown$insertion_a))
    shown$pos_b <- ifelse(is.na(shown$residue_b),"—",
      paste0(shown$residue_b,shown$insertion_b))
    shown$delta_phi <- round(shown$delta_phi,1)
    shown$delta_psi <- round(shown$delta_psi,1)
    shown$angular_displacement <- round(shown$angular_displacement,1)
    if(show_conf_a) shown$plddt_a <- round(shown$plddt_a,1)
    if(show_conf_b) shown$plddt_b <- round(shown$plddt_b,1)
    if(show_delta) shown$delta_plddt <- round(shown$delta_plddt,1)
    shown$class_changed <- ifelse(shown$class_changed,"Yes","No")
    shown$rama8000_changed <- ifelse(shown$rama8000_changed,"Yes","No")
    display <- c("chain_a","pos_a","amino_a",
      "chain_b","pos_b","amino_b","delta_phi","delta_psi",
      "angular_displacement","shift_band")
    labels <- c("Chain A","Pos A","AA A","Chain B","Pos B","AA B",
      "Δφ (°)","Δψ (°)","Backbone shift (°)","Shift band")
    if("uniprot_resi" %in% names(shown)) {
      display <- c("uniprot_resi",display)
      labels <- c("UniProt position",labels)
    }
    if(all(c("basin_a","basin_b","basin_changed") %in% names(shown))) {
      shown$basin_changed <- ifelse(shown$basin_changed,"Yes","No")
      display <- c(display,"basin_a","basin_b","basin_changed")
      labels <- c(labels,"State A","State B","State changed")
    }
    if(show_conf_a) {
      display <- c(display,"plddt_a","confidence_a")
      labels <- c(labels,"pLDDT A","Confidence A")
    }
    if(show_conf_b) {
      display <- c(display,"plddt_b","confidence_b")
      labels <- c(labels,"pLDDT B","Confidence B")
    }
    if(show_delta) {
      display <- c(display,"delta_plddt")
      labels <- c(labels,"ΔpLDDT")
    }
    display <- c(display,"class_changed","rama8000_region_a",
      "rama8000_region_b","rama8000_changed","alignment")
    labels <- c(labels,"RamplotR changed","Rama8000 A","Rama8000 B",
      "Rama8000 changed","Alignment")
    shown <- shown[,display,drop=FALSE]
    DT::datatable(shown,rownames=FALSE,colnames=labels,selection="single",
      options=list(pageLength=15,scrollX=show_conf_a||show_conf_b||show_delta ||
                      "uniprot_resi" %in% names(shown),
                   autoWidth=FALSE,dom="ftip"),
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
  group_comparison_matches <- reactive({
    value <- group_comparison_results()
    if (is.null(value)) return(NULL)
    if(identical(value$source,"atlas")) {
      if(!identical(input$groupInputMode,"atlas")) return(NULL)
      geometry <- atlas_geometry_result()
      picked <- atlas_group_transfer()
      if(is.null(geometry) || !is.null(geometry$error) ||
         is.null(picked) ||
         !identical(value$selection,picked) ||
         !identical(value$result$accession,geometry$accession))
        return(NULL)
      return(value$result)
    }
    if(identical(input$groupInputMode,"atlas")) return(NULL)
    structure <- req(loaded())
    if (!identical(value$key,structure$key) ||
        !identical(value$mode,input$validationMode) ||
        !identical(value$reference,input$bgtype) ||
        !identical(value$background,input$background)) return(NULL)
    value$result
  })

  observeEvent(loaded(), {
    structure <- loaded()
    if (is.null(structure)) return()
    updateSelectInput(session,"groupReferenceChain",
      choices=structure$chains,selected=structure$chains[[1L]])
    group_comparison_results(NULL)
  },ignoreInit=TRUE)

  observeEvent(list(input$validationMode,input$bgtype,input$background), {
    if (!is.null(group_comparison_results()))
      group_comparison_results(NULL)
  },ignoreInit=TRUE)

  observeEvent(input$runGroupComparison, {
    if(identical(input$groupInputMode,"atlas")) {
      picked <- isolate(atlas_group_transfer())
      geometry <- isolate(atlas_geometry_result())
      if(is.null(picked) || is.null(geometry) ||
         !is.null(geometry$error)) {
        showNotification("Select both verified experimental groups in Atlas first.",
          type="warning",duration=12)
        return()
      }
      result <- tryCatch(withProgress(
        message="Comparing exact UniProt-mapped Atlas groups",{
          ram_atlas_group_prepare(isolate(atlas_exact_results()),geometry,
            picked$group_a,picked$group_b,
            input$groupALabel,input$groupBLabel)
        }),error=function(e)e)
      if(inherits(result,"error")) {
        showNotification(conditionMessage(result),type="error",duration=16)
        return()
      }
      group_comparison_results(list(
        source="atlas",selection=picked,result=result))
      showNotification(sprintf(
        "Compared %d + %d verified structures over %d canonical positions.",
        result$n_a,result$n_b,result$core_positions),
        type="message",duration=10)
      return()
    }
    structure <- req(loaded())
    ref_chain <- req(input$groupReferenceChain)
    label_a <- trimws(input$groupALabel)
    label_b <- trimws(input$groupBLabel)
    if (!nzchar(label_a)) label_a <- "Group A"
    if (!nzchar(label_b)) label_b <- "Group B"
    if(identical(label_a,label_b)) {
      showNotification("Give Group A and Group B distinct labels.",
        type="warning",duration=12)
      return()
    }
    min_identity <- suppressWarnings(as.numeric(input$groupMinIdentity)/100)
    min_coverage <- suppressWarnings(as.numeric(input$groupMinCoverage)/100)
    if (!is.finite(min_identity)) min_identity <- 0.70
    if (!is.finite(min_coverage)) min_coverage <- 0.70

    reference <- classified()
    reference <- reference[reference$chain==ref_chain,,drop=FALSE]
    if (nrow(reference)<5L) {
      showNotification("The selected reference chain is too short.",
                       type="error",duration=10)
      return()
    }

    classify_uploaded <- function(upload) {
      if (is.null(upload) || !nrow(upload))
        return(list(tables=list(),labels=character()))
      tables <- vector("list",nrow(upload))
      labels <- character(nrow(upload))
      for(i in seq_len(nrow(upload))) {
        parsed <- ram_load_structure(
          path=upload$datapath[[i]],original_name=upload$name[[i]])
        torsions <- ram_extract_torsions(ram_model_at(parsed,1L))
        classified_upload <- ram_classify_torsions(
          torsions,
          reference_dir=file.path("static",input$bgtype),
          selected_reference=plot_reference(),
          mode=input$validationMode,
          threshold_fn=ram_density_thresholds
        )
        tables[[i]] <- ram_rama8000_classify(
          classified_upload,file.path("static","rama8000"))
        labels[[i]] <- tools::file_path_sans_ext(
          basename(upload$name[[i]]))
      }
      list(tables=tables,labels=labels)
    }

    withProgress(message="Comparing structure groups",value=0.05,{
      group_a <- list(); labels_a <- character()
      if (isTRUE(input$groupIncludeLoadedA)) {
        group_a <- list(classified())
        labels_a <- structure$name
      }
      upload_a <- tryCatch(classify_uploaded(input$groupAFiles),
        error=function(e) e)
      if (inherits(upload_a,"error")) {
        showNotification(conditionMessage(upload_a),type="error",duration=14)
        return()
      }
      if (length(upload_a$tables)) {
        group_a <- c(group_a,upload_a$tables)
        labels_a <- c(labels_a,upload_a$labels)
      }
      incProgress(0.25,detail=paste("Loaded",length(group_a),label_a,"structures"))

      upload_b <- tryCatch(classify_uploaded(input$groupBFiles),
        error=function(e) e)
      if (inherits(upload_b,"error")) {
        showNotification(conditionMessage(upload_b),type="error",duration=14)
        return()
      }
      group_b <- upload_b$tables
      labels_b <- upload_b$labels
      if (!length(group_a) || !length(group_b)) {
        showNotification(
          "Both groups need at least one structure. Include the loaded structure or upload files for Group A, and upload at least one Group B structure.",
          type="warning",duration=14)
        return()
      }
      incProgress(0.20,detail=paste("Loaded",length(group_b),label_b,"structures"))

      prepared_a <- tryCatch(
        ram_prepare_structure_group(reference,group_a,labels_a,
          min_identity,min_coverage),error=function(e)e)
      if (inherits(prepared_a,"error")) {
        showNotification(paste(label_a,conditionMessage(prepared_a),sep=": "),
                         type="error",duration=16)
        return()
      }
      prepared_b <- tryCatch(
        ram_prepare_structure_group(reference,group_b,labels_b,
          min_identity,min_coverage),error=function(e)e)
      if (inherits(prepared_b,"error")) {
        showNotification(paste(label_b,conditionMessage(prepared_b),sep=": "),
                         type="error",duration=16)
        return()
      }
      incProgress(0.25,detail="Calculating circular group summaries")

      comparison <- ram_group_conformation_compare(
        reference,prepared_a$models,prepared_b$models,label_a,label_b)
      members <- rbind(
        transform(prepared_a$model_summary,group=label_a),
        transform(prepared_b$model_summary,group=label_b)
      )
      fingerprint <- ram_group_fingerprint(
        prepared_a$models,prepared_b$models,label_a,label_b)
      result <- list(
        comparison=comparison,fingerprint=fingerprint,members=members,
        label_a=label_a,label_b=label_b,
        n_a=length(prepared_a$models),n_b=length(prepared_b$models),
        reference_chain=ref_chain
      )
      group_comparison_results(list(
        key=structure$key,mode=input$validationMode,
        reference=input$bgtype,background=input$background,
        result=result
      ))
      incProgress(0.25,detail="Preparing residue-level comparison")
    })
  },ignoreInit=TRUE)

  output$groupComparisonSummary <- renderUI({
    result <- group_comparison_matches()
    if (is.null(result)) return(tags$p(class="ram-field-hint",
      "Choose two structure sets and run the group analysis."))
    data <- result$comparison
    finite <- is.finite(data$angular_displacement)
    tags$div(
      if(identical(result$source,"atlas"))
        tags$div(class="ram-field-hint",
          tags$strong("Verified Atlas group comparison · UniProt ",
            result$accession),
          tags$p(sprintf(paste0(
            "%d exact shared observed canonical positions; ",
            "%d positions lack agreement on experimental monomer identity. ",
            "All structures were reused from the verified Atlas cohort."),
            result$core_positions,result$unknown_chemistry_positions)),
          if(isTRUE(result$known_chemistry_differences))
            tags$p(class="ram-confidence-warning",
              "Observed residue chemistry differs between selected structures. This is a potentially confounded structural comparison."),
          tags$p("Rama8000/native RamplotR classifications are not inferred from the Atlas backbone-only cache; local circular φ/ψ and backbone-state summaries remain available. No geometric cluster is interpreted as a validated biological state.")),
      tags$div(class="ram-confidence-metrics",
        tags$span(class="ram-confidence-metric",
          sprintf("%s: %d structures",result$label_a,result$n_a)),
        tags$span(class="ram-confidence-metric",
          sprintf("%s: %d structures",result$label_b,result$n_b)),
        tags$span(class="ram-confidence-metric",
          sprintf("%d residues compared",sum(finite))),
        tags$span(class="ram-confidence-metric",
          sprintf("%d residues with ≥30° mean shift",
            sum(finite & data$angular_displacement>=30))),
        tags$span(class="ram-confidence-metric",
          sprintf("%d high-support shifts",
            sum(data$high_support_shift,na.rm=TRUE))),
        tags$span(class="ram-confidence-metric",
          sprintf("%d low-dispersion shifts",
            sum(data$consistent_shift,na.rm=TRUE))),
        tags$span(class="ram-confidence-metric",
          sprintf("%d sparse-coverage residues",
            sum(data$evidence_profile=="Sparse coverage",na.rm=TRUE))),
        tags$span(class="ram-confidence-metric",
          sprintf("%d Rama8000 mode changes",
            sum(data$rama8000_mode_changed,na.rm=TRUE))),
        tags$span(class="ram-confidence-metric",
          sprintf("%d backbone-state changes",
            sum(data$basin_mode_changed,na.rm=TRUE)))
      ),
      tags$p(class="ram-confidence-explainer",
        "Between-group displacement compares circular mean φ/ψ values from the same complete-pair observations in each group. A high-support shift also requires at least two paired observations per group, ≥30° displacement, ≤15° within-group circular SD and ≥75% paired coverage in both groups. These are navigation criteria, not significance tests.")
    )
  })

  output$groupComparisonTrack <- renderUI({
    result <- group_comparison_matches()
    if (is.null(result)) return(NULL)
    data <- result$comparison
    finite <- which(is.finite(data$angular_displacement))
    if (!length(finite)) return(NULL)
    band_class <- function(value)
      paste0("ram-change-",tolower(gsub(" ","-",value,fixed=TRUE)))
    cells <- lapply(seq_along(finite),function(k) {
      i <- finite[[k]]
      number <- data$resi[[i]]
      insertion <- data$insertion_code[[i]]
      if (is.na(insertion)) insertion <- ""
      show_number <- k==1L || k==length(finite) ||
        (!is.na(number) && number %% 10L==0L)
      classes <- c("ram-change-cell","ram-group-cell","ram-group-pick",
                   band_class(data$shift_band[[i]]))
      if (isTRUE(data$consistent_shift[[i]]))
        classes <- c(classes,"is-consistent")
      if (isTRUE(data$high_support_shift[[i]]))
        classes <- c(classes,"is-high-support")
      if (identical(data$evidence_profile[[i]],"Sparse coverage"))
        classes <- c(classes,"is-sparse")
      if (isTRUE(data$rama8000_mode_changed[[i]]))
        classes <- c(classes,"has-standard-change")
      tags$div(class="ram-change-slot",
        tags$span(class="ram-change-position",
          if(show_number) paste0(number,insertion)
          else "\u00a0","aria-hidden"="true"),
        tags$button(type="button",class=paste(classes,collapse=" "),
          "data-chain"=data$chain[[i]],
          "data-resi"=data$resi[[i]],
          "data-insertion"=insertion,
          title=sprintf("%s %s:%s · %s → %s · Δφ %.1f° · Δψ %.1f° · shift %.1f° · max within-group SD %s · coverage %.0f%%/%.0f%% · %s",
            data$resn[[i]],data$chain[[i]],data$resi[[i]],
            result$label_a,result$label_b,
            data$delta_phi[[i]],data$delta_psi[[i]],
            data$angular_displacement[[i]],
            if(is.finite(data$max_within_group_sd[[i]]))
              sprintf("%.1f°",data$max_within_group_sd[[i]]) else "n/a",
            100*data$a_coverage[[i]],100*data$b_coverage[[i]],
            data$evidence_profile[[i]])
        )
      )
    })
    tags$section(class="ram-change-explorer ram-group-change-explorer",
      tags$div(class="ram-change-head",
        tags$div(tags$h3("Between-group backbone shift"),
          tags$p("Each cell is one reference-chain residue. Colour shows displacement between group circular means; a dark double outline marks high-support low-dispersion shifts, while faded cells indicate sparse coverage.")),
        tags$div(class="ram-change-legend",
          tags$span(class="ram-change-small","<15°"),
          tags$span(class="ram-change-moderate","15–30°"),
          tags$span(class="ram-change-large","30–60°"),
          tags$span(class="ram-change-very-large","≥60°"))
      ),
      tags$div(class="ram-change-track",role="group",
        "aria-label"="Between-group conformational-change track",cells),
      tags$p(class="ram-field-hint",
        "A secondary outline marks residues whose modal Rama8000 category differs between groups. Exact coverage, within-group SD and evidence profile remain available in the table.")
    )
  })

  output$groupComparisonSelectionInfo <- renderUI({
    result <- group_comparison_matches()
    selection <- selected_residue()
    if (is.null(result) || is.null(selection)) return(NULL)
    data <- result$comparison
    ins <- if(is.null(selection$insertion_code)) "" else
      as.character(selection$insertion_code)
    row <- data[
      as.character(data$chain)==as.character(selection$chain) &
      as.integer(data$resi)==as.integer(selection$resi) &
      ifelse(is.na(data$insertion_code),"",as.character(data$insertion_code))==ins,
      ,drop=FALSE]
    if(nrow(row)!=1L || !is.finite(row$angular_displacement[[1L]])) return(NULL)
    fmt_angle <- function(value)
      if(is.finite(value)) sprintf("%.1f°",value) else "n/a"
    fmt_pct <- function(value)
      if(is.finite(value)) sprintf("%.0f%%",100*value) else "n/a"
    tags$section(class="ram-group-evidence-card",
      tags$div(class="ram-group-evidence-head",
        tags$div(
          tags$strong(sprintf("%s %s:%s%s",row$resn[[1L]],row$chain[[1L]],
            row$resi[[1L]],
            ifelse(is.na(row$insertion_code[[1L]]),"",
              row$insertion_code[[1L]]))),
          tags$span(row$evidence_profile[[1L]])
        ),
        if(isTRUE(row$high_support_shift[[1L]]))
          tags$span(class="ram-group-evidence-badge","High-support shift")
      ),
      tags$div(class="ram-group-evidence-grid",
        tags$div(tags$small("Between-group shift"),
          tags$strong(fmt_angle(row$angular_displacement[[1L]])),
          tags$span(sprintf("Δφ %s · Δψ %s",
            fmt_angle(row$delta_phi[[1L]]),fmt_angle(row$delta_psi[[1L]])))),
        tags$div(tags$small("Within-group dispersion"),
          tags$strong(fmt_angle(row$max_within_group_sd[[1L]])),
          tags$span("maximum circular SD across both groups")),
        tags$div(tags$small("Residue coverage"),
          tags$strong(sprintf("%s / %s",
            fmt_pct(row$a_coverage[[1L]]),fmt_pct(row$b_coverage[[1L]]))),
          tags$span(sprintf("%s / %s",result$label_a,result$label_b))),
        tags$div(tags$small("Backbone state"),
          tags$strong(sprintf("%s → %s",
            ifelse(is.na(row$a_basin_mode[[1L]]),"n/a",
              row$a_basin_mode[[1L]]),
            ifelse(is.na(row$b_basin_mode[[1L]]),"n/a",
              row$b_basin_mode[[1L]]))),
          tags$span(if(isTRUE(row$basin_mode_changed[[1L]]))
            "modal state differs between groups"
            else "modal state retained")),
        tags$div(tags$small("Rama8000"),
          tags$strong(sprintf("%s → %s",
            ifelse(is.na(row$a_rama8000_mode[[1L]]),"n/a",
              row$a_rama8000_mode[[1L]]),
            ifelse(is.na(row$b_rama8000_mode[[1L]]),"n/a",
              row$b_rama8000_mode[[1L]]))),
          tags$span(if(isTRUE(row$rama8000_mode_changed[[1L]]))
            "modal category differs between groups"
            else "modal category retained"))
      ),
      tags$p(class="ram-field-hint",
        "This card explains why the residue is prioritised. It is descriptive evidence, not a statistical significance test.")
    )
  })

  group_fingerprint_selected <- reactive({
    result <- group_comparison_matches()
    if(is.null(result) || is.null(result$fingerprint) ||
       !nrow(result$comparison)) return(NULL)
    comparison <- result$comparison
    selection <- selected_residue()
    index <- integer()
    if(!is.null(selection) && length(selection$resi)==1L) {
      ins <- if(is.null(selection$insertion_code) ||
                is.na(selection$insertion_code)) "" else
        as.character(selection$insertion_code)
      index <- which(as.character(comparison$chain)==as.character(selection$chain) &
        comparison$resi==as.integer(selection$resi) &
        ifelse(is.na(comparison$insertion_code),"",
               as.character(comparison$insertion_code))==ins)
    }
    if(length(index)!=1L) {
      index <- which(is.finite(comparison$angular_displacement))
      if(!length(index)) index <- seq_len(nrow(comparison))
    }
    row <- comparison[index[[1L]],,drop=FALSE]
    selected <- ram_group_fingerprint_at(result$fingerprint,
      row$chain[[1L]],row$resi[[1L]],
      if(is.na(row$insertion_code[[1L]])) "" else row$insertion_code[[1L]])
    list(row=row,records=selected)
  })
  output$groupFingerprintPanel <- renderUI({
    result <- group_comparison_matches()
    selected <- group_fingerprint_selected()
    if(is.null(result) || is.null(selected) || !nrow(selected$records))
      return(NULL)
    row <- selected$row
    group_names <- c(result$label_a,result$label_b)
    summary <- ram_group_fingerprint_summary(selected$records,
      group_names,c(result$n_a,result$n_b))
    info <- vapply(seq_len(2L),function(i) {
      state <- if(is.na(summary$modal_state[[i]])) "no complete pair"
        else sprintf("%s (%.0f%% of classified pairs)",
          summary$modal_state[[i]],100*summary$consensus[[i]])
      sprintf("%s: %d/%d complete pairs; %s",group_names[[i]],
        summary$complete_pairs[[i]],summary$members[[i]],state)
    },character(1L))
    tags$section(class="ram-panel ram-group-fingerprint-panel",
      tags$h4(sprintf("Local conformational fingerprint · %s %s:%d%s",
        row$resn[[1L]],row$chain[[1L]],row$resi[[1L]],
        if(is.na(row$insertion_code[[1L]])) "" else
          row$insertion_code[[1L]])),
      tags$p(class="ram-field-hint",
        "One point per measured structure, using complete paired φ/ψ angles. Crosses mark circular means. Choose another residue in the track or table to update this view."),
      tags$div(class="ram-group-fingerprint-grid",
        plotOutput("groupFingerprintPlot",height="310px"),
        tags$div(
          tags$p(class="ram-field-hint",info[[1L]]),
          tags$p(class="ram-field-hint",info[[2L]]),
          DT::DTOutput("groupFingerprintMembers"))),
      uiOutput("groupFingerprintLigandPanel"),
      tags$p(class="ram-field-hint",
        "These are measured structural observations, not independent biological replicates or conformational-state probabilities. Missing angles remain missing; neither RamplotR density nor Rama8000 categories are inferred from Atlas backbone-only records.")
    )
  })
  atlas_fingerprint_contacts <- reactive({
    result <- group_comparison_matches()
    selected <- group_fingerprint_selected()
    if(is.null(result) || !identical(result$source,"atlas") ||
       is.null(result$ligand_context) || is.null(selected))
      return(NULL)
    ram_atlas_group_ligand_at(result,selected$row$resi[[1L]])
  })
  output$groupFingerprintLigandPanel <- renderUI({
    rows <- atlas_fingerprint_contacts()
    if(is.null(rows)) return(NULL)
    labels <- unique(rows$group)
    stats <- lapply(labels,function(label) {
      sub <- rows[rows$group==label,,drop=FALSE]
      sprintf("%s: %d nearby / %d examined / %d unavailable or limited",
        label,sum(sub$evidence=="Deposited proximity observed"),
        sum(sub$evidence %in% c("Deposited proximity observed",
                               "No deposited proximity reported")),
        sum(sub$evidence %in% c("Evidence unavailable",
                               "Incomplete nearest-only evidence")))
    })
    tags$div(class="ram-group-ligand-evidence",
      tags$h5("Deposited component proximity at this UniProt position"),
      tags$p(class="ram-field-hint",
        paste(unlist(stats),collapse=" · ")),
      tags$p(class="ram-field-hint",
        "4.5 Å heavy-atom proximity from model-1 deposited coordinates. A nearby component does not demonstrate functional binding. An absent or unverified record does not establish an apo state. Multiple entities in the same PDB entry are not independent experiments."),
      DT::DTOutput("groupFingerprintLigandMembers"),
      downloadButton("downloadGroupFingerprintLigand",
        "Export residue contact evidence CSV"))
  })
  output$groupFingerprintLigandMembers <- DT::renderDT({
    rows <- req(atlas_fingerprint_contacts())
    shown <- rows[,c("group","member","complete_backbone_pair",
      "evidence","component_codes","minimum_distance_A","scope"),
      drop=FALSE]
    shown$minimum_distance_A <- round(shown$minimum_distance_A,2L)
    DT::datatable(shown,rownames=FALSE,
      colnames=c("Group","PDB entity","Complete φ/ψ",
        "Deposited evidence","Components","Min distance (Å)","Coverage"),
      options=list(dom="tip",pageLength=8,scrollX=TRUE),
      class="compact stripe")
  },server=FALSE)
  output$downloadGroupFingerprintLigand <- downloadHandler(
    filename=function() "ramplotr_selected_residue_component_evidence.csv",
    content=function(file) utils::write.csv(
      req(atlas_fingerprint_contacts()),file,row.names=FALSE,na=""))

  output$groupFingerprintPlot <- renderPlot({
    selected <- req(group_fingerprint_selected())
    result <- req(group_comparison_matches())
    ram_group_fingerprint_plot(selected$records,
      c(result$label_a,result$label_b),
      contacts=atlas_fingerprint_contacts())
  })
  output$groupFingerprintMembers <- DT::renderDT({
    selected <- req(group_fingerprint_selected())
    records <- selected$records
    req(nrow(records)>0L)
    shown <- records[,c("group","member","amino_acid","phi",
      "psi","paired","backbone_state"),drop=FALSE]
    shown$phi <- round(shown$phi,1)
    shown$psi <- round(shown$psi,1)
    DT::datatable(shown,rownames=FALSE,
      colnames=c("Group","Structure","AA","φ","ψ","Complete pair","Backbone state"),
      options=list(dom="tip",pageLength=8,scrollX=TRUE),
      class="compact stripe")
  },server=FALSE)

  output$groupComparisonRows <- DT::renderDT({
    result <- group_comparison_matches()
    req(result)
    data <- result$comparison
    shown <- data[,c(
      "chain","resi","insertion_code","resn",
      "a_paired_angle_models","b_paired_angle_models",
      "a_phi_mean","b_phi_mean","delta_phi",
      "a_psi_mean","b_psi_mean","delta_psi",
      "angular_displacement","max_within_group_sd",
      "a_coverage","b_coverage","min_rama8000_consistency",
      "evidence_profile","a_basin_mode","b_basin_mode","basin_mode_changed",
      "a_rama8000_mode","b_rama8000_mode",
      "high_support_shift","consistent_shift","rama8000_mode_changed"
    ),drop=FALSE]
    for(field in c("a_phi_mean","b_phi_mean","delta_phi",
                   "a_psi_mean","b_psi_mean","delta_psi",
                   "angular_displacement","max_within_group_sd"))
      shown[[field]] <- round(shown[[field]],1L)
    shown$a_coverage <- round(100*shown$a_coverage,1L)
    shown$b_coverage <- round(100*shown$b_coverage,1L)
    shown$min_rama8000_consistency <- round(100*shown$min_rama8000_consistency,1L)
    shown$basin_mode_changed <- ifelse(shown$basin_mode_changed,"Yes","No")
    shown$high_support_shift <- ifelse(shown$high_support_shift,"Yes","No")
    shown$consistent_shift <- ifelse(shown$consistent_shift,"Yes","No")
    shown$rama8000_mode_changed <- ifelse(shown$rama8000_mode_changed,
                                           "Yes","No")
    DT::datatable(shown,rownames=FALSE,selection="single",
      colnames=c("Chain","Residue","Ins.","AA","Paired A","Paired B",
        "φ A","φ B","Δφ","ψ A","ψ B","Δψ",
        "Mean shift","Max within SD","Coverage A (%)","Coverage B (%)",
        "Min Rama8000 agreement (%)","Evidence profile",
        "State A","State B","State changed",
        "Rama8000 A","Rama8000 B","High support",
        "Consistent shift","Rama8000 changed"),
      options=list(pageLength=12,scrollX=TRUE,autoWidth=FALSE,dom="ftip"),
      class="compact stripe hover")
  },server=FALSE)

  observeEvent(input$groupComparisonRows_rows_selected, {
    result <- req(group_comparison_matches())
    selection <- input$groupComparisonRows_rows_selected
    if (is.null(selection) || !length(selection)) return()
    ix <- suppressWarnings(as.integer(selection[[1L]]))
    if(length(ix)!=1L || is.na(ix) || ix<1L ||
       ix>nrow(result$comparison)) return()
    row <- result$comparison[ix,,drop=FALSE]
    selected_residue(list(chain=as.character(row$chain[[1L]]),
      resi=as.integer(row$resi[[1L]]),
      insertion_code=if(is.na(row$insertion_code[[1L]])) ""
        else as.character(row$insertion_code[[1L]])))
  })

  observeEvent(input$ramGroupComparisonPick, {
    item <- input$ramGroupComparisonPick
    if (!is.list(item) || is.null(item$resi)) return()
    selected_residue(list(
      chain=as.character(item$chain),
      resi=as.integer(item$resi),
      insertion_code=if(is.null(item$insertion_code)) ""
        else as.character(item$insertion_code)
    ))
  },ignoreInit=TRUE)

  output$groupComparisonExports <- renderUI({
    matches <- group_comparison_matches()
    if (is.null(matches)) return(NULL)
    tags$div(class="ram-ensemble-actions ram-group-downloads",
      downloadButton("downloadGroupComparison",
        "Export residue comparison CSV"),
      downloadButton("downloadGroupMembers",
        "Export matched structure/chain CSV"),
      downloadButton("downloadGroupFingerprint",
        "Export per-model fingerprint CSV")
    )
  })

  output$downloadGroupComparison <- downloadHandler(
    filename=function() safe_filename("group-conformation-comparison.csv"),
    content=function(file) utils::write.csv(
      req(group_comparison_matches())$comparison,
      file,row.names=FALSE,na="")
  )
  output$downloadGroupFingerprint <- downloadHandler(
    filename=function() safe_filename("group-conformational-fingerprints.csv"),
    content=function(file) utils::write.csv(
      req(group_comparison_matches())$fingerprint,
      file,row.names=FALSE,na="")
  )
  output$downloadGroupMembers <- downloadHandler(
    filename=function() safe_filename("group-conformation-members.csv"),
    content=function(file) utils::write.csv(
      req(group_comparison_matches())$members,
      file,row.names=FALSE,na="")
  )

  observe({
    result <- comparison_data()
    if (!nrow(result)) return()
    main <- req(loaded()); comparison <- req(comparison_loaded())
    swapped <- isTRUE(compare_swapped())
    comparison_model <- if (is.null(input$compareModel)) 1L
      else as.integer(input$compareModel)
    reference <- plot_reference()
    palette <- active_palette()
    session$sendCustomMessage("ram-comparison", list(
      nameA=if (swapped) comparison$name else main$name,
      nameB=if (swapped) main$name else comparison$name,
      matrix=reference,
      limits=ram_density_thresholds(reference),
      backgroundColors=palette,
      backgroundName=input$background,
      phiA=result$phi_a, psiA=result$psi_a,
      phiB=result$phi_b, psiB=result$psi_b,
      rowIds=result$row_id,
      chainA=result$chain_a, posA=result$residue_a,
      insA=result$insertion_a, aminoA=result$amino_a,
      chainB=result$chain_b, posB=result$residue_b,
      insB=result$insertion_b, aminoB=result$amino_b,
      deltaPhi=result$delta_phi, deltaPsi=result$delta_psi,
      plddtA=if ("plddt_a" %in% names(result)) result$plddt_a
        else rep(NA_real_,nrow(result)),
      plddtB=if ("plddt_b" %in% names(result)) result$plddt_b
        else rep(NA_real_,nrow(result)),
      deltaPlddt=if ("delta_plddt" %in% names(result)) result$delta_plddt
        else rep(NA_real_,nrow(result)),
      confidenceA=if ("confidence_a" %in% names(result)) result$confidence_a
        else rep(NA_character_,nrow(result)),
      confidenceB=if ("confidence_b" %in% names(result)) result$confidence_b
        else rep(NA_character_,nrow(result)),
      alignment=result$alignment
    ))
    session$sendCustomMessage("ram-compare-config", list(
      chainA=input$compareChainA, chainB=input$compareChainB,
      modelA=if (swapped) comparison_model else current_model(),
      modelB=if (swapped) current_model() else comparison_model,
      multipleA=if (swapped) comparison$nmodels > 1L else main$nmodels > 1L,
      multipleB=if (swapped) main$nmodels > 1L else comparison$nmodels > 1L
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

    main <- req(loaded())
    comparison <- req(comparison_loaded())
    swapped <- isTRUE(compare_swapped())
    comparison_model <- if (is.null(input$compareModel)) 1L else
      suppressWarnings(as.integer(input$compareModel))
    if (!is.finite(comparison_model) || comparison_model < 1L ||
        comparison_model > comparison$nmodels) comparison_model <- 1L

    first_structure <- if (swapped) comparison else main
    second_structure <- if (swapped) main else comparison
    first_model <- if (swapped) comparison_model else current_model()
    second_model <- if (swapped) current_model() else comparison_model

    pair_row <- function(side) {
      number <- row[[paste0("residue_",side)]][[1L]]
      if (is.na(number)) return(NULL)
      data.frame(
        chain=as.character(row[[paste0("chain_",side)]][[1L]]),
        resi=as.integer(number),
        insertion_code=as.character(row[[paste0("insertion_",side)]][[1L]]),
        resn=as.character(row[[paste0("amino_",side)]][[1L]]),
        stringsAsFactors=FALSE
      )
    }
    local_context <- function(structure, model, side) {
      residue <- pair_row(side)
      if (is.null(residue)) return(list(available=FALSE,data=NULL,
        reason="Alignment gap"))
      if (structure$nmodels>1L && !identical(as.integer(model),1L))
        return(list(available=FALSE,data=NULL,
          reason="Hetero context is retained conservatively for model 1 only"))
      data <- ram_nearby_hetero_context(
        ram_model_at(structure$pdb,as.integer(model)),residue,
        max_distance=6,max_hits=3L)
      list(available=TRUE,data=data,reason=NULL)
    }
    context_a <- local_context(first_structure,first_model,"a")
    context_b <- local_context(second_structure,second_model,"b")
    context_label <- function(context) {
      if (!isTRUE(context$available))
        return(tags$span(class="ram-compare-context-muted",context$reason))
      data <- context$data
      if (is.null(data) || !nrow(data))
        return(tags$span(class="ram-compare-context-muted",
          "No non-water hetero residue within 6 Å"))
      labels <- vapply(seq_len(nrow(data)),function(i) {
        chain <- as.character(data$chain[[i]])
        insertion <- as.character(data$insertion_code[[i]])
        position <- paste0(if(nzchar(chain)) paste0(chain,":") else "",
          data$resi[[i]],ifelse(is.na(insertion),"",insertion))
        sprintf("%s %s · %.1f Å",data$resn[[i]],position,data$distance[[i]])
      },character(1L))
      tags$span(paste(labels,collapse="; "))
    }

    tags$div(class="ram-compare-selection",
      if("uniprot_resi" %in% names(row))
        tags$p(class="ram-field-hint",
          sprintf("Verified UniProt position %d · exact PDBe SIFTS mapping",
            row$uniprot_resi[[1L]])),
      tags$div(class="ram-compare-selection-pair",
        tags$span(class="ram-compare-primary",
          tags$small("Primary"), tags$strong(label("a")),
          tags$span(paste("φ",angle(row$phi_a[[1L]]),
                          "· ψ",angle(row$psi_a[[1L]]))),
          if ("plddt_a" %in% names(row) && is.finite(row$plddt_a[[1L]]))
            tags$span(sprintf("pLDDT %.1f%s",row$plddt_a[[1L]],
              if ("confidence_a" %in% names(row) &&
                  !is.na(row$confidence_a[[1L]]))
                paste0(" · ",row$confidence_a[[1L]]) else ""))),
        tags$span(class="ram-compare-pair-arrow", "↔", "aria-hidden"="true"),
        tags$span(class="ram-compare-secondary",
          tags$small("Comparison"), tags$strong(label("b")),
          tags$span(paste("φ",angle(row$phi_b[[1L]]),
                          "· ψ",angle(row$psi_b[[1L]]))),
          if ("plddt_b" %in% names(row) && is.finite(row$plddt_b[[1L]]))
            tags$span(sprintf("pLDDT %.1f%s",row$plddt_b[[1L]],
              if ("confidence_b" %in% names(row) &&
                  !is.na(row$confidence_b[[1L]]))
                paste0(" · ",row$confidence_b[[1L]]) else "")))
      ),
      tags$div(class="ram-compare-selection-deltas",
        tags$span(paste("Δφ",angle(row$delta_phi[[1L]]))),
        tags$span(paste("Δψ",angle(row$delta_psi[[1L]]))),
        tags$span(paste("Combined",angle(row$angular_displacement[[1L]]),
                        "·",row$shift_band[[1L]])),
        if ("basin_a" %in% names(row) &&
            !is.na(row$basin_a[[1L]]) && !is.na(row$basin_b[[1L]]))
          tags$span(sprintf("Backbone state %s → %s",
            row$basin_a[[1L]],row$basin_b[[1L]])),
        if ("delta_plddt" %in% names(row) && is.finite(row$delta_plddt[[1L]]))
          tags$span(sprintf("ΔpLDDT %+.1f",row$delta_plddt[[1L]])),
        tags$span(row$alignment[[1L]]),
        if (isTRUE(row$class_changed[[1L]])) tags$span(
          class="ram-compare-change", "RamplotR region changed"),
        if ("rama8000_region_a" %in% names(row) &&
            !is.na(row$rama8000_region_a[[1L]]))
          tags$span(sprintf("Rama8000 %s → %s",
            row$rama8000_region_a[[1L]], row$rama8000_region_b[[1L]])),
        if ("rama8000_changed" %in% names(row) &&
            isTRUE(row$rama8000_changed[[1L]]))
          tags$span(class="ram-compare-change", "Standard category changed"),
        if ("basin_changed" %in% names(row) &&
            isTRUE(row$basin_changed[[1L]]))
          tags$span(class="ram-compare-change", "Backbone state changed")
      ),
      tags$div(class="ram-compare-local-context",
        tags$div(class="ram-compare-context-side",
          tags$strong("Primary local context"),
          context_label(context_a)),
        tags$div(class="ram-compare-context-side",
          tags$strong("Comparison local context"),
          context_label(context_b)),
        tags$p(class="ram-compare-context-note",
          "Nearest heavy-atom distances to non-water hetero residues within 6 Å. Proximity is structural context, not evidence of biochemical binding.")
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
    main <- req(loaded())
    comparison <- req(comparison_loaded())
    swapped <- isTRUE(compare_swapped())
    comparison_model <- if (is.null(input$compareModel)) 1L else
      as.integer(input$compareModel)
    session$sendCustomMessage("ram-comparison-selected",
      list(rowId=id))
    session$sendCustomMessage("ram-compare-pair", list(
      a=pair("a",
        if (swapped) comparison_model else current_model(),
        if (swapped) comparison$nmodels > 1L else main$nmodels > 1L),
      b=pair("b",
        if (swapped) current_model() else comparison_model,
        if (swapped) main$nmodels > 1L else comparison$nmodels > 1L)
    ))
  })
  output$NGLCompare <- NGLVieweR::renderNGLVieweR({
    req(input$showComparison3D,input$compareChainA,input$compareChainB)
    req(comparison_data())
    main <- req(loaded()); comparison <- req(comparison_loaded())
    swapped <- isTRUE(compare_swapped())
    first <- if (swapped) comparison else main
    second <- if (swapped) main else comparison
    comparison_model <- if (is.null(input$compareModel)) 1L else
      as.integer(input$compareModel)
    first_model <- if (swapped) comparison_model else current_model()
    second_model <- if (swapped) current_model() else comparison_model
    model_a <- if (first$nmodels > 1L)
      paste0(" and /",first_model-1L) else ""
    model_b <- if (second$nmodels > 1L)
      paste0(" and /",second_model-1L) else ""
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
    is_prediction <- structure$declared_source %in%
      c("alphafold_db","alphafold2","alphafold3","esmfold","other_prediction")
    if(structure$nmodels<=1L && !is_prediction) return(NULL)

    tags$details(id="ram-ensemble-panel",class="ram-confidence-panel",
      tags$summary(
        tags$span(class="ram-confidence-title","Ensemble analysis"),
        tags$span(class="ram-confidence-subtitle",
          if(is_prediction)
            "Prediction-model agreement · circular φ/ψ variation · Rama8000 and pLDDT"
          else paste(structure$nmodels,
            "structural models · circular φ/ψ variation and region consistency"))
      ),
      tags$div(class="ram-confidence-body",
        if(structure$nmodels>1L) tagList(
          tags$h4("Models stored in this structure"),
          tags$p(class="ram-confidence-explainer",
            "Model variation is matched by chain, residue and insertion code. Circular statistics correctly handle the -180°/180° boundary; models with missing coordinates contribute only observed angles."),
          tags$div(class="ram-ensemble-actions",
            actionButton("calculateEnsemble","Analyse structural models",
                         class="btn-primary btn-sm"),
            downloadButton("downloadEnsemble","Export structural ensemble CSV")
          ),
          uiOutput("ensembleResultSummary"),
          tags$div(class="ram-residue-table",DT::DTOutput("ensembleRows"))
        ),
        if(is_prediction) tagList(
          if(structure$nmodels>1L) tags$hr(),
          tags$div(class="ram-prediction-ensemble-head",
            tags$h4("Prediction ensemble"),
            tags$p(class="ram-confidence-explainer",
              "Compare independently generated prediction models or seeds. RamplotR keeps residue-level backbone variation, pLDDT and standard validation separate from model-level ranking metrics; prediction disagreement is not experimental dynamics.")
          ),
          tags$div(class="ram-prediction-ensemble-controls",
            selectInput("predictionEnsembleSource","Prediction model type",
              choices=c("AlphaFold 2 / ColabFold"="alphafold2",
                        "AlphaFold 3 sample set"="alphafold3",
                        "ESMFold"="esmfold",
                        "Other model with pLDDT in B-factor"="other_prediction"),
              selected=if(structure$declared_source %in%
                c("alphafold3","esmfold","other_prediction"))
                  structure$declared_source else "alphafold2",
              selectize=FALSE),
            fileInput("predictionEnsembleFiles",
              "Prediction model files",
              multiple=TRUE,
              accept=c(".pdb",".ent",".cif",".mmcif",".mcif")),
            conditionalPanel(
              condition="input.predictionEnsembleSource === 'alphafold3'",
              fileInput("predictionEnsembleConfidenceFiles",
                "AF3 full confidences JSON",
                multiple=TRUE,accept=c(".json")),
              fileInput("predictionEnsembleSummaryFiles",
                "AF3 summary confidences JSON (optional)",
                multiple=TRUE,accept=c(".json")),
              tags$p(class="ram-field-hint",
                "AF3 files are paired by their official seed/sample filename stem: *_model.cif ↔ *_confidences.json ↔ optional *_summary_confidences.json. Upload order is ignored.")
            ),
            if(!identical(structure$declared_source,"alphafold3"))
              conditionalPanel(
                condition="input.predictionEnsembleSource !== 'alphafold3'",
                checkboxInput("includeLoadedPrediction",
                  paste("Include currently loaded model:",structure$name),value=TRUE)
              )
            else
              tags$p(class="ram-field-hint",
                "The currently loaded AlphaFold 3 model is not auto-included in another ensemble type. Upload it again with its matching full-confidence JSON when analysing an AF3 sample set."),
            actionButton("calculatePredictionEnsemble",
              "Analyse prediction ensemble",class="btn-primary btn-sm")
          ),
          tags$p(class="ram-field-hint",
            "AF3 requires each sample's matching full confidence JSON. pTM, ipTM and ranking score remain model-level provenance and are not folded into the residue variability measure."),
          uiOutput("predictionEnsembleSummary"),
          tags$details(class="ram-confidence-panel ram-ensemble-model-panel",
            tags$summary("Model-level confidence and provenance"),
            tags$div(class="ram-residue-table",
              DT::DTOutput("predictionEnsembleModelRows"))
          ),
          uiOutput("predictionEnsembleTrack"),
          tags$div(class="ram-residue-table",
            DT::DTOutput("predictionEnsembleRows")),
          tags$div(class="ram-ensemble-actions",
            downloadButton("downloadPredictionEnsemble",
              "Export prediction ensemble CSV"),
            downloadButton("downloadPredictionEnsembleModels",
              "Export model summary CSV"),
            downloadButton("downloadPredictionEnsembleReport",
              "HTML ensemble report")
          )
        )
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
          classifier=function(torsions) {
            classified <- ram_classify_torsions(torsions,
              reference_dir=file.path("static",input$bgtype),
              selected_reference=plot_reference(),mode=input$validationMode,
              threshold_fn=ram_density_thresholds)
            ram_rama8000_classify(classified,file.path("static","rama8000"))
          }),
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

  prediction_ensemble_input_key <- reactive({
    structure <- req(loaded())
    uploaded <- input$predictionEnsembleFiles
    signature <- function(files) if(is.null(files) || !nrow(files)) "" else
      paste(files$name,files$size,files$type,files$datapath,
            sep=":",collapse="|")
    source <- if(is.null(input$predictionEnsembleSource)) "" else
      input$predictionEnsembleSource
    include_loaded <- isTRUE(input$includeLoadedPrediction) &&
      !identical(source,"alphafold3") &&
      !identical(structure$declared_source,"alphafold3")
    paste(
      source,
      include_loaded,
      if(include_loaded) current_model() else "",
      signature(uploaded),
      if(identical(source,"alphafold3"))
        signature(input$predictionEnsembleConfidenceFiles) else "",
      if(identical(source,"alphafold3"))
        signature(input$predictionEnsembleSummaryFiles) else "",
      sep="::"
    )
  })

  prediction_ensemble_matches <- reactive({
    value <- prediction_ensemble_results()
    if(is.null(value)) return(NULL)
    structure <- req(loaded())
    if(!identical(value$key,structure$key) ||
       !identical(value$mode,input$validationMode) ||
       !identical(value$reference,input$bgtype) ||
       !identical(value$background,input$background) ||
       !identical(value$input_key,prediction_ensemble_input_key()))
      return(NULL)
    value$result
  })

  observeEvent(loaded(), {
    prediction_ensemble_results(NULL)
  }, ignoreInit=TRUE)

  observeEvent(input$calculatePredictionEnsemble, {
    structure <- req(loaded())
    source <- req(input$predictionEnsembleSource)
    permitted <- c("alphafold2","alphafold3","esmfold","other_prediction")
    if(!source %in% permitted) return()

    uploaded <- input$predictionEnsembleFiles
    include_loaded <- isTRUE(input$includeLoadedPrediction) &&
      !identical(source,"alphafold3") &&
      !identical(structure$declared_source,"alphafold3")
    source_loaded <- if(identical(structure$declared_source,"alphafold_db"))
      "alphafold2" else structure$declared_source

    if(include_loaded && !identical(source_loaded,source)) {
      showNotification(
        paste0("The loaded model is declared as ",source_loaded,
          " but the ensemble is configured as ",source,
          ". Choose the matching model type or exclude the loaded model."),
        type="error",duration=14)
      return()
    }

    af3_pairs <- NULL
    if(identical(source,"alphafold3")) {
      confidences <- input$predictionEnsembleConfidenceFiles
      summaries <- input$predictionEnsembleSummaryFiles
      af3_pairs <- tryCatch(
        ram_af3_pair_files(
          model_names=if(is.null(uploaded)) character() else uploaded$name,
          model_paths=if(is.null(uploaded)) character() else uploaded$datapath,
          confidence_names=if(is.null(confidences)) character() else confidences$name,
          confidence_paths=if(is.null(confidences)) character() else confidences$datapath,
          summary_names=if(is.null(summaries)) character() else summaries$name,
          summary_paths=if(is.null(summaries)) character() else summaries$datapath
        ),
        error=function(e) {
          showNotification(conditionMessage(e),type="error",duration=16)
          NULL
        }
      )
      if(is.null(af3_pairs)) return()
    }

    file_count <- if(identical(source,"alphafold3")) nrow(af3_pairs) else
      if(is.null(uploaded)) 0L else nrow(uploaded)
    if(file_count + as.integer(include_loaded) < 2L) {
      showNotification(
        "A prediction ensemble needs at least two models. Upload another model or include the loaded prediction.",
        type="warning",duration=12)
      return()
    }

    withProgress(message="Analysing prediction ensemble",value=0.05,{
      pdbs <- list()
      labels <- character()
      hashes <- character()
      confidence_hashes <- character()
      summary_hashes <- character()
      structure_models <- integer()
      input_roles <- character()
      sidecars <- character()
      summary_files <- character()

      if(include_loaded) {
        selected_model <- current_model()
        pdbs[[length(pdbs)+1L]] <- ram_model_at(structure$pdb,selected_model)
        labels <- c(labels,
          if(structure$nmodels>1L)
            sprintf("%s [model %s]",structure$name,selected_model)
          else structure$name)
        hashes <- c(hashes,
          if(is.character(structure$source_id) &&
             length(structure$source_id)==1L &&
             file.exists(structure$source_id))
            unname(tools::md5sum(structure$source_id)) else NA_character_)
        confidence_hashes <- c(confidence_hashes,NA_character_)
        summary_hashes <- c(summary_hashes,NA_character_)
        structure_models <- c(structure_models,selected_model)
        input_roles <- c(input_roles,"loaded")
      }

      if(file_count) {
        for(i in seq_len(file_count)) {
          model_name <- if(identical(source,"alphafold3"))
            af3_pairs$model_name[[i]] else uploaded$name[[i]]
          model_path <- if(identical(source,"alphafold3"))
            af3_pairs$model_path[[i]] else uploaded$datapath[[i]]
          incProgress(0.35/max(1L,file_count),
            detail=paste("Loading",model_name))
          model <- tryCatch(
            ram_load_structure(path=model_path,original_name=model_name),
            error=function(e) e
          )
          if(inherits(model,"error")) {
            showNotification(
              paste(model_name,conditionMessage(model),sep=": "),
              type="error",duration=14)
            return()
          }
          if(ram_model_count(model)!=1L) {
            showNotification(
              paste(model_name,
                "contains multiple structural models. Prediction-ensemble uploads must contain one model per file."),
              type="error",duration=14)
            return()
          }
          pdbs[[length(pdbs)+1L]] <- model
          labels <- c(labels,
            if(identical(source,"alphafold3")) af3_pairs$label[[i]]
            else tools::file_path_sans_ext(basename(model_name)))
          hashes <- c(hashes,unname(tools::md5sum(model_path)))
          structure_models <- c(structure_models,1L)
          input_roles <- c(input_roles,
            if(identical(source,"alphafold3")) "uploaded-af3" else "uploaded")
          if(identical(source,"alphafold3")) {
            confidence_path <- af3_pairs$confidence_path[[i]]
            summary_path <- af3_pairs$summary_path[[i]]
            sidecars <- c(sidecars,confidence_path)
            summary_files <- c(summary_files,summary_path)
            confidence_hashes <- c(confidence_hashes,
              unname(tools::md5sum(confidence_path)))
            summary_hashes <- c(summary_hashes,
              if(nzchar(summary_path)) unname(tools::md5sum(summary_path))
              else NA_character_)
          } else {
            confidence_hashes <- c(confidence_hashes,NA_character_)
            summary_hashes <- c(summary_hashes,NA_character_)
          }
        }
      }

      known_hashes <- hashes[!is.na(hashes) & nzchar(hashes)]
      if(anyDuplicated(known_hashes)) {
        showNotification(
          "The ensemble contains duplicate coordinate files. Remove duplicate seeds/models before analysing agreement.",
          type="error",duration=14)
        return()
      }

      result <- tryCatch(
        ram_prediction_ensemble_analyze(
          pdbs,source=source,labels=labels,max_models=30L,
          sidecars=if(identical(source,"alphafold3")) sidecars else NULL,
          summary_files=if(identical(source,"alphafold3")) summary_files else NULL,
          classifier=function(torsions) {
            classified <- ram_classify_torsions(
              torsions,
              reference_dir=file.path("static",input$bgtype),
              selected_reference=plot_reference(),
              mode=input$validationMode,
              threshold_fn=ram_density_thresholds
            )
            ram_rama8000_classify(
              classified,file.path("static","rama8000"))
          }
        ),
        error=function(e) {
          showNotification(conditionMessage(e),type="error",duration=14)
          NULL
        }
      )
      if(is.null(result)) return()
      n_used <- result$analyzed_models
      result$provenance <- data.frame(
        model=result$labels,
        source=result$source,
        input_role=input_roles[seq_len(n_used)],
        structure_model=structure_models[seq_len(n_used)],
        coordinate_md5=hashes[seq_len(n_used)],
        confidence_md5=confidence_hashes[seq_len(n_used)],
        summary_md5=summary_hashes[seq_len(n_used)],
        stringsAsFactors=FALSE
      )
      prediction_ensemble_results(list(
        key=structure$key,mode=input$validationMode,
        reference=input$bgtype,background=input$background,
        input_key=prediction_ensemble_input_key(),
        result=result
      ))
      incProgress(0.6,detail="Summarising model agreement")
    })
  },ignoreInit=TRUE)

  output$predictionEnsembleSummary <- renderUI({
    result <- prediction_ensemble_matches()
    if(is.null(result)) return(tags$p(class="ram-field-hint",
      "Upload at least two compatible prediction models and run the ensemble analysis."))
    data <- result$summary
    standard_changes <- if("rama8000_changes" %in% names(data))
      sum(data$rama8000_changes,na.rm=TRUE) else 0L
    angular_variable <- sum(
      pmax(data$phi_sd,data$psi_sd,na.rm=TRUE)>=20,na.rm=TRUE)
    confidence_variable <- if("plddt_sd" %in% names(data))
      sum(is.finite(data$plddt_sd) & data$plddt_sd>=10) else 0L
    tags$div(
      tags$div(class="ram-confidence-metrics",
        tags$span(class="ram-confidence-metric",
          sprintf("%s models analysed",result$analyzed_models)),
        tags$span(class="ram-confidence-metric",
          sprintf("%s residues present in every model",result$common_residues)),
        tags$span(class="ram-confidence-metric",
          sprintf("%s residues with Rama8000 disagreement",standard_changes)),
        tags$span(class="ram-confidence-metric",
          sprintf("%s residues with backbone-state disagreement",
            if("basin_changes" %in% names(data))
              sum(data$basin_changes,na.rm=TRUE) else 0L)),
        tags$span(class="ram-confidence-metric",
          sprintf("%s residues with ≥20° angular SD",angular_variable)),
        tags$span(class="ram-confidence-metric",
          sprintf("%s residues with pLDDT SD ≥10",confidence_variable))
      ),
      tags$p(class="ram-confidence-explainer",
        if(identical(result$source,"alphafold3"))
          "Residue-level disagreement is calculated from backbone geometry, Rama8000 and pLDDT. AF3 pTM, ipTM and ranking score are retained separately per sample and do not modify the residue variability map."
        else
          "These values quantify disagreement among prediction models/seeds. They do not demonstrate molecular motion or experimental conformational heterogeneity."),
      if(result$limited)
        tags$p(class="ram-confidence-warning",
          "Only the first 30 models were analysed.")
    )
  })

  output$predictionEnsembleModelRows <- DT::renderDT({
    result <- prediction_ensemble_matches()
    req(result)
    data <- result$model_summary
    if(!nrow(data)) return(DT::datatable(data,rownames=FALSE))
    fields <- c("model","residues","finite_phi_psi","rama8000_outliers",
      "plddt_mean","plddt_min","ptm","iptm","ranking_score",
      "fraction_disordered","has_clash")
    fields <- fields[fields %in% names(data)]
    shown <- data[,fields,drop=FALSE]
    for(field in intersect(c("plddt_mean","plddt_min"),names(shown)))
      shown[[field]] <- round(shown[[field]],1L)
    for(field in intersect(c("ptm","iptm","ranking_score",
                             "fraction_disordered"),names(shown)))
      shown[[field]] <- round(shown[[field]],3L)
    if("has_clash" %in% names(shown))
      shown$has_clash <- ifelse(is.na(shown$has_clash),"",
        ifelse(shown$has_clash,"Yes","No"))
    names(shown) <- c(
      model="Model/sample",residues="Residues",
      finite_phi_psi="Finite φ/ψ",rama8000_outliers="Rama8000 outliers",
      plddt_mean="pLDDT mean",plddt_min="pLDDT min",
      ptm="pTM",iptm="ipTM",ranking_score="AF3 ranking score",
      fraction_disordered="Disordered fraction",has_clash="AF3 clash flag"
    )[names(shown)]
    DT::datatable(shown,rownames=FALSE,selection="none",
      options=list(pageLength=8,scrollX=TRUE,autoWidth=FALSE,dom="tip"),
      class="compact stripe")
  },server=FALSE)

  output$predictionEnsembleTrack <- renderUI({
    result <- prediction_ensemble_matches()
    if(is.null(result) || !nrow(result$summary)) return(NULL)
    data <- result$summary
    spread <- pmax(data$phi_sd,data$psi_sd,na.rm=TRUE)
    spread[!is.finite(data$phi_sd) & !is.finite(data$psi_sd)] <- NA_real_
    band <- ifelse(!is.finite(spread),"unavailable",
      ifelse(spread<5,"stable",
        ifelse(spread<15,"moderate",
          ifelse(spread<30,"variable","high"))))
    cells <- lapply(seq_len(nrow(data)),function(i) {
      label <- paste0(data$resn[[i]]," ",data$chain[[i]],":",
        data$resi[[i]],data$insertion_code[[i]])
      standard <- if("rama8000_changes" %in% names(data) &&
                     isTRUE(data$rama8000_changes[[i]]))
        " · Rama8000 category differs across models" else ""
      basin <- if("basin_changes" %in% names(data) &&
                  isTRUE(data$basin_changes[[i]]))
        paste0(" · backbone state differs across models",
          if("basin_mode" %in% names(data) && !is.na(data$basin_mode[[i]]))
            paste0(" (mode ",data$basin_mode[[i]],")") else "") else ""
      tags$button(type="button",
        class=paste("ram-ensemble-cell",
          paste0("ram-ensemble-",band[[i]]),
          if(nzchar(standard)) "has-standard-change" else "",
          if(nzchar(basin)) "has-basin-change" else ""),
        "data-chain"=data$chain[[i]],
        "data-resi"=data$resi[[i]],
        "data-insertion"=data$insertion_code[[i]],
        title=paste0(label," · angular SD ",
          if(is.finite(spread[[i]])) sprintf("%.1f°",spread[[i]]) else "N/A",
          if("plddt_mean" %in% names(data) && is.finite(data$plddt_mean[[i]]))
            sprintf(" · mean pLDDT %.1f",data$plddt_mean[[i]]) else "",
          basin,standard),
        "aria-label"=paste("Inspect",label,"from prediction ensemble")
      )
    })
    tags$section(class="ram-ensemble-track-panel",
      tags$div(class="ram-ensemble-track-head",
        tags$strong("Prediction variability map"),
        tags$span("max circular SD of φ or ψ per residue")
      ),
      tags$div(class="ram-ensemble-track",role="group",
        "aria-label"="Prediction ensemble residue variability",cells),
      tags$div(class="ram-ensemble-track-legend",
        tags$span(class="ram-ensemble-stable","<5°"),
        tags$span(class="ram-ensemble-moderate","5–15°"),
        tags$span(class="ram-ensemble-variable","15–30°"),
        tags$span(class="ram-ensemble-high","≥30°"),
        tags$span(class="ram-ensemble-basin-mark",
          "double outline = backbone-state disagreement"),
        tags$span(class="ram-ensemble-standard-mark",
          "inner outline = Rama8000 disagreement"))
    )
  })

  output$predictionEnsembleRows <- DT::renderDT({
    result <- prediction_ensemble_matches()
    req(result)
    data <- result$summary
    if(!nrow(data)) return(DT::datatable(data,rownames=FALSE))
    fields <- c("chain","resi","insertion_code","resn","models_present",
      "phi_sd","psi_sd","basin_mode","basin_consistency",
      "rama8000_mode","rama8000_consistency",
      "plddt_mean","plddt_sd","plddt_min","plddt_max")
    fields <- fields[fields %in% names(data)]
    shown <- data[,fields,drop=FALSE]
    for(field in intersect(c("phi_sd","psi_sd","plddt_mean","plddt_sd",
                             "plddt_min","plddt_max"),names(shown)))
      shown[[field]] <- round(shown[[field]],1L)
    if("basin_consistency" %in% names(shown))
      shown$basin_consistency <- round(100*shown$basin_consistency,1L)
    if("rama8000_consistency" %in% names(shown))
      shown$rama8000_consistency <- round(100*shown$rama8000_consistency,1L)
    names(shown) <- c(
      chain="Chain",resi="Residue",insertion_code="Ins.",resn="AA",
      models_present="Models",phi_sd="φ SD (°)",psi_sd="ψ SD (°)",
      basin_mode="Backbone state",
      basin_consistency="State agreement (%)",
      rama8000_mode="Rama8000 mode",
      rama8000_consistency="Rama8000 agreement (%)",
      plddt_mean="pLDDT mean",plddt_sd="pLDDT SD",
      plddt_min="pLDDT min",plddt_max="pLDDT max"
    )[names(shown)]
    DT::datatable(shown,rownames=FALSE,selection="single",
      options=list(pageLength=12,scrollX=TRUE,autoWidth=FALSE,dom="ftip"),
      class="compact stripe hover")
  },server=FALSE)

  observeEvent(input$predictionEnsembleRows_rows_selected, {
    data <- req(prediction_ensemble_matches())$summary
    ix <- input$predictionEnsembleRows_rows_selected[[1L]]
    if(!length(ix) || !is.finite(ix) || ix<1L || ix>nrow(data)) return()
    row <- data[ix,,drop=FALSE]
    select_from(list(chain=as.character(row$chain[[1L]]),
      resi=as.integer(row$resi[[1L]]),
      insertion_code=as.character(row$insertion_code[[1L]])))
  })

  observeEvent(input$ramPredictionEnsemblePick, {
    select_from(input$ramPredictionEnsemblePick)
  },ignoreInit=TRUE)

  output$downloadPredictionEnsemble <- downloadHandler(
    filename=function() safe_filename("prediction-ensemble-residues.csv"),
    content=function(file) utils::write.csv(
      req(prediction_ensemble_matches())$summary,file,row.names=FALSE,na="")
  )
  output$downloadPredictionEnsembleModels <- downloadHandler(
    filename=function() safe_filename("prediction-ensemble-models.csv"),
    content=function(file) {
      result <- req(prediction_ensemble_matches())
      models <- result$model_summary
      if(!is.null(result$provenance)) {
        provenance <- result$provenance
        if(nrow(provenance)!=nrow(models))
          stop("Prediction ensemble provenance no longer matches model order.")
        for(field in setdiff(names(provenance),"model"))
          models[[field]] <- provenance[[field]]
      }
      utils::write.csv(models,file,row.names=FALSE,na="")
    }
  )
  output$downloadPredictionEnsembleReport <- downloadHandler(
    filename=function() safe_filename("prediction-ensemble-report.html"),
    content=function(file) {
      structure <- req(loaded())
      result <- req(prediction_ensemble_matches())
      ram_save_prediction_ensemble_report(file,result,list(
        structure=structure$name,
        reference_set=input$bgtype,
        displayed_background=input$background,
        ramplotr_classification=input$validationMode,
        loaded_prediction_source=structure$declared_source
      ))
    }
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

    standard <- data[!is.na(data$rama8000_region), , drop=FALSE]
    standard_n <- nrow(standard)
    standard_count <- function(region)
      sum(standard$rama8000_region == region, na.rm=TRUE)
    standard_pct <- function(n) if (!standard_n) "n/a" else
      sprintf("%.2f%%", 100 * n / standard_n)
    standard_row <- function(label, n) tags$tr(
      tags$td(label),
      tags$td(class="ram-numeric", format(n, big.mark=",")),
      tags$td(class="ram-numeric", standard_pct(n))
    )
    standard_outliers <- standard_count("Outlier")
    mapped_n <- if ("canonical_status" %in% names(data))
      sum(data$canonical_status=="mapped",na.rm=TRUE) else 0L
    canonical_accessions <- if ("uniprot_accession" %in% names(data))
      unique(na.omit(as.character(data$uniprot_accession))) else character()
    mapping_state <- canonical_status()

    tags$div(class = "ram-summary",
      tags$div(class = "ram-summary-metrics",
        metric("Selected residues", nrow(data), "Across selected chains"),
        metric("RamplotR not allowed", outlier, "Native density regions"),
        metric("Rama8000 outliers", standard_outliers,
               "Six-class standard validation"),
        if (!is.null(mapping_state) &&
            mapping_state$state %in% c("mapped","partial"))
          metric("UniProt mapped", mapped_n,
            if(length(canonical_accessions))
              paste(canonical_accessions,collapse=", ")
            else "Canonical coordinates")
      ),
      if (!is.null(mapping_state))
        tags$div(class="ram-canonical-summary",
          tags$strong("Canonical coordinates"),
          if (identical(mapping_state$state,"searching"))
            tags$span("Retrieving PDBe SIFTS mapping…")
          else if (identical(mapping_state$state,"mapped"))
            tags$span(sprintf("%d selected residues currently map to UniProt%s.",
              mapped_n,
              if(length(canonical_accessions))
                paste0(" ",paste(canonical_accessions,collapse=", "))
              else ""))
          else if (identical(mapping_state$state,"partial"))
            tags$span(paste0(
              "SIFTS ranges were found, but only unambiguous one-to-one ",
              "author-number ranges are expanded. Nonlinear ranges remain unresolved."))
          else if (identical(mapping_state$state,"unavailable"))
            tags$span(mapping_state$message)
          else if (identical(mapping_state$state,"error"))
            tags$span("Canonical mapping could not be retrieved; local PDB numbering remains available.")
        ),
      tags$h3("RamplotR density regions"),
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
      tags$h3("Rama8000 standard validation"),
      tags$p(class="ram-summary-note",
        "Current cctbx/Phenix-style six-class evaluation: General, Gly, cis-Pro, trans-Pro, pre-Pro and Ile/Val. This result is independent of the RamplotR display background."),
      tags$div(class="ram-summary-table-wrap",
        tags$table(class="ram-summary-table",
          tags$thead(tags$tr(tags$th("Region"),tags$th("Residues"),tags$th("Share"))),
          tags$tbody(
            standard_row("Favored", standard_count("Favored")),
            standard_row("Allowed", standard_count("Allowed")),
            standard_row("Outlier", standard_outliers),
            standard_row("Total classified", standard_n)
          ))),
      tags$div(class = "ram-summary-footnotes",
        tags$div(
          tags$strong(format(sum(is.na(data$region) &
            !data$resn %in% c("GLY", "PRO")), big.mark = ",")),
          tags$span("Missing or terminal angles (RamplotR)")
        ),
        tags$div(
          tags$strong(format(sum(is.na(data$rama8000_region)), big.mark=",")),
          tags$span("Missing or terminal angles (Rama8000)")
        )
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
  selected_prediction_ensemble_context <- reactive({
    row <- selected_row()
    result <- prediction_ensemble_matches()
    if (is.null(row) || is.null(result) || !is.data.frame(result$summary) ||
        !nrow(result$summary)) return(NULL)
    data <- result$summary
    insertion <- as.character(row$insertion_code[[1L]])
    if (is.na(insertion)) insertion <- ""
    ix <- which(
      as.character(data$chain) == as.character(row$chain[[1L]]) &
      as.integer(data$resi) == as.integer(row$resi[[1L]]) &
      as.character(data$insertion_code) == insertion &
      toupper(as.character(data$resn)) == toupper(as.character(row$resn[[1L]]))
    )
    if (!length(ix)) return(NULL)
    context <- data[ix[[1L]],,drop=FALSE]
    context$ensemble_models_total <- as.integer(result$analyzed_models)
    context$ensemble_source <- as.character(result$source)
    context
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
    structure <- loaded()
    local_context <- NULL
    # Hetero coordinates are retained separately for model 1 so normal
    # protein-only backbone/model bookkeeping is unchanged. Do not reuse
    # model-1 ligand positions for another selected structural model.
    if (!is.null(structure) &&
        (structure$nmodels==1L || identical(current_model(),1L))) {
      local_context <- ram_nearby_hetero_context(
        ram_model_at(structure$pdb,current_model()),row,
        max_distance=6,max_hits=3L
      )
    }
    ensemble_context <- selected_prediction_ensemble_context()
    evidence <- ram_residue_evidence(
      row,local_context=local_context,ensemble_context=ensemble_context)
    evidence_item <- function(item) {
      level <- as.character(item$level[[1L]])
      tags$div(class=paste("ram-evidence-item",paste0("ram-evidence-",level)),
        tags$div(class="ram-evidence-item-head",
          tags$strong(as.character(item$title[[1L]])),
          tags$span(as.character(item$source[[1L]]))),
        tags$p(as.character(item$detail[[1L]]))
      )
    }
    tags$div(class = "ram-inspector-data",
      tags$div(tags$strong(sprintf("%s %d%s · %s",
        if (nzchar(row$chain[[1L]])) paste("Chain", row$chain[[1L]]) else "Chain",
        as.integer(row$resi[[1L]]), row$insertion_code[[1L]],
        row$resn[[1L]])),
        if ("canonical_status" %in% names(row) &&
            identical(as.character(row$canonical_status[[1L]]),"mapped"))
          tags$span(class="ram-inspector-canonical",
            sprintf("UniProt %s:%d",
              row$uniprot_accession[[1L]],
              as.integer(row$uniprot_resi[[1L]]))),
        if ("canonical_status" %in% names(row) &&
            identical(as.character(row$canonical_status[[1L]]),"ambiguous"))
          tags$span(class="ram-inspector-warning",
            "UniProt mapping ambiguous"),
        tags$span(class = "ram-inspector-classification",
          if (is.na(row$region[[1L]])) "Missing angles" else row$region[[1L]])),
      tags$div(class = "ram-inspector-angles",
        tags$span(paste("φ", angle(row$phi[[1L]]))),
        tags$span(paste("ψ", angle(row$psi[[1L]]))),
        tags$span(if (is.finite(row$density[[1L]]))
          sprintf("Density percentile %.1f", row$density[[1L]]) else ""),
        if ("rama8000_region" %in% names(row) &&
            !is.na(row$rama8000_region[[1L]]))
          tags$span(class = if (identical(row$rama8000_region[[1L]], "Outlier"))
                      "ram-inspector-warning" else "ram-inspector-plddt",
            sprintf("Rama8000 %s · %s · %.2f%%",
              row$rama8000_region[[1L]], row$rama8000_group[[1L]],
              100 * row$rama8000_score[[1L]])),
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
            identical(as.character(row$rama8000_region[[1L]]), "Outlier"))
          tags$span(class = "ram-inspector-warning",
            "High model confidence with a Rama8000 outlier; inspect locally."),
        if (!is.null(ensemble_context)) {
          phi_sd <- suppressWarnings(as.numeric(ensemble_context$phi_sd[[1L]]))
          psi_sd <- suppressWarnings(as.numeric(ensemble_context$psi_sd[[1L]]))
          pmean <- if ("plddt_mean" %in% names(ensemble_context))
            suppressWarnings(as.numeric(ensemble_context$plddt_mean[[1L]]))
            else NA_real_
          psd <- if ("plddt_sd" %in% names(ensemble_context))
            suppressWarnings(as.numeric(ensemble_context$plddt_sd[[1L]]))
            else NA_real_
          basin_mode <- if ("basin_mode" %in% names(ensemble_context))
            as.character(ensemble_context$basin_mode[[1L]]) else NA_character_
          basin_consistency <- if ("basin_consistency" %in% names(ensemble_context))
            suppressWarnings(as.numeric(ensemble_context$basin_consistency[[1L]]))
            else NA_real_
          models <- suppressWarnings(as.integer(
            ensemble_context$models_present[[1L]]))
          total <- suppressWarnings(as.integer(
            ensemble_context$ensemble_models_total[[1L]]))
          tags$span(class="ram-inspector-plddt",
            paste0("Prediction ensemble · ",models,"/",total," models",
              if(is.finite(phi_sd)) sprintf(" · φ SD %.1f°",phi_sd) else "",
              if(is.finite(psi_sd)) sprintf(" · ψ SD %.1f°",psi_sd) else "",
              if(!is.na(basin_mode)) paste0(" · state ",basin_mode) else "",
              if(is.finite(basin_consistency))
                sprintf(" %.0f%% agreement",100*basin_consistency) else "",
              if(is.finite(pmean)) sprintf(" · pLDDT %.1f",pmean) else "",
              if(is.finite(psd)) sprintf(" ± %.1f",psd) else ""))
        }
      ),
      tags$details(class="ram-evidence-panel",
        open=if(nrow(evidence)>0L) "open" else NULL,
        tags$summary(
          if(nrow(evidence))
            sprintf("Why inspect this residue? · %d signal%s",
              nrow(evidence),if(nrow(evidence)==1L) "" else "s")
          else "Why inspect this residue? · no obvious issue"
        ),
        if(nrow(evidence))
          tags$div(class="ram-evidence-list",
            lapply(seq_len(nrow(evidence)),function(i)
              evidence_item(evidence[i,,drop=FALSE])))
        else
          tags$p(class="ram-evidence-none",
            "No unusual signal is present in the currently available backbone, prediction-confidence, prediction-ensemble or attached official-validation evidence.")
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
