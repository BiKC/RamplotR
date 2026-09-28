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
  "#7FC97F", "#BEAED4", "#FDC086", "#FFFF99",
  "#386CB0", "#F0027F", "#BF5B17", "#666666"
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
          class = "ram-sidebar", "aria-label" = "Analysis settings",
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
          ),
          tags$section(
            class = "ram-panel",
            tags$div(
              class = "ram-section-heading",
              tags$div(
                tags$h3("3D display"),
                tags$p("Choose which molecular features appear in the viewer.")
              )
            ),
            tags$div(
              class = "ram-toggles",
              checkboxInput("ligands", "Ligands"),
              checkboxInput("dna", "DNA"),
              checkboxInput("rna", "RNA")
            ),
            tags$div(class = "ram-panel-divider"),
            tags$div(
              class = "ram-toggles",
              checkboxInput("spinning", "Spin"),
              checkboxInput("rocking", "Rock", value = TRUE)
            )
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
                tags$span(id = "ram-current-structure",
                          class = "ram-status", "No structure loaded")
              ),
              tags$div(
                class = "ram-charts",
                tags$section(
                  class = "ram-chart-card", "aria-label" = "Ramachandran plot",
                  tags$h3(class = "ram-chart-label", "Residue distribution"),
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
                           role = "img", "aria-label" = "Interactive Ramachandran plot")
                ),
                tags$section(
                  class = "ram-chart-card", "aria-label" = "3D molecular viewer",
                  tags$h3(class = "ram-chart-label", "Molecular structure"),
                  tags$p(class = "ram-chart-help",
                         "Rotate and zoom the structure. Viewer controls are on the left."),
                  tags$div(class = "ram-ngl",
                           NGLVieweR::NGLVieweROutput("NGL")),
                  tags$p(class = "ram-ngl-note",
                         "Structure rendering is provided by NGLVieweR.")
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
  # session$onSessionEnded(stopApp)
  session$userData$previousPDB <- ""
  # reactive(bio3d::write.pdb(pdb = pdb(), file = paste0(accPDB(), '.pdb')))

  # Reference cutoffs are deterministic and calculated once per density grid.
  searchlimit <- function(matrix, percentage, x = 1) {
    ram_density_thresholds(matrix, percentage)[[1L]]
  }

  observeEvent(input$ligands,{
    if(input$ligands){
      NGLVieweR_proxy("NGL")%>%addSelection(type = "ball+stick", param=list(name="ligand",sele= "ligand"))
    }else{
      NGLVieweR_proxy("NGL")%>%removeSelection("ligand")
    }
  })
  
  observeEvent(input$dna,{
    if(input$dna){
      NGLVieweR_proxy("NGL")%>%addSelection(type = "cartoon", param=list(name="dna",sele= "dna"))
    }else{
      NGLVieweR_proxy("NGL")%>%removeSelection("dna")
    }
  })
  
  observeEvent(input$rna,{
    if(input$rna){
      NGLVieweR_proxy("NGL")%>%addSelection(type = "cartoon", param=list(name="rna",sele= "rna"))
    }else{
      NGLVieweR_proxy("NGL")%>%removeSelection("rna")
    }
  })
  observeEvent(input$rocking,{
    if(input$rocking){
      NGLVieweR_proxy("NGL")%>%updateRock()
      if(input$spinning==T){
        updateCheckboxInput(session,"spinning",value = F)
      }
    }else{
      NGLVieweR_proxy("NGL")%>%updateRock(rock = F)
      
    }
  })
  observeEvent(input$spinning,{
    if(input$spinning){
      NGLVieweR_proxy("NGL")%>%updateSpin()
      if(input$rocking==T){
        updateCheckboxInput(session,"rocking",value = F)
      }
    }else{
      NGLVieweR_proxy("NGL")%>%updateSpin(spin = F)
    }
  })
  
  updateColorInputs <- function(colors){
    updateColourInput(session, "bg1", value=colors[1])
    updateColourInput(session, "bg2", value=colors[2])
    updateColourInput(session, "bg3", value=colors[3])
    updateColourInput(session, "bg4", value=colors[4])
  }
  
  observeEvent(input$bgtype,
               {
                 files<-list.files(paste0("static/",input$bgtype))
                 choices=list("Commonly used"=files[!files %in% allAA],
                              "Per amino acid"=files[files %in% allAA])
                 if (input$background %in% files) {sel=input$background}
                 else {sel=choices[1]}
                 updateSelectInput(session,"background",choices = choices, selected = sel)
               })
  
  observeEvent(input$colorscheme,{
    if (input$colorscheme =="Rampage"){
      updateColorInputs(rampage)
    }
    else if (input$colorscheme =="PDBSum"){
      updateColorInputs(pdbsum)
      
    }
  })
  
  observeEvent(input$background,{
    if (input$background %in% allAA){
      updatePickerInput(session,"AA",selected = input$background)
    }
  })
  
  observeEvent({input$bg1
    input$bg2
    input$bg3
    input$bg4},
    {
      colors<-c(input$bg1,input$bg2,input$bg3,input$bg4)
      #print(colors == rampage)
      if (all(colors == rampage)){
        updateSelectInput(session, "colorscheme",selected = "Rampage")
        
      }
      else if (all(colors == pdbsum)){
        updateSelectInput(session, "colorscheme",selected = "PDBSum")
        
      } else {
        updateSelectInput(session, "colorscheme",selected = "custom")
      }
    })

  # Process a structure only when requested. A newly opened session starts
  # with the empty plot rather than fetching the default PDB automatically.
  observeEvent(input$submit, {
    withProgress(message = "Making plot", value = 0, {
      inputType<-""
      isolate({
        inputType <- if (identical(input$inputSource, "upload")) "file" else "code"
        if (inputType == "file" && is.null(input$structfile$datapath)) {
          showNotification("Choose a PDB or mmCIF file first.", type = "error")
          return(invisible(NULL))
        }
        accPDB <- if (inputType == "file") input$structfile$datapath else
          toupper(trimws(input$PDB))
        structure_key <- paste(inputType, accPDB, sep = ":")
        if (!identical(session$userData$previousPDB, structure_key)) {
          incProgress(1 / 4, detail = "Reading structure")
          pdb <- tryCatch(
            ram_load_structure(
              path = if (inputType == "file") accPDB else NULL,
              original_name = if (inputType == "file") input$structfile$name else NULL,
              pdb_id = if (inputType == "code") accPDB else NULL
            ),
            error = function(e) {
              showNotification(conditionMessage(e), type = "error", duration = 12)
              NULL
            }
          )
          if (is.null(pdb)) return(invisible(NULL))
          incProgress(1 / 4, detail = paste("Transforming data"))
          # Keep insertion codes and validate peptide connectivity before
          # classifying glycine, proline and pre-proline residues.
          torsion <- ram_extract_torsions(pdb)
          chains <- unique(torsion$chain)
          # Store the results in the user data
          session$userData$torsion <- torsion

          chains<-unique(session$userData$torsion$chain)

          # Shiny upload paths lack an extension; identify the format explicitly.
          viewer_format <- if (inputType == "file")
            ram_detect_format(input$structfile$name) else NULL
          nglview <- NGLVieweR(data = accPDB, format = viewer_format) %>% setRock()
          counter=1
          for (i in unique(chains)) {
            #print(i)
            nglview<-addRepresentation(NGLVieweR = nglview,type = "cartoon", param=list("sele"= paste0(":", i,"  and protein"), "color"= color_set[((counter - 1) %% length(color_set)) + 1]))
            counter=counter+1
          }
          output$NGL <- NGLVieweR::renderNGLVieweR(nglview)
          
          output$chainColors <- renderUI({
            isolate({
              widgets <- lapply(seq_along(chains), function(k) {
                colourpicker::colourInput(
                  paste0("chain", chains[[k]]),
                  label = paste("Chain", chains[[k]]),
                  value = color_set[((k - 1L) %% length(color_set)) + 1L]
                )
              })
              dropdown(widgets, label = "Chain color settings")
            })
          })

          output$chains <- renderUI({
            # pickerInput
            isolate({
              pickerInput(
                "chainselection",
                "Chain selection",
                choices = chains,
                multiple = T,
                options = list(`actions-box` = TRUE),
                selected = chains
              )
            })
          })
          
          
          
          session$userData$previousPDB <- structure_key
          incProgress(1 / 4, detail = paste("Filter data"))
        } else {
          incProgress(3 / 4, detail = paste("Filter data"))
        }
        if (input$background == "preProline"){
          matrix <- ram_read_reference(file.path("static", input$bgtype, "preProline"))
          # get a subset of only those amino acids that precede a proline
          torsionsubset <- data.frame()
          torsionsubset <- subset(session$userData$torsion,
                                 bonded_to_next & !is.na(next_resn) &
                                   next_resn == "PRO")
          torsionsubset <- subset(torsionsubset, resn %in% input$AA)
          # also subset for chains
          # since chainselection is added as uiOutput, it is not available in the beginning, so check if it exists, otherwise subset for all chains
          if (!is.null(input$chainselection)) {
            torsionsubset <- subset(torsionsubset, chain %in% input$chainselection)
          }
          
        }
        else if (!input$background %in% allAA) {
          matrix <- ram_read_reference(file.path("static", input$bgtype, "General"))
          torsionsubset <- session$userData$torsion
          torsionsubset <- subset(torsionsubset, resn %in% input$AA)
          # also subset for chains
          if (!is.null(input$chainselection)) {
            torsionsubset <- subset(torsionsubset, chain %in% input$chainselection)
          }
        } else {
          matrix <- ram_read_reference(file.path("static", input$bgtype, input$background))
          #updatePickerInput(session, "AA", selected = input$background)
          torsionsubset <- subset(session$userData$torsion, resn %in% input$AA)
          # also subset for chains
          if (!is.null(input$chainselection)) {
            torsionsubset <- subset(torsionsubset, chain %in% input$chainselection)
          }
        }
        incProgress(1 / 4, detail = paste("Creating plot"))
        # session$sendCustomMessage("updateFig", input$PDB)
        # get input chain colors from ui
        ttab <- reactive({
          ram_classify_torsions(
            session$userData$torsion,
            reference_dir = file.path("static", input$bgtype),
            selected_reference = matrix,
            mode = input$validationMode,
            threshold_fn = searchlimit
          )
        })
        output$regions <- renderDataTable({

          # filter the data frame to only include the amino acids that are in the selected region
          if (input$regionselect == "All") {
            ttabsub <- ttab()
          } else {
            ttabsub <- subset(ttab(), region == input$regionselect)
            ttabsub <- ttabsub[, c("resi", "chain", "resn", "phi", "psi", "density")]
          }
          # also filter for the selected amino acids (input$AA) and chains (input$chainselection)
          ttabsub <- subset(ttabsub, resn %in% input$AA)
          if (!is.null(input$chainselection)) {
            ttabsub <- subset(ttabsub, chain %in% input$chainselection)
          }
          ttabsub
        })

        
          # display statistics for the regions:
          # for the following, we do not include glycine and proline
          # Favoured regions (no. of residues, %)
          # Allowed regions (no. of residues, %)
          # Generously allowed regions (no. of residues, %)
          # Not allowed regions (no. of residues, %)
          # Total no. of residues (no. of residues, %)
          # --------------------------
          # End-residues (Excl. Gly and Pro)
          # --------------------------
          # Glycine residues (no. of residues)
          # Proline residues (no. of residues)
          # --------------------------
          # Total no. of residues (no. of residues)

          # All statistics use the same amino-acid and chain selection as the plot.
          stats_table <- subset(ttab(), resn %in% input$AA)
          if (!is.null(input$chainselection)) {
            stats_table <- subset(stats_table, chain %in% input$chainselection)
          }
          # define a function to negate the %in% operator
          `%nin%` <- Negate(`%in%`)
          # exclude glycine and proline
          exclude <- c("GLY", "PRO")
          
          missing_angles <- is.na(stats_table$region)
          end_count <- sum(missing_angles & !stats_table$resn %in% exclude)
          eligible <- stats_table[!missing_angles & !stats_table$resn %in% exclude, , drop = FALSE]
          count_no_gly_pro <- nrow(eligible)
          region_count <- function(name) sum(eligible$region == name, na.rm = TRUE)
          region_percent <- function(n) {
            if (count_no_gly_pro == 0L) return(NA_real_)
            round(100 * n / count_no_gly_pro, 2)
          }
          fr_count <- region_count("Favoured")
          ar_count <- region_count("Allowed")
          gar_count <- region_count("Generously allowed")
          nar_count <- region_count("Not allowed")
          fr_percent <- region_percent(fr_count)
          ar_percent <- region_percent(ar_count)
          gar_percent <- region_percent(gar_count)
          nar_percent <- region_percent(nar_count)
          total_count <- count_no_gly_pro
          total_percent <- region_percent(total_count)
          gly_count <- sum(stats_table$resn == "GLY")
          pro_count <- sum(stats_table$resn == "PRO")
          total_count2 <- nrow(stats_table)

          # create html output to display the statistics
        
            output$summary <- renderUI({
  percent_label <- function(p) {
    if (is.na(p)) "n/a" else sprintf("%.2f%%", p)
  }
  metric <- function(label, value, hint) {
    tags$div(
      class = "ram-summary-metric",
      tags$span(class = "ram-summary-metric-label", label),
      tags$strong(as.character(value)),
      tags$small(hint)
    )
  }
  region_row <- function(label, count, percent) {
    tags$tr(
      tags$td(label),
      tags$td(class = "ram-numeric", format(count, big.mark = ",")),
      tags$td(class = "ram-numeric", percent_label(percent))
    )
  }
  tags$div(
    class = "ram-summary",
    tags$div(
      class = "ram-summary-metrics",
      metric("Selected residues", total_count2, "Across selected chains"),
      metric("Classified, excluding Gly/Pro", total_count,
             "Residues with defined backbone angles"),
      metric("Outliers", nar_count, "Outside the selected reference regions")
    ),
    tags$h3("Region breakdown"),
    tags$p(class = "ram-summary-note",
           "Percentages use classified residues other than glycine and proline as the denominator."),
    tags$div(
      class = "ram-summary-table-wrap",
      tags$table(
        class = "ram-summary-table",
        tags$thead(
          tags$tr(tags$th("Region"), tags$th("Residues"), tags$th("Share"))
        ),
        tags$tbody(
          region_row("Favoured", fr_count, fr_percent),
          region_row("Allowed", ar_count, ar_percent),
          region_row("Generously allowed", gar_count, gar_percent),
          region_row("Not allowed", nar_count, nar_percent),
          region_row("Total classified", total_count, total_percent)
        )
      )
    ),
    tags$div(
      class = "ram-summary-footnotes",
      tags$div(
        tags$strong(format(end_count, big.mark = ",")),
        tags$span("Missing or terminal angles (excluding Gly/Pro)")
      ),
      tags$div(
        tags$strong(format(gly_count, big.mark = ",")),
        tags$span("Glycine residues")
      ),
      tags$div(
        tags$strong(format(pro_count, big.mark = ",")),
        tags$span("Proline residues")
      )
    )
  )
})


        name<-ifelse(inputType=="file",tools::file_path_sans_ext(basename(input$structfile$name)),accPDB)

        session$sendCustomMessage(
          "process",
          list(
            df = torsionsubset,
            matrix = matrix,
            name=name,
            pdb = accPDB,
            backgroundColors = c(input$bg1, input$bg2, input$bg3, input$bg4),
            # get the colors from the chain ui color settings and add them to a vector
            chainColors = unlist(lapply(unique(torsionsubset$chain), function(x) {
              input[[paste0("chain", x)]]
            })),
            limits = ram_density_thresholds(matrix)
          )
        )
      })
    })
  }, ignoreInit = TRUE)
}


# Run the application
shinyApp(ui = ui, server = server)