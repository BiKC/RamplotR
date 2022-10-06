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
library(shinycssloaders)
library(shinyWidgets)
library(colourpicker)

# used for pdb files
library(bio3d)
library(NGLVieweR)

# Used for processing data
library(plyr)

color_set <- c(
  "#7FC97F", "#BEAED4", "#FDC086", "#FFFF99",
  "#386CB0", "#F0027F", "#BF5B17", "#666666"
)

# load bg image
fig <- readRDS("static/bgfig.RDS")
densMatGeneral <- readRDS("static/densMatGeneral")
#densMatGeneral$z<-t(densMatGeneral$z)
densMatGly <- readRDS("static/densMatGly")
#densMatGly$z<-t(densMatGly$z)

densMatPrepro <- readRDS("static/densMatprePro")
#densMatPrepro$z<-t(densMatPrepro$z)

densMatPro <- readRDS("static/densMatPro")
#densMatPro$z<-t(densMatPro$z)

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

# Define UI for application that draws a histogram
ui <- fluidPage(
  # Application title
  titlePanel("ShinyRam"),
  tags$head(
    tags$script(src = "https://cdn.rawgit.com/arose/ngl/v0.10.4-1/dist/ngl.js"),
    tags$script(src = "https://cdn.plot.ly/plotly-2.14.0.min.js")
  ),
  # Sidebar with a slider input for number of bins
  fluidPage(fluidRow(
    column(
      4,
      textInput(
        inputId = "PDB",
        label = "PDB-code:",
        value = "1BBB"
      ),
      actionButton("submit", "Apply Changes", class = "btn-primary btn-lg"),
      hr(),
      selectInput(
        "background", "Choose background plot",
        c("General", "Glycine", "Preproline", "Proline")
      ),
      pickerInput(
        "AA",
        "Amino acid selection",
        allAA,
        multiple = T,
        options = list(`actions-box` = TRUE),
        selected = allAA
      ),
      hr(),
      dropdown(
        list(
          colourpicker::colourInput("bg1", "Not allowed region",
            value = "#f1eef6"
          ),
          colourpicker::colourInput("bg2", "Generously allowed regions",
            value = "#bdc9e1"
          ),
          colourpicker::colourInput("bg3", "Allowed regions",
            value = "#74a9cf"
          ),
          colourpicker::colourInput("bg4", "Favoured regions",
            value = "#0570b0"
          )
        ),
        label = "Background density settings"
      ),
      hr(),
      uiOutput("chainColors"),
      hr(),
      fluidRow(
        NGLVieweR::NGLVieweROutput("NGL"),
        # tags$div(id = "viewport", style = "width:600px; height:600px;"),
        # tags$script(
        #  "var stage = new NGL.Stage('viewport',{backgroundColor:'white'})"
        # ),
        # tags$script(
        #  "
        # Shiny.addCustomMessageHandler('updateFig', function(pdbcode) {
        #  stage.removeAllComponents();
        #  stage.loadFile('rcsb://'+pdbcode, {defaultRepresentation: true});
        # });
        # "
        # ),
        tags$script(src = "custom.js")
      )
    ),
    column(
      8,
      # Output panel (tabsetpanel)
      mainPanel(
        tabsetPanel(
          tabPanel("Ramachandran plot", tags$div(id = "plotly", style = "height:800px;width:891px")),
          tabPanel(
            "Residues list",
            selectInput("regionselect", "Show residues for region", c("Not allowed", "Generously allowed", "Allowed", "Favoured")),
            dataTableOutput("regions")
          )
        )
      )
    )
  ))
)

# Define server logic required to draw a histogram
server <- function(input, output, session) {
  output$dummy <- reactive(FALSE)
  outputOptions(output, "dummy", suspendWhenHidden = FALSE)
  # session$onSessionEnded(stopApp)
  options(warn = -1)
  session$userData$previousPDB <- ""
  session$userData$background <- ""
  # reactive(bio3d::write.pdb(pdb = pdb(), file = paste0(accPDB(), '.pdb')))

  searchlimit <- function(matrix, percentage, x = 1) {
    totalGen <- 100 * sum(matrix$z[matrix$z > x]) / sum(matrix$z)
    error <- abs(totalGen - percentage)
    if (error < 0.1) {
      return(x)
    } else {
      if (totalGen > percentage) {
        searchlimit(matrix, percentage, x = x + x / 2)
      } else {
        searchlimit(matrix, percentage, x = x - x / 2)
      }
    }
  }
  densityToPercent <- function(matrix, x) {
    100 * sum(matrix$z[matrix$z > x]) / sum(matrix$z)
  }

  output$dummy <- reactive({
    input$submit
    withProgress(message = "Making plot", value = 0, {
      isolate({
        accPDB <- input$PDB
        if (session$userData$previousPDB != accPDB) {
          incProgress(1 / 4, detail = paste("Fetching sequence"))
          pdb <- tryCatch(expr = {
            bio3d::read.cif(accPDB)
          }, error = function(e) {
            print(e)
            bio3d::read.pdb(accPDB)
          })
          incProgress(1 / 4, detail = paste("Transforming data"))
          torsion <- torsion.pdb(pdb)
          tortab <- torsion[["tbl"]][, c("phi", "psi")]
          spltor <- strsplit(rownames(tortab), split = ".", fixed = T)
          torsion <- cbind(tortab, as.data.frame(do.call(rbind, spltor)))
          session$userData$torsion <-
            rename(torsion, c(
              "V1" = "resi",
              "V2" = "chain",
              "V3" = "resn"
            ))
          output$chainColors <- renderUI({
            isolate({
              x <- vector("list", length(unique(torsionsubset$chain)))
              c <- 1
              for (i in unique(torsionsubset$chain)) {
                x[[i]] <-
                  list(colourpicker::colourInput(
                    paste0("chain", i),
                    label = paste("Chain", i),
                    value = color_set[((c - 1) %% length(color_set)) + 1]
                  ))
                c <- c + 1
              }

              dropdown(x, label = "Chain color settings")
            })
          })
          session$userData$previousPDB <- accPDB
          incProgress(1 / 4, detail = paste("Filter data"))
        } else {
          incProgress(3 / 4, detail = paste("Filter data"))
        }
        AAselection <- input$AA
        if (session$userData$background != input$background) {
          session$userData$background <- input$background
          updatePickerInput(session, "AA",
            selected =
              AAselection <- allAA
          )
        }

        if (input$background == "General") {
          matrix <- densMatGeneral
          torsionsubset <- session$userData$torsion
          torsionsubset <- subset(torsionsubset, resn %in% AAselection)
        } else if (input$background == "Glycine") {
          matrix <- densMatGly
          # update input$AA to only include glycine
          updatePickerInput(session, "AA", selected = c("GLY"))
          torsionsubset <- subset(session$userData$torsion, resn == "GLY")
        } else if (input$background == "Preproline") {
          matrix <- densMatPrepro
          # get a subset of only those amino acids that precede a proline
          torsionsubset <- data.frame()
          for (i in 1:length(session$userData$torsion$resn) - 1) {
            if (session$userData$torsion$resn[i + 1] == "PRO") {
              torsionsubset <- rbind(torsionsubset, session$userData$torsion[i, ])
            }
          }
          torsionsubset <- subset(torsionsubset, resn %in% AAselection)
        } else if (input$background == "Proline") {
          matrix <- densMatPro
          updatePickerInput(session, "AA", selected = c("PRO"))
          torsionsubset <- subset(session$userData$torsion, resn == "PRO")
        }
        incProgress(1 / 4, detail = paste("Creating plot"))
        # session$sendCustomMessage("updateFig", input$PDB)
        output$NGL <- NGLVieweR::renderNGLVieweR(NGLVieweR(data = input$PDB) %>% addRepresentation("cartoon"))
        # get input chain colors from ui
        chain_colors <- vector("list", length(unique(torsionsubset$chain)))
        for (i in unique(torsionsubset$chain)) {
          chain_colors[[i]] <- input[[paste0("chain", i)]]
        }

        ttab <- reactive({
          limits <- c(
            searchlimit(matrix, 85), # smaller than this value is favoured
            searchlimit(matrix, 98), # smaller than this value is allowed
            searchlimit(matrix, 99.95) # smaller than this value is generously allowed and larger than this value is not allowed
          )
          ttab <- session$userData$torsion
          # add a column to the tortab data frame that contains the region of the amino acid
          # based on the limits
          # tortab has two columns containing the phi and psi angles but they should be rounded to 0 decimals
          # matrix has the density matrix where x and y are the phi and psi angles and z is the density
          # z is a two dimensional matrix with the density for each phi and psi angle
          ttab$region <- NA
          ttab$density <- NA
          for (i in 1:nrow(ttab)) {
            # check if phi and psi are not na
            if (!is.na(ttab$phi[i]) && !is.na(ttab$psi[i])) {
              # round phi and psi to 0 decimals
              phi <- round(ttab$phi[i], 0)
              psi <- round(ttab$psi[i], 0)
              # get the density for the phi and psi angles
              # index of phi in matrix
              phi_index <- which(matrix$x == phi)
              # index of psi in matrix
              psi_index <- which(matrix$y == psi)
              # get the density
              ttab$density[i] <- matrix$z[psi_index, phi_index]
              # get the region based on the limits
              if (ttab$density[i] > limits[1]) {
                ttab$region[i] <- "Favoured"
              } else if (ttab$density[i] > limits[2]) {
                ttab$region[i] <- "Allowed"
              } else if (ttab$density[i] > limits[3]) {
                ttab$region[i] <- "Generously allowed"
              } else {
                ttab$region[i] <- "Not allowed"
              }
              # change density to percentage using densityToPercent function
              ttab$density[i] <- densityToPercent(matrix, ttab$density[i])
            }
          }
          


          ttab
        })
        output$regions <- renderDataTable({

          # filter the data frame to only include the amino acids that are in the selected region
          ttabsub <- subset(ttab(), region == input$regionselect)
          ttabsub <- ttabsub[, c("resi", "chain", "resn", "phi", "psi", "density")]
          ttabsub
        })

        session$sendCustomMessage(
          "process",
          list(
            df = torsionsubset,
            matrix = matrix,
            pdb = accPDB,
            backgroundColors = c(input$bg1, input$bg2, input$bg3, input$bg4),
            # get the colors from the chain ui color settings and add them to a vector
            chainColors = unlist(lapply(unique(torsionsubset$chain), function(x) {
              input[[paste0("chain", x)]]
            })),
            limits = c(
              searchlimit(matrix, 85),
              searchlimit(matrix, 98),
              searchlimit(matrix, 99.95)
            )
          )
        )
      })
    })
  })
}


# Run the application
shinyApp(ui = ui, server = server)