#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

#shiny related packages
library(shiny)
library(shinycssloaders)
library(shinyWidgets)
library(shinyjs)
library(colourpicker)

#used for pdb files
library(bio3d)
library(Rpdb)
#used for graphs
library(plotly)
library(RColorBrewer)

#Used for processing data
library(plyr)

colorSet <- RColorBrewer::brewer.pal(8, "Accent")
previousPDB<-""

#load bg image
fig <- readRDS("static/bgfig.RDS")
densMatGeneral <- readRDS("static/densMatGeneral")
densMatGly <- readRDS("static/densMatGly")
densMatprePro <- readRDS("static/densMatprePro")
densMatPro <- readRDS("static/densMatPro")

allAA<-c(
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
      actionButton("submit","Apply Changes",class="btn-primary btn-lg"),
      hr(),
      
      selectInput("background", "Choose background plot", c("General","Glycine","Preproline","Proline")),
      
      pickerInput(
        "AA",
        "Amino acid selection",
        allAA,
        multiple = T,
        options = list(`actions-box` = TRUE),
        selected = allAA
      )
      ,
      hr(),
      dropdown(
        list(
          colourpicker::colourInput("bg1", "Not allowed region", value = "#f1eef6"),
          colourpicker::colourInput("bg2", "Generously allowed regions", value =
                                      '#bdc9e1'),
          colourpicker::colourInput("bg3", "Allowed regions", value = '#74a9cf'),
          colourpicker::colourInput("bg4", "Favoured regions", value = '#0570b0')
        ),
        label = "Background density settings"
      ),
      hr(),
      uiOutput("chainColors"),
      hr(),
      
      fluidRow(
        tags$div(id = "viewport", style = "width:600px; height:600px;"),
        tags$script(
          "var stage = new NGL.Stage('viewport',{backgroundColor:'white'})"
        ),
        tags$script(
          "
      Shiny.addCustomMessageHandler('updateFig', function(pdbcode) {
        stage.removeAllComponents();
        stage.loadFile('rcsb://'+pdbcode, {defaultRepresentation: true});
      });
    "
        ),
        tags$script(src = "custom.js")
        
      )
      
    )
    ,
    column(8,
           # Output panel
           mainPanel(
             #plotlyOutput("ramplot"), type = 6
             tags$div(id = "plotly", style = "height:800px;width:891px")

             
           ))
  ))
)

# Define server logic required to draw a histogram
server <- function(input, output, session) {
  output$dummy <- reactive(FALSE)
  outputOptions(output, "dummy", suspendWhenHidden = FALSE)
  #session$onSessionEnded(stopApp)
  options(warn = -1)
  session$userData$previousPDB<-""
  session$userData$background<-""
  #reactive(bio3d::write.pdb(pdb = pdb(), file = paste0(accPDB(), '.pdb')))
  
  searchlimit <- function(matrix,percentage, x=1 ) {
    totalGen <- 100 * sum(matrix$z[matrix$z > x]) / sum(matrix$z)
    error <- abs(totalGen - percentage)
    if (error < 0.1) {
      return(x)
    }
    else {
      if (totalGen > percentage) {
        searchlimit(matrix,percentage, x = x + x/2)
      }
      else {
        searchlimit(matrix,percentage, x = x - x/2)
      }
    }
  }
  
  output$dummy <- reactive({
    input$submit
    withProgress(message = 'Making plot', value = 0, {
      isolate({
        accPDB <- input$PDB
        if(session$userData$previousPDB!=accPDB){
          incProgress(1 / 4, detail = paste("Fetching sequence"))
          pdb <- bio3d::read.pdb(accPDB)
          
          incProgress(1 / 4, detail = paste("Transforming data"))
          torsion <- torsion.pdb(pdb)
          print("step 1")
          tortab <- torsion[["tbl"]][, c("phi", "psi")]
          print("step 2")
          spltor <- strsplit(rownames(tortab), split = ".", fixed = T)
          print("step 3")
          torsion <- cbind(tortab, as.data.frame(do.call(rbind, spltor)))
          print("step 4")
          session$userData$torsion <-
            rename(torsion, c(
              'V1' = "resi",
              "V2" = "chain",
              "V3" = "resn"
            ))
          output$chainColors <- renderUI({
            isolate({
              x <- vector('list', length(unique(torsionsubset$chain)))
              c <- 1
              for (i in unique(torsionsubset$chain)) {
                x[[i]] <-
                  list(colourpicker::colourInput(
                    paste0("chain", i),
                    label = paste("Chain", i),
                    value = colorSet[((c-1)%%length(colorSet))+1]
                  ))
                c <- c + 1
              }
              
              dropdown(x, label = "Chain color settings")
            })
          })
          session$userData$previousPDB<-accPDB
          incProgress(1 / 4, detail = paste("Filter data"))
          
        }
        else {
          incProgress(3 / 4, detail = paste("Filter data"))
        }
        AAselection<-input$AA
        if (session$userData$background != input$background){
          session$userData$background <- input$background 
          updatePickerInput(session, "AA", selected = 
          AAselection<-allAA)}

        if (input$background == "General") {
          matrix <- densMatGeneral
          torsionsubset <- session$userData$torsion
          torsionsubset <- subset(torsionsubset, resn %in% AAselection)
          
        } else if (input$background == "Glycine") {
          matrix <- densMatGly
          # update input$AA to only include glycine
          updatePickerInput(session, "AA", selected = c("GLY")) 
          torsionsubset <- subset(session$userData$torsion, resn =="GLY")
        } else if (input$background == "Preproline") {
          matrix <- densMatprePro
          # get a subset of only those amino acids that precede a proline
          torsionsubset <- data.frame()
          for (i in 1:length(session$userData$torsion$resn)-1) {
            if (session$userData$torsion$resn[i+1] == "PRO") {
              torsionsubset <- rbind(torsionsubset, session$userData$torsion[i,])
            }
          }
          torsionsubset <- subset(torsionsubset, resn %in% AAselection)
          
        } else if (input$background == "Proline") {
          matrix <- densMatPro
          updatePickerInput(session, "AA", selected = c("PRO")) 
          torsionsubset <- subset(session$userData$torsion, resn =="PRO")
        }
        incProgress(1 / 4, detail = paste("Creating plot"))
        session$sendCustomMessage("updateFig", input$PDB)
        #get input chain colors from ui
        chainColors <- vector('list', length(unique(torsionsubset$chain)))
        for (i in unique(torsionsubset$chain)) {
          chainColors[[i]] <- input[[paste0("chain", i)]]
        }

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
            limits=c(searchlimit(matrix,85),
                     searchlimit(matrix,98),
                     searchlimit(matrix,99.95)
                     )
            
          )
        )
        #use output to update figures
        #!ouput$dummy
      })
    })
  })
  
}


# Run the application
shinyApp(ui = ui, server = server)
