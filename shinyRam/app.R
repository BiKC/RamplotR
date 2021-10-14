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

#load bg image
fig <- readRDS("bgfig.RDS")
densMatgeneral <- readRDS("densMatGeneral")

# Define UI for application that draws a histogram
ui <- fluidPage(# Application title
    titlePanel("ShinyRam"),
    tags$head(
        tags$script(src = "https://cdn.rawgit.com/arose/ngl/v0.10.4-1/dist/ngl.js"),
        
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
            submitButton(),
            hr(),
            
              pickerInput("AA",
                          "Amino acid selection",
                          c("ALA","ARG","ASN","ASP","CYS","GLU","GLN","GLY","HIS","ILE","LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL"),
                          multiple = T,options = list(`actions-box` = TRUE),
                          selected = c("ALA","ARG","ASN","ASP","CYS","GLU","GLN","GLY","HIS","ILE","LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL"))
            ,
            hr(),
            dropdown(
                list(
                    colourpicker::colourInput("bg1", "Not allowed region", value = "#f1eef6"),
                    colourpicker::colourInput("bg2", "Generously allowed regions", value ='#bdc9e1'),
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
                tags$script("var stage = new NGL.Stage('viewport',{backgroundColor:'white'})"),
                tags$script("
      Shiny.addCustomMessageHandler('updateFig', function(pdbcode) {
        stage.removeAllComponents();
        stage.loadFile('rcsb://'+pdbcode, {defaultRepresentation: true});
      });
    ")
            )
            
        )
        ,
        column(8,
               # Output panel
               mainPanel(withSpinner(
                   plotlyOutput("ramplot"), type = 6
               )))
    )))

# Define server logic required to draw a histogram
server <- function(input, output, session) {
    session$onSessionEnded(stopApp)
    options(warn = -1)
    print("hi")
    accPDB <- reactive(input$PDB)
    pdb <- reactive(bio3d::read.pdb(accPDB()))
    
    reactive(bio3d::write.pdb(pdb = pdb(), file = paste0(accPDB(), '.pdb')))
    
    tor <- reactive({
        x <- torsion.pdb(pdb())
        tortab <- x[["tbl"]]
        spltor <-
            strsplit(rownames(tortab), split = ".", fixed = T)
        x <- cbind(tortab,as.data.frame(do.call(rbind,spltor)))
        x<-rename(x,c('V1'="resi","V2"="chain","V3"="resn"))
        x
    })
    
    torsubset <- reactive({
      subset(tor(),resn %in% input$AA)
    })
    print("there")

    
    output$chainColors <- renderUI({
        x <- vector('list', length(unique(torsubset()$chain)))
        c <- 1
        for (i in unique(torsubset()$chain)) {
            x[[i]] <-
                list(colourpicker::colourInput(
                    paste0("chain", i),
                    label = paste("Chain", i),
                    value = colorSet[c]
                ))
            c <- c + 1
        }
        dropdown(x, label = "Chain color settings")
    })
    
    output$ramplot <- renderPlotly({
        session$sendCustomMessage("updateFig", input$PDB)
        
        densMat <- densMatgeneral
        fig <-
            plot_ly(
                x = densMat$x,
                y = densMat$y,
                z = t(densMat$z),
                type = "contour",
                colorscale = "hot",
                showscale = F,
                showlegend = F,
                hoverinfo = "none",
                contours = list(
                    type = "constraint",
                    operation = "<",
                    value = 0.00000004293,
                    showlines = F
                ),
                fillcolor = input$bg2
            ) %>% add_contour(
                z = t(densMat$z),
                contours = list(
                    type = "constraint",
                    operation = "<",
                    value = 0.0000012354
                ),
                fillcolor = input$bg3,
                opacity = 0.8
            ) %>% add_contour(
                z = t(densMat$z),
                contours = list(
                    type = "constraint",
                    operation = "<",
                    value = 0.00001265
                ),
                fillcolor = input$bg4,
                opacity = 0.8
            ) %>% add_contour(
                z = t(densMat$z),
                contours = list(
                    type = "constraint",
                    operation = ">",
                    value = 0.00000004293
                ),
                fillcolor = input$bg1,
                opacity = 0.8
            ) 
        counter <- 1
        t <- torsubset()
        for (i in unique(t$chain)) {
            x <- input[[paste0("chain", i)]]
            if (!is.null(x)) {
                c = input[[paste0("chain", i)]]
            }
            else{
                c = colorSet[counter]
            }
            fig <- add_trace(
                p = fig,
                x = t$phi[which(t$chain == i)],
                y = t$psi[which(t$chain == i)],
                type = "scatter",
                mode = "markers",
                showlegend = T,
                marker = list(color = c),
                name = i,
                text = paste0(
                    i,
                    ": ",
                    t$resi[which(t$chain == i)],
                    " ",
                    t$resn[which(t$chain == i)],
                    ": (",
                    round(t$phi[which(t$chain == i)], 2),
                    ",",
                    round(t$psi[which(t$chain == i)], 2),
                    ")"
                ),
                hoverinfo = 'text'
            )
            
            counter <- counter + 1
            
        }
        
        fig <- fig %>% layout(
            autosize = F,
            width = 891,
            height = 800,
            legend = list(x = 1, y = 0.5),
            title = accPDB(),
            xaxis = list(title = "Phi \U03D5 (degrees)"),
            yaxis = list(title = "Psi ψ (degrees)"),
            annotations = list(
                text = 'Made using shinyRam',
                font = list(size = 12),
                showarrow = FALSE,
                xref = 'paper',
                x = 1,
                yref = 'paper',
                y = -0.065
            )
        )
        fig
    })
    
    
    
    
}


# Run the application
shinyApp(ui = ui, server = server)
