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

# Define UI for application that draws a histogram
ui <- fluidPage(
  # Application title
  titlePanel("RamplotR"),
  tags$head(
    tags$script(src = "https://cdn.rawgit.com/arose/ngl/v0.10.4-1/dist/ngl.js"),
    tags$script(src = "https://cdn.plot.ly/plotly-2.14.0.min.js")
  ),
  # Sidebar with a slider input for number of bins
  fluidPage(fluidRow(
    column(
      4,
      fluidRow(
        column(6,
          textInput(
            inputId = "PDB",
            label = "PDB-code:",
            value = "1BBB"
          ),
          actionButton("submit", "Apply Changes", class = "btn-primary btn-lg")
        ),
        column(6,
          fileInput("structfile",label = "Or upload a custom file")
        )
      ),
      
      hr(),
      h4("Background settings"),
      fluidRow(
        column(6,selectInput(
          "bgtype","Reference dataset for background density",
          c("original","alphafold","alphafold_filtered","astral2.08","custom_high_resolution")
        )),
        column(6,selectInput(
          "background", "Choose background plot",
          c("General", "Glycine", "Preproline", "Proline")
        ))
      ),
      dropdown(
        list(
          colourpicker::colourInput("bg1", "Not allowed region",
                                    value = "#F1EEF6"
          ),
          colourpicker::colourInput("bg2", "Generously allowed regions",
                                    value = "#BDC9E1"
          ),
          colourpicker::colourInput("bg3", "Allowed regions",
                                    value = "#74A9CF"
          ),
          colourpicker::colourInput("bg4", "Favoured regions",
                                    value = "#0570B0"
          ),
          selectInput("colorscheme","Default colorscheme", choices = c("Rampage","PDBSum","custom"),selected = "Rampage")
        ),
        label = "Background colors"
      ),
      hr(),
      fluidRow(
        column(6,
          pickerInput(
            "AA",
            "Amino acid selection",
            allAA,
            multiple = T,
            options = list(`actions-box` = TRUE),
            selected = allAA
          )
        ),
        column(6,
          uiOutput("chains")
        )
      ),
      
      hr(),
      tags$label("Chain colors",class="control-label"),
      uiOutput("chainColors"),
      hr(),
      h4("3D view options"),
      fluidRow(
        column(4,
               checkboxInput("ligands",label = "Show ligands"),
        ),
        column(4,
               checkboxInput("dna",label = "Show DNA"),
        ),
        column(4,
               checkboxInput("rna",label = "Show RNA"),
        )
      ),
      fluidRow(
        column(4,
               checkboxInput("spinning",label = "Spinning"),
        ),
        column(4,
               checkboxInput("rocking",label = "Rocking",value=T),
        )
      ),
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
            selectInput("regionselect", "Show residues for region", c("Not allowed", "Generously allowed", "Allowed", "Favoured", "All")),
            dataTableOutput("regions")
          ),
          tabPanel(
            "Summary statistics",
            htmlOutput("summary")
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

  output$dummy <- reactive({
    input$submit
    withProgress(message = "Making plot", value = 0, {
      inputType<-""
      isolate({
        if (is.null(input$structfile$datapath)){
          inputType<-"code"
          accPDB <- input$PDB
        }
        else {
          inputType<-"file"
          accPDB <- input$structfile$datapath
        }
        if (session$userData$previousPDB != accPDB) {
          #print("yup")
          incProgress(1 / 4, detail = paste("Fetching sequence"))
          if (inputType=="file") {
            pdb<-bio3d::read.pdb(accPDB)
          }
          else {
            pdb <- tryCatch(expr = {
              bio3d::read.cif(accPDB)
            }, error = function(e) {
              #print(e)
              bio3d::read.pdb(accPDB)
            })
          }
          incProgress(1 / 4, detail = paste("Transforming data"))
          pdb<<-pdb
          # Get the list of unique chains
          chains <- unique(pdb$atom[,"chain"])

          # Initialize an empty data frame to store the results
          torsion <- data.frame()

          # Loop over each chain
          for (chain in chains) {

            print(chain)
            # Subset the pdb data for the current chain
            pdb_chain <- trim.pdb(pdb, chain = chain)
            # Identify the residue numbers that correspond to the residue names in allAA
            resno_allAA <- unique(pdb_chain$atom[pdb_chain$atom[,"resid"] %in% allAA, "resno"])

            # Subset the pdb_chain object to only include those residue numbers
            pdb_chain <- trim.pdb(pdb_chain, resno = resno_allAA)
            
            # if there is no data for the current chain, skip it
            if (nrow(pdb_chain$atom) == 0) {
              next
            }
            # Calculate the torsion angles for the current chain
            torsion_chain <- torsion.pdb(pdb_chain)
            
            # Add the chain information to the torsion angles data frame
            tortab_chain <- torsion_chain[["tbl"]][, c("phi", "psi")]
            spltor_chain <- strsplit(rownames(tortab_chain), split = ".", fixed = T)
            torsion_chain <- cbind(tortab_chain, as.data.frame(do.call(rbind, spltor_chain)))
            torsion_chain <- rename(torsion_chain, c("V1" = "resi", "V2" = "chain", "V3" = "resn"))
            
            # Combine the results row wise
            torsion <- rbind(torsion, torsion_chain)
          }
          print(torsion)
          # Store the results in the user data
          session$userData$torsion <- torsion

          chains<-unique(session$userData$torsion$chain)

          nglview<-NGLVieweR(data = accPDB)%>%setRock()
          counter=1
          for (i in unique(chains)) {
            #print(i)
            nglview<-addRepresentation(NGLVieweR = nglview,type = "cartoon", param=list("sele"= paste0(":", i,"  and protein"), "color"= color_set[((counter - 1) %% length(color_set)) + 1]))
            counter=counter+1
          }
          output$NGL <- NGLVieweR::renderNGLVieweR(nglview)
          
          output$chainColors <- renderUI({
            isolate({
              x <- vector("list", length(chains))
              c <- 1
              for (i in chains) {
                col<-color_set[((c - 1) %% length(color_set)) + 1]
                x[[i]] <-
                  list(colourpicker::colourInput(
                    paste0("chain", i),
                    label = paste("Chain", i),
                    value = col
                  ))
                c <- c + 1
              }

              dropdown(x, label = "Chain color settings")
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
          
          
          
          session$userData$previousPDB <- accPDB
          incProgress(1 / 4, detail = paste("Filter data"))
        } else {
          incProgress(3 / 4, detail = paste("Filter data"))
        }
        if (input$background == "preProline"){
          matrix <- readRDS(paste0("static/",input$bgtype,"/preProline"))
          # get a subset of only those amino acids that precede a proline
          torsionsubset <- data.frame()
          for (i in 1:length(session$userData$torsion$resn) - 1) {
            if (session$userData$torsion$resn[i + 1] == "PRO") {
              torsionsubset <- rbind(torsionsubset, session$userData$torsion[i, ])
            }
          }
          torsionsubset <- subset(torsionsubset, resn %in% input$AA)
          # also subset for chains
          # since chainselection is added as uiOutput, it is not available in the beginning, so check if it exists, otherwise subset for all chains
          if (!is.null(input$chainselection)) {
            torsionsubset <- subset(torsionsubset, chain %in% input$chainselection)
          }
          
        }
        else if (!input$background %in% allAA) {
          matrix <- readRDS(paste0("static/",input$bgtype,"/General"))
          torsionsubset <- session$userData$torsion
          torsionsubset <- subset(torsionsubset, resn %in% input$AA)
          # also subset for chains
          if (!is.null(input$chainselection)) {
            torsionsubset <- subset(torsionsubset, chain %in% input$chainselection)
          }
        } else {
          matrix <- readRDS(paste0("static/",input$bgtype,"/",input$background))
          #updatePickerInput(session, "AA", selected = input$background)
          torsionsubset <- subset(session$userData$torsion, resn %in% input$AA)
          # also subset for chains
          print(input$chainselection)
          if (!is.null(input$chainselection)) {
            torsionsubset <- subset(torsionsubset, chain %in% input$chainselection)
          }
        }
        incProgress(1 / 4, detail = paste("Creating plot"))
        # session$sendCustomMessage("updateFig", input$PDB)
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

          # define a function to negate the %in% operator
          `%nin%` <- Negate(`%in%`)
          # exclude glycine and proline
          exclude <- c("GLY", "PRO")
          
          # get the number of end-residues
          # end-residues are the residues that are at the beginning or end of a chain, so they don't have a region (na) and they are not glycine or proline
          ttabsub2 <- ttab()
          ttabsub2 <- ttabsub2[is.na(ttabsub2$region), ]
          print(ttabsub2)
          end_count <- nrow(subset(ttabsub2, resn %nin% exclude & resn %in% allAA))

          # get the number of residues for each region
          # for the following, we do not include glycine and proline
          count_no_gly_pro <- nrow(subset(ttab(), resn %nin% exclude & resn %in% allAA)) - end_count
          fr_count <- nrow(subset(ttab(), region == "Favoured" & resn %nin% exclude & resn %in% allAA))
          fr_percent <- round(100 * fr_count / count_no_gly_pro, 2) # round to 2 decimals
          ar_count <- nrow(subset(ttab(), region == "Allowed" & resn %nin% exclude & resn %in% allAA))
          ar_percent <- round(100 * ar_count / count_no_gly_pro, 2)
          gar_count <- nrow(subset(ttab(), region == "Generously allowed" & resn %nin% exclude & resn %in% allAA))
          gar_percent <- round(100 * gar_count / count_no_gly_pro, 2)
          nar_count <- nrow(subset(ttab(), region == "Not allowed" & resn %nin% exclude & resn %in% allAA))
          nar_percent <- round(100 * nar_count / count_no_gly_pro, 2)
          total_count <- count_no_gly_pro
          total_percent <- round(100 * total_count / count_no_gly_pro, 2)

          
          # exclude glycine and proline
          ttabsub2 <- subset(ttabsub2, resn %nin% exclude & resn %in% allAA)
          end_count <- nrow(ttabsub2)
          # get the number of glycine and proline residues
          gly_count <- nrow(subset(ttab(), resn == "GLY"))
          pro_count <- nrow(subset(ttab(), resn == "PRO"))
          # get the total number of residues
          total_count2 <- nrow(subset(ttab(), resn %in% allAA))

          # create html output to display the statistics
        
            output$summary <- renderUI({
  HTML(paste0(
    "<div style='font-size:16px; line-height:1.6;'>",
    "<p><b>Statistics for the regions:</b></p>",
    "<table>",
    "<tr><th style='padding: 0 1em;'>Region</th><th style='padding: 0 1em;'>No. of residues</th><th style='padding: 0 1em;'>%</th></tr>",
    "<tr><td style='padding: 0 1em;'>Favoured regions:</td><td style='padding: 0 1em;'>", fr_count, "</td><td style='padding: 0 1em;'>(", fr_percent, "%)</td></tr>",
    "<tr><td style='padding: 0 1em;'>Allowed regions:</td><td style='padding: 0 1em;'>", ar_count, "</td><td style='padding: 0 1em;'>(", ar_percent, "%)</td></tr>",
    "<tr><td style='padding: 0 1em;'>Generously allowed regions:</td><td style='padding: 0 1em;'>", gar_count, "</td><td style='padding: 0 1em;'>(", gar_percent, "%)</td></tr>",
    "<tr><td style='padding: 0 1em;'>Not allowed regions:</td><td style='padding: 0 1em;'>", nar_count, "</td><td style='padding: 0 1em;'>(", nar_percent, "%)</td></tr>",
    "<tr><td style='padding: 0 1em;'>Non-glycine and non-proline residues:</td><td style='padding: 0 1em;'>", total_count, "</td><td style='padding: 0 1em;'>(", total_percent, "%)</td></tr>",
    "</table>",
    "<hr>",
    "<p><b>End-residues (Excl. Gly and Pro)</b></p>",
    "<table>",
    "<tr><td style='padding: 0 1em;'>Total no. of end-residues:</td><td style='padding: 0 1em;'>", end_count, "</td></tr>",
    "</table>",
    "<hr>",
    "<p><b>Glycine and proline residues</b></p>",
    "<table>",
    "<tr><td style='padding: 0 1em;'>Glycine residues:</td><td style='padding: 0 1em;'>", gly_count, "</td></tr>",
    "<tr><td style='padding: 0 1em;'>Proline residues:</td><td style='padding: 0 1em;'>", pro_count, "</td></tr>",
    "</table>",
    "<hr>",
    "<p><b>Total no. of residues</b></p>",
    "<table>",
    "<tr><td style='padding: 0 1em;'>Total no. of residues:</td><td style='padding: 0 1em;'>", total_count2, "</td></tr>",
    "</table>",
    "</div>"
  ))
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