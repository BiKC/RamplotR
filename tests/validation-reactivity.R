library(shiny)
source('shinyRam/R/ramachandran.R'); source('shinyRam/R/rama8000.R')
source('shinyRam/R/predictions.R'); source('shinyRam/R/geometry.R')
source('shinyRam/R/canonical.R')
expressions <- parse('shinyRam/app.R')
assignment <- function(x,name) is.call(x)&&identical(x[[1]],as.name('<-'))&&identical(x[[2]],as.name(name))
server <- Filter(function(x)assignment(x,'server'),as.list(expressions))[[1]][[3]]
statements <- as.list(server[[3]])[-1]
definition <- function(name) Filter(function(x)assignment(x,name),statements)[[1]]
fixture <- data.frame(chain=c('A','B'),resi=1:2,insertion_code='',resn=c('ALA','GLY'),
 phi=c(-60,-70),psi=c(-40,30),omega_prev=NA_real_,next_resn=NA_character_,bonded_to_next=FALSE)
count <- new.env();count$standard <- 0L
classify_standard <- ram_rama8000_classify
ram_rama8000_classify <- function(...) { count$standard<-count$standard+1L;classify_standard(...) }
shiny::testServer(function(input,output,session) {
 loaded <- reactive(list(key='fixture',prediction=NULL))
 model_torsions <- reactive(fixture)
 model_geometry <- reactive(data.frame(omega=c(180,180)))
 current_model <- reactive(1L)
 plot_reference <- reactive(ram_read_reference(file.path('shinyRam/static/original',input$background)))
 external_validation <- reactive(NULL); canonical_mapping <- reactive(NULL)
 # Use the actual reference loader with a path relative to the test working dir.
 classify_native <- ram_classify_torsions
 ram_classify_torsions <- function(...,reference_dir) classify_native(...,reference_dir=file.path('shinyRam',reference_dir))
 ram_rama8000_classify <- function(x,directory) {
  count$standard<-count$standard+1L;classify_standard(x,file.path('shinyRam',directory))
 }
 for(name in c('native_classified','standard_classified','scientific_classified','classified'))
  eval(definition(name),envir=environment())
}, {
 session$setInputs(validationMode='residue',bgtype='original',background='General')
 original <- scientific_classified()
 expected <- classify_native(fixture,'shinyRam/static/original',NULL,'residue',ram_density_thresholds)
 expected <- classify_standard(expected,'shinyRam/static/rama8000')
 stopifnot(identical(original$region,expected$region),identical(original$rama8000_score,expected$rama8000_score))
 count_before <- count$standard
 session$setInputs(background='GLY')
 stopifnot(identical(original,scientific_classified()),count$standard==count_before)
 session$setInputs(validationMode='legacy'); scientific_classified()
 stopifnot(count$standard==count_before)
 session$setInputs(validationMode='residue'); scientific_classified()
 stopifnot(count$standard==count_before)
})
message('Validation caches preserve scores and reuse standard validation across display changes.')
