# Exercise the actual app reactive definitions, without downloads or a viewer.
stopifnot(requireNamespace('shiny', quietly=TRUE))
library(shiny)
library(magrittr)
source('shinyRam/R/inspection.R')
expressions <- parse('shinyRam/app.R')
assignment <- function(x, name) is.call(x) && identical(x[[1L]], as.name('<-')) &&
  identical(x[[2L]], as.name(name))
server <- Filter(function(x) assignment(x,'server'),as.list(expressions))[[1L]][[3L]]
statements <- as.list(server[[3L]])[-1L]
definition <- function(name) Filter(function(x) assignment(x,name),statements)[[1L]]
counts <- new.env(); counts$plot <- 0L; counts$sequence <- 0L
fixture <- data.frame(chain=c('A','B'),resn=c('ALA','GLY'),resi=1:2,
  insertion_code='',phi=c(-60,-70),psi=c(-40,30),region='Favoured')
shiny::testServer(function(input,output,session) {
  classified <- reactive(fixture)
  scientific_classified <- classified
  eval(definition("residue_key"), envir = environment())
  loaded <- reactive(list(name='fixture',chains=c('A','B')))
  plot_reference <- reactive(list(z=matrix(1,2,2)))
  active_palette <- reactive(c('white','blue','green','red'))
  current_chain_colors <- reactive(c('#000000','#ffffff'))
  ram_density_thresholds <- function(x) c(.1,.2,.3)
  for (name in c('plot_residues','displayed','plot_limits','plot_state','sequence_groups'))
    eval(definition(name),envir=environment())
  observe({ plot_state(); counts$plot <- counts$plot+1L })
  observe({ sequence_groups(); counts$sequence <- counts$sequence+1L })
}, {
  session$setInputs(AA=c('ALA','GLY'),background='General',chainselection=c('A','B'))
  session$elapse(150)
  before <- c(counts$plot,counts$sequence)
  stopifnot(nrow(displayed())==2L)
  session$setInputs(chainselection='B'); session$elapse(150)
  stopifnot(identical(displayed()$chain,'B'),identical(before,c(counts$plot,counts$sequence)))
  session$setInputs(chainselection=character()); session$elapse(150)
  stopifnot(nrow(displayed())==0L,identical(before,c(counts$plot,counts$sequence)))
})
message('Actual app reactives: chain filtering updates data with zero plot or sequence invalidations.')
