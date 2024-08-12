
# Packages ----
library(shiny)
library(dplyr, warn.conflicts = FALSE)
library(dbplyr, warn.conflicts = FALSE)
library(tidyr)

library(DT, warn.conflicts = FALSE)

library(ggplot2)
library(plotly, warn.conflicts = FALSE)

library(shinyFeedback)

# Data ----
dpsidb <- DBI::dbConnect(duckdb::duckdb(),
                         dbdir = "data_use/240726_dpsi.duckdb",
                         read_only = TRUE)
onStop(function() {
  DBI::dbDisconnect(dpsidb)
})

tbl_dpsidb <- tbl(dpsidb, "dpsi")


# hardcoded parameters
LIMIT_NB_EVENTS_TO_PLOT <- 100
VERBOSE <- 0


measured_neurons <- readLines("data_use/measured_neurons.txt")


neurons_table <- readr::read_csv("data_use/neuron_properties.csv",
                                 col_types = "cccc") %>%
  filter(Neuron_type %in% measured_neurons) %>%
  mutate(across(Modality:Neurotransmitter, stringr::str_to_lower))
readr::stop_for_problems(neurons_table)


# UI ----
ui <- fluidPage(
  
  if(file.exists("google-analytics.html")){
    tags$head(includeHTML(("google-analytics.html")))
  },
  
  shinyFeedback::useShinyFeedback(),
  
  titlePanel("Local quantification analyses"),
  
  tabsetPanel(
    
    #~ Pairs ----
    tabPanel(
      "Pair of neurons",
      fluidRow(
        textInput(inputId = "pair_selected_neurA",
                  label = "First neuron:",
                  value = "OLL"),
        textInput(inputId = "pair_selected_neurB",
                  label = "Second neuron:",
                  value = "OLQ"),
        
        numericInput("selected_p20", "Probability that the event is DAS",
                     value = 0.5, min = 0, max = 1, step = .1),
        
        numericInput("selected_p05", "Probability that the event is not DAS",
                     value = 0.05, min = 0, max = 1, step = .1),
        
        numericInput("selected_dpsi", "Minimal deltaPSI",
                     value = 0.5, min = 0, max = 1, step = .1),
        
        actionButton("pair_submit",
                     "Run this pair of neurons")
      ),
      
      fluidRow(
        
        div(
          htmlOutput("pair_summary",
                     style = "width:600px; border: 3px dashed grey; padding: 50px; margin: 20px;")
        ),
        
        dataTableOutput("pair_table_das_sjs",
                        width = "80%"),
        
        plotly::plotlyOutput("pair_gg_das_genes",
                             width = "60%")
        
      )
    ),
    
    #~ Sets ----
    tabPanel(
      "Sets of neurons",
      fluidRow(
        textInput(inputId = "sets_selected_neursA",
                  label = "First set of neuron(s):",
                  value = "OLL,OLQ"),
        textInput(inputId = "sets_selected_neursB",
                  label = "Second set of neuron(s):",
                  value = "AVM,PVM"),
        
        numericInput("sets_selected_fdr", "FDR threshold",
                     value = 0.05, min = 0, max = 1, step = .1),
        
        numericInput("sets_selected_deltapsi", "DeltaPSI threshold",
                     value = 0.5, min = 0, max = 1, step = .1),
        
        
        actionButton("sets_submit",
                     "Run these sets of neurons")
      ),
      
      fluidRow(
        
        div(
          htmlOutput("sets_summary"),
          style = "width:600px; border: 3px dashed grey; padding: 50px; margin: 20px;"
        ),
        
        dataTableOutput("sets_table_das_sjs",
                        width = "80%"),
        
        plotly::plotlyOutput("sets_gg_das_genes",
                             width = "60%")
        
      )
    )
    
  ),
  # Footer
  hr(),
  tags$footer(
    "See app documentation at",
    tags$a("splicing.cengen.org", href="http://splicing.cengen.org"),
    " and source code ",
    tags$a("on Github", href="https://github.com/cengenproject/das_by_neuron")
  )
)




# SERVER ----
server <- function(input, output) {
  
  
  
  # ** Pairs ** ----
  
  #~ check neur names ----
  r_pair_selected_neurA <- reactive({
    showNotification("Processing neuron...", duration = 2)
    
    neur <- input$pair_selected_neurA |>
      toupper() |>
      split_text_to_vector() |>
      validate_neurons(neurons_table)
    
    is_valid <- (length(neur) > 0)
    shinyFeedback::feedbackDanger("pair_selected_neurA", !is_valid, "Invalid neuron name")
    req(is_valid, cancelOutput = TRUE)
    
    if(VERBOSE) message("   Pair neur A ", neur)
    neur
  })
  
  r_pair_selected_neurB <- reactive({
    showNotification("Processing neuron...", duration = 2)
    
    neur <- input$pair_selected_neurB |>
      toupper() |>
      split_text_to_vector() |>
      validate_neurons(neurons_table)
    
    is_valid <- (length(neur) > 0)
    shinyFeedback::feedbackDanger("pair_selected_neurB", !is_valid, "Invalid neuron name")
    req(is_valid, cancelOutput = TRUE)
    
    if(VERBOSE) message("   Pair neur B ", neur)
    neur
  })
  
  
  
  #~ get DAS events ----
  r_das_events <- eventReactive(
    eventExpr = input$pair_submit,
    valueExpr = {
      
      message("   Pair: get DAS events")
      tbl_dpsidb |>
        filter(
          (
            (neurA %in% !!r_pair_selected_neurA() ) & 
              (neurB %in% !!r_pair_selected_neurB() )
          ) | (
            (neurA %in% !!r_pair_selected_neurB() ) & 
              (neurB %in% !!r_pair_selected_neurA() )
          ),
          p20 >= input$selected_p20,
          p05 <= input$selected_p05,
          abs(dpsi) >= input$selected_dpsi
        ) |>
        collect()
    })
  
  r_lsvs_das <- reactive( r_das_events() |> pull(lsv_id) |> unique() )
  
  
  #~ text summary ----
  r_pair_summary <- eventReactive(
    eventExpr = input$pair_submit,
    valueExpr = {
      
      
      message("   Pair: text summary")
      
      nb_events_das <- r_das_events() |>
        pull(lsv_id) |>
        unique() |>
        length()
      
      nb_genes_das <- r_das_events() |>
        pull(gene_id) |>
        unique() |>
        length()
      
      text_experiment <- "Using statistics from MAJIQ."
      
      text_res <- paste0("Found ", nb_events_das, " DAS events in ", nb_genes_das, " genes.")
      
      
      
      if(length(r_pair_selected_neurA()) == 0 |
         length(r_pair_selected_neurB()) == 0 ){
        
        
        text_err <- "Neuron not recognized."
        text <- list(p(text_err, style = "color:red"))
        return(text)
      }
      
      
      if(length(r_lsvs_das()) > LIMIT_NB_EVENTS_TO_PLOT){
        
        text_err <- paste0("Too many DAS events (",length(r_lsvs_das()),"), will not plot.
                                   Consider using a more stringent filter (e.g. increase minimal deltaPSI).")
        
      } else{
        text_err <- ""
      }
      
      list(p(text_experiment),
           p(text_res),
           p(text_err, style = "color:red"))
    })
  
  output$pair_summary <- renderUI(r_pair_summary())
  
  
  #~ table ----
  r_pair_table_das_sjs <- eventReactive(
    eventExpr = input$pair_submit,
    valueExpr = {
      
      message("   Pair: get table")
      
      r_das_events() |>
        mutate(junction_name = paste0(event_name, "-", junction_id)) |>
        mutate(gene_id = lapply(gene_id, \(gid){
          paste0('<a href="http://splicingapps.cengen.org/voila/gene/',gid,'/">',gid,'</a>')
        } )) |>
        select(`Gene Name` = gene_name, `Junction Name`=junction_name,
               `Gene ID` = gene_id, `Event ID` = lsv_id, `LSV type` = lsv_type,
               `Junction Coordinates` = sj_coords,
               `DeltaPSI` = dpsi,
               `Probability non DAS` = p05, `Probability DAS` = p20)
      
    })
  
  
  output$pair_table_das_sjs <- renderDT(r_pair_table_das_sjs())
  
  #~ plot ----
  r_pair_gg_das_genes <- eventReactive(
    eventExpr = input$pair_submit,
    valueExpr = {
      
      message("   Pair: plot")
      
      if(length(r_lsvs_das()) > LIMIT_NB_EVENTS_TO_PLOT){
        return(invisible(NULL))
      }
      
      toplot <- tbl_dpsidb |>
        filter((neurA %in% !!r_pair_selected_neurA() & neurB %in% !!r_pair_selected_neurB() ) |
                 (neurA %in% !!r_pair_selected_neurB() & neurB %in% !!r_pair_selected_neurA() ),
               lsv_id %in% !!(r_lsvs_das())) |>
        select(event_name, gene_name, junction_id, neurA,neurB,psiA,psiB) |>
        as_tibble() |>
        pivot_longer(-c(event_name, gene_name, junction_id),
                     names_to = c(".value",NA),
                     names_pattern = "(neur|psi)(A|B)") |>
        rename(Neuron = neur, PSI = psi) |>
        arrange(gene_name, event_name, junction_id) |>
        mutate(event_name_annot = paste0(gene_name, " - ", event_name))
      
      
      
      gg_das_genes <- toplot |>
        ggplot() +
        theme_minimal() +
        theme(legend.position = 'none') +
        theme(strip.text.y = element_text(angle = 0),
              strip.clip = "off") +
        ylab(NULL) +
        facet_grid(rows = vars(event_name_annot)) +
        geom_col(aes(x = PSI, y = Neuron, fill = junction_id),
                 position = position_stack())
      
      
      ggplotly(gg_das_genes,
               width = 640,
               height = 15*nrow(toplot)) |>
        layout(margin = list(r=200))
      
    })
  
  output$pair_gg_das_genes <- renderPlotly(r_pair_gg_das_genes())
  
  
  
  # ** Sets ** ----
  
  
  
  #~ check neur names ----
  r_sets_selected_neursA <- reactive({
    showNotification("Processing neurons...", duration = 2)
    
    neurs <- input$sets_selected_neursA |>
      toupper() |>
      split_text_to_vector() |>
      validate_neurons(neurons_table)
    
    is_valid <- (length(neurs) > 0)
    shinyFeedback::feedbackDanger("sets_selected_neursA", !is_valid, "No valid neuron name")
    req(is_valid, cancelOutput = TRUE)
    
    if(VERBOSE) message("   Set A: ", neurs)
    neurs
  })
  
  
  r_sets_selected_neursB <- reactive({
    showNotification("Processing neurons...", duration = 2)
    
    neurs <- input$sets_selected_neursB |>
      toupper() |>
      split_text_to_vector() |>
      validate_neurons(neurons_table)
    
    is_valid <- (length(neurs) > 0)
    shinyFeedback::feedbackDanger("sets_selected_neursB", !is_valid, "No valid neuron name")
    req(is_valid, cancelOutput = TRUE)
    
    if(VERBOSE) message("   Set B: ", neurs)
    neurs
  })
  
  
  #~ get DAS events ----
  r_sets_psis <- eventReactive(
    eventExpr = input$sets_submit,
    valueExpr = {
      
      message("   Set: get psis subset")
      
      sub <- tbl_dpsidb |>
        filter((neurA %in% !!r_sets_selected_neursA() & neurB %in% !!r_sets_selected_neursB()) |
                 (neurA %in% !!r_sets_selected_neursB() & neurB %in% !!r_sets_selected_neursA())) |>
        mutate(psi_setA = if_else(neurA %in% !!r_sets_selected_neursA(),
                                  psiA, psiB),
               psi_setB = if_else(neurA %in% !!r_sets_selected_neursB(),
                                  psiA, psiB),
               junction_name = paste0(event_name, "-", junction_id))
      
      if(VERBOSE) message("got sub: ", length(sub |> pull(junction_name)))
      
      sub
    })
  
  
  
  
  
  
  r_sets_das_events <- eventReactive(
    eventExpr = input$sets_submit,
    valueExpr = {
      
      message("   Set: das events")
      
      
      junction_names <- r_sets_psis() |> pull(junction_name)
      sets_psiA <- r_sets_psis() |> pull(psi_setA) |> split(f = junction_names)
      sets_psiB <- r_sets_psis() |> pull(psi_setB) |> split(f = junction_names)
      group_names <- names(sets_psiA)
      
      n <- length(sets_psiA)
      stopifnot(length(sets_psiB) == n)
      stopifnot(length(unique(junction_names)) == n)
      
      if(VERBOSE) message("Nb of psis: ", length(junction_names))
      if(VERBOSE) message("Nb of tests: ", n)
      
      all_pvals <- numeric(n) |> setNames(group_names)
      all_mean_deltapsi <- numeric(n) |> setNames(group_names)
      for(i in seq_len(n)){
        all_pvals[[i]] <- tryCatch(t.test(sets_psiA[[i]], sets_psiB[[i]])[["p.value"]],
                                   error = \(x) NA_real_)
        all_mean_deltapsi[[i]] <- mean(sets_psiA[[i]]) - mean(sets_psiB[[i]])
      }
      
      all_fdr <- p.adjust(all_pvals, method = "BH") |> setNames(group_names)
      
      keep <- which(
        all_fdr <= input$sets_selected_fdr &
          abs(all_mean_deltapsi) >= input$sets_selected_deltapsi
      )
      
      if(VERBOSE) message("filtered: ", length(keep)," pass")
      
      jcts_to_keep <- names(keep)
      
      set_of_das_evs <- r_sets_psis() |>
        filter(junction_name %in% jcts_to_keep) |>
        select(gene_id, gene_name, junction_name, lsv_id, lsv_type,
               sj_coords) |>
        distinct() |>
        collect()
      
      
      set_of_das_evs$p_t <- all_pvals[ set_of_das_evs$junction_name ]
      set_of_das_evs$fdr <- all_fdr[ set_of_das_evs$junction_name ]
      set_of_das_evs$mean_deltapsi <- all_mean_deltapsi[ set_of_das_evs$junction_name ]
      
      
      if(VERBOSE) message("got it")
      set_of_das_evs
      
    })
  
  
  r_sets_lsvs_das <- reactive( r_sets_das_events() |> pull(lsv_id) |> unique() )
  
  #~ text summary ----
  r_sets_summary <- eventReactive(
    eventExpr = input$sets_submit,
    valueExpr = {
      
      message("   Set: text summary")
      
      nb_events_das <- r_sets_das_events() |>
        pull(lsv_id) |>
        unique() |>
        length()
      
      nb_genes_das <- r_sets_das_events() |>
        pull(gene_id) |>
        unique() |>
        length()
      
      text_experiment <- paste("Performing a t-test between sets of PSI,",
                               "followed by Benjamini-Hochberg FDR estimate.",
                               "Please note: depending on the sets of neurons,",
                               "the test's assumptions may not be fullfilled,",
                               "and the FDR may not be controlled.",
                               "These results are provided as a help to select",
                               "interesting events, but not are not appropriate",
                               "for comparison of sets of DAS genes.")
      
      text_res <- paste0("Found ", nb_events_das, " DAS events in ", nb_genes_das, " genes.")
      
      
      
      if(length(r_sets_selected_neursA()) == 0 |
         length(r_sets_selected_neursB()) == 0 ){
        
        
        text_err <- "No valid neuron name."
        text <- list(p(text_err, style = "color:red"))
        return(text)
      }
      
      
      if(length(r_sets_lsvs_das()) > LIMIT_NB_EVENTS_TO_PLOT){
        
        text_err <- paste0("Too many DAS events (",length(r_sets_lsvs_das()),"), will not plot.
                                   Consider using a more stringent filter (e.g. increase minimal deltaPSI).")
        
      } else{
        text_err <- ""
      }
      
      if(VERBOSE) message("Collected text")
      
      list(p(text_experiment),
           p(text_res),
           p(text_err, style = "color:red"))
    })
  
  output$sets_summary <- renderUI(r_sets_summary())
  
  
  #~ table ----
  r_sets_table_das_sjs <- eventReactive(
    eventExpr = input$sets_submit,
    valueExpr = {
      
      message("   Set: table")
      
      tab <- r_sets_das_events() |>
        mutate(gene_id = lapply(gene_id, \(gid){
          paste0('<a href="http://splicingapps.cengen.org/voila/gene/',gid,'/">',gid,'</a>')
        } )) |>
        mutate(p_t = round(p_t, 3),
               fdr = round(fdr, 3),
               mean_deltapsi = round(mean_deltapsi, 2)) |>
        select(`Gene Name` = gene_name, `Junction Name`=junction_name,
               `Gene ID` = gene_id, `Event ID` = lsv_id, `LSV type` = lsv_type,
               `Junction Coordinates` = sj_coords,
               `Mean DeltaPSI` = mean_deltapsi,
               `p value` = p_t, `FDR` = fdr)
      
      if(VERBOSE) message("Collected table, ", nrow(tab))
      
      tab
    })
  
  
  output$sets_table_das_sjs <- renderDT(r_sets_table_das_sjs())
  
  #~ plot ----
  r_sets_gg_das_genes <- eventReactive(
    eventExpr = input$sets_submit,
    valueExpr = {
      
      message("   Set: plot")
      
      if(length(r_sets_lsvs_das()) > LIMIT_NB_EVENTS_TO_PLOT){
        return(invisible(NULL))
      }
      
      if(VERBOSE) message("  plotting!")
      
      toplot <- r_sets_psis() |>
        filter(lsv_id %in% !!( r_sets_lsvs_das() )) |>
        select(event_name, gene_name, junction_id, neurA,neurB,psiA,psiB) |>
        as_tibble() |>
        pivot_longer(-c(event_name, gene_name, junction_id),
                     names_to = c(".value",NA),
                     names_pattern = "(neur|psi)(A|B)") |>
        rename(Neuron = neur, PSI = psi) |>
        arrange(gene_name, event_name, junction_id) |>
        mutate(event_name_annot = paste0(gene_name, " - ", event_name)) |>
        distinct()
      
      if(VERBOSE) message("Collected toplot: ", length(toplot |> pull(event_name_annot)))
      
      
      gg_das_genes <- toplot |>
        ggplot() +
        theme_minimal() +
        theme(legend.position = 'none') +
        theme(strip.text.y = element_text(angle = 0),
              strip.clip = "off") +
        ylab(NULL) +
        facet_grid(rows = vars(event_name_annot)) +
        geom_col(aes(x = PSI, y = Neuron, fill = junction_id),
                 position = position_stack())
      
      if(VERBOSE) message("Generated ggplot")
      
      ggplotly(gg_das_genes,
               width = 640,
               height = 10*nrow(toplot)) |>
        layout(margin = list(r=200))
      
    })
  
  output$sets_gg_das_genes <- renderPlotly(r_sets_gg_das_genes())
  
  
}

# Run the application 
shinyApp(ui = ui, server = server)
