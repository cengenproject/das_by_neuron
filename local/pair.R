input <- list(pair_selected_neurA = "PVM",
              pair_selected_neurB = "NSM",
              selected_p20 = 0.5,
              selected_p05 = 0.05,
              selected_dpsi = 0.5)

source("R/utils.R")




r_pair_selected_neurA <- {

  neur <- input$pair_selected_neurA |>
    toupper() |>
    split_text_to_vector() |>
    validate_neurons(neurons_table)
  
  is_valid <- (length(neur) > 0)
  
  stopifnot(is_valid)
  
  if(VERBOSE) message("   Pair neur A ", neur)
  neur
}

r_pair_selected_neurB <- {
 
  neur <- input$pair_selected_neurB |>
    toupper() |>
    split_text_to_vector() |>
    validate_neurons(neurons_table)
  
  is_valid <- (length(neur) > 0)
  stopifnot(is_valid)
  
  if(VERBOSE) message("   Pair neur B ", neur)
  neur
}



#~ get DAS events ----
r_das_events <- {
    
    message("   Pair: get DAS events")
    tbl_dpsidb |>
      filter(
        (
          (neurA %in% r_pair_selected_neurA ) & 
            (neurB %in% r_pair_selected_neurB )
        ) | (
          (neurA %in% r_pair_selected_neurA ) & 
            (neurB %in% r_pair_selected_neurB )
        ),
        p20 >= input$selected_p20,
        p05 <= input$selected_p05,
        abs(dpsi) >= input$selected_dpsi
      ) |>
      collect()
  }

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