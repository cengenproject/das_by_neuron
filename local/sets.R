input <- list(sets_selected_neursA = "ASK,ADF,ASG,ASEL,ASER,AFD,BAG,AWA,AIY,AWB",
               sets_selected_neursB = "ADL,AIN,ASI,AVA,AVE,AVG,AVH,AVK,AVM,AWC,CAN,DA,I5,IL1,IL2,NSM,OLL,OLQ,PHA,PVC,PVD,PVM,RIA,RIC,RIM,RIS,RMD,SMD,VB,VC,DD,VD",
               sets_selected_fdr = 0.05,
               sets_selected_deltapsi = 0.5)

input <- list(sets_selected_neursA = "OLL,OLQ",
              sets_selected_neursB = "AVM,PVM",
              sets_selected_fdr = 0.05,
              sets_selected_deltapsi = 0.5)

source("R/utils.R")

#~ check neur names ----
r_sets_selected_neursA <-  input$sets_selected_neursA |>
    toupper() |>
    split_text_to_vector() |>
    validate_neurons(neurons_table)
  


r_sets_selected_neursB <- input$sets_selected_neursB |>
    toupper() |>
    split_text_to_vector() |>
    validate_neurons(neurons_table)




#~ get DAS events ----
r_sets_psis <- {
    
    message("   Set: get psis subset")
    
  neurs_selected <- union(r_sets_selected_neursA,
                          r_sets_selected_neursB)
  
    sub <- tbl_psidb |>
      filter( neur %in% neurs_selected ) |>
      mutate(set = if_else(neur %in% r_sets_selected_neursA,
                           "A", "B"),
             junction_name = paste0(event_name, "-", junction_id))
    
    if(VERBOSE) message("got sub: ", length(sub |> pull(junction_name)))
    
    sub
  }



r_sets_das_events <- {
    
    message("   Set: das events")
    
    junction_names_A <- r_sets_psis |>
      filter(set == "A") |>
      arrange(junction_name) |>
      pull(junction_name)
    
    psi_setA <- r_sets_psis |>
      filter(set == "A") |>
      arrange(junction_name) |>
      pull(mean_psi_per_lsv_junction) |>
      split(f = junction_names_A)
    
    
    junction_names_B <- r_sets_psis |>
      filter(set == "B") |>
      arrange(junction_name) |>
      pull(junction_name)
    
    psi_setB <- r_sets_psis |>
      filter(set == "B") |>
      arrange(junction_name) |>
      pull(mean_psi_per_lsv_junction) |>
      split(f = junction_names_B)
    
    # first filter: need to be measured in both groups
    group_names <- intersect(names(psi_setA), names(psi_setB))
    
    psi_setA <- psi_setA[group_names]
    psi_setB <- psi_setB[group_names]
    
    # second filter: need to have several measures in both groups
    nb_values <- pmin(sapply(psi_setA, length),
                      sapply(psi_setB, length))
    
    group_names <- group_names[nb_values > 1]
    psi_setA <- psi_setA[group_names]
    psi_setB <- psi_setB[group_names]
    
    
    n <- length(psi_setA)
    
    if(VERBOSE) message("Nb of tests: ", n)
    
    all_pvals <- numeric(n) |> setNames(group_names)
    all_mean_deltapsi <- numeric(n) |> setNames(group_names)
    for(i in seq_len(n)){
      all_pvals[[i]] <- tryCatch(t.test(psi_setA[[i]], psi_setB[[i]])[["p.value"]],
                                 error = \(x) NA_real_)
      all_mean_deltapsi[[i]] <- mean(psi_setA[[i]]) - mean(psi_setB[[i]])
    }
    
    # hist(all_pvals, breaks = 50)
    # plot(all_mean_deltapsi, -log10(all_fdr)); abline(h = -log10(0.05))
    # hist(all_fdr, breaks = 50)
    # table(all_fdr < .05)
    
    all_fdr <- p.adjust(all_pvals, method = "BH") |> setNames(group_names)
    
    
    
    keep <- which(
      all_fdr <= input$sets_selected_fdr &
        abs(all_mean_deltapsi) >= input$sets_selected_deltapsi
    )
    
    if(VERBOSE) message("filtered: ", length(keep)," pass")
    
    jcts_to_keep <- names(keep)
    
    set_of_das_evs <- r_sets_psis |>
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
    
  }


r_sets_lsvs_das <-  r_sets_das_events |> pull(lsv_id) |> unique() 



#~ table ----
r_sets_table_das_sjs <- {
    
    message("   Set: table")
    
    tab <- r_sets_das_events |>
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
  }




#~ plot ----
r_sets_gg_das_genes <- {
    
    message("   Set: plot")
    
    if(length(r_sets_lsvs_das) > LIMIT_NB_EVENTS_TO_PLOT){
      return(invisible(NULL))
    }
    
    if(VERBOSE) message("  plotting!")
    
    toplot <- r_sets_psis |>
      filter(lsv_id %in% !!( r_sets_lsvs_das )) |>
      select(event_name, gene_name, junction_id, neur, mean_psi_per_lsv_junction, set) |>
      as_tibble() |>
      rename(Neuron = neur, PSI = mean_psi_per_lsv_junction) |>
      arrange(set, Neuron, gene_name, event_name, junction_id) |>
      mutate(event_name_annot = paste0(gene_name, " - ", event_name),
             Neuron = forcats::fct_inorder(Neuron)) |>
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
    
  }

