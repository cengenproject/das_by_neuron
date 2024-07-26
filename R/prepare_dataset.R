# 
# # Only run for new dataset, results saved as duckdb for the app
# 
# # Note: takes a while to run
# #           dpsi read file: ~ 1 min
# #           dpsi split rows: maybe 20 min?
# #           psi: a few seconds
# 
# 
# library(tidyverse)
# library(wbData)
# 
# 
# gids <- wb_load_gene_ids(289)
# 
# 
# #~ deltapsi ----
# data_dir <- "data_raw/2024-03-04_outs/deltapsi/"
# 
# 
# files_dpsi <- list.files(data_dir,
#                          pattern = "\\.tsv$",
#                          full.names = FALSE)
# 
# dpsidta <- map_dfr(files_dpsi,
#                    ~read_tsv(file.path(data_dir, .x),
#                              col_names = c("gene_id", "lsv_id", "lsv_type",
#                                            "dpsi", "p20", "p05", "psiA", "psiB",
#                                            "nb_sj", "nb_exons",
#                                            "sj_coords", "ir_coords"),
#                              col_types = cols(
#                                gene_id = col_character(),
#                                lsv_id = col_character(),
#                                lsv_type = col_character(),
#                                dpsi = col_character(),
#                                p20 = col_character(),
#                                p05 = col_character(),
#                                psiA = col_character(),
#                                psiB = col_character(),
#                                nb_sj = col_double(),
#                                nb_exons = col_double(),
#                                sj_coords = col_character(),
#                                ir_coords = col_character()
#                              ),
#                              skip = 1,
#                              na = "na",
#                              progress = FALSE) |>
#                      add_column(neurA = str_split(.x,"[-\\.]")[[1]][1],
#                                 neurB = str_split(.x,"[-\\.]")[[1]][2]),
#                    .progress = TRUE)
# 
# 
# 
# # separate data in E(PSI) and Std(PSI) fields
# # We get one row per junction
# dpsi <- dpsidta |>
#   separate_rows(dpsi, p20, p05, psiA, psiB,
#                 sep = ";") |>
#   group_by(neurA, neurB, lsv_id) |>
#   mutate(junction_id = row_number()) |>
#   ungroup()  |>
#   mutate(across(c(dpsi,p20, p05, psiA, psiB), as.double),
#          junction_id = factor(junction_id))
# 
# 
# dpsi <- dpsi |>
#   mutate(gene_name = i2s(gene_id, gids, warn_missing = TRUE))
# 
# # save as DuckDB
# con <- DBI::dbConnect(duckdb::duckdb(), "data_use/240726_dpsi.duckdb")
# DBI::dbWriteTable(con, name = "dpsi",value = dpsi)
# 
# DBI::dbDisconnect(con)
# 
# 
# 
# #~ Load PSI ----
# data_dir <- "data_raw/2024-03-04_outs/psi/"
# 
# files_psi <- list.files(data_dir,
#                         pattern = "\\.tsv$",
#                         full.names = FALSE)
# 
# 
# psidta <- map_dfr(files_psi,
#                   ~read_tsv(file.path(data_dir, .x),
#                             col_names = c("gene_id", "lsv_id", "lsv_type",
#                                           "mean_psi_per_lsv_junction", "stdev_psi_per_lsv_junction",
#                                           "nb_sj", "nb_exons",
#                                           "sj_coords", "ir_coords"),
#                             col_types = cols(
#                               gene_id = col_character(),
#                               lsv_id = col_character(),
#                               lsv_type = col_character(),
#                               mean_psi_per_lsv_junction = col_character(),
#                               stdev_psi_per_lsv_junction = col_character(),
#                               nb_sj = col_double(),
#                               nb_exons = col_double(),
#                               sj_coords = col_character(),
#                               ir_coords = col_character()
#                             ),
#                             skip = 1,
#                             na = "na",
#                             progress = FALSE) |>
#                     add_column(neuron = str_split_1(.x, "\\.")[[1]]),
#                   .progress = TRUE)
# 
# psi <- psidta |>
#   separate_rows(mean_psi_per_lsv_junction, stdev_psi_per_lsv_junction,
#                 sep = ";") |>
#   group_by(neuron, lsv_id) |>
#   mutate(junction_id = row_number()) |>
#   ungroup()  |>
#   mutate(across(c(mean_psi_per_lsv_junction,stdev_psi_per_lsv_junction), as.double),
#          junction_id = factor(junction_id)) |>
#   mutate(gene_name = i2s(gene_id, gids, warn_missing = TRUE))
# 
# 
# # save as DuckDB
# con <- DBI::dbConnect(duckdb::duckdb(), "data_use/240726_psi.duckdb")
# DBI::dbWriteTable(con, name = "psi",value = psi)
# 
# DBI::dbDisconnect(con)
# 
# 
# 
