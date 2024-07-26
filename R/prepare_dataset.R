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
#   separate_rows(dpsi, p20, p05, psiA, psiB, sj_coords,
#                 sep = ";") |>
#   group_by(neurA, neurB, lsv_id) |>
#   mutate(junction_id = row_number()) |>
#   ungroup()  |>
#   mutate(across(c(dpsi,p20, p05, psiA, psiB), as.double),
#          junction_id = factor(junction_id)) |>
#   mutate(gene_name = i2s(gene_id, gids, warn_missing = TRUE),
#          .after = gene_id)
# 
# 
# 
# 
# # set.seed(123)
# # table_event_names <- tibble(lsv_id = unique(dpsi$lsv_id),
# #                             event_name = babynames::babynames$name |> unique() |> sample(length(unique(dpsi$lsv_id))))
# # write_csv(table_event_names, "data_use/240726_table_event_names.csv")
# 
# table_event_names <- read_csv("data_use/240726_table_event_names.csv")
# dpsi <- dpsi |>
#   left_join(table_event_names,
#             by = "lsv_id")
# 
# # save as DuckDB
# con <- DBI::dbConnect(duckdb::duckdb(), "data_use/240726_dpsi.duckdb")
# DBI::dbWriteTable(con, name = "dpsi",value = dpsi)
# 
# DBI::dbDisconnect(con)

