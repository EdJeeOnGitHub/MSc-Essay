size_sim_func <- function(draw,
                          n_covariates = 4,
                          n_tiles = 2,
                          N = 1000){
  sim_df <- matrix(rnorm(N*n_covariates), N, n_covariates) %>% 
    as_tibble() %>% 
    mutate(D = rnorm(N),
           Z = rbinom(N, 1, 0.5))
  
  
  margins_output <- sim_df %>% 
    run_first_stage_interactions_fast(dataset = .,
                                      dependent_variable = "D",
                                      instrument = "Z",
                                      n_tiles = n_tiles)
  pval_output <- margins_output %>%
    select(pval_one_neg) %>% 
    distinct() %>% 
    pull()
  n_distinct_subgroup <- margins_output %>% 
    select(unique_subgroup) %>% 
    n_distinct()
  
  sim_output <- tibble(pvals = pval_output,
                       n_groups = n_distinct_subgroup,
                       draw = draw)
  return(sim_output)
}



######### Size Sims SFS


n_cov_list <- matrix(1:10, 200, 10) %>% 
  as.list()
# Covariates
run_sim <- FALSE
if (run_sim){
  plan(multisession)
  new_size_sims_cov <-  n_cov_list %>% 
    future_imap_dfr(~size_sim_func(draw = .y,
                                   n_covariates = .x) %>% 
                      mutate(n_cov = .x,
                             n_tiles = 2),
                    .options = future_options(globals = c("run_first_stage_interactions_fast",
                                                          "discretise_df",
                                                          "find_SEs",
                                                          "create_fake_data",
                                                          "simulation_func",
                                                          "power_test",
                                                          "sim_first_stage_forest",
                                                          "size_sim_func"),
                                              packages = c("dplyr",
                                                           "margins",
                                                           "broom",
                                                           "grf",
                                                           "tidyr")),
                    .progress = TRUE)
  
  plan(sequential)
  write.csv(new_size_sims_cov,
            file = "simulations/size/new_size_sims_cov.csv",
            row.names = FALSE)
} else{
  new_size_sims_cov <- readr::read_csv("simulations/size/new_size_sims_cov.csv")
}

# tiles


n_tile_list <- matrix(2:9, 100, 8) %>% 
  as.list()
run_sim <- FALSE
if (run_sim){
  plan(multisession)
  new_size_sims_tile <-  n_tile_list %>% 
    future_imap_dfr(~size_sim_func(draw = .y,
                                   n_covariates = 4,
                                   n_tiles = .x) %>% 
                      mutate(n_cov = 4,
                             n_tiles = .x),
                    .options = future_options(globals = c("run_first_stage_interactions_fast",
                                                          "discretise_df",
                                                          "find_SEs",
                                                          "create_fake_data",
                                                          "simulation_func",
                                                          "power_test",
                                                          "sim_first_stage_forest",
                                                          "size_sim_func"),
                                              packages = c("dplyr",
                                                           "margins",
                                                           "broom",
                                                           "grf",
                                                           "tidyr")),
                    .progress = TRUE)
  
  
  plan(sequential) 
  write.csv(new_size_sims_tile,
            file = "simulations/size/new_size_sims_tile.csv",
            row.names = FALSE)
} else {
  new_size_sims_tile <- readr::read_csv("simulations/size/new_size_sims_tile.csv")
}

rm(
  size_sim_func,
  n_tile_list,
  run_sim,
  n_cov_list
)
