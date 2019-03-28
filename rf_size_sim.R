## Random forest size sims





size_sim_func <- function(draw,
                          n_covariates = 4,
                          n_tiles = 2,
                          N = 1000){
  sim_df <- matrix(rnorm(N*n_covariates), N, n_covariates) %>% 
    as_tibble() %>% 
    mutate(D = rnorm(N),
           Z = rbinom(N, 1, 0.5))
  
  
  forest_output <- sim_df %>% 
    sim_first_stage_forest()


  return(forest_output)
}




plan(multisession)
rf_size_sims <- 1:1000 %>% 
  future_map_dfr(size_sim_func,
                .options = future_options(globals = c(                                                      "create_fake_data",
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
