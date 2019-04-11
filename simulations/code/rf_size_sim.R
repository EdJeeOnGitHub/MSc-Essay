
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
    sim_first_stage_forest() %>% 
    mutate(draw = draw)


  return(forest_output)
}



run_sim <- FALSE
if (run_sim){
  library(furrr)
  plan(multisession)
  rf_size_sims <- 1:1000 %>% 
    future_map_dfr(~size_sim_func(draw = .) %>% 
                     filter(pval_one_neg == min(pval_one_neg)) %>% 
                              summarise_all(first),
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
  write.csv(rf_size_sims,
            file = "simulations/size/rf_size_sims.csv",
            row.names = FALSE)
} else {
  rf_size_sims <- readr::read_csv("simulations/size/rf_size_sims.csv")
}

rm(
  run_sim,
  size_sim_func
)



