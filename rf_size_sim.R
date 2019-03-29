

#### Copied from essay chnk  on size

p.adjust_df <- function(p, method, n){
  df <- p.adjust.methods %>% 
    map_dfc(~p.adjust(p = p,
                      method = .x,
                      n = n))
  colnames(df) <- p.adjust.methods
  return(df)
}


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
  rf_size_sims <- 1:100 %>% 
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


rf_size_sims %>% 
  mutate(reject_null = ifelse(pval_one_neg < 0.05, TRUE, FALSE)) %>% 
  summarise(n_rejected = sum(reject_null),
            n_draws = n(),
            proportion_rejected = n_rejected/n_draws)


  
holm_adjusted_pvals <- seq(from = 1000, to = 500000, by = 1000) %>% 
  map_df(~p.adjust(n = ., method = "holm", p = rf_size_sims$pval_one_neg) %>% 
           enframe("row", "pval") %>% 
           mutate(n = .x))

holm_adjustment_summary <- holm_adjusted_pvals %>% 
  mutate(reject_null = ifelse(pval < 0.05, TRUE, FALSE)) %>% 
  group_by(n) %>% 
  summarise(empirical_alpha = sum(ifelse(reject_null, 1, 0)) / n())



holm_adjustment_summary %>%  
  ggplot(aes(x = n, y = empirical_alpha)) +
  geom_line() + 
  geom_hline(yintercept = 0.05, linetype = "longdash") +
  scale_y_log10() +
  theme_minimal()




