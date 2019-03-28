

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

adjusted_pvals_rf_size <- rf_size_sims %>% 
  mutate(bonferroni = pval_one_neg*200000,
         none = pval_one_neg)


rf_size_sims %>% 
  nest(-draw) %>% 
  mutate(adjusted_pvals = map(.x = data,
                              ~p.adjust_df(p = .x$pval_one_neg,
                                           method = "bonferroni",
                                           n = 1)))


adjusted_pvals_rf_size %>% 
  gather(method, pval, bonferroni:none) %>% 
  mutate(reject_null = ifelse(pval < 0.05, TRUE, FALSE)) %>% 
  group_by(method) %>% 
  summarise(empirical_alpha = sum(ifelse(reject_null, 1, 0)) / n())
##

adjusted_pvals_cov %>% 
  unnest() %>% 
  select(-pvals) %>% 
  gather(method, pval, holm:none) %>% 
  mutate(reject_null = ifelse(pval < 0.05, TRUE, FALSE)) %>% 
  group_by(draw, n_cov, n_tiles,method) %>% 
  summarise(null_rejected_in_test = ifelse(sum(reject_null) > 0, TRUE, FALSE),
            n_groups = mean(n_groups)) %>% 
  ungroup() %>% 
  group_by(n_cov,n_tiles, method) %>% 
  summarise(empirical_alpha = sum(null_rejected_in_test) / n(),
            n_groups = first(n_groups)) %>% 
  arrange(n_cov) %>% 
  mutate(n_comp = n_cov^n_tiles,
         theoretical_alpha = 1 - 0.95^n_comp)
##

