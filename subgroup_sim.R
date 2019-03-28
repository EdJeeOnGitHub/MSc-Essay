library(dplyr)
library(Rfast)
library(margins)
library(tibble)
library(ggplot2)
library(broom)
library(purrr)
library(grf)
library(tidyr)
set.seed(1234)
##### Fake Data Generation #####


# 2 + 2k = n
coefs_positive <- runif(10, 0, 1)

create_fake_data <- function(N, model_betas, force_positive = FALSE){
  mu <- rep(0, 4)
  
  A <- matrix(runif(4^2)*2-1, ncol=4) 
  Sigma <- t(A) %*% A
  
  if (force_positive == TRUE){
    X_data <- abs(MASS::mvrnorm(n = N, mu = mu, Sigma = Sigma)) %>% 
      as_tibble() %>% 
      mutate(Z = rbinom(n = N, 1, 0.5),
             V4 = ifelse(V4 > 0.5, 1, 0))
  } else {
    X_data <- MASS::mvrnorm(n = N, mu = mu, Sigma = Sigma) %>% 
      as_tibble() %>% 
      mutate(Z = rbinom(n = N, 1, 0.5),
             V4 = ifelse(V4 > 0, 1, 0))
  }
  
  X_data_interactions <- model.matrix(~ Z*., data = X_data)
  D <- X_data_interactions %*% model_betas
  sim_data <- X_data_interactions %>% 
    as_tibble() %>% 
    mutate(D = D[, 1] + rnorm(N)*5)
  return(sim_data)
}

##### Power Simulations #####
## Split up data generation from power test so we can ensure enough data in each band.
power_test <- function(coefs, sim_data){
  names(coefs) <- names(sim_data)[1:10]
  
  
  first_stage_model <- sim_data %>% 
    select_at(vars(-contains(":"), -`(Intercept)`)) %>% 
    run_first_stage_interactions_fast(.,
                                      "D",
                                      "Z") %>% 
    mutate(row_id = row_number(),
           pval_one_pos = pnorm(-t_stat),
           model = "saturated first stage")
  
  interaction_coefs <- c(coefs[2],
                         coefs[7:10])
  interaction_X <- sim_data %>% 
    select_at(vars(Z, contains("V"), -contains(":"))) %>% 
    mutate(Z = 1)
  interaction_X <- model.matrix(~Z*., data = interaction_X) %>% 
    as_tibble() %>% 
    select_at(vars(Z, contains(":"))) %>% 
    as.matrix()
  
  true_dydx <- interaction_X %*% interaction_coefs
  first_stage_model$true_dydx <- true_dydx[,1]
  forest_model <- sim_first_stage_forest(sim_data)
  
  forest_full_wide <- left_join(forest_model,
                                first_stage_model %>% 
                                  select_at(vars(row_id,
                                                 contains("V", ignore.case = FALSE),
                                                 instrument,
                                                 true_dydx,
                                                 dependent_variable)),
                                by = "row_id") %>% 
    rename("dydx_instrument" = predictions,
           "SE_dydx_instrument" = sigma_hat,
           "dydx_lo" = prediction_lo,
           "dydx_hi" = prediction_hi) %>% 
    select(-variance.estimates,
           -debiased.error) %>% 
    mutate(model = "forest")
  
  models_df_long <- suppressWarnings(bind_rows(first_stage_model %>% 
                                                 select(-pval_holm),
                                               forest_full_wide))
  
  return(models_df_long)
}



coef_matrix <- matrix(runif(100*10, -1, 1), 100, 10)
coef_matrix_df <- coef_matrix %>% 
  as_tibble()
factor_pos <- ifelse(rowMeans(coef_matrix) < 0, -1, 1) 

coef_list <- lapply(seq_len(nrow(coef_matrix)), function(i) coef_matrix[i,]*factor_pos[i])


simulation_func <- function(x, force_positive = FALSE){
  sim_data <- create_fake_data(N = 10000,
                               model_betas = x,
                               force_positive = force_positive)
  
  pw_test_data <- power_test(x, sim_data)
  
  estimated_coef <- lm(dependent_variable ~ instrument + V1 + V2 + V3 + V4, data = pw_test_data) %>% 
    tidy() %>% 
    filter(term == "instrument") %>% 
    select(estimate)
  
  pw_test_data$estimate <- estimated_coef[[1]]
  
  power_sim_df <- pw_test_data %>%   
    mutate(defier_true = !(sign(true_dydx) == (sign(estimate))),
           defier_estimated = !(sign(dydx_instrument) == (sign(estimate)))) %>% 
    group_by(model) %>% 
    summarise(min_tstat = min(t_stat),
              max_tstat = max(t_stat),
              pval_one_neg = min(pval_one_neg),
              mean_dydx = mean(dydx_instrument),
              pval_one_pos = min(pval_one_pos),
              n_defiers_true = sum(defier_true),
              pct_defiers_true = n_defiers_true/nrow(pw_test_data),
              n_defiers_estimated = sum(defier_estimated),
              pct_defiers_estimated = n_defiers_estimated/nrow(pw_test_data),
              estimate = first(estimate)) %>% 
    mutate(pval_defier = ifelse(estimate > 0,
                                pval_one_neg,
                                pval_one_pos))
  return(power_sim_df)
}

run_sim <- FALSE

if (run_sim) {
  library(furrr)
  plan(multisession)
  
  simulations_power <- coef_list %>% 
    future_map_dfr(simulation_func,
                   .options = future_options(globals = c("sim_data",
                                                         "run_first_stage_interactions_fast",
                                                         "find_SEs",
                                                         "create_fake_data",
                                                         "simulation_func",
                                                         "power_test",
                                                         "sim_first_stage_forest"),
                                             packages = c("dplyr",
                                                          "modelr",
                                                          "margins",
                                                          "broom",
                                                          "grf")),
                   .progress = TRUE)
  
  
  
  plan(sequential)
  # write.csv(simulations_power, file = "power_simulations.csv", row.names = FALSE)
  # write.csv(simulations_power, file = "power_simulations_both.csv", row.names = FALSE)  
} else {
  simulations_power <- readr::read_csv("power_simulations_both.csv")
}


##### Test Size Simulations ####
## DO this after while loops done.


run_sim <- FALSE
if (run_sim) {
  library(furrr)
  plan(multisession)
  
  simulations_size <- coef_list %>% 
    map(abs) %>% 
    future_map_dfr(simulation_func,
                   force_positive = TRUE,
                   .options = future_options(globals = c("sim_data",
                                                         "run_first_stage_interactions_fast",
                                                         "find_SEs",
                                                         "create_fake_data",
                                                         "simulation_func",
                                                         "power_test",
                                                         "sim_first_stage_forest"),
                                             packages = c("dplyr",
                                                          "modelr",
                                                          "margins",
                                                          "broom",
                                                          "grf")),
                   .progress = TRUE)
  
  
  
  plan(sequential)
  write.csv(simulations_size, file = "size_simulations_both.csv", row.names = FALSE)
  # write.csv(simulations_size, file = "size_simulations_both.csv", row.names = FALSE)  
} else {
  simulations_size <- readr::read_csv("size_simulations_both.csv")
}
