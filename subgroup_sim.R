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
           model = "Saturated First Stage")
  
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
    mutate(model = "Forest")
  
  models_df_long <- suppressWarnings(bind_rows(first_stage_model %>% 
                                                 select(-pval_holm),
                                               forest_full_wide))
  
  return(models_df_long)
}

coef_matrix <- matrix(runif(5000*10, -1, 1), 5000, 10)
coef_matrix_df <- coef_matrix %>% 
  as_tibble()
factor_pos <- ifelse(rowMeans(coef_matrix) < 0, -1, 1) 

coef_list <- lapply(seq_len(nrow(coef_matrix)), function(i) coef_matrix[i,]*factor_pos[i])


simulation_func <- function(x, N, force_positive = FALSE){
  sim_data <- create_fake_data(N = N,
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


##### Actual Running of Simulations

# Small
run_sim <- FALSE
if (run_sim) {
  library(furrr)
  plan(multisession)
  
  simulations_power_subgroup_small <- coef_list %>% 
    future_map_dfr(simulation_func,
                   N = 1000,
                   .options = future_options(globals = c("sim_data",
                                                         "run_first_stage_interactions_fast",
                                                         "find_SEs",
                                                         "create_fake_data",
                                                         "simulation_func",
                                                         "power_test",
                                                         "sim_first_stage_forest"),
                                             packages = c("dplyr",
                                                          "margins",
                                                          "broom",
                                                          "grf")),
                   .progress = TRUE) %>% 
    mutate(N = "Small")
  
  
  
  plan(sequential)
  write.csv(simulations_power_subgroup_small, file = "simulations/power/simulations_power_subgroup_small.csv", row.names = FALSE)
} else {
  simulations_power_subgroup_small <- readr::read_csv("simulations/power/simulations_power_subgroup_small.csv")
}

# Medium
run_sim <- FALSE
if (run_sim) {
  library(furrr)
  plan(multisession)
  
  simulations_power_subgroup_medium <- coef_list %>% 
    future_map_dfr(simulation_func,
                   N = 5000,
                   .options = future_options(globals = c("sim_data",
                                                         "run_first_stage_interactions_fast",
                                                         "find_SEs",
                                                         "create_fake_data",
                                                         "simulation_func",
                                                         "power_test",
                                                         "sim_first_stage_forest"),
                                             packages = c("dplyr",
                                                          "margins",
                                                          "broom",
                                                          "grf")),
                   .progress = TRUE) %>% 
    mutate(N = "Medium")
  
  
  
  plan(sequential)
  write.csv(simulations_power_subgroup_medium, file = "simulations/power/simulations_power_subgroup_medium.csv", row.names = FALSE)
} else {
  simulations_power_subgroup_medium <- readr::read_csv("simulations/power/simulations_power_subgroup_medium.csv")
}

# Large
run_sim <- FALSE
if (run_sim) {
  library(furrr)
  plan(multisession)
  
  simulations_power_subgroup_large <- coef_list %>% 
    future_map_dfr(simulation_func,
                   N = 10000,
                   .options = future_options(globals = c("sim_data",
                                                         "run_first_stage_interactions_fast",
                                                         "find_SEs",
                                                         "create_fake_data",
                                                         "simulation_func",
                                                         "power_test",
                                                         "sim_first_stage_forest"),
                                             packages = c("dplyr",
                                                          "margins",
                                                          "broom",
                                                          "grf")),
                   .progress = TRUE) %>% 
    mutate(N = "Large")
  
  
  
  plan(sequential)
  write.csv(simulations_power_subgroup_large, file = "simulations/power/simulations_power_subgroup_large.csv", row.names = FALSE)
} else {
  simulations_power_subgroup_large <- readr::read_csv("simulations/power/simulations_power_subgroup_large.csv")
}

##### Test Size Simulations ####


# Small
run_sim <- FALSE
if (run_sim) {
  library(furrr)
  plan(multisession)
  
  simulations_size_subgroup_small <- coef_list %>% 
    map(abs) %>% 
    future_map_dfr(simulation_func,
                   N = 1000,
                   force_positive = TRUE,
                   .options = future_options(globals = c("sim_data",
                                                         "run_first_stage_interactions_fast",
                                                         "find_SEs",
                                                         "create_fake_data",
                                                         "simulation_func",
                                                         "power_test",
                                                         "sim_first_stage_forest"),
                                             packages = c("dplyr",
                                                          "margins",
                                                          "broom",
                                                          "grf")),
                   .progress = TRUE) %>% 
    mutate(N = "Small")
  
  
  
  plan(sequential)
  write.csv(simulations_size_subgroup_small, file = "simulations/size/simulations_size_subgroup_small.csv", row.names = FALSE)
} else {
  simulations_size_subgroup_small <- readr::read_csv("simulations/size/simulations_size_subgroup_small.csv")
}


# Medium
run_sim <- FALSE
if (run_sim) {
  library(furrr)
  plan(multisession)
  
  simulations_size_subgroup_medium <- coef_list %>% 
    map(abs) %>% 
    future_map_dfr(simulation_func,
                   N = 5000,
                   force_positive = TRUE,
                   .options = future_options(globals = c("sim_data",
                                                         "run_first_stage_interactions_fast",
                                                         "find_SEs",
                                                         "create_fake_data",
                                                         "simulation_func",
                                                         "power_test",
                                                         "sim_first_stage_forest"),
                                             packages = c("dplyr",
                                                          "margins",
                                                          "broom",
                                                          "grf")),
                   .progress = TRUE) %>% 
    mutate(N = "Medium")
  
  
  
  plan(sequential)
  write.csv(simulations_size_subgroup_medium, file = "simulations/size/simulations_size_subgroup_medium.csv", row.names = FALSE)
} else {
  simulations_size_subgroup_medium <- readr::read_csv("simulations/size/simulations_size_subgroup_medium.csv")
}

# Large
run_sim <- FALSE
if (run_sim) {
  library(furrr)
  plan(multisession)
  
  simulations_size_subgroup_large <- coef_list %>% 
    map(abs) %>% 
    future_map_dfr(simulation_func,
                   N = 10000,
                   force_positive = TRUE,
                   .options = future_options(globals = c("sim_data",
                                                         "run_first_stage_interactions_fast",
                                                         "find_SEs",
                                                         "create_fake_data",
                                                         "simulation_func",
                                                         "power_test",
                                                         "sim_first_stage_forest"),
                                             packages = c("dplyr",
                                                          "margins",
                                                          "broom",
                                                          "grf")),
                   .progress = TRUE) %>% 
    mutate(N = "Large")
  
  
  
  plan(sequential)
  write.csv(simulations_size_subgroup_large, file = "simulations/size/simulations_size_subgroup_large.csv", row.names = FALSE)
} else {
  simulations_size_subgroup_large <- readr::read_csv("simulations/size/simulations_size_subgroup_large.csv")
}


rm(
  coef_list,
  coef_matrix,
  coef_matrix_df,
  factor_pos,
  run_sim,
  create_fake_data,
  power_test,
  simulation_func
)

