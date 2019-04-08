
create_fake_heterogeneous_data <- function(N,
                                           model_betas,
                                           force_positive = FALSE){
  ## X covariate means and variances for MV normal
  mu <- rep(0, 4)
  A <- matrix(runif(4^2)*2-1, ncol = 4)
  Sigma <- t(A) %*% A
  
  ## Whether or not to force all to be +ve for size simulations
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
  
  # Expanding model to include interactions of main effects with Z
  X_data_interactions <- model.matrix(~ Z*., data = X_data)
  
  # Takes NxK  xdata  and multiplies by transpose of NxK coef matrix then takes diag to give Nx1
  D <- diag(X_data_interactions %*% t(model_betas))
  
  # Population error
  errors_N <- rnorm(N) 
  
  sim_data <- X_data_interactions %>% 
    as_tibble() %>% 
    mutate(D = D + 5*errors_N)
  return(sim_data)
}

# TODO Fix this
create_hetero_coefs <- function(N){ # hardcoded 10 variables (constant + instrument + 4 main + interacts)
  # Multivariate normal with means on standard unif
  mu <- runif(10, min = -5, max = 5)
  # independent coefs with variance 2
  Sigma <- diag(1, 10)
  # return MV normal
  
  coef_matrix <- MASS::mvrnorm(n = N, mu = mu, Sigma)
  
  # factor_pos <- ifelse(rowMeans(coefs) < 0, -1, 1) 
  # 
  # coef_list <- lapply(seq_len(nrow(coefs)), function(i) coefs[i,]*factor_pos[i])
  # 
  # coef_matrix <- matrix(unlist(coef_list), N, 10, byrow = TRUE)
  return(coef_matrix)
}


find_true_dydx <- function(true_coefs,
                           true_data){
  
  # only data for dydx
  interaction_data <- true_data %>% 
    select_at(vars(contains("Z"))) %>% 
    mutate(Z = 1) %>% 
    as.matrix()
  
  # adding names so easier to debug
  colnames(true_coefs) <- colnames(true_data %>% 
                                     select(-D))
  # grab main effect and all interactions
  interaction_coefs <- matrix(NA, nrow(true_data), 5)
  interaction_coefs[, 1] <- true_coefs[, 2]
  interaction_coefs[, 2:5] <- true_coefs[, 7:10]
  
  # find true dydx
  dydx_true <- diag(interaction_data %*% t(interaction_coefs))
  return(dydx_true)
}




## Count n estimated defiers, n true defiers etc.
hetero_simul_func <- function(N){
  hetero_coefs <- create_hetero_coefs(N)
  hetero_fake_data <- create_fake_heterogeneous_data(N = N,
                                                     model_betas = hetero_coefs,
                                                     force_positive = FALSE)
  true_dydx <- find_true_dydx(hetero_coefs,
                              hetero_fake_data)
  first_stage_sign <- sign(sum(sign(true_dydx)))
  
  saturated_first_stage <- hetero_fake_data %>% 
    select_at(vars(-contains(":"), -`(Intercept)`)) %>% 
    run_first_stage_interactions_fast(.,
                                      "D",
                                      "Z") %>% 
    mutate(row_id = row_number(),
           pval_one_pos = pnorm(-t_stat),
           model = "Saturated First Stage",
           true_dydx = true_dydx,
           first_stage_sign = first_stage_sign) %>% 
    rename(dependent_variable = D)
  
  forest_model <- sim_first_stage_forest(hetero_fake_data)
  
  forest_full_wide <- left_join(forest_model,
                                saturated_first_stage %>% 
                                  select_at(vars(row_id,
                                                 contains("V", ignore.case = FALSE),
                                                 instrument,
                                                 true_dydx,
                                                 dependent_variable,
                                                 first_stage_sign)),
                                by = "row_id") %>% 
    rename("dydx_instrument" = predictions,
           "SE_dydx_instrument" = sigma_hat,
           "dydx_lo" = prediction_lo,
           "dydx_hi" = prediction_hi) %>% 
    select(-variance.estimates,
           -debiased.error) %>% 
    mutate(model = "Forest")
  
  models_df_long <- suppressWarnings(bind_rows(saturated_first_stage,
                                               forest_full_wide))
  return(models_df_long)
}



extract_hetero_draw_stats <- function(hetero_draw_df){
  
  summ_stats <- hetero_draw_df %>% 
    mutate(defier_true = !(first_stage_sign == sign(true_dydx)),
           defier_estimated = !(first_stage_sign == sign(dydx_instrument))) %>% 
    group_by(model) %>% 
    summarise(
      min_tstat = min(t_stat),
      max_tstat = max(t_stat),
      pval_one_neg = min(pval_one_neg),
      mean_dydx = mean(dydx_instrument),
      pval_one_pos = min(pval_one_pos),
      n_defiers_true = sum(defier_true),
      pct_defiers_true = n_defiers_true/(nrow(hetero_draw_df)/2),
      n_defiers_estimated = sum(defier_estimated),
      pct_defiers_estimated = n_defiers_estimated/(nrow(hetero_draw_df)/2),
      first_stage_sign = unique(first_stage_sign)
    ) %>% 
    mutate(pval_defier = ifelse(first_stage_sign > 0,
                                pval_one_neg,
                                pval_one_pos))
  return(summ_stats)
}


anon_sim_func <- function(dummy_arg, N){
  df <- hetero_simul_func(N) %>% 
    extract_hetero_draw_stats()
  return(df)
}

### Running sims

# Small
run_sim <- FALSE

if (run_sim) {
  library(furrr)
  plan(multisession)
  
  simulations_power_hetero_small <- 1:2000 %>% 
    future_map_dfr(anon_sim_func,
                   N = 1000,
                   .options = future_options(globals = c("sim_first_stage_forest",
                                                         "run_first_stage_interactions_fast",
                                                         "find_SEs",
                                                         "hetero_simul_func",
                                                         "find_true_dydx",
                                                         "extract_hetero_draw_stats",
                                                         "create_fake_heterogeneous_data",
                                                         "create_hetero_coefs",
                                                         "discretise_df"),
                                             packages = c("dplyr",
                                                          "margins",
                                                          "grf",
                                                          "tidyr")),
                   .progress = TRUE) %>% 
    mutate(N = "Small")
  
  
  
  plan(sequential)
  write.csv(simulations_power_hetero_small, file = "simulations/power/simulations_power_hetero_small.csv", row.names = FALSE)
} else {
  simulations_power_hetero_small <- readr::read_csv("simulations/power/simulations_power_hetero_small.csv")
}


# Medium
run_sim <- FALSE

if (run_sim) {
  library(furrr)
  plan(multisession)
  
  simulations_power_hetero_medium <- 1:2000 %>% 
    future_map_dfr(anon_sim_func,
                   N = 5000,
                   .options = future_options(globals = c("sim_first_stage_forest",
                                                         "run_first_stage_interactions_fast",
                                                         "find_SEs",
                                                         "hetero_simul_func",
                                                         "find_true_dydx",
                                                         "extract_hetero_draw_stats",
                                                         "create_fake_heterogeneous_data",
                                                         "create_hetero_coefs",
                                                         "discretise_df"),
                                             packages = c("dplyr",
                                                          "margins",
                                                          "grf",
                                                          "tidyr")),
                   .progress = TRUE) %>% 
    mutate(N = "Medium")
  
  
  
  plan(sequential)
  write.csv(simulations_power_hetero_medium, file = "simulations/power/simulations_power_hetero_medium.csv", row.names = FALSE)
} else {
  simulations_power_hetero_medium <- readr::read_csv("simulations/power/simulations_power_hetero_medium.csv")
}


# Large
run_sim <- FALSE

if (run_sim) {
  library(furrr)
  plan(multisession)
  
  simulations_power_hetero_large <- 1:2000 %>% 
    future_map_dfr(anon_sim_func,
                   N = 10000,
                   .options = future_options(globals = c("sim_first_stage_forest",
                                                         "run_first_stage_interactions_fast",
                                                         "find_SEs",
                                                         "hetero_simul_func",
                                                         "find_true_dydx",
                                                         "extract_hetero_draw_stats",
                                                         "create_fake_heterogeneous_data",
                                                         "create_hetero_coefs",
                                                         "discretise_df"),
                                             packages = c("dplyr",
                                                          "margins",
                                                          "grf",
                                                          "tidyr")),
                   .progress = TRUE) %>% 
    mutate(N = "Large")
  
  
  
  plan(sequential)
  write.csv(simulations_power_hetero_large, file = "simulations/power/simulations_power_hetero_large.csv", row.names = FALSE)
} else {
  simulations_power_hetero_large <- readr::read_csv("simulations/power/simulations_power_hetero_large.csv")
}



rm(
  run_sim,
  anon_sim_func,
  create_fake_heterogeneous_data,
  create_hetero_coefs,
  extract_hetero_draw_stats,
  find_true_dydx,
  hetero_simul_func
)
