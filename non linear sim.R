library(dplyr)
library(margins)
library(grf)
library(ggplot2)
library(furrr)


set.seed(1234)



sim_first_stage_forest <- function(dataset){
  forest_sim <- dataset %>% 
    select_at(vars(contains("V"), -contains(":"))) %>% 
    as.matrix() %>% 
    causal_forest(X = .,
                  Y = dataset$D,
                  W = dataset$Z,
                  num.trees = 4000) %>% 
    predict(., estimate.variance = TRUE) %>% 
    as_tibble() %>% 
    mutate(sigma_hat = sqrt(variance.estimates),
           prediction_lo = predictions - 1.96*sigma_hat,
           prediction_hi = predictions + 1.96*sigma_hat,
           t_stat = predictions / sigma_hat,
           pval_one_neg = pnorm(t_stat),
           pval_one_pos = pnorm(-t_stat),
           row_id = row_number())
  return(forest_sim)
}
find_SEs <- function(model_data_no_y, model_vcov, instrument){
  model_data_no_y <- model_data_no_y %>% 
    rename("instrument" = instrument)
  model_formula <- as.formula(paste0("~", "instrument", "*."))
  
  model_matrix <- model_data_no_y %>% 
    mutate(instrument = 1) %>% 
    model.matrix(model_formula, data = .)
  
  gradient_matrix_kinda <- model_matrix %>%
    as_tibble() %>% 
    mutate_at(vars(-contains(":")), ~(. = 0)) %>% 
    mutate(instrument = 1) %>% 
    as.matrix()
  
  split_indices <- seq(from = 1, to = nrow(gradient_matrix_kinda), by = 10000) %>% 
    c(nrow(gradient_matrix_kinda))
  
  if (length(split_indices) < 3){
    vcov_dydx_intermediate <- Rfast::mat.mult(gradient_matrix_kinda, model_vcov)
    vcov_dydx <- Rfast::mat.mult(vcov_dydx_intermediate,
                                 Rfast::transpose(gradient_matrix_kinda))
    SE_dydx <- sqrt(diag(vcov_dydx))
    
  } else {
    baby_SE <- matrix(nrow = nrow(gradient_matrix_kinda), ncol = 1)
    for (i in 1:(length(split_indices)-1)){
      baby_matrix <- gradient_matrix_kinda[split_indices[[i]]:split_indices[[i+1]], ]
      # print(paste0(split_indices[i],"to", split_indices[i+1]))
      baby_vcov <- Rfast::mat.mult(baby_matrix, model_vcov)
      baby_vcov <- Rfast::mat.mult(baby_vcov, Rfast::transpose(baby_matrix))
      baby_SE[split_indices[[i]]:split_indices[[i+1]], ] <- sqrt(diag(baby_vcov))
      i <- i + 1
    }
    SE_dydx <- baby_SE[, 1]
  }
  
  return(SE_dydx)
  
} 
run_first_stage_interactions_fast <- function(dataset,
                                              dependent_variable,
                                              instrument,
                                              weights = NULL,
                                              vcov_func = vcov,
                                              ...){
  dataset <- dataset %>% 
    rename("instrument" = instrument)
  model_formula <- as.formula(paste0(dependent_variable,  "~ ", "instrument", "*."))
  first_stage_fit <- lm(data = dataset, formula = model_formula, weights = weights)
  degrees_freedom <- first_stage_fit$df.residual
  instrument_dummy_val <- dataset$instrument[1]
  
  dataset_unique <- dataset %>% 
    select(-dependent_variable, -instrument) %>% 
    mutate(instrument = instrument_dummy_val)
  first_stage_margins <- dydx(model = first_stage_fit, data = dataset_unique, variable = "instrument") 
  df_margins <- bind_cols(first_stage_margins %>% 
                            select(contains("dydx")),
                          dataset_unique) %>% 
    as_tibble()
  df_margins$instrument <- dataset$instrument
  df_margins$SE_dydx_instrument <- find_SEs(dataset_unique, vcov_func(first_stage_fit, ...), instrument)
  
  df_margins <- df_margins %>% 
    mutate(dydx_lo = dydx_instrument - qt(0.975, degrees_freedom) * SE_dydx_instrument,
           dydx_hi = dydx_instrument + qt(0.975, degrees_freedom) * SE_dydx_instrument,
           t_stat = dydx_instrument/SE_dydx_instrument,
           pval_one_neg = pt(t_stat, degrees_freedom),
           pval_holm = p.adjust(pval_one_neg, method = "holm"))
  df_margins$dependent_variable <- dataset %>% 
    select(dependent_variable) %>% 
    pull()
  return(df_margins)
}


## create the X data according to thingy so can tell true defiers.
create_fake_non_linear_data <- function(N,
                                           force_positive = FALSE){
  ## X covariate means and variances for MV normal
  mu <- rep(1, 4)
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
  # Population error
  errors_N <- rnorm(N) 

  
  partial_effects <- X_data_interactions %>% 
    as_tibble() %>% 
    transmute(beta_intercept = rnorm(1),
              beta_Z = pmax(V1, -1.5),
              beta_1 = rnorm(1),
              beta_2 = rnorm(1),
              beta_3 = rnorm(1),
              beta_4 = rnorm(1),
              beta_Zi1 = 0*`Z:V1`^3*rnorm(1),
              beta_Zi2 = 0*abs(`Z:V2`)*rnorm(1),
              beta_Zi3 = 0*pmax(`Z:V3`, 0),
              beta_Zi4 = 0*`Z:V4`*rnorm(1)) 

    
    
    # transmute(beta_intercept = rnorm(1),
    #           beta_Z = pmax(V1, 0)^2,
    #           beta_1 = rnorm(1),
    #           beta_2 = pmin(V2, 0),
    #           beta_3 = pmax(V3, 1)*0.5,
    #           beta_4 = V4*0.2,
    #           beta_Zi1 = `Z:V1`^3*rnorm(1),
    #           beta_Zi2 = abs(`Z:V2`)*rnorm(1),
    #           beta_Zi3 = pmax(`Z:V3`, 0),
    #           beta_Zi4 = `Z:V4`*rnorm(1)) 
  
  
  D <- diag(X_data_interactions %*% t(as.matrix(partial_effects)))
  
  sim_data <- X_data_interactions %>% 
    as_tibble() %>% 
    mutate(D = D + errors_N*5)
  
  
  return(list("sim_data" = sim_data,
              "partial_effects" = partial_effects))
}


## Count n estimated defiers, n true defiers etc.
non_linear_simul_func <- function(N){
  non_linear_output <- create_fake_non_linear_data(N = N,
                                                   force_positive = FALSE)
  fake_data <- non_linear_output$sim_data
  
  partial_effects_coefs <- non_linear_output$partial_effects %>% 
    select_at(vars(contains("Z"))) %>% 
    as.matrix()
  
  partial_effects_data <- fake_data %>% 
    select_at(vars(contains("Z"))) %>% 
    mutate(Z = 1) %>% 
    as.matrix()
  
  true_dydx <- diag(partial_effects_data %*% t(partial_effects_coefs))
  
  first_stage_sign <- sign(sum(sign(true_dydx)))
  
  saturated_first_stage <- fake_data %>% 
    select_at(vars(-contains(":"), -`(Intercept)`)) %>% 
    run_first_stage_interactions_fast(.,
                                      "D",
                                      "Z") %>% 
    mutate(row_id = row_number(),
           pval_one_pos = pnorm(-t_stat),
           model = "saturated first stage",
           true_dydx = true_dydx,
           first_stage_sign = first_stage_sign)
  
  forest_model <- sim_first_stage_forest(fake_data)
  
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
    mutate(model = "forest")
  
  models_df_long <- suppressWarnings(bind_rows(saturated_first_stage %>% 
                                                 select(-pval_holm),
                                               forest_full_wide))
  return(models_df_long)
}



extract_draw_stats <- function(hetero_draw_df){
  
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
  df <- non_linear_simul_func(N) %>% 
    extract_draw_stats()
  return(df)
}


anon_sim_func(100, 1000)
### Running sims


run_sim <- TRUE

if (run_sim) {
  library(furrr)
  plan(multisession)
  # small
  simulations_power_non_linear <- 1:2000 %>% 
    future_map_dfr(anon_sim_func,
                   N = 1000,
                   .options = future_options(globals = c("sim_first_stage_forest",
                                                         "run_first_stage_interactions_fast",
                                                         "find_SEs",
                                                         "non_linear_simul_func",
                                                         "extract_draw_stats",
                                                         "create_fake_non_linear_data"),
                                             packages = c("dplyr",
                                                          "margins",
                                                          "grf")),
                   .progress = TRUE)
  
  
  
  plan(sequential)
  write.csv(simulations_power_non_linear, file = "C:/Users/ed/Dropbox/Ed/power_simulations_non_linear_small.csv", row.names = FALSE)
} 


run_sim <- TRUE

if (run_sim) {
  library(furrr)
  plan(multisession)
  # medium
  simulations_power_non_linear <- 1:2000 %>% 
    future_map_dfr(anon_sim_func,
                   N = 5000,
                   .options = future_options(globals = c("sim_first_stage_forest",
                                                         "run_first_stage_interactions_fast",
                                                         "find_SEs",
                                                         "non_linear_simul_func",
                                                         "extract_draw_stats",
                                                         "create_fake_non_linear_data"),
                                             packages = c("dplyr",
                                                          "margins",
                                                          "grf")),
                   .progress = TRUE)
  
  
  
  plan(sequential)
  write.csv(simulations_power_non_linear, file = "C:/Users/ed/Dropbox/Ed/power_simulations_non_linear_medium.csv", row.names = FALSE)
} 



run_sim <- TRUE

if (run_sim) {
  library(furrr)
  plan(multisession)
  # large
  simulations_power_non_linear <- 1:2000 %>% 
    future_map_dfr(anon_sim_func,
                   N = 10000,
                   .options = future_options(globals = c("sim_first_stage_forest",
                                                         "run_first_stage_interactions_fast",
                                                         "find_SEs",
                                                         "non_linear_simul_func",
                                                         "extract_draw_stats",
                                                         "create_fake_non_linear_data"),
                                             packages = c("dplyr",
                                                          "margins",
                                                          "grf")),
                   .progress = TRUE)
  
  
  
  plan(sequential)
  write.csv(simulations_power_non_linear, file = "C:/Users/ed/Dropbox/Ed/power_simulations_non_linear_large.csv", row.names = FALSE)
} 