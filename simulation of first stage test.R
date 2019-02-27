library(dplyr)
library(Rfast)
library(margins)
library(tibble)
library(ggplot2)
library(broom)
library(purrr)
##### Functions

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
      print(paste0(split_indices[i],"to", split_indices[i+1]))
      
      baby_vcov <- Rfast::mat.mult(baby_matrix, model_vcov)
      baby_vcov <- Rfast::mat.mult(baby_vcov, Rfast::transpose(baby_matrix))
      baby_SE[split_indices[[i]]:split_indices[[i+1]], ] <- sqrt(diag(baby_vcov))
      i <- i + 1
    }
    SE_dydx <- baby_SE[, 1]
  }
  
  return(SE_dydx)
  
} 
run_first_stage_interactions_fast <- function(dataset, dependent_variable, instrument, weights = NULL){
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
  df_margins$SE_dydx_instrument <- find_SEs(dataset_unique, vcov(first_stage_fit), instrument)
  
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



##### Fake Data Generation


# 2 + 2k = n
coefs_positive <- runif(10, 0, 1)

create_fake_data <- function(N, model_betas){
  X_data <- matrix(abs(rnorm(N*4)), N, 4) %>% 
    as_tibble() %>% 
    mutate(Z = rbinom(n = N, 1, 0.5),
           V4 = ifelse(V4 > 0.5, 1, 0))
  X_data_interactions <- model.matrix(~ Z*., data = X_data)
  D <- X_data_interactions %*% model_betas
  sim_data <- X_data_interactions %>% 
    as_tibble() %>% 
    mutate(D = D[, 1] + rnorm(N)*10)
  return(sim_data)
}
sim_data <- create_fake_data(1000, coefs_positive)


##### Running models


sim_positive_ols <- lm(D ~ ., data = sim_data) %>% 
  tidy() %>% 
  mutate(true_beta = coefs_positive)
sim_positive_ols


sim_positive <- sim_data %>% 
  select_at(vars(-contains(":"), -`(Intercept)`)) %>% 
  run_first_stage_interactions_fast(.,
                                    dependent_variable = "D",
                                    instrument = "Z") %>% 
  mutate(pval_twosided = 2*pnorm(-abs(t_stat)))

p_adj <- p.adjust(sim_positive$pval_twosided, method = "holm")

sim_positive$p_2_adjust <- p_adj


library(tidyr)


sim_positive %>% 
  select(p_2_adjust, pval_twosided)

##### Romano Wolf
# library(StepwiseTest)
# library(modelr)
# library(furrr)
# plan(multisession)
# 
# anon_function <- function(dummy_arg,data){
#   data %>% 
#     select_at(vars(-contains(":"), -`(Intercept)`)) %>% 
#     resample_bootstrap() %>% 
#     as.data.frame() %>% 
#     run_first_stage_interactions_fast(.,
#                                       dependent_variable = "D",
#                                       instrument = "Z") %>% 
#     select(t_stat) %>% 
#     pull()
# }
# 
# model_boot <- 1:1000 %>% 
#   future_map_dfc(., anon_function, data = sim_data, .progress = TRUE,
#              .options = future_options(globals = c("sim_data",
#                                                    "run_first_stage_interactions_fast",
#                                                    "find_SEs"),
#                                        packages = c("dplyr",
#                                                     "modelr",
#                                                     "margins"))) %>% 
#   as.matrix()
# 
# adjusted_val <- FWERkControl(test_stat = sim_positive$t_stat,
#              boot_stat = model_boot,
#              k = 1,
#              alpha = 0.05)
# 
# adjusted_val$Reject %>% 
#   c() %>% 
#   enframe(.) %>% 
#   filter(value == 1)
# 
# adjusted_val$CV
# 
# 
# 
# 
# fdp <- FDPControl(sim_positive$t_stat,
#                   boot_stat = model_boot,
#                   gamma = 0.50,
#                   0.05)
# 


##### Power stuff

coefs_negative <- runif(10, -1, 1)
coefs_negative[1:2] <- abs(coefs_negative[1:2])


sim_negative <- coefs_negative %>% 
  create_fake_data(N = 1000,
                   model_betas = .) %>% 
  select_at(vars(-contains(":"), -`(Intercept)`)) %>%
  run_first_stage_interactions_fast(.,
                                    "D",
                                    "Z")

fixed_negative_coefs <- c(0.5, 0.5, -0.5, -0.5, 0, 1, 0.25, 0.3, -0.2, -0.2)

power_neg <- 1:10 %>% 
  map_df(., ~(create_fake_data(
    N = 1000,
    model_betas = fixed_negative_coefs) %>% 
      select_at(vars(-contains(":"), - `(Intercept)`)) %>% 
      run_first_stage_interactions_fast(.,
                                        "D",
                                        "Z") %>% 
      mutate(draw = .x)))
power_neg


power_test <- function(coefs){
  sim_data <- create_fake_data(N = 1000,
                               model_betas = coefs)
  names(coefs) <- names(sim_data)[1:10]
  
  
  first_stage_model <- sim_data %>% 
    select_at(vars(-contains(":"), -`(Intercept)`)) %>% 
    run_first_stage_interactions_fast(.,
                                      "D",
                                      "Z")
  
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
  return(first_stage_model)
}

power_test(fixed_negative_coefs) %>% 
  ggplot(aes(x = true_dydx, y = dydx_instrument)) +
  geom_point()



power_test(fixed_negative_coefs) %>% 
  filter(true_dydx < 0 | dydx_instrument < 0)



## MAKE A TEST
# power_test(coefs_negative) %>% 
#   select(dydx_instrument, V1, V2, V3, V4, instrument, true_dydx) %>% 
#   mutate(diff = true_dydx - dydx_instrument)



coef_matrix <- matrix(rnorm(1000*10, mean = -1, sd = 3), 1000, 10)
coef_matrix[, 2] <- abs(coef_matrix[, 2])
coef_matrix_df <- coef_matrix %>% 
  as_tibble()
factor_pos <- ifelse(rowMeans(coef_matrix) < 0, -1, 1) 

coef_list <- lapply(seq_len(nrow(coef_matrix)), function(i) coef_matrix[i,]*factor_pos[i])


simulation_func <- function(x){
  pw_test_data <- power_test(x)
  
  summ_stats <- pw_test_data %>%  
    mutate(defier = ifelse(true_dydx < 0, 1, 0)) %>% 
    summarise(min_tstat = min(t_stat),
              pval_one_neg = min(pval_one_neg),
              pct_defier = sum(defier)/nrow(.),
              mean_dydx = mean(dydx_instrument))
  model_coefs <- lm(dependent_variable ~ instrument + V1 + V2 + V3 + V4, data = pw_test_data) %>% 
    tidy() %>% 
    filter(term == "instrument") %>% 
    select(estimate) %>% 
    bind_cols(summ_stats)
  return(model_coefs)
}

library(furrr)
plan(multisession)

simulations <- coef_list %>% 
  future_map_dfr(simulation_func,
                 .options = future_options(globals = c("sim_data",
                                                       "run_first_stage_interactions_fast",
                                                       "find_SEs",
                                                       "create_fake_data",
                                                       "simulation_func",
                                                       "power_test"),
                                           packages = c("dplyr",
                                                        "modelr",
                                                        "margins",
                                                        "broom")),
                 .progress = TRUE)

## TODO discretize and bin up

simulations %>% 
  filter(estimate > 0) %>%
  filter(pct_defier > 0) %>% 
  ggplot(aes(x = pct_defier, y = pval_one_neg)) +
  geom_point() +
  geom_hline(yintercept = 0.05) +
  scale_y_reverse() +
  geom_smooth()

simulations %>% 
  filter(estimate > 0) %>% 
  filter(pct_defier > 0) %>% 
  ggplot(aes(x = mean_dydx, y = pval_one_neg)) +
  geom_point() +
  geom_smooth()


sim_bin <- simulations %>% 
  group_by(gr=cut(pct_defier, breaks= seq(0, 1, by = 0.05)) ) %>% 
  mutate(n= n()) %>%
  arrange(as.numeric(gr))
plan(sequential)


sim_bin %>% 
  mutate(pct_rejected = sum(ifelse(pct_defier < 0.05, 1, 0))/n) %>% 
  ggplot(aes(x = gr, y = pct_rejected)) +
  geom_point()


sim_bin %>% 
  mutate(n_rejected = ifelse(pval_one_neg < 0.05, 1, 0),
         pct_rejected = sum(n_rejected)/n) %>% 
  na.omit() %>% 
  ggplot(aes(x = gr, y = pct_rejected)) +
  geom_point() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

library(ggridges)


sim_bin %>% 
  na.omit() %>% 
  ggplot(aes(x = pval_one_neg, y = gr, fill = gr)) +
  geom_density_ridges() +
  theme_ridges() +
  guides(fill = "none")


sim_bin %>% 
  na.omit() %>% 
  ggplot(aes(x = pval_one_neg, fill = gr)) +
  geom_histogram() +
  facet_wrap(~gr, ncol = 3) +
  guides(fill = "none")

##### Power checks

sim_negative %>% 
  filter(pval_one_neg < 0.05)

sim_negative %>% 
  filter(pval_one_neg == min(pval_one_neg))


##### Graphics


sim_positive %>%
  arrange(dydx_instrument) %>% 
  mutate(rank = row_number()) %>% 
  ggplot(aes(x = rank,
             y = dydx_instrument,
             ymin = dydx_lo,
             ymax = dydx_hi)) +
  geom_ribbon(alpha = 0.1) +
  theme_minimal() +
  geom_point() +
  labs(title = "Heterogeneous treatment effects partial derivative")



sim_negative %>%
  arrange(dydx_instrument) %>% 
  mutate(rank = row_number()) %>% 
  ggplot(aes(x = rank,
             y = dydx_instrument,
             ymin = dydx_lo,
             ymax = dydx_hi)) +
  geom_ribbon(alpha = 0.1) +
  theme_minimal() +
  geom_point() +
  labs(title = "Heterogeneous treatment effects partial derivative")






