library(haven)
library(dplyr)
library(AER)

df_china <- read_dta("Autor Dorn Hansen/dta/workfile_china.dta")


iv_model <- df_china %>%  
  ivreg(d_sh_empl_mfg ~ d_tradeusch_pw + t2 | d_tradeotch_pw_lag + t2, weights = timepwt48, data = .)
stargazer::stargazer(iv_model)
iv_model %>% broom::tidy()



full_form <- d_sh_empl_mfg ~ . -timepwt48 - d_tradeotch_pw_lag | . - timepwt48
iv_model_full <- df_china %>% 
  select(d_sh_empl_mfg,
         d_tradeusch_pw,
         d_tradeotch_pw_lag,
         l_shind_manuf_cbp,
         starts_with("reg"),
         l_sh_popedu_c,
         l_sh_popfborn,
         l_sh_empl_f,
         l_sh_routine33,
         l_task_outsource,
         t2,
         timepwt48) %>% 
  ivreg(formula = full_form, data = ., weights = timepwt48)
iv_model_full %>% broom::tidy()

iv_model_full_manual <- df_china %>% 
  ivreg(d_sh_empl_mfg ~ d_tradeusch_pw + l_shind_manuf_cbp +
          reg_midatl + reg_encen + reg_wncen + reg_satl + reg_escen +
          reg_wscen + reg_mount + reg_pacif + l_sh_popedu_c + l_sh_popfborn + l_sh_empl_f + l_sh_routine33 + l_task_outsource + t2 |
          d_tradeotch_pw_lag +
          l_shind_manuf_cbp +
          reg_midatl + reg_encen + reg_wncen + reg_satl + reg_escen +
          reg_wscen + reg_mount + reg_pacif + l_sh_popedu_c + l_sh_popfborn + l_sh_empl_f + l_sh_routine33 + l_task_outsource + t2,
        data = ., weights = timepwt48)

iv_model_full_manual %>% broom::tidy()
