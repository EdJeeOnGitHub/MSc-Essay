library(haven)
library(dplyr)
library(tidylog, warn.conflicts = FALSE)
library(tidyr)

df_cwhsc <- read_dta("dataverse_files Angrist Veteran/cwhsc_new.dta")
df_cpi <- read_dta("dataverse_files Angrist Veteran/cpi_angrist1990.dta")

df_merged <- left_join(df_cwhsc,
                       df_cpi,
                       by = "year") 
cpi_78 <- df_merged %>% 
  filter(year == 78) %>% 
  summarise(mean(cpi)) %>% 
  pull()
df_transformed <- df_merged %>% 
  mutate(cpi = round(cpi / cpi_78, 3),
         cpi_squared = cpi^2,
         smplsz = nj-nj0) %>% 
  filter(year >= 81) %>% 
  filter(byr <= 52 & byr >= 50) %>% 
  mutate(rlvar = (1/iweight)*smplsz,
         var = rlvar*cpi_squared,
         nomearn = earnings*cpi,
         eligible = 0,
         eligible = ifelse(byr == 50 & interval == 1, 1, eligible),
         eligible = ifelse(byr == 51 & interval <= 25, 1, eligible),
         eligible = ifelse((byr == 52 | byr == 53) & interval <= 19, 1, eligible),
         white = ifelse(race == 1, 1, 0)) %>% 
  group_by(white, byr, year, eligible, type) %>% 
  mutate(sumwt = sum(smplsz),
         wtmult = 1/sumwt,
         var_cm = wtmult*var,
         numtype = ifelse(type == "TAXAB", 1, NA),
         numtype = ifelse(type == "ADJ", 2, numtype),
         numtype = ifelse(type == "TOTAL", 3, numtype)) %>% 
  ungroup()

final_df <- df_transformed %>%
  group_by(white, byr, year, eligible, numtype) %>% 
  mutate(var_cm = weighted.mean(var_cm, smplsz),
         nomearn = weighted.mean(nomearn, smplsz),
         earnings = weighted.mean(earnings, smplsz),
         cpi = weighted.mean(cpi, smplsz),
         id = paste0(white, year, numtype, byr))



# 
# df %>% 
#   gather(variable, value, -(month:student)) %>%
#   unite(temp, student, variable) %>%
#   spread(temp, value)