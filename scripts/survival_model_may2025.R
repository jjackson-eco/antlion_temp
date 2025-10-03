#####################################################
##                                                 ##
##             Antlion temperature study           ##
##                                                 ##
##         Basic survival analysis pipeline        ##
##                                                 ##
##                JJ - 13/05/2025                  ##
##                                                 ##
#####################################################
rm(list = ls())
options(width = 100)

library(tidyverse)
library(survival)
library(ggfortify)
library(survminer)

#_______________________________________________________________________________
#### 1. Loading + wrangling data ####

load("data/antlion_lh_raw.RData")

## converting to survival type data
antlion_survdat <- 
  antlion_lh_raw %>% 
  mutate(time = experiment_length,
         time = if_else(is.na(death_days) == F, death_days, experiment_length),
         status = if_else(is.na(death_days) == F, 1, 0)) %>% 
  dplyr::select(ID, temperature, time, status)

#_______________________________________________________________________________
#### 2. Simple Kaplan Meier ####

km_fit <- survfit(Surv(time, status) ~ 1, data = antlion_survdat)

autoplot(km_fit) +
  labs(x = "Experiment days", y = "Larval survival") +
  theme_test(base_size = 12)


km_fit_temp <- survfit(Surv(time, status) ~ temperature, data = antlion_survdat)

viridis::inferno(n = 5, begin = 0.2, end = 0.8)

jpeg("output/exploration_may25/survplot_km_test.jpeg", 
     width = 20, height = 17, units = "cm", res = 400)
ggsurvplot(km_fit_temp,
           palette = viridis::magma(n = 5, begin = 0.2, end = 0.8),
           conf.int = T,
           ggtheme = theme_test(base_size = 12))
dev.off()


#_______________________________________________________________________________
#### 3. Cox proportional hazards ####

