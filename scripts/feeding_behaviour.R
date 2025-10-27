#####################################################
##                                                 ##
##             Antlion temperature study           ##
##                                                 ##
##                Feeding behaviour                ##
##                                                 ##
##                 JJ - 06/06/2025                 ##
##                                                 ##
#####################################################

rm(list = ls())
options(width = 100)

library(tidyverse)
library(brms)
library(posterior)     
library(bayesplot)   
library(tidybayes)
library(mgcv)
library(patchwork)
library(flextable)
library(see) 

#_______________________________________________________________________________
#### 1. Loading and wrangling data ####

feeding_behaviour <- read_csv(file = "data/feeding_behaviour.csv",  locale = locale(encoding = "Latin1")) %>% 
  pivot_longer(`Day 1`:`Day 61`) %>% 
  na.omit() %>% 
  mutate(fed = if_else(value == "F", 1, 0),
         day = as.numeric(gsub("Day ", "", name)),
         day_z = scale(day)) %>% 
  dplyr::select(id = ID, spp = sp, stage = Stage, temperature = `t (¡C)`,
                day, day_z, fed)

#_______________________________________________________________________________
#### 2. Workflow ####

## Models
base_model <- brm(fed ~ s(day_z, k = 20) + (1|id), 
                  data = feeding_behaviour, 
                  family = bernoulli(link = "logit"), 
                  prior = c(prior(normal(0, 0.5), class = "b"), 
                            prior(normal(0, 0.5), class = "Intercept"), 
                            prior(normal(0, 0.25), class = "sd", group = "id", lb = 0),
                            prior(normal(0, 1), class = "sds", lb = 0)), 
                  control = list(adapt_delta = 0.95, max_treedepth = 12),
                  chains = 4, cores = 4, seed = 420,
                  iter = 4000, warmup = 2000, backend = "cmdstanr")

temp_model <- brm(fed ~ s(day_z, by = temperature, k = 20) + (1|id), 
                  data = feeding_behaviour, 
                  family = bernoulli(link = "logit"), 
                  prior = c(prior(normal(0, 0.5), class = "b"), 
                            prior(normal(0, 0.5), class = "Intercept"), 
                            prior(normal(0, 0.25), class = "sd", group = "id", lb = 0),
                            prior(normal(0, 1), class = "sds", lb = 0)), 
                  control = list(adapt_delta = 0.95, max_treedepth = 12),
                  chains = 4, cores = 4, seed = 420,
                  iter = 4000, warmup = 2000, backend = "cmdstanr")

## Model comparison
base_model <- add_criterion(base_model, criterion = c("loo","waic"))
temp_model <- add_criterion(temp_model, criterion = c("loo","waic"))

as.data.frame(loo_compare(base_model, temp_model, criterion = "loo"))

## Predictions
# Pull out the scaling attributes that scale() stored
day_mean <- attr(feeding_behaviour$day_z, "scaled:center")
day_sd   <- attr(feeding_behaviour$day_z, "scaled:scale")

# Make newdata grid over the full range of days you observed
newdat <- expand_grid(day = 1:max(feeding_behaviour$day),
                      temperature = unique(feeding_behaviour$temperature)) %>% 
  mutate(day_z = (day - day_mean) / day_sd)

## Daily probability of feeding
haz_draws <- posterior_epred(
  temp_model,
  newdata    = newdat,
  re_formula = NA    # marginalise over random effects
)     

## 6.1 Fitted data plot
feeding_prediciton <- newdat %>% 
  mutate(fit = apply(haz_draws, 2, median),
         lwr = apply(haz_draws, 2, quantile, prob = c(0.025)),
         upr = apply(haz_draws, 2, quantile, prob = c(0.975))) %>% 
  mutate(temperature_lab = paste0(temperature, "\u00B0C")) 

feeding_plot <- 
  ggplot(feeding_prediciton, 
         aes(x = day, y = fit, fill = temperature_lab, colour = temperature_lab)) +
  geom_ribbon(aes(ymax = upr, ymin = lwr, colour = NULL), alpha = 0.3, show.legend = F) +
  geom_line(size = 1) +
  guides(colour = guide_legend(override.aes = list(linewidth = 2))) +
  scale_colour_viridis_d(option = "A", begin = 0.4, 
                         end = 0.6, aesthetics = c("colour", "fill")) +
  labs(x = "Experiment day", y = "Predicted probability of feeding", 
       colour = "Temperature\ntreatment") +
  theme_test(base_size = 11)

ggsave(feeding_plot, filename = "../feeding_behaviour.tiff", width = 14, 
       height = 10, units = "cm", dpi = 600)



