#####################################################
##                                                 ##
##             Antlion temperature study           ##
##                                                 ##
##                 Mortality models                ##
##                                                 ##
##                 JJ - 14/05/2025                 ##
##                                                 ##
#####################################################

## Sex determination based on emergence and so ignored at this stage, 
## will do a separate analysis for sex differences with a reduced sample

rm(list = ls())
options(width = 100)

library(tidyverse)
library(brms)
library(posterior)     # tidy access to Stan draws
library(bayesplot)     # visual helpers
library(tidybayes)
library(mgcv)
library(patchwork)
library(flextable)
library(see) ## half violin test

#_______________________________________________________________________________
#### 1. Loading data ####

load("data/antlion_lh.RData", verbose = TRUE)

## exploration plot
antlion_mortality %>%
  filter(temperature != "25") %>% 
  group_by(day, temperature) %>% 
  summarise(mne = mean(dead), sem = (sd(dead)/sqrt(n()))) %>% 
  ggplot(aes(x = day, y = mne, colour = temperature)) +
  geom_errorbar(aes(ymax = mne + sem, ymin = mne - sem), width = 0.01) +
  geom_point()  

#_______________________________________________________________________________
#### 2. PPS ####

# Note on regularising + basis dimension smoothers:
#   1. basis dimension smoothers s(), are controlled with the weighting parameters (as in mgcv),
#   but in brms (also in mgcv but more obvious here) these weighting parameters are treated as random effect
#   variances through "sds" priors, and regularising these helps a lot.
#   
#   2. prior predictive simulation is crucial for balancing regularised priors and smoother wiggliness
#

# model but only sample from the priors and not posterior
model_prior_only <- 
  brm(dead ~ s(day_z, by = temperature, k = 20), 
                  data = antlion_mortality, 
                  family = bernoulli(link = "cloglog"), 
                  prior = c(prior(normal(0, 0.5), class = "b"), 
                            prior(normal(0, 0.5), class = "Intercept"),
                            prior(normal(0, 1), class = "sds", lb = 0)),
                  sample_prior = "only", 
                  chains = 2, cores = 2, seed = 420,
                  iter = 2000, warmup = 1000)

## Predict from the priors
grid <- expand_grid(
  temperature = unique(antlion_mortality$temperature),
  day_z       = seq(-1.1, 3.6, length = 120)
)

haz_prior <- add_epred_draws(model_prior_only, newdata = grid, re_formula = NA)  # hazard draws

haz_prior %>% 
  mutate(day = day_z * attr(antlion_mortality$day_z, "scaled:scale") +
           attr(antlion_mortality$day_z, "scaled:center")) %>% 
  ggplot(aes(day, .epred, group = .draw)) +
  geom_line(alpha = .05) +
  facet_wrap(~temperature) +
  labs(y = "Daily mortality hazard (prior)", title = "Prior-predictive curves")

#_______________________________________________________________________________
#### 3. Core models ####

amd <- antlion_mortality %>%
  filter(temperature != "25")

base_model <- brm(dead ~ s(day_z, k = 20) + spp + (1|id), 
                  data = amd, 
                  family = bernoulli(link = "cloglog"), 
                  prior = c(prior(normal(0, 0.5), class = "b"), 
                            prior(normal(0, 0.5), class = "Intercept"), 
                            ## Regularise more strongly on the random effect variances
                            prior(normal(0, 0.25), class = "sd", group = "id", lb = 0),
                            prior(normal(0, 1), class = "sds", lb = 0)), ## basis weights get treated like a random EFFECT!!! This was the key
                  control = list(adapt_delta = 0.95, max_treedepth = 12),
                  chains = 4, cores = 4, seed = 420,
                  iter = 4000, warmup = 2000, backend = "cmdstanr")

spp_model <- brm(dead ~ s(day_z, by = spp, k = 20) + (1|id), 
                  data = amd, 
                  family = bernoulli(link = "cloglog"), 
                  prior = c(prior(normal(0, 0.5), class = "b"), 
                            prior(normal(0, 0.5), class = "Intercept"), 
                            ## Regularise more strongly on the random effect variances
                            prior(normal(0, 0.25), class = "sd", group = "id", lb = 0),
                            prior(normal(0, 1), class = "sds", lb = 0)), ## basis weights get treated like a random EFFECT!!!
                  control = list(adapt_delta = 0.95, max_treedepth = 12),
                  chains = 4, cores = 4, seed = 420,
                  iter = 4000, warmup = 2000, backend = "cmdstanr")

temp_model <- brm(dead ~ s(day_z, by = temperature, k = 20) + spp + (1|id), 
                  data = amd, 
                  family = bernoulli(link = "cloglog"), 
                  prior = c(prior(normal(0, 0.5), class = "b"), 
                            prior(normal(0, 0.5), class = "Intercept"),
                            ## Regularise more strongly on the random effect variances
                            prior(normal(0, 0.25), class = "sd", group = "id", lb = 0),
                            prior(normal(0, 1), class = "sds", lb = 0)), ## basis weights get treated like a random EFFECT!!!
                  control = list(adapt_delta = 0.95, max_treedepth = 12),
                  chains = 4, cores = 4, seed = 420,
                  iter = 4000, warmup = 2000, backend = "cmdstanr")

temp_spp_model <- brm(dead ~ s(day_z, by = temp_spp, k = 20) + (1|id), 
                  data = amd, 
                  family = bernoulli(link = "cloglog"), 
                  prior = c(prior(normal(0, 0.5), class = "b"), 
                            prior(normal(0, 0.5), class = "Intercept"),
                            ## Regularise more strongly on the random effect variances
                            prior(normal(0, 0.25), class = "sd", group = "id", lb = 0),
                            prior(normal(0, 1), class = "sds", lb = 0)), ## basis weights get treated like a random EFFECT!!!
                  control = list(adapt_delta = 0.95, max_treedepth = 12),
                  chains = 4, cores = 4, seed = 420,
                  iter = 4000, warmup = 2000, backend = "cmdstanr")
#_______________________________________________________________________________
#### 4. Comparison ####

base_model <- add_criterion(base_model, criterion = c("loo","waic"))
spp_model <- add_criterion(spp_model, criterion = c("loo","waic"))
temp_model <- add_criterion(temp_model, criterion = c("loo","waic"))
temp_spp_model <- add_criterion(temp_spp_model, criterion = c("loo","waic"))

## Too big for github but for offline
# save(base_model, spp_model, temp_model, temp_spp_model,
#      file = "data/mortality_models.RData")

load("data/mortality_models.RData")

as.data.frame(loo_compare(base_model, spp_model, temp_model, temp_spp_model, criterion = "loo")) 

as.data.frame(loo_compare(base_model, spp_model, temp_model, temp_spp_model, criterion = "loo")) %>% 
  mutate(Model = c("Species only", "Temperature and species model", "Temperature", "Base")) %>%
  dplyr::select(Model, elpd = elpd_loo, elpd_diff,
                `elpd se` = se_elpd_loo, `LOOIC` = looic) %>%
  mutate(across(elpd:LOOIC, \(x) round(x, 1))) %>%
  flextable(cwidth = 1) %>%
  width(j = 1, width = 3) %>%
  compose(part = "header", j = "elpd_diff", value = as_paragraph("\U0394", "elpd")) %>%
  add_header_row(values = "Mortality", colwidths = 5) %>% 
  italic(i = 1, part = "header") %>% 
  theme_zebra(odd_header = alpha("grey90", alpha = 0.3)) %>%
  save_as_image("output/model_comparison_mortality.png",res = 600)

#_______________________________________________________________________________
#### 5. Predictions ####

# Pull out the scaling attributes that scale() stored
day_mean <- attr(amd$day_z, "scaled:center")
day_sd   <- attr(amd$day_z, "scaled:scale")

# Make newdata grid over the full range of days you observed
newdat <- expand_grid(day = 1:max(amd$day),
                       temp_spp = unique(amd$temp_spp)) %>% 
  mutate(day_z = (day - day_mean) / day_sd,
         spp = paste0(sapply(strsplit(temp_spp, "_"), `[`, 2), "_",
                      sapply(strsplit(temp_spp, "_"), `[`, 3)),
         temperature = sapply(strsplit(temp_spp, "_"), `[`, 1))

# Posterior expected probability of mortality *on that day*
# (for cloglog–Bernoulli this is already the discrete hazard)
haz_draws <- posterior_epred(
  temp_spp_model,
  newdata    = newdat,
  re_formula = NA    # marginalise over random effects
)     

## 6.1 Fitted data plot
mortality_prediciton <- newdat %>% 
  mutate(fit = apply(haz_draws, 2, median),
         lwr = apply(haz_draws, 2, quantile, prob = c(0.025)),
         upr = apply(haz_draws, 2, quantile, prob = c(0.975))) %>% 
  mutate(temperature_lab = paste0(temperature, "\u00B0C"),
         spp_lab = gsub("_", " ", spp)) 

## raw_summary
mort_raw <- antlion_mortality %>%
  mutate(temperature_lab = paste0(temperature, "\u00B0C"),
         day_5 = day - (day %% 5)) %>% 
  filter(temperature != "25") %>% 
  group_by(day, spp, temperature_lab) %>% 
  summarise(mne = mean(dead), sem = (sd(dead)/sqrt(n()))) 

mortality_pred_plot <- 
  ggplot(mortality_prediciton, aes(x = day, y = fit, group = temp_spp)) +
  geom_ribbon(aes(ymax = upr, ymin = lwr, colour = NULL, fill = temperature_lab), 
              alpha = 0.2, show.legend = F) +
  geom_line(aes(colour = temperature_lab), size = 1) +
  facet_wrap(~ spp_lab, ncol = 2) +
  guides(colour = guide_legend(override.aes = list(linewidth = 2))) +
  scale_x_continuous(limits = c(1,76)) +
  scale_colour_viridis_d(option = "A", begin = 0.2, 
                         end = 0.8, aesthetics = c("colour", "fill")) +
  labs(x = "Experiment day", y = "Expected probability\nof mortality", 
       colour = "Temperature\ntreatment") +
  scale_y_continuous(limits = c(0,0.03), breaks= seq(0,0.03, by = 0.005)) +
  theme_test(base_size = 13) +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "italic"))

#_______________________________________________________________________________
#### 7. Combining plots ####

ggsave(mortality_pred_plot,
       filename = "output/mortality_predictions_may25.jpeg",
       width = 27, height = 15, units = "cm", dpi = 600)

# ggsave(mortality_pred_plot,
#        filename = "../mortality_predictions_may25.tiff",
#        width = 27, height = 15, units = "cm", dpi = 600)
