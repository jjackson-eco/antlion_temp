#####################################################
##                                                 ##
##             Antlion temperature study           ##
##                                                 ##
##                 Pupation models                 ##
##                                                 ##
##                 JJ - 22/05/2025                 ##
##                                                 ##
#####################################################

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

## exploration
antlion_pupation %>%
  filter(temperature != "25") %>% 
  group_by(day) %>% 
  summarise(mnp = mean(pupated), sem = (sd(pupated)/sqrt(n()))) %>% 
  ggplot(aes(x = day, y = mnp)) +
  geom_errorbar(aes(ymax = mnp + sem, ymin = mnp - sem), width = 0.01) +
  geom_point()

#_______________________________________________________________________________
#### 2. Core models ####

apd <- antlion_pupation %>%
  filter(temperature != "25")

base_model <- brm(pupated ~ s(day_z, k = 30) + spp + (1|id), 
                  data = apd, 
                  family = bernoulli(link = "cloglog"), 
                  prior = c(prior(normal(0, 0.5), class = "b"), 
                            prior(normal(0, 0.5), class = "Intercept"), 
                            ## Regularise more strongly on the random effect variances
                            prior(normal(0, 0.25), class = "sd", group = "id", lb = 0),
                            prior(normal(0, 1), class = "sds", lb = 0)), ## basis weights get treated like a random EFFECT!!! This was the key
                  control = list(adapt_delta = 0.95, max_treedepth = 12),
                  chains = 4, cores = 4, seed = 420,
                  iter = 4000, warmup = 2000, backend = "cmdstanr")

spp_model <- brm(pupated ~ s(day_z, by = spp, k = 20) + (1|id), 
                  data = apd, 
                  family = bernoulli(link = "cloglog"), 
                  prior = c(prior(normal(0, 0.5), class = "b"), 
                            prior(normal(0, 0.5), class = "Intercept"), 
                            ## Regularise more strongly on the random effect variances
                            prior(normal(0, 0.25), class = "sd", group = "id", lb = 0),
                            prior(normal(0, 1), class = "sds", lb = 0)), ## basis weights get treated like a random EFFECT!!!
                  control = list(adapt_delta = 0.95, max_treedepth = 12),
                  chains = 4, cores = 4, seed = 420,
                  iter = 4000, warmup = 2000, backend = "cmdstanr")

temp_model <- brm(pupated ~ s(day_z, by = temperature, k = 20) + spp + (1|id), 
                  data = apd, 
                  family = bernoulli(link = "cloglog"), 
                  prior = c(prior(normal(0, 0.5), class = "b"), 
                            prior(normal(0, 0.5), class = "Intercept"),
                            ## Regularise more strongly on the random effect variances
                            prior(normal(0, 0.25), class = "sd", group = "id", lb = 0),
                            prior(normal(0, 1), class = "sds", lb = 0)), ## basis weights get treated like a random EFFECT!!!
                  control = list(adapt_delta = 0.95, max_treedepth = 12),
                  chains = 4, cores = 4, seed = 420,
                  iter = 4000, warmup = 2000, backend = "cmdstanr")

temp_spp_model <- brm(pupated ~ s(day_z, by = temp_spp, k = 20) + (1|id), 
                  data = apd, 
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

## Too big for github but for offline files
# save(base_model, spp_model, temp_model, temp_spp_model,
#      file = "data/pupation_models.RData")

load("data/pupation_models.RData")

### REMEMBER TO ALWAYS COMPARE TO THIS WHEN MAKING A TABLE WITH MANUAL LABELS
as.data.frame(loo_compare(base_model, spp_model, temp_model, temp_spp_model, criterion = "loo")) 

as.data.frame(loo_compare(base_model, spp_model, temp_model, temp_spp_model, criterion = "loo")) %>% 
  mutate(Model = c("Base","Species only", "Temperature and species model", "Temperature")) %>%
  dplyr::select(Model, elpd = elpd_loo, elpd_diff,
                `elpd se` = se_elpd_loo, `LOOIC` = looic) %>%
  mutate(across(elpd:LOOIC, \(x) round(x, 1))) %>%
  flextable(cwidth = 1) %>%
  width(j = 1, width = 3) %>%
  compose(part = "header", j = "elpd_diff", value = as_paragraph("\U0394", "elpd")) %>%
  add_header_row(values = "Pupation", colwidths = 5) %>% 
  italic(i = 1, part = "header") %>% 
  theme_zebra(odd_header = alpha("grey90", alpha = 0.3)) %>%
  save_as_image("output/model_comparison_pupation.png", res = 600)

#_______________________________________________________________________________
#### 5. Predictions ####

# Pull out the scaling attributes that scale() stored
day_mean <- attr(apd$day_z, "scaled:center")
day_sd   <- attr(apd$day_z, "scaled:scale")

# Make newdata grid over the full range of days you observed
newdat <- expand_grid(day = 1:max(apd$day),
                       spp = unique(apd$spp)) %>% 
  mutate(day_z = (day - day_mean) / day_sd)

# Posterior expected probability of pupation *on that day*
# (for cloglog–Bernoulli this is already the discrete hazard)
haz_draws <- posterior_epred(
  base_model,
  newdata    = newdat,
  re_formula = NA    # marginalise over random effects
)     

## 6.1 Fitted data plot
pupation_prediciton <- newdat %>% 
  mutate(fit = apply(haz_draws, 2, median),
         lwr = apply(haz_draws, 2, quantile, prob = c(0.025)),
         upr = apply(haz_draws, 2, quantile, prob = c(0.975))) %>% 
  filter(spp == "Myrmeleon_almohadarum") 

pupation_pred_plot <- 
  ggplot(pupation_prediciton, aes(x = day, y = fit)) +
  geom_ribbon(aes(ymax = upr, ymin = lwr), alpha = 0.2, show.legend = F) +
  geom_line(size = 1) +
  scale_x_continuous(limits = c(1,61)) +
  labs(x = "Experiment day", y = "Expected probability\nof pupation") +
  scale_y_continuous(limits = c(0,0.4)) +
  theme_test(base_size = 13)

ggsave(pupation_pred_plot,
       filename = "output/pupation_predictions_may25.jpeg",
       width = 15, height = 13, units = "cm", dpi = 600)
#_______________________________________________________________________________
#### 6. Reporting ####

newdat %>% 
  slice(rep(1:n(), times = 8000)) %>% 
  mutate(post = c(t(haz_draws))) %>% 
  group_by(spp) %>% 
  summarise(mn = mean(post), 
            upr = quantile(post, 0.975), 
            lwr = quantile(post, 0.025))
