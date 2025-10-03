#####################################################
##                                                 ##
##             Antlion temperature study           ##
##                                                 ##
##                 Emergence models                ##
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
library(AutoScore)

#_______________________________________________________________________________
#### 1. Loading data ####

load("data/antlion_lh.RData", verbose = TRUE)

## exploration plot
antlion_emergence %>%
  filter(temperature != "25") %>% 
  group_by(day, temperature) %>% 
  summarise(mne = mean(emerged), sem = (sd(emerged)/sqrt(n()))) %>% 
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
  brm(emerged ~ s(day_z, by = temperature, k = 20), 
                  data = antlion_emergence, 
                  family = bernoulli(link = "cloglog"), 
                  prior = c(prior(normal(0, 0.5), class = "b"), 
                            prior(normal(0, 0.5), class = "Intercept"),
                            prior(normal(0, 1), class = "sds", lb = 0)),
                  sample_prior = "only", 
                  chains = 2, cores = 2, seed = 420,
                  iter = 2000, warmup = 1000)

## Predict from the priors
grid <- expand_grid(
  temperature = unique(antlion_emergence$temperature),
  day_z       = seq(-1.1, 3.6, length = 120)
)

haz_prior <- add_epred_draws(model_prior_only, newdata = grid, re_formula = NA)  # hazard draws

haz_prior %>% 
  mutate(day = day_z * attr(antlion_emergence$day_z, "scaled:scale") +
           attr(antlion_emergence$day_z, "scaled:center")) %>% 
  ggplot(aes(day, .epred, group = .draw)) +
  geom_line(alpha = .05) +
  facet_wrap(~temperature) +
  labs(y = "Daily emergence hazard (prior)", title = "Prior-predictive curves")

#_______________________________________________________________________________
#### 3. Core models ####

aed <- antlion_emergence %>%
  filter(temperature != "25")

base_model <- brm(emerged ~ s(day_z, k = 20) + spp + (1|id), 
                  data = aed, 
                  family = bernoulli(link = "cloglog"), 
                  prior = c(prior(normal(0, 0.5), class = "b"), 
                            prior(normal(0, 0.5), class = "Intercept"), 
                            ## Regularise more strongly on the random effect variances
                            prior(normal(0, 0.25), class = "sd", group = "id", lb = 0),
                            prior(normal(0, 1), class = "sds", lb = 0)), ## basis weights get treated like a random EFFECT!!! This was the key
                  control = list(adapt_delta = 0.95, max_treedepth = 12),
                  chains = 4, cores = 4, seed = 420,
                  iter = 4000, warmup = 2000, backend = "cmdstanr")

spp_model <- brm(emerged ~ s(day_z, by = spp, k = 20) + (1|id), 
                  data = aed, 
                  family = bernoulli(link = "cloglog"), 
                  prior = c(prior(normal(0, 0.5), class = "b"), 
                            prior(normal(0, 0.5), class = "Intercept"), 
                            ## Regularise more strongly on the random effect variances
                            prior(normal(0, 0.25), class = "sd", group = "id", lb = 0),
                            prior(normal(0, 1), class = "sds", lb = 0)), ## basis weights get treated like a random EFFECT!!!
                  control = list(adapt_delta = 0.95, max_treedepth = 12),
                  chains = 4, cores = 4, seed = 420,
                  iter = 4000, warmup = 2000, backend = "cmdstanr")

temp_model <- brm(emerged ~ s(day_z, by = temperature, k = 20) + spp + (1|id), 
                  data = aed, 
                  family = bernoulli(link = "cloglog"), 
                  prior = c(prior(normal(0, 0.5), class = "b"), 
                            prior(normal(0, 0.5), class = "Intercept"),
                            ## Regularise more strongly on the random effect variances
                            prior(normal(0, 0.25), class = "sd", group = "id", lb = 0),
                            prior(normal(0, 1), class = "sds", lb = 0)), ## basis weights get treated like a random EFFECT!!!
                  control = list(adapt_delta = 0.95, max_treedepth = 12),
                  chains = 4, cores = 4, seed = 420,
                  iter = 4000, warmup = 2000, backend = "cmdstanr")

temp_spp_model <- brm(emerged ~ s(day_z, by = temp_spp, k = 20) + (1|id), 
                  data = aed, 
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

# save(base_model, spp_model, temp_model, temp_spp_model,
#      file = "data/emergence_models.RData")

load("data/emergence_models.RData")

as.data.frame(loo_compare(base_model, spp_model, temp_model, temp_spp_model, criterion = "loo")) 

as.data.frame(loo_compare(base_model, spp_model, temp_model, temp_spp_model, criterion = "loo")) %>% 
  mutate(Model = c("Temperature", "Base", "Species only", "Temperature and species model")) %>%
  dplyr::select(Model, elpd = elpd_loo, elpd_diff,
                `elpd se` = se_elpd_loo, `LOOIC` = looic) %>%
  mutate(across(elpd:LOOIC, \(x) round(x, 1))) %>%
  flextable(cwidth = 1) %>%
  width(j = 1, width = 3) %>%
  compose(part = "header", j = "elpd_diff", value = as_paragraph("\U0394", "elpd")) %>%
  add_header_row(values = "Emergence", colwidths = 5) %>% 
  italic(i = 1, part = "header") %>% 
  theme_zebra(odd_header = alpha("grey90", alpha = 0.3)) %>%
  save_as_image("output/model_comparison_emergence.png",res = 600)

#_______________________________________________________________________________
#### 5. Divergent transition exploration - with help from CHAD ####

# ## temp model
# draws_temp <- as_draws_df(temp_model)          # Draws from posterior in tidy form
# nuts_temp  <- bayesplot::nuts_params(temp_model$fit)      # status of all of the draws
# 
# div_temp   <- nuts_temp %>% 
#   filter(Parameter == "divergent__") %>%          # keep just that diagnostic
#   arrange(Chain, Iteration) %>%                   # ensure draw order
#   pull(Value) == 1                
# 
# div_score <- sapply(draws_temp, function(x) {
#   if (!is.numeric(x)) return(NA_real_)
#   # share of divergent draws that fall in the extreme 1 % tail of the parameter
#   mean(div_temp & x > quantile(x, .99), na.rm = TRUE) +
#     mean(div_temp & x < quantile(x, .01), na.rm = TRUE)
# })
# 
# tail(sort(div_score, decreasing = TRUE), 20)   # top-20 culprits

# ### Which parameters are implicated by
# bayesplot::mcmc_parcoord(temp_model$fit, 
#                          pars = names(temp_model$fit)[1:10],    # first few pars
#                          np = nuts_params(temp_model$fit))  

#_______________________________________________________________________________
#### 6. Predictions ####

# Pull out the scaling attributes that scale() stored
day_mean <- attr(aed$day_z, "scaled:center")
day_sd   <- attr(aed$day_z, "scaled:scale")

# Make newdata grid over the full range of days you observed
newdat <- expand_grid(day = 1:max(aed$day),
                       temperature = unique(aed$temperature),
                       spp = unique(aed$spp)) %>% 
  mutate(day_z = (day - day_mean) / day_sd)

# Posterior expected probability of emergence *on that day*
# (for cloglog–Bernoulli this is already the discrete hazard)
haz_draws <- posterior_epred(
  temp_model,
  newdata    = newdat,
  re_formula = NA    # marginalise over random effects
)     

## 6.1 Fitted data plot
emergence_prediciton <- newdat %>% 
  mutate(fit = apply(haz_draws, 2, median),
         lwr = apply(haz_draws, 2, quantile, prob = c(0.025)),
         upr = apply(haz_draws, 2, quantile, prob = c(0.975))) %>% 
  mutate(temperature_lab = paste0(temperature, "\u00B0C")) %>% 
  filter(spp == "Myrmeleon_almohadarum") 

## raw_summary
emerge_raw <- antlion_emergence %>%
  mutate(temperature_lab = paste0(temperature, "\u00B0C"),
         day_5 = day - (day %% 5)) %>% 
  filter(temperature != "25") %>% 
  group_by(day, temperature_lab) %>% 
  summarise(mne = mean(emerged), sem = (sd(emerged)/sqrt(n()))) 

emergence_pred_plot <- 
  ggplot(emergence_prediciton, aes(x = day, y = fit, colour = temperature_lab)) +
  geom_ribbon(aes(ymax = upr, ymin = lwr, colour = NULL, fill = temperature_lab), alpha = 0.2, show.legend = F) +
  geom_line(size = 1) +
  guides(colour = guide_legend(override.aes = list(linewidth = 2))) +
  scale_x_continuous(limits = c(1,76)) +
  scale_colour_viridis_d(option = "A", begin = 0.2, 
                         end = 0.8, aesthetics = c("colour", "fill")) +
  labs(x = "Experiment day", y = "Expected probability\nof emergence", 
       colour = "Temperature\ntreatment", tag = "a)") +
  scale_y_continuous(limits = c(0,0.4)) +
  theme_test(base_size = 13)

### Peak of emergence for each draw of posterior summarised across temperatures
emergence_peak_dat <- 
  as_tibble(haz_draws) %>% 
  mutate(draw = 1:8000) %>% 
  pivot_longer(-draw) %>% 
  bind_cols(slice(newdat, rep(1:n(), times = 8000))) %>% 
  group_by(draw, temperature, spp) %>% 
  filter(value == max(value)) %>% 
  ungroup() %>% 
  mutate(temperature_lab = paste0(temperature, "\u00B0C")) 

emergence_peak_day_plot <- emergence_peak_dat %>% 
  group_by(temperature_lab, spp) %>% 
  summarise(peak_day = mean(day),
            upr_day = quantile(day, 0.975),
            lwr_day = quantile(day, 0.025)) %>% 
  filter(spp == "Myrmeleon_almohadarum") %>% 
  ggplot(aes(x = temperature_lab, y = peak_day, colour = temperature_lab, fill = temperature_lab)) +
  geom_violinhalf(data = filter(emergence_peak_dat, spp == "Myrmeleon_almohadarum"),
                  aes(y = day), colour = "white", alpha = 0.2, show.legend = F) +
  geom_errorbar(aes(ymin = lwr_day, ymax = upr_day), width = 0.05, show.legend = F) +
  geom_point(size = 3, show.legend = F) +
  scale_colour_viridis_d(option = "A", begin = 0.2, end = 0.8, 
                         aesthetics = c("fill", "colour")) +
  scale_y_continuous(limits = c(1,76)) +
  labs(x = "Temperature treatment", y = "Day of peak emergence", tag = "c)") +
  coord_flip() +
  theme_test(base_size = 13)

emergence_pred_summary <- emergence_prediciton %>% 
  mutate(temperature_lab = paste0(temperature, "\u00B0C")) %>% 
  filter(spp == "Myrmeleon_almohadarum") %>% 
  group_by(temperature_lab) %>% 
  summarise(fit = mean(fit),
           upr = mean(upr),
           lwr = mean(lwr))

## Emergence peak values plot
emergence_peak_value_plot <- 
  emergence_peak_dat %>% 
  group_by(temperature_lab, spp) %>% 
  summarise(peak_value = mean(value),
            upr_value = quantile(value, 0.975),
            lwr_value = quantile(value, 0.025)) %>% 
  filter(spp == "Myrmeleon_almohadarum") %>% 
  ggplot(aes(x = temperature_lab, y = peak_value, colour = temperature_lab, fill = temperature_lab)) +
  geom_violinhalf(data = filter(emergence_peak_dat, spp == "Myrmeleon_almohadarum"),
                  aes(y = value), colour = "white", alpha = 0.2, show.legend = F) +
  geom_errorbar(aes(ymin = lwr_value, ymax = upr_value), width = 0.05, show.legend = F) +
  geom_point(aes(shape = "Peak"), size = 3) +
  geom_errorbar(data = emergence_pred_summary,
                aes(y = NULL, ymin = lwr, ymax = upr), width = 0.05, show.legend = F) +
  geom_point(data = emergence_pred_summary,
             aes(y = fit, shape = "Average"), size = 3) +
  scale_colour_viridis_d(option = "A", begin = 0.2, end = 0.8, 
                         aesthetics = c("fill", "colour"), guide = "none") +
  scale_y_continuous(limits = c(0,0.4)) +
  labs(x = "Temperature treatment", y = "Daily probability\nof emergence", shape = "", tag = "b)") +
  theme_test(base_size = 13)

## Supplementary plot with raw data too
emergence_pred_plot_raw <- 
  ggplot(emergence_prediciton, aes(x = day, y = fit, colour = temperature_lab)) +
  geom_ribbon(aes(ymax = upr, ymin = lwr, colour = NULL, fill = temperature_lab), alpha = 0.2, show.legend = F) +
  geom_point(data = emerge_raw, aes(x = day, y = mne, size = sem), 
             show.legend = F, alpha = 0.6, shape = 16) +
  geom_line(size = 1) +
  guides(colour = guide_legend(override.aes = list(linewidth = 2))) +
  scale_x_continuous(limits = c(1,76)) +
  scale_size_continuous(range = c(1,3)) +
  scale_colour_viridis_d(option = "A", begin = 0.2, 
                         end = 0.8, aesthetics = c("colour", "fill")) +
  labs(x = "Experiment day", y = "Expected probability\nof emergence", 
       colour = "Temperature\ntreatment") +
  scale_y_continuous(limits = c(0,0.4)) +
  theme_test(base_size = 13)

#_______________________________________________________________________________
#### 7. Combining plots ####

(emergence_pred_plot + emergence_peak_value_plot) /
  (emergence_peak_day_plot)

ggsave(emergence_pred_plot + emergence_peak_value_plot +
         emergence_peak_day_plot +
         plot_layout(ncol = 2, widths = c(3,1), heights = c(2.3,1)),
       filename = "output/emergence_predictions_may25.jpeg",
       width = 33, height = 15, units = "cm", dpi = 600)

### For SI
ggsave(emergence_pred_plot_raw,
       filename = "output/emergence_prediction_raw_may25.jpeg",
       width = 19, height = 10, units = "cm", dpi = 600)


#_______________________________________________________________________________
#### 8. Reporting ####

as.data.frame(loo_compare(base_model, spp_model, temp_model, temp_spp_model, criterion = "loo")) 

temp_model


newdat %>% 
  slice(rep(1:n(), times = 8000)) %>% 
  mutate(post = c(t(haz_draws))) %>% 
  group_by(spp) %>% 
  summarise(mn = mean(post), 
            upr = quantile(post, 0.975), 
            lwr = quantile(post, 0.025))



