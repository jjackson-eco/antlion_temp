#####################################################
##                                                 ##
##             Antlion temperature study           ##
##                                                 ##
##                Pit diameter change              ##
##                                                 ##
##                  JJ - 28/05/2025                ##
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
#### 1. Loading data ####

load("data/pit_diameter_timeseries.RData")

pd_rawplot <- pit_diameter_timeseries %>% 
  mutate(spp_lab = gsub("_", " ", spp),
         temperature_lab = paste0(temperature, "\u00B0C")) %>%
  ggplot(aes(x = day, y = pit_diameter, group = id)) +
  geom_smooth(aes(colour = temperature_lab), linewidth = 0.8, se = F, show.legend = F) +
  facet_grid(temperature_lab ~ spp_lab) +
  coord_cartesian(ylim = c(0,100)) + 
  scale_colour_viridis_d(option = "A", begin = 0.2, end = 0.8) +
  labs(x = "Experiment day", y = "Pit diameter (mm)") +
  theme_test(base_size = 13) +
  theme(strip.text.x = element_text(face = "italic"),
        strip.background = element_blank())

# pit_diameter_timeseries %>% 
#   filter(pit_diameter != 0) %>% 
#   ggplot(aes(x = temperature, y = pit_diameter)) +
#   geom_violin(aes(fill = temperature), show.legend = F) 

ggsave(pd_rawplot, filename = "output/pit_diameter_raw_smooth.jpeg",
       width = 21, height = 22, units = "cm", dpi = 600)

#_______________________________________________________________________________
#### 2. Core models ####

apd <- pit_diameter_timeseries %>%
 mutate(day_z = scale(day),
        pd = scale(pit_diameter),
        id_temp = paste0(id, "_", temperature))

base_model <- brm(pd ~ s(day_z, k = 20) + spp + (1|id), 
                  data = apd, 
                  family = gaussian(), 
                  prior = c(prior(normal(0, 0.5), class = "b"), 
                            prior(normal(0, 0.5), class = "Intercept"), 
                            prior(normal(0, 1), class = "sds", lb = 0)),
                  control = list(adapt_delta = 0.95, max_treedepth = 12),
                  chains = 4, cores = 4, seed = 420,
                  iter = 4000, warmup = 2000, backend = "cmdstanr")


temp_model  <- brm(pd ~ s(day_z, by = temperature, k = 20) + spp + (1|id), 
                         data = apd, family = gaussian(), 
                         prior = c(prior(normal(0, 0.5), class = "b"), 
                                   prior(normal(0, 0.5), class = "Intercept"), 
                                   prior(normal(0, 0.5), class = "sds", lb = 0)),
                         control = list(adapt_delta = 0.95, max_treedepth = 12),
                         chains = 4, cores = 4, seed = 420,
                         iter = 4000, warmup = 2000, backend = "cmdstanr")

temp_model_rs <- brm(pd ~ s(day_z, by = temperature, k = 20) + spp + (1 + temperature|id), 
                     data = apd, family = gaussian(), 
                     prior = c(prior(normal(0, 0.5), class = "b"), 
                               prior(normal(0, 0.5), class = "Intercept"), 
                               prior(normal(0, 0.5), class = "sds", lb = 0)),
                     control = list(adapt_delta = 0.97, max_treedepth = 12),
                     chains = 4, cores = 4, seed = 420,
                     iter = 4000, warmup = 2000, backend = "cmdstanr")

id_model <- brm(pd ~ s(day_z, by = id, k = 20) + spp, 
                         data = apd, family = gaussian(), 
                         prior = c(prior(normal(0, 0.5), class = "b"), 
                                   prior(normal(0, 0.5), class = "Intercept"), 
                                   prior(normal(0, 0.5), class = "sds", lb = 0)),
                         control = list(adapt_delta = 0.95, max_treedepth = 12),
                         chains = 4, cores = 4, seed = 420,
                         iter = 4000, warmup = 2000, backend = "cmdstanr")


base_model <- add_criterion(base_model, criterion = c("loo","waic"))
temp_model <- add_criterion(temp_model, criterion = c("loo","waic"))
temp_model_rs <- add_criterion(temp_model_rs, criterion = c("loo","waic"))
id_model <- add_criterion(id_model, criterion = c("loo","waic"))

save(base_model, temp_model, temp_model_rs, id_model, file = "data/pitdiameter_models.RData")

#_______________________________________________________________________________
#### 3. Model comparison table ####

load("data/pitdiameter_models.RData", verbose = T)

as.data.frame(loo_compare(base_model, temp_model, temp_model_rs, id_model, criterion = "loo"))

as.data.frame(loo_compare(base_model, temp_model, temp_model_rs, id_model, criterion = "loo")) %>% 
  mutate(Model = c("Temperature", "Temperature random slope", "Base model", "ID-level smooths")) %>%
  dplyr::select(Model, elpd = elpd_loo, elpd_diff,
                `elpd se` = se_elpd_loo, `LOOIC` = looic) %>%
  mutate(across(elpd:LOOIC, \(x) round(x, 1))) %>%
  flextable(cwidth = 1) %>%
  width(j = 1, width = 3) %>%
  compose(part = "header", j = "elpd_diff", value = as_paragraph("\U0394", "elpd")) %>%
  add_header_row(values = "Pit diameter", colwidths = 5) %>% 
  italic(i = 1, part = "header") %>% 
  theme_zebra(odd_header = alpha("grey90", alpha = 0.3)) %>%
  save_as_image("output/model_comparison_pitdiameter.png",res = 600)

#_______________________________________________________________________________
#### 4. Predictions ####

# Pull out the scaling attributes that scale() stored
day_mean <- attr(apd$day_z, "scaled:center")
day_sd   <- attr(apd$day_z, "scaled:scale")

pd_mean <- attr(apd$pd, "scaled:center")
pd_sd   <- attr(apd$pd, "scaled:scale")

## maximum observation day by temperature
max_days <- apd %>% 
  group_by(temperature) %>% 
  summarise(max_day = max(day))

# Make newdata grid over the full range of days you observed
newdat <- expand_grid(day = 1:max(apd$day),
                      temperature = unique(apd$temperature),
                      spp = unique(apd$spp)) %>% 
  mutate(day_z = (day - day_mean) / day_sd)

# Posterior expected probability of emergence *on that day*
# (for cloglog–Bernoulli this is already the discrete hazard)
haz_draws <- posterior_epred(
  temp_model,
  newdata    = newdat,
  re_formula = NA    # marginalise over random effects
)     

## 6.1 Fitted data plot
pitdiameter_prediction <- newdat %>% 
  mutate(fit = apply(haz_draws, 2, median),
         lwr = apply(haz_draws, 2, quantile, prob = c(0.025)),
         upr = apply(haz_draws, 2, quantile, prob = c(0.975))) %>% 
  ## Back-transform response variable
  mutate(across(fit:upr, ~ (.x*pd_sd) + pd_mean)) %>% 
  mutate(temperature_lab = paste0(temperature, "\u00B0C")) %>% 
  left_join(x = ., y = max_days) %>% 
  filter(spp == "Myrmeleon_almohadarum" & day <= max_day & fit >= 0)
  
pitdiameter_pred_plot <- 
  ggplot(pitdiameter_prediction, aes(x = day, y = fit, colour = temperature_lab)) +
  geom_ribbon(aes(ymax = upr, ymin = lwr, colour = NULL, fill = temperature_lab), alpha = 0.1, show.legend = F) +
  geom_line(linewidth = 1) +
  guides(colour = guide_legend(override.aes = list(linewidth = 2))) +
  scale_colour_viridis_d(option = "A", begin = 0.2, 
                         end = 0.8, aesthetics = c("colour", "fill")) +
  labs(x = "Experiment day", y = "Pit diameter (mm)", 
       colour = "Temperature\ntreatment", tag = "a)") +
  coord_cartesian(ylim = c(0,50)) +
  theme_test(base_size = 13)

ggsave(pitdiameter_pred_plot, filename = "output/pitdiameter_prediction.jpeg",
       width = 22, height = 11, units = "cm", dpi = 600)

#_______________________________________________________________________________
#### 5. Time with a pit trap  ####

time_pit_ind <- apd %>% 
  group_by(id) %>% 
  filter(pit_diameter != 0) %>% 
  summarise(spp = spp[1], temperature = temperature[1],
            days_pit = n())

base_model_d <- brm(days_pit ~ spp + (1|id),
                  data = time_pit_ind, family = poisson(), 
                  prior = c(prior(normal(0, 0.5), class = "b"), 
                            prior(normal(0, 0.5), class = "Intercept")),
                  control = list(adapt_delta = 0.97, max_treedepth = 12),
                  chains = 4, cores = 4, seed = 420,
                  iter = 4000, warmup = 2000, backend = "cmdstanr")

temp_model_d <- brm(days_pit ~ spp + temperature + (1|id),
                  data = time_pit_ind, family = poisson(), 
                  prior = c(prior(normal(0, 0.5), class = "b"), 
                            prior(normal(0, 0.5), class = "Intercept")),
                  control = list(adapt_delta = 0.97, max_treedepth = 12),
                  chains = 4, cores = 4, seed = 420,
                  iter = 4000, warmup = 2000, backend = "cmdstanr")

temp_spp_model_d <- brm(days_pit ~ spp*temperature + (1|id),
                  data = time_pit_ind, family = poisson(), 
                  prior = c(prior(normal(0, 0.5), class = "b"), 
                            prior(normal(0, 0.5), class = "Intercept")),
                  control = list(adapt_delta = 0.97, max_treedepth = 12),
                  chains = 4, cores = 4, seed = 420,
                  iter = 4000, warmup = 2000, backend = "cmdstanr")


base_model_d <- add_criterion(base_model_d, criterion = c("loo","waic"))
temp_model_d <- add_criterion(temp_model_d, criterion = c("loo","waic"))
temp_spp_model_d <- add_criterion(temp_spp_model_d, criterion = c("loo","waic"))

as.data.frame(loo_compare(base_model_d, temp_model_d, temp_spp_model_d, criterion = "loo"))

as.data.frame(loo_compare(base_model_d, temp_model_d, temp_spp_model_d, criterion = "loo")) %>% 
  mutate(Model = c("Temperature and species model", "Temperature model", "Base model")) %>%
  dplyr::select(Model, elpd = elpd_loo, elpd_diff,
                `elpd se` = se_elpd_loo, `LOOIC` = looic) %>%
  mutate(across(elpd:LOOIC, \(x) round(x, 1))) %>%
  flextable(cwidth = 1) %>%
  width(j = 1, width = 3) %>%
  compose(part = "header", j = "elpd_diff", value = as_paragraph("\U0394", "elpd")) %>%
  add_header_row(values = "Pit occurence", colwidths = 5) %>% 
  italic(i = 1, part = "header") %>% 
  theme_zebra(odd_header = alpha("grey90", alpha = 0.3)) %>%
  save_as_image("output/model_comparison_pitoccurence.png",res = 600)

# Make newdata grid over the full range of days you observed
newdat <- expand_grid(temperature = unique(apd$temperature),
                      spp = unique(apd$spp))

# Posterior expected probability of emergence *on that day*
# (for cloglog–Bernoulli this is already the discrete hazard)
haz_draws <- posterior_epred(
  temp_spp_model_d,
  newdata    = newdat,
  re_formula = NA    # marginalise over random effects
)     


post_full <- as.data.frame(haz_draws) %>% 
  mutate(draw = 1:n()) %>% 
  pivot_longer(-draw) %>% 
  bind_cols(slice(newdat, rep(1:n(), 8000))) %>% 
  mutate(temperature_lab = paste0(temperature, "\u00B0C"),
         spp = gsub("_", " ", spp))

## 6.1 Fitted data plot
pitdays_prediction <- newdat %>% 
  mutate(fit = apply(haz_draws, 2, median),
         lwr = apply(haz_draws, 2, quantile, prob = c(0.025)),
         upr = apply(haz_draws, 2, quantile, prob = c(0.975))) %>% 
  mutate(temperature_lab = paste0(temperature, "\u00B0C"),
         spp = gsub("_", " ", spp))

pitdays_predplot <- ggplot(pitdays_prediction, aes(x = temperature_lab, colour = temperature_lab)) + 
  geom_violin(data = post_full, aes(y = value, fill = temperature_lab),
              alpha = 0.2, colour = "white", show.legend = F) +
  geom_errorbar(aes(ymax = upr, ymin = lwr), width = 0.05, show.legend = F) +
  geom_point(aes(y = fit), size = 3, show.legend = F) +
  facet_wrap(~spp) +
  scale_colour_viridis_d(option = "A", begin = 0.2, 
                         end = 0.8, aesthetics = c("colour", "fill")) +
  labs(x = "Temperature treatment", y = "Days with a pit trap", tag = "b)") +
  theme_test(base_size = 13) +
  theme(strip.text = element_text(face = "italic"),
        strip.background = element_blank())

ggsave(pitdiameter_pred_plot / pitdays_predplot, filename = "output/pit_trap_predictions.jpeg",
       width = 22, height = 20, units = "cm", dpi = 600)

