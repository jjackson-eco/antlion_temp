#####################################################
##                                                 ##
##             Antlion temperature study           ##
##                                                 ##
##                 Morphology models               ##
##                                                 ##
##                  JJ - 26/05/2025                ##
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

imago_morph_data <- antlion_lh_raw %>% 
  dplyr::select(ID, spp = sp, sex, stage_capture = stage, temperature,
                body_length_I:condition_I) %>% 
  na.omit() %>% 
  filter(condition_I != "undev") %>% 
  mutate(across(body_length_I:weight_I, .fns = scale)) %>% 
  rename(body_length = body_length_I,
         wing_length = length_front_wing_I,
         weight = weight_I)

#_______________________________________________________________________________
#### 2. Exploratory plots ####

imago_morph_data %>% 
  pivot_longer(body_length:weight) %>% 
  ggplot(aes(x = value)) +
  geom_density(fill = "grey80") +
  facet_grid(name ~ spp) +
  theme_test() +
  theme(strip.background = element_rect(fill = "white"))

imago_morph_data %>% 
  pivot_longer(body_length:weight) %>% 
  ggplot(aes(x = value)) +
  geom_density(fill = "grey80") +
  facet_grid(name ~ sex) +
  theme_test() +
  theme(strip.background = element_rect(fill = "white"))

#_______________________________________________________________________________
#### 3. Core models ####

## 3a. Regression components - Base
weight_base <- bf(weight ~ spp*sex)
wing_base <- bf(wing_length ~ spp*sex)
body_base <- bf(body_length ~ spp*sex)

base_formula <- weight_base + wing_base + body_base + set_rescor(TRUE)

## 3b. Regression components - Temperature simple
weight_temp <- bf(weight ~ spp*sex + temperature)
wing_temp <- bf(wing_length ~ spp*sex + temperature)
body_temp <- bf(body_length ~ spp*sex + temperature)

temp_formula <- weight_temp + wing_temp + body_temp + set_rescor(TRUE)

## 3c. Regression components - Temperature + sex
weight_temp_sex <- bf(weight ~ spp*sex + temperature*sex)
wing_temp_sex <- bf(wing_length ~ spp*sex + temperature*sex)
body_temp_sex <- bf(body_length ~ spp*sex + temperature*sex)

temp_sex_formula <- weight_temp_sex + wing_temp_sex + body_temp_sex + set_rescor(TRUE)

## 3d. Regression components - Temperature + species
weight_temp_spp <- bf(weight ~ spp*sex + temperature*spp)
wing_temp_spp <- bf(wing_length ~ spp*sex + temperature*spp)
body_temp_spp <- bf(body_length ~ spp*sex + temperature*spp)

temp_spp_formula <- weight_temp_spp + wing_temp_spp + body_temp_spp + set_rescor(TRUE)

## 3e. Priors - general: weakly informative, same scale for all fixed effects
morph_priors <- c(
  prior(normal(0, 0.5), class = "b"),          # temperature effect
  prior(normal(0, 0.5), class = "Intercept"),
  prior(lkj(2), class = "rescor")      # residual correlation LKJ(2) ~ mild shrink
)

## 3f. Models
base_model <- brm(base_formula, data = imago_morph_data, prior  = morph_priors,
                  chains = 4, cores = 4, iter   = 4000, warmup = 2000, seed = 420,
                  control = list(adapt_delta = 0.95, max_treedepth = 12), backend = "cmdstanr")

temp_model <- brm(temp_formula, data = imago_morph_data, prior  = morph_priors,
                  chains = 4, cores = 4, iter   = 4000, warmup = 2000, seed = 420,
                  control = list(adapt_delta = 0.95, max_treedepth = 12), backend = "cmdstanr")

temp_sex_model <- brm(temp_sex_formula, data = imago_morph_data, prior  = morph_priors,
                      chains = 4, cores = 4, iter   = 4000, warmup = 2000, seed = 420,
                      control = list(adapt_delta = 0.95, max_treedepth = 12), backend = "cmdstanr")

temp_spp_model <- brm(temp_spp_formula, data = imago_morph_data, prior  = morph_priors,
                      chains = 4, cores = 4, iter   = 4000, warmup = 2000, seed = 420,
                      control = list(adapt_delta = 0.95, max_treedepth = 12), backend = "cmdstanr")

#_______________________________________________________________________________
#### 3. Comparison models ####

base_model <- add_criterion(base_model, criterion = c("loo","waic"))
temp_model <- add_criterion(temp_model, criterion = c("loo","waic"))
temp_sex_model <- add_criterion(temp_sex_model, criterion = c("loo","waic"))
temp_spp_model <- add_criterion(temp_spp_model, criterion = c("loo","waic"))

## Save (in zip but too big for github repo)
# save(base_model, temp_model, temp_spp_model, temp_sex_model,
#      file = "data/morphology_models.RData")

as.data.frame(loo_compare(base_model, temp_model, temp_spp_model, temp_sex_model, criterion = "loo")) 

as.data.frame(loo_compare(base_model, temp_model, temp_spp_model, temp_sex_model, criterion = "loo")) %>% 
  mutate(Model = c("Base", "Temperature and sex", "Temperature", "Temperature and species")) %>%
  dplyr::select(Model, elpd = elpd_loo, elpd_diff,
                `elpd se` = se_elpd_loo, `LOOIC` = looic) %>%
  mutate(across(elpd:LOOIC, \(x) round(x, 1))) %>%
  flextable(cwidth = 1) %>%
  width(j = 1, width = 3) %>%
  compose(part = "header", j = "elpd_diff", value = as_paragraph("\U0394", "elpd")) %>%
  add_header_row(values = "Imago morphology multivariate", colwidths = 5) %>% 
  italic(i = 1, part = "header") %>% 
  theme_zebra(odd_header = alpha("grey90", alpha = 0.3)) %>%
  save_as_image("output/model_comparison_morphology.png",res = 600)

#_______________________________________________________________________________
#### 4. Summary plot ####

morph_plot <- imago_morph_data %>% 
  pivot_longer(body_length:weight) %>% 
  mutate(name = str_to_sentence(gsub("_", " ", name)),
         sex = str_to_sentence(sex),
         temperature_lab = paste0(temperature, "\u00B0C"),
         spp_lab = gsub("_", " ", spp)) %>% 
  ggplot(aes(x = temperature_lab, y = value)) +
  geom_violin(aes(fill = temperature_lab),
              alpha = 0.2, colour = "white", show.legend = F) +
  geom_jitter(aes(colour = temperature_lab), 
              shape = 16, width = 0.1, alpha = 0.6, show.legend = F) +
  facet_grid(name ~ spp_lab + sex) +
  scale_colour_viridis_d(option = "A", begin = 0.2, 
                         end = 0.8, aesthetics = c("colour", "fill")) +
  labs(x = "Temperature treatment", y = "Scaled morphology value") +
  theme_test() +
  theme(strip.background = element_blank(),
        strip.text.x = element_text(face = "italic"))

ggsave(morph_plot, filename = "output/morphology_summary_plot_may25.jpeg",
       width = 21, height = 14, units = "cm", dpi = 600)
