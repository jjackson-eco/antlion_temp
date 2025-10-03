#####################################################
##                                                 ##
##             Antlion temperature study           ##
##                                                 ##
##         Cleaning + Merging baseline data        ##
##                                                 ##
##                JJ - 12/05/2025                  ##
##                                                 ##
#####################################################
rm(list = ls())
options(width = 100)

library(tidyverse)

#_______________________________________________________________________________
#### 1. Loading data ####

antlion_lh_raw_1 <- read_csv("data/base_lh_1_antlion.csv",  locale = locale(encoding = "Latin1"))
antlion_lh_raw_2 <- read_csv("data/base_lh_2_antlion.csv",  locale = locale(encoding = "Latin1"))

#_______________________________________________________________________________
#### 2. Merging data 2 ####

## renaming columns in data 1

antlion_lh_raw_1 <- rename(antlion_lh_raw_1,
                           stage = Stage, weight = `Weight  (g)`,
                           length = `Length (mm)`, temperature = `t (°C)`,
                           death_days = Numb_days_before_death,
                           weight_dead = `Weight_dead_L_P (g)`,
                           stage_death = `Stage of death`,
                           days_P = Number_days_P, days_I = Number_days_I,
                           sex = Gender, body_length_I = `Body_length_I\n(mm)`,
                           length_front_wing_I = `Length_front_\nwing_I (mm)`,
                           weight_I = `Weight_I (g)`, condition_I = Condition_I)


antlion_lh_raw_2 <- rename(antlion_lh_raw_2,
                           length = length_larva, temperature = `t (°C)`, 
                           days_P = days_to_pupa, days_I = days_to_imago,
                           body_length_I = length_body, length_front_wing_I = length_front_wing) %>% 
  dplyr::select(ID, sp, temperature, length, days_P, days_I, body_length_I, length_front_wing_I, experiment_length) %>% 
  mutate(temperature = as.character(temperature))


antlion_lh_raw <- bind_rows(antlion_lh_raw_1, antlion_lh_raw_2) %>% 
  mutate(sex = case_when(sex == "m" ~ "male",
                         sex == "w" ~ "female")) %>% 
  ## Fixing IDs - repeats from the datasets :o
  mutate(ID = case_when(experiment_length == 99 ~ paste0(ID, "-1"),
                        experiment_length == 76 ~ paste0(ID, "-2"),
                        experiment_length == 64 ~ paste0(ID, "-3")),
         ID = gsub("-", "_", ID))

#_______________________________________________________________________________
#### 3. Exploration figures ####

## 3a. adult development time
imago_exp <- 
  ggplot(antlion_lh_raw, aes(x = temperature, y = days_I, colour = temperature, fill = temperature)) +
  geom_violin(colour = "white", alpha = 0.15, show.legend = F) +
  geom_jitter(alpha = 0.8, width = 0.1, size = 2, show.legend = F) +
  scale_colour_viridis_d(option = "A", begin = 0.2, 
                         end = 0.8, aesthetics = c("colour", "fill")) +
  labs(x = "Temperature (\u00B0C)", y = "Days to Imago (adult)") +
  theme_test(base_size = 14)

ggsave(imago_exp, filename = "output/exploration_may25/imago_temperature_exploration.jpeg",
       width = 16, height = 14, units = "cm", dpi = 300)

## 3b. morphological traits
imago_morph <- antlion_lh_raw %>% 
  dplyr::select(ID, temperature, body_length_I:weight_I) %>% 
  pivot_longer(body_length_I:weight_I) %>% 
  mutate(name = str_to_sentence(gsub("_", " ", gsub("_I", "", name)))) %>% 
  na.omit() %>% 
  ggplot(aes(x = temperature, y = value, colour = temperature, fill = temperature)) +
  geom_violin(colour = "white", alpha = 0.15, show.legend = F) +
  geom_jitter(alpha = 0.8, width = 0.1, size = 2, show.legend = F) +
  facet_wrap(~name, ncol = 3, scales = "free_y") +
  scale_colour_viridis_d(option = "A", begin = 0.2, 
                               end = 0.8, aesthetics = c("colour", "fill")) +
  labs(x = "Temperature (\u00B0C)", y = "Morphological trait (Imago)") +
  theme_test(base_size = 12) +
  theme(strip.background = element_blank())

ggsave(imago_morph, filename = "output/exploration_may25/imago_morphology_temperature_exploration.jpeg",
       width = 30, height = 14, units = "cm", dpi = 400)

jpeg(filename = "output/exploration_may25/morphology_covariance.jpeg",width = 12, height = 12, units = "cm", res = 400)
pairs(antlion_lh_raw[,13:15])
dev.off()

#_______________________________________________________________________________
#### 4. Time-event data ####

## Emergence time-event
antlion_emergence <- 
  antlion_lh_raw %>% 
  # filter(is.na(sex) == F) %>% 
  dplyr::select(id = ID, spp = sp, sex, temperature, death_days, days_I, experiment_length) %>% 
  ## First get correct departure day in one column and status in the other
  mutate(departure = case_when(
                    is.na(days_I) == F ~ days_I,
                    is.na(death_days) == F ~ death_days,
                    TRUE ~ experiment_length),
         status = if_else(is.na(days_I) == T, 0L, 1L)) %>% 
  ## Now do rowwise days up to departure for each
  rowwise() %>% 
  mutate(day = list(seq_len(departure))) %>% # creates a list for each row with your new var of interest
  ungroup() %>% 
  unnest(day) %>% # unnest to expand out by the listed column(s)
  mutate(emerged = as.integer(status == 1L & day == departure)) %>% # get an event based on the day being the day of departure
  mutate(day_z = scale(day, center = TRUE, scale = TRUE)) %>% 
  dplyr::select(id, spp, sex, temperature, day, day_z, emerged) %>% 
  mutate(temp_sex = paste0(temperature,"_", sex),
         temp_spp = paste0(temperature,"_", spp),
         temp_sex_spp = paste0(temperature,"_", sex, "_", spp))

## checking some random ids
filter(antlion_emergence, id == "E-1") %>% 
  as.data.frame()

## Pupation time-event
antlion_pupation <- 
  antlion_lh_raw %>% 
  # filter(is.na(sex) == F) %>% 
  dplyr::select(id = ID, spp = sp, sex, temperature, death_days, days_P, experiment_length) %>% 
  ## First get correct departure day in one column and status in the other
  mutate(departure = case_when(
    is.na(days_P) == F ~ days_P,
    is.na(death_days) == F ~ death_days,
    TRUE ~ experiment_length),
    status = if_else(is.na(days_P) == T, 0L, 1L)) %>% 
  ## Now do rowwise days up to departure for each
  rowwise() %>% 
  mutate(day = list(seq_len(departure))) %>% # creates a list for each row with your new var of interest
  ungroup() %>% 
  unnest(day) %>% # unnest to expand out by the listed column(s)
  mutate(pupated = as.integer(status == 1L & day == departure)) %>% # get an event based on the day being the day of departure
  mutate(day_z = scale(day, center = TRUE, scale = TRUE)) %>% 
  dplyr::select(id, spp, sex, temperature, day, day_z, pupated) %>% 
  mutate(temp_sex = paste0(temperature,"_", sex),
         temp_spp = paste0(temperature,"_", spp),
         temp_sex_spp = paste0(temperature,"_", sex, "_", spp))

## Mortality time-event
antlion_mortality <- 
  antlion_lh_raw %>% 
  # filter(is.na(sex) == F) %>% 
  dplyr::select(id = ID, spp = sp, sex, temperature, death_days, experiment_length) %>% 
  ## First get correct departure day in one column and status in the other
  mutate(departure = if_else(is.na(death_days) == T, experiment_length, death_days),
         status = if_else(is.na(death_days) == T, 0L, 1L)) %>% 
  ## Now do rowwise days up to departure for each
  rowwise() %>% 
  mutate(day = list(seq_len(departure))) %>% # creates a list for each row with your new var of interest
  ungroup() %>% 
  unnest(day) %>% # unnest to expand out by the listed column(s)
  mutate(dead = as.integer(status == 1L & day == departure)) %>% # get an event based on the day being the day of departure
  mutate(day_z = scale(day, center = TRUE, scale = TRUE)) %>% 
  dplyr::select(id, spp, sex, temperature, day, day_z, dead) %>% 
  mutate(temp_sex = paste0(temperature,"_", sex),
         temp_spp = paste0(temperature,"_", spp),
         temp_sex_spp = paste0(temperature,"_", sex, "_", spp))

#_______________________________________________________________________________
#### 5. Extra checks ####

## checking days to pupa and to adult
antlion_lh_check <- antlion_lh_raw %>% 
  mutate(p_na = is.na(days_P),
         i_na = is.na(days_I))

## there are number differences in presence of ecolsion data - treat pupa and imago as separate datasets
table(antlion_lh_check$p_na, antlion_lh_check$i_na)

antlion_pupated_test <- antlion_lh_raw %>% filter(is.na(stage) == F & is.na(days_P) == F) 
antlion_emerged_test <- antlion_lh_raw %>% filter(is.na(stage) == F & is.na(days_I) == F) 

antlion_lh_raw %>% 
  filter(is.na(stage) == F & is.na(days_I) == F) %>% 
  mutate(stage = paste0(stage, "L")) %>% 
  ggplot(aes(x = stage, y = days_I)) +
  geom_violin(aes(fill = stage), alpha = 0.2, show.legend = F, colour = "white") +
  geom_jitter(aes(colour = stage), width = 0.06, alpha = 0.5, size = 2, show.legend = F) +
  scale_color_viridis_d(option = "D", begin = 0.2, end = 0.8, aesthetics = c("colour", "fill")) +
  labs(x = "Larval stage at capture", y = "Days to emergence") +
  theme_test(base_size = 14)


#_______________________________________________________________________________
#### 6. Save ####

save(antlion_lh_raw, antlion_emergence, antlion_pupation, antlion_mortality,
     file = "data/antlion_lh.RData")

#_______________________________________________________________________________
#### 7. Pit diameter timeseries ####

pup_dead <- antlion_lh_raw %>% 
  dplyr::select(ID, days_P, death_days)

pit_1_raw <- read_csv("data/pit_diameter_1_2023.csv",  locale = locale(encoding = "Latin1")) %>% 
  mutate(ID = gsub("-", "_", paste0(ID, "-1")))
pit_2_raw <- read_csv("data/pit_diameter_2_2024.csv",  locale = locale(encoding = "Latin1")) %>% 
  mutate(ID = gsub("-", "_", paste0(ID, "-2")))
pit_3_raw <- read_csv("data/pit_diameter_3_2024.csv",  locale = locale(encoding = "Latin1")) %>% 
  mutate(ID = gsub("-", "_", paste0(ID, "-3")))

pit_1_clean <- pit_1_raw %>% 
  pivot_longer(day1:day37) %>% 
  mutate(day = as.numeric(gsub("day", "", name))) %>% 
  left_join(x = ., y = pup_dead, by = "ID") %>% 
  mutate(day_filter = case_when(is.na(days_P) == F ~ days_P,
                                is.na(death_days) == F ~ death_days,
                                TRUE ~ max(.$day))) %>% 
  filter(day < day_filter, # has to be before the day of death
         is.na(value) == F) %>% 
  dplyr::select(id = ID, spp = sp, temperature = `t (¡C)`, day, pit_diameter = value) %>% 
  mutate(temperature = as.character(temperature)) 

pit_2_clean <- pit_2_raw %>% 
  mutate(across(`Day 1`:`Day 76`, as.character)) %>% 
  pivot_longer(`Day 1`:`Day 76`) %>% 
  mutate(day = as.numeric(gsub("Day ", "", name))) %>% 
  left_join(x = ., y = pup_dead, by = "ID") %>% 
  mutate(day_filter = case_when(is.na(days_P) == F ~ days_P, # for some reason P is always day before pupation
                                is.na(death_days) == F ~ death_days - 1, # for some reason D is always day before death
                                TRUE ~ max(.$day))) %>% 
  filter(day < day_filter,
         is.na(value) == F) %>% 
  mutate(value = as.numeric(value)) %>% 
  filter(is.na(value) == F) %>% 
  dplyr::select(id = ID, spp = sp, temperature = `t (¡C)`, day, pit_diameter = value) %>% 
  mutate(temperature = as.character(temperature))


pit_3_clean <- pit_3_raw %>% 
  mutate(across(`Day 1`:`Day 64`, as.character)) %>% 
  pivot_longer(`Day 1`:`Day 64`) %>% 
  mutate(day = as.numeric(gsub("Day ", "", name))) %>% 
  left_join(x = ., y = pup_dead, by = "ID") %>% 
  mutate(day_filter = case_when(is.na(days_P) == F ~ days_P - 1, # for some reason P is always day before pupation
                                is.na(death_days) == F ~ death_days - 1, # for some reason D is always day before death
                                TRUE ~ max(.$day))) %>% 
  filter(day < day_filter,
         is.na(value) == F) %>% 
  mutate(value = as.numeric(value)) %>% 
  filter(is.na(value) == F) %>% 
  dplyr::select(id = ID, spp = sp, temperature = `t (¡C)`, day, pit_diameter = value) %>% 
  mutate(temperature = as.character(temperature))

pit_diameter_timeseries <- bind_rows(pit_1_clean, pit_2_clean, pit_3_clean) %>% 
  ## data selection - select individuals with at least 10 days of data
  group_by(id) %>% 
  mutate(n = n()) %>% 
  ungroup() %>% 
  filter(n >= 10) %>% 
  dplyr::select(-n)

save(pit_diameter_timeseries, file = "data/pit_diameter_timeseries.RData")








