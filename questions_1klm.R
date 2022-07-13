
# Load pkgs and data (from orig Rmd) --------------------------------------

library(sf)
library(readxl)
library(dplyr)
library(tidyr)
library(data.table)
library(stringr)
library(zoo)
library(ggplot2)
library(lme4)
library(emmeans)
library(brms)
library(tidybayes)
library(brmsmargins)
library(gt)
library(glue)

theme_set(
  theme_bw() +
    theme(strip.background = element_blank(),
          panel.grid = element_blank())
)

options(mc.cores = 4, brms.backend = 'cmdstanr', brms.file_refit = 'on_change')

w_antibiotic_conc <- fread('project/csv/water_sample_antibiotic_conc.csv')
w_ar <- fread('project/csv/water_sample_ar.csv')
w_arg_copynum <- fread('project/csv/water_sample_arg_copynum.csv')
w_bl <- fread('project/csv/water_sample_bl.csv')
w_env <- fread('project/csv/water_sample_environmental.csv')
ww_antibiotic_conc <- fread('project/csv/wastewater_sample_antibiotic_conc.csv')
ww_ar <- fread('project/csv/wastewater_sample_ar.csv')
ww_arg_copynum <- fread('project/csv/wastewater_sample_arg_copynum.csv')
ww_bl <- fread('project/csv/wastewater_sample_bl.csv')
ww_env <- fread('project/csv/wastewater_sample_ecolicounts.csv')

w_id_cols <- c('sample_no', 'season', 'site', 'area')
ww_id_cols <- c('sample_no', 'season', 'site', 'type')

sites <- st_read('project/sites.gpkg')
watersheds <- st_read('project/study_region_huc12s.gpkg')

w_env[temperature >= 32, temperature := NA]


# new data processing and analysis ----------------------------------------

# 1k. Use environmental factors to predict beta-lactamase species
# 1l. Use environmental factors to predict AR gene copy numbers
# 1m. Use environmental factors to predict antibiotic concentration in the water


# 1K. beta-lactamase species ----------------------------------------------

# Visualize distribution of outcomes overall, and by each predictor variable, with tables

w_env_bl <- w_env[w_bl, on = w_id_cols]
season_order <- c('fall', 'winter', 'spring', 'summer')
w_env_bl[, season := factor(season, levels = season_order)]

table(w_env_bl$Bl_species, useNA = "always") # Some have a - and some are empty. The empty ones should be removed.

w_env_bl <- w_env_bl[!Bl_species %in% '']
w_env_bl[Bl_species == '-', Bl_species := 'none']
w_env_bl[, Bl_species := relevel(factor(Bl_species), ref = 'none')]

# Put together the data for the present absent regression
env_cols <- c('temperature', 'pH', 'conductivity', 'turbidity', 'Ecoli_counts')
w_env_bl_presabs <- w_env_bl[, .(Bl_present = any(Bl_species != 'none')), by = c(w_id_cols, env_cols)]

# Scale the predictors as done previously
bl_presabs_scaled_predictors <- w_env_bl_presabs[, .(temperature, pH, conductivity, turbidity, Ecoli_counts)]
bl_presabs_scaled_predictors[, Ecoli_counts := log1p(Ecoli_counts)]
bl_presabs_scaled_predictors[, conductivity := log1p(conductivity)]
bl_presabs_scaled_predictors[, turbidity := log1p(turbidity)]
bl_presabs_scaled_predictors <- scale(bl_presabs_scaled_predictors)

bl_presabs_model_data <- cbind(
  w_env_bl_presabs[, .(sample_no, season, site, Bl_present)],
  bl_presabs_scaled_predictors)

# Put together the data for the categorical regression
bl_cat_scaled_predictors <- w_env_bl[, .(temperature, pH, conductivity, turbidity, Ecoli_counts)]
bl_cat_scaled_predictors[, Ecoli_counts := log1p(Ecoli_counts)]
bl_cat_scaled_predictors[, conductivity := log1p(conductivity)]
bl_cat_scaled_predictors[, turbidity := log1p(turbidity)]
bl_cat_scaled_predictors <- scale(bl_cat_scaled_predictors)

bl_cat_model_data <- cbind(
  w_env_bl[, .(sample_no, season, site, Bl_species)],
  bl_cat_scaled_predictors)
bl_cat_model_data[, sample_no := factor(sample_no)]

# One-hot encoding for the Bl_species column for the multinomial model
Bl_species_matrix <- data.table(ID = 1:nrow(w_env_bl), sp = w_env_bl$Bl_species) |>
  dcast(ID ~ sp, length)
Bl_species_matrix[, ID := NULL]
w_env_bl_model_data[, Bl_species := do.call(cbind, Bl_species_matrix)]

# It will be very difficult to predict Klebsiella, Enterobacter, and probably E.coli, because there are few of each.
# First, let's simplify to a binomial or logistic regression. (any present vs. any absent). Bernoulli is used.
# But otherwise, as done by the UGA statistician, a categorical or softmax regression is appropriate
# For the categorical, we will sometimes have more than one taxon present from the same sample number
# This can either be ignored or we can add a random effect for sample number if the data supports it.

bl_bern_fit <- brm(
  bf(Bl_present ~ temperature + pH + conductivity + turbidity + Ecoli_counts + season + (1|site)), 
  data = bl_presabs_model_data, family = bernoulli(link = 'logit'),
  prior = c(
    prior(normal(0, 3), class = b)
  ),
  chains = 4, iter = 3000, warmup = 2000, seed = 70924,
  file = 'project/bl_bern_fit'
)

bl_cat_fit <- brm(
  bf(Bl_species | 1 ~ temperature + pH + conductivity + turbidity + Ecoli_counts + season + (1|site) + (1|sample_no)), 
  data = bl_cat_model_data, family = categorical(link = 'logit'),
  prior = c(
    prior(normal(0, 3), class = b)
  ),
  chains = 4, iter = 3000, warmup = 2000, seed = 71322,
  file = 'project/bl_cat_fit'
)

# FIXME Extract some results.

# 1M. antibiotic concentration --------------------------------------------

w_env_ant_conc <- w_env[w_antibiotic_conc, on = w_id_cols]
season_order <- c('fall', 'winter', 'spring', 'summer')
w_env_ant_conc[, season := factor(season, levels = season_order)]

# Visualize distributions of the different antibiotics, ignoring predictors.

antibiotic_names <- setdiff(names(w_antibiotic_conc), w_id_cols)
w_env_ant_conc_long <- melt(w_env_ant_conc, id.vars = w_id_cols, measure.vars = antibiotic_names, variable.name = 'antibiotic', value.name = 'concentration')

# Pseudolog is used to keep zeroes in
ggplot(w_env_ant_conc_long, aes(x = concentration)) +
  geom_density() +
  facet_wrap(~ antibiotic, scales = 'free_y') +
  scale_x_continuous(trans = 'pseudo_log', breaks = c(0, 1, 10, 100, 1000))

# Look at only total. A hurdle model (gamma or lognormal) appears to be indicated.
ggplot(w_env_ant_conc, aes(x = total)) +
  geom_density() 

# Scale predictors
antconc_scaled_predictors <- w_env_ant_conc[, .(temperature, pH, conductivity, turbidity, Ecoli_counts)]
antconc_scaled_predictors[, Ecoli_counts := log1p(Ecoli_counts)]
antconc_scaled_predictors[, conductivity := log1p(conductivity)]
antconc_scaled_predictors[, turbidity := log1p(turbidity)]
antconc_scaled_predictors <- scale(antconc_scaled_predictors)

antconc_model_data <- cbind(
  w_env_ant_conc[, .(sample_no, season, site, total)],
  antconc_scaled_predictors)

ant_conc_hugamma_fit <- brm(
  bf(total ~ temperature + pH + conductivity + turbidity + Ecoli_counts + season + (1|site)), 
  data = antconc_model_data, family = hurdle_gamma(link = 'log', link_shape = 'log', link_hu = 'logit'),
  prior = c(
    prior(normal(0, 3), class = b)
  ),
  chains = 4, iter = 3000, warmup = 2000, seed = 70922,
  file = 'project/ant_conc_hugamma_fit'
)

ant_conc_hulognorm_fit <- brm(
  bf(total ~ temperature + pH + conductivity + turbidity + Ecoli_counts + season + (1|site)), 
  data = antconc_model_data, family = hurdle_lognormal(link = 'identity', link_sigma = 'log', link_hu = 'logit'),
  prior = c(
    prior(normal(0, 3), class = b)
  ),
  chains = 4, iter = 3000, warmup = 2000, seed = 70923,
  file = 'project/ant_conc_hulognorm_fit'
)
