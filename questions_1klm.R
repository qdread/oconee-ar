
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

scale_predictors <- function(DT) {
  DTscale <- DT[, .(temperature, pH, conductivity, turbidity, Ecoli_counts)]
  DTscale[, Ecoli_counts := log1p(Ecoli_counts)]
  DTscale[, conductivity := log1p(conductivity)]
  DTscale[, turbidity := log1p(turbidity)]
  scale(DTscale)
}

##### new data processing and analysis 

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

# Plot parameter estimates from Bernoulli fit. We don't have separate ones by taxon so this is simpler.
bl_bern_slopes <- bl_bern_fit %>%
  gather_draws(`b_.*`, regex = TRUE) %>%
  filter(!.variable %in% 'b_Intercept') %>%
  mutate(.variable = gsub('b_', '', .variable))

# Weak evidence for a positive pH trend, where there is more likelihood to be a beta-lac bacterium where it is more basic. 
# Also varies by season 
ggplot(bl_bern_slopes, aes(y = .variable, x = .value)) +
  stat_interval(.width = c(0.66, 0.9, 0.95)) +
  stat_pointinterval(geom = "point", size = 2) +
  scale_color_brewer(palette = 'Blues', name = 'credible interval') +
  geom_vline(xintercept = 0, linetype = 'dotted', color = 'gray50') +
  labs(x = 'parameter estimate') +
  theme(axis.title.y = element_blank())

# Plot marginal effects for each of the continuous predictors
scaled_centers <- attr(bl_presabs_scaled_predictors, 'scaled:center')
scaled_scales <- attr(bl_presabs_scaled_predictors, 'scaled:scale')

x_break_list <- list(
  temperature = c(5, 10, 15, 20, 25),
  pH = c(6.5, 7, 7.5),
  turbidity = c(0, 3, 10, 30),
  conductivity = c(50, 100, 500),
  Ecoli_counts = c(3, 10, 30, 100, 300, 3000)
)

variables <- c('temperature', 'pH', 'turbidity', 'conductivity', 'Ecoli_counts')
blue_cols <- RColorBrewer::brewer.pal(3, "Blues")[1:2]

plot_marginal_effect <- function(variable, fit, y_label, n_x = 50, y_scale_trans = 'identity') {
  x_range <- range(bl_presabs_model_data[[variable]], na.rm = TRUE)
  x_r_seq <- seq(x_range[1], x_range[2], length.out = n_x)
  x_r_seq_backtransformed <- scaled_scales[variable] * x_r_seq + scaled_centers[variable]
  if (variable %in% c('turbidity', 'conductivity', 'Ecoli_counts')) {
    x_r_seq_backtransformed <- expm1(x_r_seq_backtransformed)
    x_scale_trans <- 'log1p'
  } else {
    x_scale_trans <- 'identity'
  }
  
  pred_dat <- setNames(data.frame(x = x_r_seq), variable)
  pred_dat_backtransformed <- setNames(data.frame(x = x_r_seq_backtransformed), variable)
  
  marginal_draws <- brmsmargins(fit, at = pred_dat)
  marginal_summ <- Rutilitybelt::pred_quantile(x_pred = pred_dat_backtransformed, y_pred = marginal_draws$Posterior)
  
  ggplot(marginal_summ, aes(x = !!ensym(variable), y = q0.5)) + 
    geom_ribbon(aes(ymin = q0.025, ymax = q0.975), color = NA, fill = blue_cols[1]) + 
    geom_ribbon(aes(ymin = q0.17, ymax = q0.83), color = NA, fill = blue_cols[2]) + 
    geom_line(size = 0.8) + 
    theme(legend.position = 'none') +
    scale_fill_brewer(palette = 'Dark2') + scale_color_brewer(palette = 'Dark2') +
    scale_x_continuous(name = variable, trans = x_scale_trans, breaks = x_break_list[[variable]]) +
    scale_y_continuous(name = y_label, trans = y_scale_trans)
}

bl_bern_marginal_plots <- lapply(variables, plot_marginal_effect, fit = bl_bern_fit, y_label = 'modeled probability of any beta-lactamase presence')

# Plot marginal means for each season
bl_bern_season_means <- emmeans(bl_bern_fit, ~ season) %>%
  gather_emmeans_draws() %>%
  mutate(.value = plogis(.value))

ggplot(bl_bern_season_means, aes(x = season, y = .value)) +
  stat_interval(.width = c(0.66, 0.95)) +
  stat_pointinterval(geom = "point", size = 2) +
  scale_color_brewer(palette = 'Blues', name = 'credible interval') +
  labs(y = 'modeled probability of any beta-lactamase presence') +
  theme(legend.position = c(0.1, 0.8))

#### Results from categorical fit

# Plot parameter estimates from softmax or categorical fit. We have a slope for each combination of taxon x predictor.
bl_cat_slopes <- bl_cat_fit %>%
  gather_draws(`b_.*`, regex = TRUE) %>%
  mutate(.variable = gsub('b_mu', '', .variable)) %>%
  filter(!grepl('Intercept', .variable)) %>%
  separate(.variable, into = c('taxon', '.variable'), sep = '_', extra = 'merge')

# Looks like the Klebsiella species are more likely to be found in basic water samples, same with Serratia to some extent.
# The rarer ones do not have much of a pattern.
# Serratia is less common in the fall relative to the other three seasons.
ggplot(bl_cat_slopes, aes(y = .variable, x = .value)) +
  stat_interval(.width = c(0.66, 0.9, 0.95)) +
  stat_summary(fun = median, geom = 'point', size = 2) +
  facet_wrap(~ taxon) +
  scale_color_brewer(name = 'credible interval', palette = 'Blues') +
  geom_vline(xintercept = 0, linetype = 'dotted', color = 'gray50') +
  labs(x = 'parameter estimate') +
  theme(axis.title.y = element_blank(), legend.position = c(0.8, 0.2))

# Plot marginal effects for each of the continuous predictors from softmax fit
# brmsmargins doesn't currently support models with multivariate outcome. emmeans is not supported either.
# FIXME later I will try to add something here.

# 1L. antibiotic resistance gene copy numbers -----------------------------

w_env_arg_cn <- w_env[w_arg_copynum, on = w_id_cols]
season_order <- c('fall', 'winter', 'spring', 'summer')
w_env_arg_cn[, season := factor(season, levels = season_order)]

# Visualize distributions of the different ARGs, ignoring predictors.
arg_abs_names <- grep('^abs', names(w_arg_copynum), value = TRUE)
arg_rlt_names <- grep('^rlt', names(w_arg_copynum), value = TRUE)
w_env_arg_long <- melt(w_env_arg_cn, id.vars = w_id_cols, measure.vars = c(arg_abs_names, arg_rlt_names), variable.name = 'ARG', value.name = 'copy_number') %>%
  separate(ARG, into = c('type', 'ARG'), sep = '_') %>%
  dcast(sample_no + season + site + area + ARG ~ type)

ggplot(w_env_arg_long, aes(x = abs)) +
  geom_density() +
  facet_wrap(~ ARG, scales = 'free') +
  scale_x_continuous(trans = 'pseudo_log', breaks = c(0, 1, 10, 100, 1000))

ggplot(w_env_arg_long, aes(x = rlt)) +
  geom_density() +
  facet_wrap(~ ARG, scales = 'free') +
  scale_x_continuous(trans = 'pseudo_log')

# As the UGA statistician noted, abs and rlt have the same distribution.
# A hurdle model (gamma or lognormal) is indicated.
ggplot(w_env_arg_long[ARG == 'total'], aes(x = abs)) +
  geom_density() 

ggplot(w_env_arg_long[ARG == 'total'], aes(x = rlt)) +
  geom_density() 

# Scale predictors
arg_scaled_predictors <- scale_predictors(w_env_arg_cn)

arg_model_data <- cbind(
  w_env_arg_cn[, .(sample_no, season, site, abs_total)],
  arg_scaled_predictors)

arg_cn_hugamma_fit <- brm(
  bf(abs_total ~ temperature + pH + conductivity + turbidity + Ecoli_counts + season + (1|site)), 
  data = arg_model_data, family = hurdle_gamma(link = 'log', link_shape = 'log', link_hu = 'logit'),
  prior = c(
    prior(normal(0, 3), class = b)
  ),
  chains = 4, iter = 3000, warmup = 2000, seed = 71422,
  file = 'project/arg_cn_hugamma_fit'
)

arg_cn_hulognorm_fit <- brm(
  bf(abs_total ~ temperature + pH + conductivity + turbidity + Ecoli_counts + season + (1|site)), 
  data = arg_model_data, family = hurdle_lognormal(link = 'identity', link_sigma = 'log', link_hu = 'logit'),
  prior = c(
    prior(normal(0, 3), class = b)
  ),
  chains = 4, iter = 3000, warmup = 2000, seed = 71423,
  file = 'project/arg_cn_hulognorm_fit'
)

# Compare the models
arg_cn_hugamma_fit <- add_criterion(arg_cn_hugamma_fit, 'loo')
arg_cn_hulognorm_fit <- add_criterion(arg_cn_hulognorm_fit, 'loo')

loo_compare(arg_cn_hugamma_fit, arg_cn_hulognorm_fit) # Hurdle gamma is somewhat better.

pp_check(arg_cn_hugamma_fit) + scale_x_log10() # This looks better.
pp_check(arg_cn_hulognorm_fit) + scale_x_log10()

# Plot parameter estimates from hurdle gamma fit. We don't have separate ones by taxon so this is simpler.
arg_cn_slopes <- arg_cn_hugamma_fit %>%
  gather_draws(`b_.*`, regex = TRUE) %>%
  filter(!.variable %in% 'b_Intercept') %>%
  mutate(.variable = gsub('b_', '', .variable))

# We see decent evidence for a negative pH trend (more ARG copy numbers in acidic water)
# and a positive Ecoli trend (more ARG copy numbers in samples with more Ecoli)
# conductivity may also show a positive trend.
ggplot(arg_cn_slopes, aes(y = .variable, x = .value)) +
  stat_pointinterval() +
  geom_vline(xintercept = 0, linetype = 'dotted', color = 'gray50') +
  labs(x = 'parameter estimate') +
  theme(axis.title.y = element_blank())

# Plot marginal effects for each of the continuous predictors
scaled_centers <- attr(arg_scaled_predictors, 'scaled:center')
scaled_scales <- attr(arg_scaled_predictors, 'scaled:scale')

arg_cn_marginal_plots <- lapply(variables, plot_marginal_effect, fit = arg_cn_hugamma_fit, y_label = 'modeled total AR gene copy numbers', y_scale_trans = 'log10')

# Plot marginal means for each season
arg_cn_season_means <- emmeans(arg_cn_hugamma_fit, ~ season) %>%
  gather_emmeans_draws() %>%
  mutate(.value = exp(.value))

ggplot(arg_cn_season_means, aes(x = season, y = .value)) +
  stat_interval(.width = c(0.66, 0.95)) +
  stat_pointinterval(geom = "point", size = 2) +
  scale_color_brewer(palette = 'Blues', name = 'credible interval') +
  scale_y_log10(name = 'modeled total AR gene copy numbers') +
  theme(legend.position = c(0.1, 0.8)) 


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
antconc_scaled_predictors <- scale_predictors(w_env_ant_conc)

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

# Compare the models
ant_conc_hugamma_fit <- add_criterion(ant_conc_hugamma_fit, 'loo')
ant_conc_hulognorm_fit <- add_criterion(ant_conc_hulognorm_fit, 'loo')

loo_compare(ant_conc_hugamma_fit, ant_conc_hulognorm_fit) # Hurdle gamma is somewhat better.

pp_check(ant_conc_hugamma_fit) + scale_x_log10() # This looks a tiny bit better.
pp_check(ant_conc_hulognorm_fit) + scale_x_log10()

# Plot parameter estimates from hurdle gamma fit. We don't have separate ones by taxon so this is simpler.
ant_conc_slopes <- ant_conc_hugamma_fit %>%
  gather_draws(`b_.*`, regex = TRUE) %>%
  filter(!.variable %in% 'b_Intercept') %>%
  mutate(.variable = gsub('b_', '', .variable))

# Essentially no trend except for differences between seasons.
ggplot(ant_conc_slopes, aes(y = .variable, x = .value)) +
  stat_pointinterval() +
  geom_vline(xintercept = 0, linetype = 'dotted', color = 'gray50') +
  labs(x = 'parameter estimate') +
  theme(axis.title.y = element_blank())

# Plot marginal effects for each of the continuous predictors
scaled_centers <- attr(antconc_scaled_predictors, 'scaled:center')
scaled_scales <- attr(antconc_scaled_predictors, 'scaled:scale')

ant_conc_marginal_plots <- lapply(variables, plot_marginal_effect, fit = ant_conc_hugamma_fit, y_label = 'modeled total antibiotic concentration', y_scale_trans = 'log10')

# Plot marginal means for each season
ant_conc_season_means <- emmeans(ant_conc_hugamma_fit, ~ season) %>%
  gather_emmeans_draws() %>%
  mutate(.value = exp(.value))

ggplot(ant_conc_season_means, aes(x = season, y = .value)) +
  stat_interval(.width = c(0.66, 0.95)) +
  stat_pointinterval(geom = "point", size = 2) +
  scale_color_brewer(palette = 'Blues', name = 'credible interval') +
  scale_y_log10(name = 'modeled total antibiotic concentration') +
  theme(legend.position = c(0.1, 0.8))
