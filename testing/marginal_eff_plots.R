effect <- 'temperature:taxon'
scaled_centers <- attr(scaled_predictors, 'scaled:center')
scaled_scales <- attr(scaled_predictors, 'scaled:scale')

# Custom function to plot the conditional effects by season and also overall, with two kinds of credible interval

CE_season_95 <- conditional_effects(mv_ar_fit, effects = effect, conditions = data.frame(season = c('fall','winter','spring','summer')), prob = 0.95, plot = FALSE)

CE_overall_95 <- conditional_effects(mv_ar_fit, effects = effect, prob = 0.95, plot = FALSE)

CE_season_66 <- conditional_effects(mv_ar_fit, effects = effect, conditions = data.frame(season = c('fall','winter','spring','summer')), prob = 0.66, plot = FALSE)

CE_overall_66 <- conditional_effects(mv_ar_fit, effects = effect, prob = 0.66, plot = FALSE)


#### Try this with the brmsmargins package
# This should average across all the other predictors.
library(brmsmargins)
n_x <- 10
temp_r <- range(w_env_ar_long$temperature, na.rm = TRUE)
temp_r_seq <- seq(temp_r[1], temp_r[2], length.out = n_x)
temp_r_seq_backtransformed <- scaled_scales['temperature'] * temp_r_seq + scaled_centers['temperature']
pred_dat <- expand.grid(temperature = temp_r_seq, taxon = unique(w_env_ar_long$taxon))
pred_dat_backtransformed <- expand.grid(temperature = temp_r_seq_backtransformed, taxon = unique(w_env_ar_long$taxon))

temp_marginal_draws <- brmsmargins(mv_ar_fit, at = pred_dat)

temp_marginal_summ <- Rutilitybelt::pred_quantile(x_pred = pred_dat_backtransformed, y_pred = temp_marginal_draws$Posterior)

ggplot(temp_marginal_summ, aes(x=temperature, y=q0.5, color = taxon, fill = taxon)) + geom_ribbon(aes(ymin = q0.025, ymax = q0.975), alpha = 0.3) + geom_line() + facet_wrap(~ taxon)
