library("tidyverse")
library("flexsurv")
library("here")
library("gridExtra")
#library("githubinstall")
#gh_install_packages("cuRe")  
library("cuRe")   # For GenFlexCureModel
library("rstpm2")  # For stpm2

df_sims = read_rds(here("df_sims.rds"))
my_times = c(seq(from=0, to = 5, by = 0.5), seq(from=6, to = 30, by = 1))
df2b = read_csv(here("Surv_tru.csv"))

# https://github.com/LasseHjort/LossOfLifetimeEstimation/blob/master/R/Analysis.R
# Function to get knot locations
get.knots <- function(data, k){
  quantile(data$Obs_Surv[data$Event ==1], seq(0, 1, length.out = k))
}

# Fit a standard log-logistic
set.seed(675)
df_LL = mutate(df_sims, mod = map(df_IPD, function(x) flexsurvreg(Surv(Obs_Surv, Event) ~ 1, data = x, bhazard = x$Pop_Haz, dist = "llogis")),
               mod_surv = map(mod, function(x) summary(x, t=my_times, type="survival", ci = FALSE, tidy = TRUE)),
               mod_haz = map(mod, function(x) summary(x, t=my_times, type="hazard", ci = FALSE, tidy = TRUE)), Model = "LL")

# Fit flexible models: https://doi.org/10.1186/s12874-019-0661-8
set.seed(675)
df_FCM = mutate(df_sims, knots = map(df_IPD, function(x) log(get.knots(x, 6))),
                mod = map2(knots, df_IPD, function(k, df) GenFlexCureModel(Surv(Obs_Surv, Event) ~ -1, data = df, bhazard = df$Pop_Haz,
                                                              smooth.formula = ~cb(x = log(Obs_Surv), knots = k), ini.types = "cure")),
                mod_surv = map(mod, function(x) predict(x, type = "surv", time = my_times, var.type = "n")),
                mod_haz = map(mod, function(x) predict(x, type = "hazard", time = my_times, var.type = "n")), Model = "FCM")

set.seed(675)
df_NRS = mutate(df_sims, knots = map(df_IPD, function(x) log(get.knots(x, 6))),
                mod = map2(knots, df_IPD, function(k, df) stpm2(Surv(Obs_Surv, Event) ~ -1, data = df, bhazard = df$Pop_Haz,
                                                                smooth.formula = ~cb(x = log(Obs_Surv), knots = k))),
                mod_surv = map(mod, function(x) predict(x, newdata = tibble(Obs_Surv = my_times), type = "surv")),
                mod_haz = map(mod, function(x) predict(x, newdata = tibble(Obs_Surv = my_times), type = "hazard")), Model = "NRS")

# Save estimates
df_mods = rbind(df_LL, select(df_FCM,-knots), select(df_NRS,-knots))
write_rds(df_mods, here("df_mods.rds"))
