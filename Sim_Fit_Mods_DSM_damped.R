library("tidyverse")
library("here")
library("rstan")
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

df_sims = read_rds(here("df_sims.rds"))
my_times = c(seq(from=0, to = 5, by = 0.5), seq(from=6, to = 30, by = 1))
df2b = read_csv(here("Surv_tru.csv"))

###########################################
# DRSM
haz_rate = function(x, t, out = "prob"){ # Function to convert between rates and probability
  tmp  = t - lag(t, default = 0)
  if (out == "rate"){
    y = - (log(1 - x)) / tmp
  } else if (out == "prob") {
    y = 1 - exp(- x * tmp)
  } else {
    "error!"
  }
  return (y)
}

# First need some extra data: convert 'my_times' to a data frame, and get general population hazards
new_df = tibble(Time = my_times)
Eng_HMD = read.table(here("HMD", "UK_HMDv2.txt"), header=TRUE)
Eng_2016 = filter(Eng_HMD, Year == 2016) %>% mutate(Age2 = fct_recode(Age, `110` = "110+"),
                                                    Years = as.numeric(levels(Age2))[Age2] - 63,
                                                    Hazard = haz_rate(qx, Years + 1, "rate")) %>%
  filter(Years >= 0) # Median age = 63, treat as mean

Gen_pop = as_tibble(approx(x = Eng_2016$Years, y = Eng_2016$Hazard, xout = new_df$Time, rule = 2))
names(Gen_pop) = c("Time", "Hazard")

# Function to fit models
mod_DRSM = function(df, my_file){ # Fit DRSMs
  my_data = list(
    y = tail(df$nevents, -1),
    T = length(tail(df$nevents, -1)),
    n = tail(df$AtRisk, -1),
    tau = tail(df$Ln_Tau, -2),
    Pop = tail(df$Pop_Haz, -1)  # General population hazard
  )
  init_list = list(list(Z = runif(1, 0, 0.2)),
                   list(Z = runif(1, 0, 0.2)),
                   list(Z = runif(1, 0, 0.2)),
                   list(Z = runif(1, 0, 0.2)))
  fit = stan(
    file = here("DSM", my_file),          # Stan program
    data = my_data,          # named list of data
    chains = 4,              # number of Markov chains
    warmup = 1000,           # number of warmup iterations per chain
    iter = 2000,             # total number of iterations per chain
    cores = 4,               # number of cores
    init = init_list,        # Initial values for DSM models
    refresh = 0,              # show progress every 'refresh' iterations
    seed = 4827,
    verbose = FALSE,
    control = list(adapt_delta = 0.999, max_treedepth = 15)      # http://mc-stan.org/misc/warnings.html#maximum-treedepth-exceeded
    )
  return(fit)
}

# Function to get estimates from model
fun_DSM_est = function(mod, df, mod_name, g_pop, n_boot){
  mod_df = as.data.frame(mod)
  betas = select(mod_df, grep("beta1", names(mod_df), fixed = TRUE))
  # Within-sample estimates and CI - note these are for observed time
  tmp = tibble(time = df$time, est=apply(betas, 2, mean),
            lcl = apply(betas, 2, quantile, 0.025),
            ucl = apply(betas, 2, quantile, 1-0.025), Model = mod_name)
  # Now get estimates for required times (assume linear changes)
  max_fu = max(tmp$time)
  time_in = filter(new_df, Time <= max_fu)
  df_in = tibble(time = time_in$Time, est=exp(approx(x=tmp$time, y=tmp$est, xout=time_in$Time, rule=1)$y),
            lcl = exp(approx(x=tmp$time, y=tmp$lcl, xout=time_in$Time, rule=2)$y),
            ucl = exp(approx(x=tmp$time, y=tmp$ucl, xout=time_in$Time, rule=2)$y), Model = mod_name)
  # PSA for within-sample
  params = sample_n(betas, size = n_boot, replace = TRUE)
  params2 = split(params, seq(nrow(params))) # Turn into list, to use with map.
  PSA_in = map(params2, function(par) approx(x=tmp$time, y=c(par), xout=time_in$Time, rule=1)$y)

  # Now for extrapolations - for now use bootstrapping
  time_ext = filter(new_df, Time > max_fu)

  tmp = tibble(Level = sample(mod_df$level, size = n_boot, replace = TRUE),  # PSA samples
           Trend = sample(mod_df$trend, size = n_boot, replace = TRUE), Phi = sample(mod_df$phi, size = n_boot, replace = TRUE))
      #### Full posterior 
    full_level = select(mod_df, grep("level", names(mod_df), fixed = TRUE))  # This and below - 2 diff ways for same output
    full_trend = as_tibble(extract(mod, pars = c("trend")))
    full_phi = as_tibble(extract(mod, pars = c("phi")))
    tmp_full = tibble(Level = full_level$level, Trend = full_trend$trend, Phi= full_phi$phi) 

    fu_time = filter(new_df, Time > max_fu) %>% select(Time) %>% mutate(Time = Time)   # Times we want extraps for
    PSA_out = bind_cols(Time = fu_time$Time, 
                        pmap_dfc(tmp, function(Level, Trend, Phi) tibble(val = Level + Trend * Phi * (1 - Phi^(log1p(fu_time$Time) - log1p(max_fu)))/(1 - Phi)))) %>%
      gather(-Time, key = "ID", value = "Est") %>% nest(Vals = c(Time, Est)) %>%
      mutate(yout = map(Vals, function(df) approx(x = df$Time, y = df$Est, xout = fu_time$Time)$y)) %>% select(yout) 
    PSA_out = PSA_out$yout[[1]]
    
    tmp_df = tibble(Time = log1p(fu_time$Time) - log1p(max_fu)) %>% mutate(level = extract(mod, pars = c("level")),
                                                                           trend = extract(mod, pars = c("trend")), phi = extract(mod, pars = c("phi")),
                                                                           ext_est = pmap(list(level, trend, phi, Time), function(level, trend, phi, Time) level + trend * phi * (1 - phi^Time)/(1 - phi)),
                                                                           Pred = map_dbl(ext_est, mean), Time = fu_time$Time,  # Note over-writing time here
                                                                           Low = map_dbl(ext_est, function(x) quantile(x, probs = 0.025)),
                                                                           Upp = map_dbl(ext_est, function(x) quantile(x, probs = 1 - 0.025))) %>% select(Time, Pred, Low, Upp)
    df_out = tibble(time = time_ext$Time, est = exp(approx(x = tmp_df$Time, y = tmp_df$Pred, xout = fu_time$Time)$y),
                    lcl = exp(approx(x = tmp_df$Time, y = tmp_df$Low, xout = fu_time$Time)$y),
                    ucl = exp(approx(x = tmp_df$Time, y = tmp_df$Upp, xout = fu_time$Time)$y), Model = mod_name)
    
  df_full = bind_rows(df_in, df_out) %>% mutate(est = est + g_pop$Hazard)
  PSA_df = map2(PSA_in, PSA_out, function(x, y) exp(c(x, y)) + g_pop$Hazard)
  #return(list(df_full = df_full, PSA_df = PSA_df))
  return(df_full)
}

# Fit model in global environment once to avoid recompiling
my_data <- list(
  y = rep(1,10), T = 10, n = seq(10,1,-1), tau = rep(1,9), Pop = rep(0.01,10)
)

temp_mod <- stan(
  file = here("DSM", "DampedTrend_GenPop.stan"), # Stan program
  data = my_data,         # named list of data
  chains = 1
)
set.seed(675)
df_DRSM = mutate(df_sims, mod = map(df_Agg, function(x) mod_DRSM(x, "DampedTrend_GenPop.stan")),
                 mod_haz = map2(df_Agg, mod, function(x, y) fun_DSM_est(y, tail(x, -1), "Damped", Gen_pop, 1)), Model = "Damped")

write_rds(df_DRSM, here("df_Damped.rds"))

# df_fig = select(df_DRSM, ID, mod_haz) %>% unnest(cols = c(mod_haz)) %>% replace_na(list(est=0))
# ggplot() + geom_line(data = df_fig, aes(x = time, y = est, group = ID), alpha = 0.8, colour = "purple") + 
#   geom_line(data = df2b, aes(x = time, y = haz_tru), colour = "black", size = 1) + theme_light() +
#   labs(x = "Years", y = "Hazard") + coord_cartesian(ylim=c(0,1))

