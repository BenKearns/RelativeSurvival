# Convert IPD to aggregate data for DSMs
library("tidyverse")
library("rstan")
library("here")
library("gridExtra")
theme_set(theme_light())  # GGplot theme
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
set.seed(4361)
small_num = 1*10^-5 

g_legend = function(fig){  # Function used later to get legend from figures.
  tmp = ggplot_gtable(ggplot_build(fig))
  leg = which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend = tmp$grobs[[leg]]
  return(legend)
}

max_h = 40
my_time = 12*seq(from=0.05, to=max_h, by=0.05) # Time values for which we want estimates.
new_df = tibble(Time = my_time, AtRisk = 1)

# Read in the OS data.
OS_doc_IPD = read.table(here("Prep", "IPDdata_OS_doc.txt" ), header = TRUE) %>% arrange(time)
OS_niv_IPD = read.table(here("Prep", "IPDdata_OS_niv.txt" ), header = TRUE) %>% arrange(time)

# Get grouped estimates for plotting
list_OS = readRDS(here("Files", "list_OS.rds"))
  OS_Agg = list_OS$full_df %>% mutate(AtRisk = At_risk * Tau) %>% filter(is.finite(AtRisk))

# Aggregrate so have 1 observation per event.
fun_IPD = function(x){  # If name of event, time need changing this needs changing
  tmp_time = arrange(subset(x, event==1), time)$time
  tmp_time[length(tmp_time)] = max(x$time) # For if last observation is censored
  tmp = mutate(x, Index = findInterval(x$time, vec=tmp_time, left.open = TRUE)) %>% #Get intervals for events (events at t = 0 can be issue)
    group_by(Index) %>%  # Per interval, summary stats
    summarise(Count = n(),
              Events = sum(event),
              Censor = Count - Events)
    tmp = tmp %>% mutate(Time = tmp_time[Index+1],  # New variables
                       Alive = sum(Count) - cumsum(lag(Count, default=0)),
                       Tau = Time - lag(Time, default = 0),
                       Ln_Tau = log(Time+1) - lag(log(Time+1), default = 0),
                       AtRisk = (Alive - Censor/2) * Tau,
                       haz_t = Events / AtRisk,
                       p_t = 1 - exp(-haz_t * (lead(Time) - Time)))# %>% na.omit %>% select(-Index)
    return(tmp)
}

IPD_doc = fun_IPD(OS_doc_IPD)
IPD_niv = fun_IPD(OS_niv_IPD) 

IPD_DSM = list(Doc = IPD_doc, Niv = IPD_niv)
saveRDS(IPD_DSM, here("Files", "IPD_DSM.rds"))

mod_DSM = function(df, my_file){ # Fit DSMs
  my_data = list(
    y = df$Events,
    T = length(df$Events),
    n = df$AtRisk,
    tau = tail(df$Ln_Tau, -1)  # Or tail(df$Tau, -1)
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
    control = list(adapt_delta = 0.999, max_treedepth = 15)      # http://mc-stan.org/misc/warnings.html#maximum-treedepth-exceeded
    )
  return(fit)
}

DSM_LT_niv = mod_DSM(IPD_niv, "LT_IG.stan")
#   DSM_AT_niv = mod_DSM(IPD_niv, "AdaptedTrend_IG.stan") # Not run
 DSM_DT_niv = mod_DSM(IPD_niv, "DT_IG.stan") # 1 divergent (tau)

DSM_LT_doc = mod_DSM(IPD_doc, "LT_IG.stan") 
# traceplot(DSM_LT_doc, pars = c("Z", "beta_01", "beta_02"), inc_warmup = FALSE)
# summary(DSM_LT_doc, pars = c("Z", "beta_01", "beta_02", "level", "trend"))$summary
#  DSM_AT_doc = mod_DSM(IPD_doc, "AdaptedTrend_IG.stan") # Not run
  DSM_DT_doc = mod_DSM(IPD_doc, "DT_IG.stan") 

# Save models for later use
DSM_IPD = list(DSM_LT_niv = DSM_LT_niv,
                #DSM_AT_niv = DSM_AT_niv,
                DSM_DT_niv = DSM_DT_niv,
                DSM_LT_doc = DSM_LT_doc,
                #DSM_AT_doc = DSM_AT_doc,
                DSM_DT_doc = DSM_DT_doc)
saveRDS(DSM_IPD, here("Files", "Mods_DSM_IPD_IG_Ln.rds"))  # Or "Mods_DSM_IPD_IG_Ln.rds"

#### Below is currently for ln_tau.

# Plots of interest.
fun_DSM_plot = function(mod, df, my_str, mod_type){
# Graph of trend overtime
  fig_trend = plot(mod, pars = grep("beta2", names(mod), value = TRUE)) + theme_bw() + scale_y_continuous(trans = "reverse") +
    coord_flip() + theme(axis.text.x = element_blank()) + geom_vline(xintercept = 0, size = 1.2, colour = "blue") +
    labs(y = "Time", x = "Trend", subtitle = my_str)
# Alternative graph
  tmp_df = as_tibble(rstan::summary(mod, pars=grep("beta2", names(mod), value = TRUE),
                                    probs = c(0.025, 0.975))$summary, rownames = NA) %>%
    rownames_to_column(var = "Stat") %>% bind_cols(Time = head(df$Time, -1))
  fig_trend2 = ggplot(tmp_df, aes(x=Time/12, y = mean)) + geom_line(color="darkblue", size = 1) + 
    geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`), fill="skyblue", alpha = 0.25) +
    geom_hline(yintercept = 0, size = 1, colour = "black") + labs(y = "Trend", x="Time", subtitle = my_str)
  
# Trace plots for key parameters
  fig_trace = plot(mod, plotfun = "trace", pars = c("level", "trend", "Z")) + labs(subtitle = my_str)
  
# Plot observed vs modelled hazard
  mod_df = as.data.frame(mod)
   params = select(mod_df, grep("beta1", names(mod_df), fixed = TRUE)) %>% colMeans() %>% data.frame()
   tmp_df = data.frame(modHaz = exp(params[,1]), time = df$Time)
  fig_haz = ggplot(data = tmp_df, aes(x = time, y=modHaz)) + geom_line(color="purple") + 
     geom_line(data=filter(list_OS$group_df, Arm== my_str), aes(x=Time, y=Haz), color="black", size=1) + 
     guides(color = "none") + coord_cartesian(xlim = c(0, 25), ylim=c(0,0.25)) + labs(x="", y="", subtitle= my_str)

  my_list = list(fig_trend = fig_trend,
                 fig_trend2 = fig_trend2,
                 fig_trace = fig_trace,
                 fig_haz = fig_haz)
  return(my_list)
} 
  
# LT_Niv = fun_DSM_plot(DSM_LT_niv, IPD_niv, "Nivolumab", "Local")
#   AT_Niv = fun_DSM_plot(DSM_AT_niv, IPD_niv, "Nivolumab", "Adapted")
#   DT_Niv = fun_DSM_plot(DSM_DT_niv, IPD_niv, "Nivolumab", "Damped")
LT_Doc = fun_DSM_plot(DSM_LT_doc, IPD_doc, "Docetaxel", "Local")
  #AT_Doc = fun_DSM_plot(DSM_AT_doc, IPD_doc, "Docetaxel", "Adapted")
  DT_Doc = fun_DSM_plot(DSM_DT_doc, IPD_doc, "Docetaxel", "Damped")
  
pdf(here("Graphs", "Trend_DSM3.pdf"), width=10, height=6)
  grid.arrange(LT_Doc$fig_trend2 + labs(subtitle = "Local trend") + coord_cartesian(ylim=c(-2, 1.5)),
               DT_Doc$fig_trend2 + labs(subtitle = "Damped trend") + coord_cartesian(ylim=c(-2, 1.5)), ncol=2)
dev.off()
  #+ coord_cartesian(ylim=c(-0.7, 0.6)

# Get estimates and extrapolations - only doing for local and damped trend
fun_DSM_est = function(mod, df, mod_type){  # THIS IS FOR TAU. Inputs = model, df used to fit, model type
  beta1 = extract(mod, pars = c("beta1"))
  int_haz = map_dfr(beta1, function(x) colMeans(x))
  tmp2 = tibble(Time = df$Time, mean = int_haz$beta1)
    max_fu = max(tmp2$Time) # = tmp2$Time[length(tmp2$Time)] as ordered
  time_int = filter(new_df, Time <= max_fu) 
    int_est = approx(x=tmp2$Time, y=tmp2$mean, xout=time_int$Time, rule=2)$y
  
if (mod_type == "Local") {
  tmp = as_tibble(rstan::summary(mod, pars=c("level", "trend"), probs = 0.5)$summary, rownames = NA) %>% rownames_to_column(var = "Stat")
  tmp_df = filter(new_df, Time > max_fu) %>% select(Time) %>% mutate(Time = Time - max_fu,
           #level = filter(tmp, Stat=="level")$mean, trend = filter(tmp, Stat=="trend")$mean,  # Det mean
           level = extract(mod, pars = c("level")), trend = extract(mod, pars = c("trend")),   # Prob mean
           ext_est = pmap(list(level, trend, Time), function(level, trend, Time) level + trend * Time),
           Pred = map_dbl(ext_est, mean),
           ext_est2 = pmap(list(level, trend, Time), function(level, trend, Time) exp(level + trend * Time)),
           Pred2 = map_dbl(ext_est2, mean),
           Low = map_dbl(ext_est, function(x) quantile(x, probs = 0.025)),
           Upp = map_dbl(ext_est, function(x) quantile(x, probs = 1 - 0.025)),
           Trend_m = map_dbl(trend, mean), Level_m = map_dbl(level, mean)) %>%
            select(Time, Pred, Pred2, Low, Upp, Trend_m, Level_m)
  
 mod_est = tibble(Time = my_time, Pred = exp(c(int_est, tmp_df$Pred)), Pred2 = c(exp(int_est), tmp_df$Pred2),
                 Low = exp(c(int_est, tmp_df$Low)),
                 Upp = exp(c(int_est, tmp_df$Upp))) 
} else if (mod_type == "Damped") {
# Forecast future widths
  my_data = list(
    T = length(df$Events),
    tau = log(df$Tau)
  )
  # init_list = list(list(Z = runif(1, 0, 0.2), beta_01 = log(df$Tau)[1],
  #                list(Z = runif(1, 0, 0.2), beta_01 = log(df$Tau)[1]),
  #                list(Z = runif(1, 0, 0.2), beta_01 = log(df$Tau)[1]),
  #                list(Z = runif(1, 0, 0.2)), beta_01 = log(df$Tau)[1])
  fit = stan(
    file = here("DSM", "Widths.stan"),          # Stan program
    data = my_data,          # named list of data
    chains = 4,              # number of Markov chains
    warmup = 1000,           # number of warmup iterations per chain
    iter = 2000,             # total number of iterations per chain
    #init = init_list,        # Initial values for DSM models
    cores = 4,               # number of cores
    refresh = 0,              # show progress every 'refresh' iterations
    control = list(adapt_delta = 0.999, max_treedepth = 15)      # http://mc-stan.org/misc/warnings.html#maximum-treedepth-exceeded
    )

  tmp_df = tibble(Time = seq(1:100), level = extract(fit, pars = c("level")), trend = extract(fit, pars = c("trend")),
           ext_est = pmap(list(level, trend, Time), function(level, trend, Time) level + trend * Time),
           Pred = exp(map_dbl(ext_est, mean))) %>% select(Time, Pred) %>% filter(Pred < max_h)  # So cumsum doesn't get too big
  w_df = tibble(time_new = cumsum(c(0,tmp_df$Pred)))

  # Estimates
  tmp_df = tibble(Time = w_df$time_new) %>% mutate(level = extract(mod, pars = c("level")),
                     trend = extract(mod, pars = c("trend")), phi = extract(mod, pars = c("phi")),
   ext_est = pmap(list(level, trend, phi, Time), function(level, trend, phi, Time) level + trend * phi * (1 - phi^Time)/(1 - phi)),
           Pred = map_dbl(ext_est, mean),
           Low = map_dbl(ext_est, function(x) quantile(x, probs = 0.025)),
           Upp = map_dbl(ext_est, function(x) quantile(x, probs = 1 - 0.025)),
           Phi_m = map_dbl(phi, mean), Trend_m = map_dbl(trend, mean), Level_m = map_dbl(level, mean)) %>%
    select(Time, Pred, Low, Upp, Phi_m, Trend_m, Level_m)
    
  fu_time = filter(new_df, Time > max_fu) %>% select(Time) %>% mutate(Time = Time - max_fu)
    
 mod_est = tibble(Time = my_time,
                  Pred = exp(c(int_est, approx(x = tmp_df$Time, y = tmp_df$Pred, xout = fu_time$Time)$y)),
                  Low = exp(c(int_est, approx(x = tmp_df$Time, y = tmp_df$Low, xout = fu_time$Time)$y)),
                  Upp = exp(c(int_est, approx(x = tmp_df$Time, y = tmp_df$Upp, xout = fu_time$Time)$y))) 
} else stop("Wrong model")  
  return(mod_est)
}

fun_DSM_est_Ln = function(mod, df, mod_type){  # Ln_Tau. Inputs = model, df used to fit, model type
  beta1 = extract(mod, pars = c("beta1"))
  int_haz = map_dfr(beta1, function(x) colMeans(x))
  tmp2 = tibble(Time = df$Time, mean = int_haz$beta1)
    max_fu = max(tmp2$Time) # = tmp2$Time[length(tmp2$Time)] as ordered
  time_int = filter(new_df, Time <= max_fu) 
    int_est = approx(x=tmp2$Time, y=tmp2$mean, xout=time_int$Time, rule=2)$y
  
if (mod_type == "Local") {
  tmp_df = filter(new_df, Time > max_fu) %>% select(Time) %>% mutate(Time = log1p(Time) - log1p(max_fu),
           level = extract(mod, pars = c("level")), trend = extract(mod, pars = c("trend")),   # Prob mean
           ext_est = pmap(list(level, trend, Time), function(level, trend, Time) level + trend * Time),
           Pred = map_dbl(ext_est, mean),
           ext_est2 = pmap(list(level, trend, Time), function(level, trend, Time) exp(level + trend * Time)),
           Pred2 = map_dbl(ext_est2, mean),
           Low = map_dbl(ext_est, function(x) quantile(x, probs = 0.025)),
           Upp = map_dbl(ext_est, function(x) quantile(x, probs = 1 - 0.025)),
           Trend_m = map_dbl(trend, mean), Level_m = map_dbl(level, mean)) %>%
            select(Time, Pred, Pred2, Low, Upp, Trend_m, Level_m)
  
 mod_est = tibble(Time = my_time, Pred = exp(c(int_est, tmp_df$Pred)),
                 Low = exp(c(int_est, tmp_df$Low)),
                 Upp = exp(c(int_est, tmp_df$Upp))) 
} else if (mod_type == "Damped") {
# Forecast future widths
  my_data = list(
    T = length(df$Events),
    tau = log(df$Ln_Tau)
  )
  # init_list = list(list(Z = runif(1, 0, 0.2), beta_01 = log(df$Tau)[1],
  #                list(Z = runif(1, 0, 0.2), beta_01 = log(df$Tau)[1]),
  #                list(Z = runif(1, 0, 0.2), beta_01 = log(df$Tau)[1]),
  #                list(Z = runif(1, 0, 0.2)), beta_01 = log(df$Tau)[1])
  fit = stan(
    file = here("DSM", "Widths.stan"),          # Stan program
    data = my_data,          # named list of data
    chains = 4,              # number of Markov chains
    warmup = 1000,           # number of warmup iterations per chain
    iter = 2000,             # total number of iterations per chain
    #init = init_list,        # Initial values for DSM models
    cores = 4,               # number of cores
    refresh = 0,              # show progress every 'refresh' iterations
    control = list(adapt_delta = 0.999, max_treedepth = 15)      # http://mc-stan.org/misc/warnings.html#maximum-treedepth-exceeded
    )

  tmp_df = tibble(Time = seq(1:100), level = extract(fit, pars = c("level")), trend = extract(fit, pars = c("trend")),
           ext_est = pmap(list(level, trend, Time), function(level, trend, Time) level + trend * Time),
           Pred = exp(map_dbl(ext_est, mean))) %>% select(Time, Pred) %>% filter(Pred < max_h)  # So cumsum doesn't get too big
  w_df = tibble(time_new = cumsum(c(max_fu, tmp_df$Pred)))

  # Estimates
  tmp_df = tibble(Time = log1p(w_df$time_new) - log1p(max_fu)) %>% mutate(level = extract(mod, pars = c("level")),
                     trend = extract(mod, pars = c("trend")), phi = extract(mod, pars = c("phi")),
   ext_est = pmap(list(level, trend, phi, Time), function(level, trend, phi, Time) level + trend * phi * (1 - phi^Time)/(1 - phi)),
           Pred = map_dbl(ext_est, mean), Time = w_df$time_new,  # Note over-writing time here
           Low = map_dbl(ext_est, function(x) quantile(x, probs = 0.025)),
           Upp = map_dbl(ext_est, function(x) quantile(x, probs = 1 - 0.025)),
           Phi_m = map_dbl(phi, mean), Trend_m = map_dbl(trend, mean), Level_m = map_dbl(level, mean)) %>%
    select(Time, Pred, Low, Upp, Phi_m, Trend_m, Level_m)
    
  fu_time = filter(new_df, Time > max_fu) %>% select(Time) %>% mutate(Time = Time)
    
 mod_est = tibble(Time = my_time,
                  Pred = exp(c(int_est, approx(x = tmp_df$Time, y = tmp_df$Pred, xout = fu_time$Time)$y)),
                  Low = exp(c(int_est, approx(x = tmp_df$Time, y = tmp_df$Low, xout = fu_time$Time)$y)),
                  Upp = exp(c(int_est, approx(x = tmp_df$Time, y = tmp_df$Upp, xout = fu_time$Time)$y))) 
} else stop("Wrong model")  
  return(mod_est)
}
  
LT_Niv_est = fun_DSM_est_Ln(DSM_LT_niv, IPD_niv, "Local")
  DT_Niv_est = fun_DSM_est_Ln(DSM_DT_niv, IPD_niv, "Damped")  # 1 divergent
LT_Doc_est = fun_DSM_est_Ln(DSM_LT_doc, IPD_doc, "Local")
  DT_Doc_est = fun_DSM_est_Ln(DSM_DT_doc, IPD_doc, "Damped")  # 5 divergent






