# Get transition probabilities to use in CEA model
library("tidyverse")
library("flexsurv")
library("mgcv")
library("here")
library("rstan")
library("gridExtra")
library("numDeriv")
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
set.seed(4361)

theme_set(theme_light())
max_h = 20*13
my_time = seq(from=(1/4), to=max_h, by=(1/4)) # Time values for which we want estimates.
new_df = tibble(Time = my_time, AtRisk = 1)
num_PSA = 2000

g_legend = function(fig){  # Function used later to get legend from figures.
  tmp = ggplot_gtable(ggplot_build(fig))
  leg = which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend = tmp$grobs[[leg]]
  return(legend)
}

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

# Load data - concentrating on OS estimates and extrapolations (both with 95% CI)
IPD_DSM = readRDS(here("Files", "IPD_DSM.rds"))
    df_Doc = IPD_DSM$Doc
    df_Niv = IPD_DSM$Niv
OS_IPD_Niv = read.table(here("Prep", "IPDdata_OS_niv.txt" ), header = TRUE)
OS_IPD_Doc = read.table(here("Prep", "IPDdata_OS_doc.txt" ), header = TRUE)

###########################################
# Log-logistic
tmp = flexsurvreg(Surv(time, event) ~ 1, data = OS_IPD_Doc, dist = "llogis")
  ll_OS_Doc = summary(tmp, t=my_time, type="hazard", tidy=TRUE) %>% bind_cols(Model = rep("Log-logistic", length(my_time)))
  # PSA bit
  params = normboot.flexsurvreg(tmp, num_PSA)
  #PSA_ll_OS_Doc = data.table::transpose(map2_dfc(params[,1], params[,2], function(x, y) var = hllogis(my_time, x, y)))
  PSA_ll_OS_Doc = map2(params[,1], params[,2], function(x, y) var = hllogis(my_time, x, y))
  
tmp = flexsurvreg(Surv(time, event) ~ 1, data = OS_IPD_Niv, dist = "llogis")
  ll_OS_Niv = summary(tmp, t=my_time, type="hazard", tidy=TRUE) %>% bind_cols(Model = rep("Log-logistic", length(my_time)))
  # PSA bit
  params = normboot.flexsurvreg(tmp, num_PSA)
  PSA_ll_OS_Niv =  map2(params[,1], params[,2], function(x, y) hllogis(my_time, x, y))

###########################################
# RP, 2 knots, odds
tmp = flexsurvspline(Surv(time, event) ~ 1, k = 2, data = OS_IPD_Doc, scale = "odds")
  RP_OS_Doc = summary(tmp, t=my_time, type="hazard", tidy=TRUE) %>% bind_cols(Model = rep("Royston-Parmar", length(my_time)))
  params = normboot.flexsurvreg(tmp, num_PSA)
  params2 = split(params, seq(nrow(params))) # Turn into list, to use with map
  PSA_RP_OS_Doc = map(params2, function(x) hsurvspline(my_time, gamma = x, knots = tmp$knots, scale = "odds", timescale="log"))

tmp = flexsurvspline(Surv(time, event) ~ 1, k = 2, data = OS_IPD_Niv, scale = "odds")
  RP_OS_Niv = summary(tmp, t=my_time, type="hazard", tidy=TRUE) %>% bind_cols(Model = rep("Royston-Parmar", length(my_time)))
  params = normboot.flexsurvreg(tmp, num_PSA)
  params2 = split(params, seq(nrow(params))) # Turn into list, to use with map
  PSA_RP_OS_Niv = map(params2, function(x) hsurvspline(my_time, gamma = x, knots = tmp$knots, scale = "odds", timescale="log"))

###########################################
# GAM, logit link, log-time
df_Doc = df_Doc %>% mutate(NonEvent = pmax(AtRisk - Events, 0))
tmp = gam(cbind(Events, NonEvent) ~ s(log(Time)), data=df_Doc, family=binomial(link="logit"))
  my_res = predict(object=tmp, newdata=new_df, type="link", se.fit=TRUE)
  ilink = family(tmp)$linkinv  # https://www.fromthebottomoftheheap.net/2018/12/10/confidence-intervals-for-glms/
  GAM_OS_Doc = tibble(time = my_time, est=ilink(my_res$`fit`),
            lcl = ilink(my_res$`fit` + qnorm(0.025) * my_res$se.fit),
            ucl = ilink(my_res$`fit` + qnorm(1 - 0.025) * my_res$se.fit), Model = "GAM")
  # PSA - see p342-343 of 2nd edition of Wood GAM book.
  Xp = predict(object=tmp, newdata=new_df, type="lpmatrix")
  params = rmvn(n = num_PSA, coef(tmp), vcov(tmp))
  params2 = split(params, seq(nrow(params)))
  PSA_GAM_OS_Doc = map(params2, function(x) ilink(Xp %*% x))
  
df_Niv = df_Niv %>% mutate(NonEvent = pmax(AtRisk - Events, 0))
tmp = gam(cbind(Events, NonEvent) ~ s(log(Time)), data=df_Niv, family=binomial(link="logit"))
  my_res = predict(object=tmp, newdata=new_df, type="link", se.fit=TRUE)
  GAM_OS_Niv = tibble(time = my_time, est=ilink(my_res$`fit`),
            lcl = ilink(my_res$`fit` + qnorm(0.025) * my_res$se.fit),
            ucl = ilink(my_res$`fit` + qnorm(1 - 0.025) * my_res$se.fit), Model = "GAM")  
  Xp = predict(object=tmp, newdata=new_df, type="lpmatrix")
  params = rmvn(n = num_PSA, coef(tmp), vcov(tmp))
  params2 = split(params, seq(nrow(params)))
  PSA_GAM_OS_Niv = map(params2, function(x) ilink(Xp %*% x))
  
###########################################
# FP(2), Powers = c(1,1)
# Not doing PSA for now (there is a simulate option, but this simulates the actual events, so is of limited use)
# If I were to include this in PSA, would probs use vcov(tmp), feed into rmvn, then manually generate predictions.
tmp = glm(Events ~ I(Time^1) + I(Time^1 * log(Time)) + offset(log(AtRisk)), data=df_Doc, family=poisson(link="log"))
  my_res = predict(object=tmp, newdata=new_df, type="link", se.fit=TRUE)
  ilink = family(tmp)$linkinv 
  FP_OS_Doc = tibble(time = my_time, est=ilink(my_res$`fit`),
            lcl = ilink(my_res$`fit` + qnorm(0.025) * my_res$se.fit),
            ucl = ilink(my_res$`fit` + qnorm(1 - 0.025) * my_res$se.fit), Model = "FP2")
tmp = glm(Events ~ I(Time^1) + I(Time^1 * log(Time)) + offset(log(AtRisk)), data=df_Niv, family=poisson(link="log"))
  my_res = predict(object=tmp, newdata=new_df, type="link", se.fit=TRUE)
  FP_OS_Niv = tibble(time = my_time, est=ilink(my_res$`fit`),
            lcl = ilink(my_res$`fit` + qnorm(0.025) * my_res$se.fit),
            ucl = ilink(my_res$`fit` + qnorm(1 - 0.025) * my_res$se.fit), Model = "FP2")

###########################################
# DSM: local trend and damped trend, load from DSM_IPDv2.R
DSM_IPD = readRDS(here("Files", "Mods_DSM_IPD_IG_Ln.rds"))  # Or "Mods_DSM_IPD_IG.rds"
LT_mod_Doc = DSM_IPD$DSM_LT_doc
  DT_mod_Doc = DSM_IPD$DSM_DT_doc
LT_mod_Niv = DSM_IPD$DSM_LT_niv
  DT_mod_Niv = DSM_IPD$DSM_DT_niv
  
fun_DSM_est = function(mod, df, mod_name, n_boot){
  mod_df = as.data.frame(mod)
  betas = select(mod_df, grep("beta1", names(mod_df), fixed = TRUE))
  # Within-sample estimates and CI - note these are for observed time
  tmp = tibble(time = df$Time, est=apply(exp(betas), 2, mean),
            lcl = apply(exp(betas), 2, quantile, 0.025),
            ucl = apply(exp(betas), 2, quantile, 1-0.025), Model = mod_name)
  # Now get estimates for required times (assume linear changes)
  max_fu = max(tmp$time)
  time_in = filter(new_df, Time <= max_fu)
  df_in = tibble(time = time_in$Time, est=approx(x=tmp$time, y=tmp$est, xout=time_in$Time, rule=2)$y,
            lcl = approx(x=tmp$time, y=tmp$lcl, xout=time_in$Time, rule=2)$y,
            ucl = approx(x=tmp$time, y=tmp$ucl, xout=time_in$Time, rule=2)$y, Model = mod_name)
  # PSA for within-sample
  params = sample_n(betas, size = n_boot, replace = TRUE)
  params2 = split(params, seq(nrow(params))) # Turn into list, to use with map.
  PSA_in = map(params2, function(par) approx(x=tmp$time, y=par, xout=time_in$Time, rule=2)$y)

  # Now for extrapolations - for now use bootstrapping
  time_ext = filter(new_df, Time > max_fu)
  tmp = tibble(Level = sample(mod_df$level, size = n_boot, replace = TRUE),  # PSA samples
               Trend = sample(mod_df$trend, size = n_boot, replace = TRUE))
  #### Full posterior 
  full_level = select(mod_df, grep("level", names(mod_df), fixed = TRUE))  # This and below - 2 diff ways for same output
  full_trend = as_tibble(extract(mod, pars = c("trend")))
  tmp_full = tibble(Level = full_level$level, Trend = full_trend$trend)   
  
  if (mod_name == "Local") {
    PSA_out = pmap(tmp, function(Level, Trend) Level + Trend * (log1p(time_ext$Time) - log1p(max_fu)))
    df_tmp = pmap_dfc(tmp_full, function(Level, Trend) exp(Level + Trend * (log1p(time_ext$Time) - log1p(max_fu))))
    df_out = tibble(time = time_ext$Time, est = apply(df_tmp, 1, mean),
                    lcl = apply(df_tmp, 1, quantile, 0.025),
                    ucl = apply(df_tmp, 1, quantile, 1-0.025), Model = mod_name)
    # PSA_out = pmap(tmp, function(Level, Trend) Level + Trend * (time_ext$Time - max_fu))
  # } else if (mod_name == "Adapted"){
  #   tmp = bind_cols(tmp, Phi = sample(mod_df$phi, size = n_boot, replace = TRUE))
  #   ext_est = pmap_dfc(tmp, function(Level, Trend, Phi) Level + Trend * Phi * (time_ext$Time - max_fu))
  #   PSA_out = pmap(tmp, function(Level, Trend, Phi) Level + Trend * Phi * (time_ext$Time - max_fu))
  } else { # Damped
    tmp = bind_cols(tmp, Phi = sample(mod_df$phi, size = n_boot, replace = TRUE))  # Add info on phi
    full_phi = as_tibble(extract(mod, pars = c("phi")))
    tmp_full =  bind_cols(tmp_full, Phi= full_phi$phi)  
######## Forecast times of future dampenings.
  my_data = list(
    T = length(df$Events),
    tau = log(df$Ln_Tau)
  )
  fit = stan(
    file = here("DSM", "Widths.stan"),      # Also have Widths_Level, but lot slower and struggles with convergence
    data = my_data,          # named list of data
    chains = 4,              # number of Markov chains
    warmup = 1000,           # number of warmup iterations per chain
    iter = 2000,             # total number of iterations per chain
    #init = init_list,        # Initial values for DSM models
    cores = 4,               # number of cores
    refresh = 0,              # show progress every 'refresh' iterations
    control = list(adapt_delta = 0.999, max_treedepth = 15)      # http://mc-stan.org/misc/warnings.html#maximum-treedepth-exceeded
    )
  
tmp_df = tibble(Time = seq(1:1000), level = extract(fit, pars = c("level")), trend = extract(fit, pars = c("trend")),
           ext_est = pmap(list(level, trend, Time), function(level, trend, Time) exp(level + trend * Time)),
         Pred = map_dbl(ext_est, mean)) %>% select(Time, Pred) #%>% filter(Pred < max_h) # So cumsum doesn't get too big
w_df = tibble(time_new = cumsum(c(max_fu, tmp_df$Pred)))

######## Use forecasted times to get extraps.
  fu_time = filter(new_df, Time > max_fu) %>% select(Time) %>% mutate(Time = Time)   # Times we want extraps for
  PSA_out = bind_cols(Time = w_df$time_new, 
              pmap_dfc(tmp, function(Level, Trend, Phi) Level + Trend * Phi * (1 - Phi^(log1p(w_df$time_new) - log1p(max_fu)))/(1 - Phi))) %>%
    gather(-Time, key = "ID", value = "Est") %>% nest(-ID, .key = "Vals") %>%
    mutate(yout = map(Vals, function(df) approx(x = df$Time, y = df$Est, xout = fu_time$Time)$y)) %>% select(yout) 
  PSA_out = PSA_out$yout

  tmp_df = tibble(Time = log1p(w_df$time_new) - log1p(max_fu)) %>% mutate(level = extract(mod, pars = c("level")),
                     trend = extract(mod, pars = c("trend")), phi = extract(mod, pars = c("phi")),
   tmp_est = pmap(list(level, trend, phi, Time), function(level, trend, phi, Time) level + trend * phi * (1 - phi^Time)/(1 - phi)),
   ext_est = map(tmp_est, exp),
           Pred = map_dbl(ext_est, mean), Time = w_df$time_new,  # Note over-writing time here
           Low = map_dbl(ext_est, function(x) quantile(x, probs = 0.025)),
           Upp = map_dbl(ext_est, function(x) quantile(x, probs = 1 - 0.025))) %>% select(Time, Pred, Low, Upp)
  df_out = tibble(time = time_ext$Time, est = approx(x = tmp_df$Time, y = tmp_df$Pred, xout = fu_time$Time)$y,
                lcl = approx(x = tmp_df$Time, y = tmp_df$Low, xout = fu_time$Time)$y,
                ucl = approx(x = tmp_df$Time, y = tmp_df$Upp, xout = fu_time$Time)$y, Model = mod_name)
  }
  df_full = bind_rows(df_in, df_out)
  PSA_df = map2(PSA_in, PSA_out, function(x, y) exp(c(x, y)))
  return(list(df_full = df_full, PSA_df = PSA_df))
}

num_b = 2000
tmp = fun_DSM_est(LT_mod_Doc, df_Doc, "Local", num_b)
  LT_OS_Doc = tmp$df_full
  PSA_LT_OS_Doc = tmp$PSA_df
tmp = fun_DSM_est(DT_mod_Doc, df_Doc, "Damped", num_b)
  DT_OS_Doc = tmp$df_full
  PSA_DT_OS_Doc = tmp$PSA_df  
tmp = fun_DSM_est(LT_mod_Niv, df_Niv, "Local", num_b)
  LT_OS_Niv = tmp$df_full
  PSA_LT_OS_Niv = tmp$PSA_df
tmp = fun_DSM_est(DT_mod_Niv, df_Niv, "Damped", num_b)
  DT_OS_Niv = tmp$df_full
  PSA_DT_OS_Niv = tmp$PSA_df  

################# Consolidate into 1 df by arm and save ###########################
GAM_OS_Doc = GAM_OS_Doc %>% mutate(est = exp(qlogis(est)),   # Convert GAM into a rate
                                   lcl = exp(qlogis(lcl)),
                                   ucl = exp(qlogis(ucl)))
GAM_OS_Niv = GAM_OS_Niv %>% mutate(est = exp(qlogis(est)),
                                   lcl = exp(qlogis(lcl)),
                                   ucl = exp(qlogis(ucl)))
  
df_OS_Doc = bind_rows(ll_OS_Doc, RP_OS_Doc, GAM_OS_Doc, FP_OS_Doc, LT_OS_Doc, DT_OS_Doc)
df_OS_Niv = bind_rows(ll_OS_Niv, RP_OS_Niv, GAM_OS_Niv, FP_OS_Niv, LT_OS_Niv, DT_OS_Niv)
saveRDS(df_OS_Doc, here("Files", "TP_Doc.rds"))
saveRDS(df_OS_Niv, here("Files", "TP_Niv.rds"))

PSA_OS_Doc = tibble(PSA_ID = seq(1:num_PSA), ll = PSA_ll_OS_Doc, RP = PSA_RP_OS_Doc, GAM = PSA_GAM_OS_Doc, LT = PSA_LT_OS_Doc, DT = PSA_DT_OS_Doc)
PSA_OS_Niv = tibble(PSA_ID = seq(1:num_PSA), ll = PSA_ll_OS_Niv, RP = PSA_RP_OS_Niv, GAM = PSA_GAM_OS_Niv, LT = PSA_LT_OS_Niv, DT = PSA_DT_OS_Niv)
saveRDS(PSA_OS_Doc, here("Files", "TP_PSA_Doc.rds"))
saveRDS(PSA_OS_Niv, here("Files", "TP_PSA_Niv.rds"))

#############################################

# Graphs of hazards with CI over time
# df = df_OS_Doc %>% mutate(Model = fct_relevel(Model, c("Log-logistic","Royston-Parmar","GAM","FP2","Local","Damped")),
#                           est = case_when(Model %in% c("GAMv2", "Localv2", "Localv3") ~ Prob, TRUE ~ est)) #%>%
# Add DRSMs from GenPopV2.R
# PSA_OS_Doc2 = readRDS(here("Files", "TP_PSA_Doc_GPDSM.rds"))
df_OS_Doc2 = readRDS(here("Files", "TP_Doc_GPDSM.rds"))
df_OS_Doc$Model = fct_recode(df_OS_Doc$Model, "Trend (DSM)" = "Local", "Damped (DSM)" = "Damped")
df_OS_Doc2$Model = fct_recode(df_OS_Doc2$Model, "Trend (DRSM)" = "Trend", "Damped (DRSM)" = "Damped", "Level (DRSM)" = "Level")
df_OS_Doc = rbind(df_OS_Doc, df_OS_Doc2)

df = df_OS_Doc %>% mutate(Model = fct_relevel(Model, c("Log-logistic","Trend (DSM)","Damped (DSM)",
                                                       "Trend (DRSM)","Damped (DRSM)","Level (DRSM)",
                                                       "Royston-Parmar","GAM","FP2")))
ggplot(df, aes(x = time/12, y = est)) + geom_line() + geom_ribbon(aes(ymin = lcl, ymax = ucl), fill = "skyblue", alpha = 0.4) +
  labs(x = "Time", y = "Hazard (95% confidence interval)") + facet_wrap(~Model) + coord_cartesian(ylim = c(0, 0.5))
ggsave(here("Graphs", "Haz_OS_Doc.pdf"), width=9, height=6)
ggplot(df, aes(x = time/12, y = log(est))) + geom_line() + geom_ribbon(aes(ymin = log(lcl), ymax = log(ucl)), fill = "skyblue", alpha = 0.4) +
  labs(x = "Time", y = "Log hazard (95% confidence interval)") + facet_wrap(vars(Model)) + coord_cartesian(ylim = c(-10,0))
ggsave(here("Graphs", "LnHaz_OS_Doc.pdf"), width=9, height=6)

df = df_OS_Niv %>% mutate(Model = fct_relevel(Model, c("Log-logistic","Royston-Parmar","GAM","FP2","Local","Damped")))
df$Model = fct_recode(df$Model, "Local trend" = "Local", "Damped trend" = "Damped") 
ggplot(filter(df, !Model %in% c("Royston-Parmar", "FP2")), aes(x = time, y = est)) + geom_line() + 
  geom_ribbon(aes(ymin = lcl, ymax = ucl), fill = "skyblue", alpha = 0.4) +
  labs(x = "Time", y = "Hazard (95% confidence interval)") + facet_wrap(vars(Model)) + coord_cartesian(ylim = c(0,0.2))
ggsave(here("Graphs", "Haz_OS_Niv.pdf"), width=9, height=6)
ggplot(filter(df, !Model %in% c("Royston-Parmar", "FP2")), aes(x = time, y = log(est))) + geom_line() + 
  geom_ribbon(aes(ymin = log(lcl), ymax = log(ucl)), fill = "skyblue", alpha = 0.4) +
  labs(x = "Time", y = "Log-Hazard (95% confidence interval)") + facet_wrap(vars(Model)) + coord_cartesian(ylim = c(-10,0))
ggsave(here("Graphs", "LnHaz_OS_Niv.pdf"), width=9, height=6)
