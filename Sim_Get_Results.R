library("tidyverse")
library("here")
library("gridExtra")
theme_set(theme_light())  # GGplot theme

df_mods = read_rds(here("df_mods.rds"))
df_DRSM = read_rds(here("df_DRSM.rds"))
df_damped = read_rds(here("df_damped.rds"))
my_times = c(seq(from=0, to = 5, by = 0.5), seq(from=6, to = 30, by = 1))
df2b = read_csv(here("Surv_tru.csv"))
small_n = 1*10^-100

# For FCM and NRS need to add gen pop hazards
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

new_df = tibble(Time = my_times)
Eng_HMD = read.table(here("HMD", "UK_HMDv2.txt"), header=TRUE)
Eng_2016 = filter(Eng_HMD, Year == 2016) %>% mutate(Age2 = fct_recode(Age, `110` = "110+"),
                                                    Years = as.numeric(levels(Age2))[Age2] - 63,
                                                    Hazard = haz_rate(qx, Years + 1, "rate")) %>%
  filter(Years >= 0) # Median age = 63, treat as mean

Gen_pop = as_tibble(approx(x = Eng_2016$Years, y = Eng_2016$Hazard, xout = new_df$Time, rule = 2))
names(Gen_pop) = c("time", "Gen_Haz")

# Current results are as lists (with slightly different structures) - add info on true hazard
df_LL = select(df_mods, c(ID, mod_haz, Model)) %>% filter(Model == "LL") %>% unnest(cols = c(mod_haz)) %>% 
  left_join(df2b, by="time") %>% relocate(time, .after = Model) %>% left_join(Gen_pop, by="time") %>% mutate(est = est + Gen_Haz)
df_NRS = select(df_mods, c(ID, mod_haz, Model)) %>% filter(Model == "NRS") %>% unnest(cols = c(mod_haz))
  df_NRS = df_NRS %>% bind_cols(tibble(time = rep(my_times, dim(table(df_NRS$ID))))) %>% left_join(df2b, by="time") %>% rename(est = mod_haz) %>%
    left_join(Gen_pop, by="time") %>% mutate(est = est + Gen_Haz)
df_FCM = select(df_mods, c(ID, mod_haz, Model)) %>% filter(Model == "FCM") %>% unnest(cols = c(mod_haz)) %>% unnest(cols = c(mod_haz)) 
  df_FCM = df_FCM %>% bind_cols(tibble(time = rep(my_times, dim(table(df_FCM$ID))))) %>% left_join(df2b, by="time") %>% rename(est = Estimate) %>%
    left_join(Gen_pop, by="time") %>% mutate(est = est + Gen_Haz)
df_DSM = select(df_DRSM, c(ID, mod_haz)) %>% unnest(cols = c(mod_haz)) %>% left_join(df2b, by="time") %>% select(-c(lcl, ucl)) %>%
  relocate(time, .after = Model) %>% left_join(Gen_pop, by="time")
df_damp = select(df_damped, c(ID, mod_haz)) %>% unnest(cols = c(mod_haz)) %>% left_join(df2b, by="time") %>% select(-c(lcl, ucl)) %>%
  relocate(time, .after = Model) %>% left_join(Gen_pop, by="time")

df_all = rbind(df_LL, df_NRS, df_FCM, df_DSM, df_damp) %>% 
  mutate(est = case_when(est>1 ~ 1, TRUE ~ est),
         Err = log(est + small_n) - log(haz_tru + small_n),
         Err2 = (log(est + small_n) - log(haz_tru + small_n))^2,
         Model = fct_relevel(Model, c("LL", "NRS", "FCM", "Trend", "Damped")),
         Model = fct_recode(Model, "Log-logistic" = "LL"))


df_all2 = group_by(df_all, Model, time) %>% na.omit %>% summarise(MSE = mean(Err2), Bias = mean(Err), my_n = n(),
                                                                  MSE_MC = (mean((Err2 - MSE)^2) / (my_n-1))^0.5,
                                                                  Bias_MC = (MSE * my_n / (my_n-1))^0.5) %>%
  filter(time != 0) # Errors over time (not applicable at time zero)

rm(df_mods, df_DRSM, df_damped, df_LL, df_NRS, df_FCM, df_DSM, df_damp)
gc()
gc()

# Calculate MSE, Bias weighted over time-points
n_sim = 1000
df_all3 = left_join(df_all2, df2b, by = "time") %>% ungroup() %>%  # Check if need left_join
  mutate(MSE_new = map2(MSE, MSE_MC, function(x, y) rgamma(n_sim, shape = x^2 / y^2, rate = x / y^2)),
         Bias_new = map2(Bias, Bias_MC, function(x, y) rnorm(n_sim, mean = x, sd = y))) %>%
  unnest(cols = c(MSE_new, Bias_new)) %>% group_by(time) %>% mutate(id = row_number()) %>% ungroup() %>%
  group_by(Model, id) %>% nest() %>% # This gies n_sim rows per model-sim
  mutate(MSE_full = map_dbl(data, function(x) weighted.mean(x = x$MSE, w = x$surv_tru, na.rm=TRUE)),
         MSE_full2 = map_dbl(data, function(x) weighted.mean(x = x$MSE_new, w = x$surv_tru, na.rm=TRUE)),
         Bias_full = map_dbl(data, function(x) weighted.mean(x = abs(x$Bias), w = x$surv_tru, na.rm=TRUE)),
         Bias_full2 = map_dbl(data, function(x) weighted.mean(x = abs(x$Bias_new), w = x$surv_tru, na.rm=TRUE))) %>%
  group_by(Model) %>%  # (check if needed)
  summarise(MSE_full = mean(MSE_full), MSE_SE = sd(MSE_full2), MSE_full2 = mean(MSE_full2),
            Bias_full = mean(Bias_full), Bias_SE = sd(Bias_full2), Bias_full2 = mean(Bias_full2),
            MSE_low = MSE_full - 1.96 * MSE_SE, MSE_upp = MSE_full + 1.96 * MSE_SE,
            Bias_low = Bias_full - 1.96 * Bias_SE, Bias_upp = Bias_full + 1.96 * Bias_SE) %>% ungroup()
write_csv(df_all3, here("Res.csv"))

# Graphs
g_legend = function(fig){  # Function used to get legend from figures.
  tmp = ggplot_gtable(ggplot_build(fig))
  leg = which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend = tmp$grobs[[leg]]
  return(legend)
}

ggplot(df_all, aes(x=time, colour = Model)) + geom_line(aes(y=est, group=ID), alpha = 0.4) + geom_line(aes(y=haz_tru), colour="black") + 
  labs(x = "Time", y = "Hazard") + scale_color_brewer(palette="Set1") + facet_grid(cols = vars(Model)) + theme(legend.position = "none")
ggsave(here("Figs", "Fig_ests.pdf"), width=9, height=6)

f1 = ggplot(df_all2, aes(x=time, y=Bias, colour=Model)) + geom_line() + geom_hline(yintercept=0) + 
  scale_color_brewer(palette="Set1") + labs(x="Time") + theme(legend.position = "bottom")
f2 = ggplot(df_all2, aes(x=time, y=MSE, colour=Model)) + geom_line() + geom_hline(yintercept=0) + 
  scale_color_brewer(palette="Set1") + labs(x="Time") + theme(legend.position = "none")
f3 = g_legend(f1)

pdf(here("Figs", "Stats.pdf"), width=9, height=6)
  grid.arrange(arrangeGrob(f2, f1 + theme(legend.position = "none"), ncol=2), f3, nrow=2, heights = c(12, 1))
dev.off()

