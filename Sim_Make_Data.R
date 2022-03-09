library("tidyverse")
library("flexsurv")
library("here")

set.seed(4368)
my_n = 10^7  #^7 = 10million
small_n = 10^-4
my_times = c(seq(from=0, to = 6, by = 0.25), seq(from=7, to = 30, by = 1))
my_rand1 = runif(n=my_n)
my_rand2 = runif(n=my_n)
my_rand3 = runif(n=my_n)

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

OS_doc_IPD = read.table(here("Prep", "IPDdata_OS_doc.txt" ), header = TRUE)

mod = flexsurvreg(Surv(time, event) ~ 1, data = OS_doc_IPD, dist = "llogis")
# Simulate survival data from this model
Surv_dis= (1/my_rand1 - 1)^(-mod$coefficients[1]) * mod$coefficients[2]

# Add general population mortality
Eng_HMD = read.table(here("HMD", "UK_HMDv2.txt"), header=TRUE)
Eng_2016 = filter(Eng_HMD, Year == 2016) %>% mutate(Age2 = fct_recode(Age, `110` = "110+"),
                                                    Years = as.numeric(levels(Age2))[Age2] - 63) %>% filter(Years >= 0) %>%
  mutate(Hazard = haz_rate(qx, Years + 1, "rate"), CumHaz = cumsum(Hazard), CumProb = haz_rate(CumHaz, Years + 1, "prob")) # Median age = 63
  
# Sample year of death
Surv_gen = findInterval(my_rand2, Eng_2016$CumProb) + my_rand3

# Add censoring - uniform random 5 to 6 years
Censor = runif(n=my_n, min=5, max=6)

df = tibble(Obs_Surv = pmin(Surv_dis, Surv_gen, Censor), Event = case_when(Censor == Obs_Surv ~ 0, TRUE ~ 1),
            Tru_Surv = pmin(Surv_dis, Surv_gen), Event_tru=1, Dis_Surv = Surv_dis, Pop_Surv = Surv_gen, Censor = Censor,
            Pop_Haz = approx(x=Eng_2016$Years, y=Eng_2016$Hazard, xout=Obs_Surv)$y)

Gen_dth = tibble(Dth = case_when(df$Tru_Surv == df$Pop_Surv ~ 1, TRUE ~ 0),
                 Dth_Cens = case_when(df$Obs_Surv == df$Pop_Surv ~ 1, TRUE ~ 0))
mean(Gen_dth$Dth)  # 3.765% die of other causes... 
mean(Gen_dth$Dth_Cens)  # Of which 2.237 are observed 

write_csv(df, here("Surv_full.csv"))

fit = survfit(Surv(Tru_Surv, Event_tru) ~ 1, data = df)

df2 = summary(fit, times=my_times)
df2b = tibble(time = df2$time, n=df2$n, nevents = df2$n.event, nrisk = df2$n.risk, surv_tru = df2$surv) %>%
  mutate(haz_tru = nevents / (lag(nrisk, default = 0) * (time - lag(time, default = 0))))
df2b$haz_tru[1] = 0
write_csv(df2b, here("Surv_tru.csv"))

set.seed(1501)

n_sim = 200 # Number of simulations
n_obs = 300 # Observations per simulation
df_sims = tibble(ID = seq(1:n_sim)) %>% mutate(df_IPD = map(ID, function(x) slice_sample(df, n=n_obs)),
                                               df_tmp1 = map(df_IPD, function(x) survfit(Surv(Obs_Surv, Event) ~ 1, data = x)),
                                               df_tmp2 = map(df_tmp1, function(x) summary(x, times=my_times)),
                                               df_Agg = map(df_tmp2, function(x) tibble(time = x$time, nevents = x$n.event, nrisk = x$n.risk,
                                                                                        Tau = x$time - lag(x$time), AtRisk = nrisk * Tau,
                                                                                        Ln_Tau = log1p(x$time) - log1p(lag(x$time)),
                                                                                        Pop_Haz = approx(x=Eng_2016$Years, y=Eng_2016$Hazard, xout=time)$y))) %>%
  select(-c(df_tmp1, df_tmp2))
write_rds(df_sims, here("df_sims.rds"))

