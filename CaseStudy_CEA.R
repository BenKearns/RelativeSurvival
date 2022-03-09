library("tidyverse")
library("flexsurv")
library("here")
set.seed(4361)

theme_set(theme_light())
my_time = seq(from=(1/4), to=(20*13), by=(1/4)) # Time values for which we want estimates.
new_df = tibble(Time = my_time, AtRisk = 1)
num_PSA = 2000
time_cycles = length(my_time) + 1  # Counting time = 0
cycle_times = seq(from=1, to=time_cycles, by=1) - 1 # Could half-cycle correct these
cycle_length = max(my_time)/ length(my_time)

disc_cost = 0.035
disc_QALY = 0.035

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

to_list = function(x){ # Convert data-frame or matrix to list
  y = split(x, seq(nrow(x)))
}

########## Load data for PSA (OS estimates)
PSA_OS_Doc = readRDS(here("Files", "TP_PSA_Doc.rds"))
PSA_Niv_HR = rlnorm(num_PSA, -0.52756051, 0.099751345)

########## PSA - PFS: RP 2 knots
PFS_doc_IPD = read.table(here("Prep", "IPDdata_PFS_doc.txt" ), header = TRUE)  # Digitised data
PFS_niv_IPD = read.table(here("Prep", "IPDdata_PFS_niv.txt" ), header = TRUE)

PFS_RP_PSA = function(df){
  tmp = flexsurvspline(Surv(time, event) ~ 1, k = 2, data = df, scale = "hazard")
  params = normboot.flexsurvreg(tmp, num_PSA)
  params2 = to_list(params) # Turn into list, to use with map
  tmp2 = map(params2, function(x) hsurvspline(my_time, gamma = x, knots = tmp$knots, scale = "hazard", timescale="log"))
  return(tmp2)
}

PSA_Niv_PFS = PFS_RP_PSA(PFS_niv_IPD)
PSA_Doc_PFS = PFS_RP_PSA(PFS_doc_IPD)

########## PSA - costs
PSA_costs = matrix(nrow = num_PSA, ncol = 9)
colnames(PSA_costs) = c("Cost_dis_manage_SD", "Cost_dis_manage_PD", "Cost_drug_acq_SD_niv", "Cost_drug_acq_SD_doc",
                        "Cost_drug_acq_PD", "Cost_admin_niv", "Cost_admin_doc", "Cost_monitor", "Cost_death")
PSA_costs[,"Cost_dis_manage_SD"] = rgamma(num_PSA, shape = 10*313.55, scale = 1/10) / 4
PSA_costs[,"Cost_dis_manage_PD"] = rgamma(num_PSA, shape = 10*766.62, scale = 1/10) / 4
PSA_costs[,"Cost_drug_acq_SD_niv"] = rgamma(num_PSA, shape = 10*2634, scale = 1/10) / 2
PSA_costs[,"Cost_drug_acq_SD_doc"] = rgamma(num_PSA, shape = 10*900, scale = 1/10) / 3
PSA_costs[,"Cost_drug_acq_PD"] = rep(0, num_PSA)
PSA_costs[,"Cost_admin_niv"] = rgamma(num_PSA, shape = 10*269.92, scale = 1/10)
PSA_costs[,"Cost_admin_doc"] = rgamma(num_PSA, shape = 10*167.34, scale = 1/10)
PSA_costs[,"Cost_monitor"] = rgamma(num_PSA, shape = 10*151.89, scale = 1/10) / 4
PSA_costs[,"Cost_death"] = rgamma(num_PSA, shape = 10*3628.70, scale = 1/10)

PSA_Niv_cost = matrix(nrow = num_PSA, ncol = 3)
PSA_Niv_cost[,1] = PSA_costs[,"Cost_dis_manage_SD"] + PSA_costs[,"Cost_drug_acq_SD_niv"] + PSA_costs[,"Cost_admin_niv"] + PSA_costs[,"Cost_monitor"]
PSA_Niv_cost[,2] = PSA_costs[,"Cost_dis_manage_PD"] + PSA_costs[,"Cost_drug_acq_PD"]
PSA_Niv_cost[,3] = PSA_costs[,"Cost_death"]
PSA_Doc_cost = PSA_Niv_cost
PSA_Doc_cost[,1] = PSA_costs[,"Cost_dis_manage_SD"] + PSA_costs[,"Cost_drug_acq_SD_doc"] + PSA_costs[,"Cost_admin_doc"] + PSA_costs[,"Cost_monitor"]

########## PSA - QALYs
fun_beta_param = function(my_mean, my_stdev){
  alpha = (((my_mean * (1- my_mean)) / my_stdev^2) - 1) * my_mean
  beta = alpha * (1 - my_mean) / my_mean
  return(c(alpha, beta))
}
PSA_Niv_qalys = matrix(nrow = num_PSA, ncol = 3)
PSA_Niv_qalys = cbind(rbeta(num_PSA, fun_beta_param(0.75, 0.236)[1], fun_beta_param(0.75, 0.236)[2]),
                      rbeta(num_PSA, fun_beta_param(0.592, 0.315)[1], fun_beta_param(0.592, 0.315)[2]), rep(0, num_PSA))
PSA_Doc_qalys = PSA_Niv_qalys

########## Combine PSA inputs into a single object
PSA_Doc = PSA_OS_Doc %>% mutate(PFS = PSA_Doc_PFS, Cost = to_list(PSA_Doc_cost), QALY = to_list(PSA_Doc_qalys))
PSA_Niv = PSA_OS_Doc %>% mutate(HR = PSA_Niv_HR, PFS = PSA_Niv_PFS, Cost = to_list(PSA_Niv_cost), QALY = to_list(PSA_Niv_qalys))

########## Code for DAM
fun_PartSA = function(OS, PFS, Cost, Util, Name="Tx", Total = FALSE, Gen_pop = 0){
  # Function to generate the Markov trace and derived values from OS and PFS estimates.
  # If Total = TRUE, only return totals
  
  # Inputs are hazards, convert to probabilities.
  tmp = my_time - lag(my_time, default=0)
  PFS = 1 - exp(-PFS * tmp)
  OS = 1 - exp(-OS * tmp)

  # Constrain OS hazard to never fall below Gen_pop value
  OS = pmax(OS, Gen_pop)
  
  # Check constraint that hazard for PFS >= hazard for OS
  PFS = case_when(PFS < OS ~ OS,
                  TRUE ~ PFS)
    
  # Define an array for each health state
  vec_stable = vector(mode = "numeric", length = time_cycles)
  vec_prog = vector(mode = "numeric", length = time_cycles)
  vec_death = vector(mode = "numeric", length = time_cycles)
  
  # Initial values
  vec_stable[1] = 1
  vec_prog[1] = 0
  vec_death[1] = 0
  
  for (i in 2:time_cycles) {
    vec_death[i] = vec_death[i-1] + (vec_stable[i-1] + vec_prog[i-1]) * OS[i-1]
    vec_stable[i] = min(vec_stable[i-1] * (1 - PFS[i-1]), 1 - vec_death[i])
    vec_prog[i] = 1 - (vec_death[i] + vec_stable[i])
  }
  
  trace_mat = bind_cols("Stable" = vec_stable,
                        "Progressed" = vec_prog,
                        "Dead" = vec_death)
  
  trace_costs = sweep(trace_mat, MARGIN=2, Cost, "*")
  # Replace death cost with proper cost
  trace_costs$Dead = (trace_mat$Dead - lag(trace_mat$Dead, default = 0)) * Cost[3]
  trace_costs_disc = apply(trace_costs, 2, function(x) x / (1 + disc_cost)^c(0,floor(my_time/12)))
  trace_lys = (trace_mat * cycle_length)/12
  trace_lys$Dead = 0
  trace_qalys = sweep(trace_lys, MARGIN=2, Util, "*")
  trace_qalys_disc = apply(trace_qalys, 2, function(x) x / (1 + disc_QALY)^c(0,floor(my_time/12)))
  
  my_list1 = list(trace_mat = trace_mat,
                  trace_lys = trace_lys,
                  trace_costs = trace_costs,
                  trace_qalys = trace_qalys,
                  trace_costs_disc = trace_costs_disc,
                  trace_qalys_disc = trace_qalys_disc,
                  Name = Name)
  
  my_list2 = list(total_lys = sum(trace_lys),
                  total_costs = sum(trace_costs),
                  total_qalys = sum(trace_qalys),
                  total_costs_disc = sum(trace_costs_disc),
                  total_qalys_disc = sum(trace_qalys_disc))
  
  my_list = if(Total) {my_list2} else {my_list1}
  
  return(my_list)
}

PSA_Doc = PSA_Doc %>% gather(key = "Model", value = "OS", -c(PSA_ID, PFS, Cost, QALY)) %>%
  mutate(Res = pmap(list(OS=OS, PFS=PFS, Cost=Cost, Util=QALY),
                    function(OS, PFS, Cost, Util) fun_PartSA(OS, PFS, Cost, Util, Total=TRUE)))
PSA_Doc_res = do.call(rbind.data.frame, PSA_Doc$Res) %>% bind_cols(Model = PSA_Doc$Model)
tmp = PSA_Doc_res %>% group_by(Model) %>% summarise(Cost = mean(total_costs_disc),
                                              QALY = mean(total_qalys_disc))
write.csv(PSA_Doc_res, here("Tables", "PSA_Doc.csv"))
write.csv(tmp, here("Tables", "PSA_Doc2.csv"))

PSA_Niv = PSA_Niv %>% gather(key = "Model", value = "OS", -c(HR, PSA_ID, PFS, Cost, QALY)) %>%
  mutate(Res = pmap(list(OS=OS, HR=HR, PFS=PFS, Cost=Cost, Util=QALY),
                    function(OS, HR, PFS, Cost, Util) fun_PartSA(OS*HR, PFS, Cost, Util, Total=TRUE)))
PSA_Niv_res = do.call(rbind.data.frame, PSA_Niv$Res) %>% bind_cols(Model = PSA_Niv$Model)
tmp = PSA_Niv_res %>% group_by(Model) %>% summarise(Cost = mean(total_costs_disc),
                                                    QALY = mean(total_qalys_disc))
write.csv(PSA_Niv_res, here("Tables", "PSA_Niv.csv"))
write.csv(tmp, here("Tables", "PSA_Niv2.csv"))

# Above, but incorporate general population data
Eng_HMD = read.table(here("HMD", "UK_HMDv2.txt"), header=TRUE)
Eng_2016 = filter(Eng_HMD, Year == 2016) %>% mutate(Age2 = fct_recode(Age, `110` = "110+"),
                                                    Years = as.numeric(levels(Age2))[Age2] - 63,
                                                    Hazard = haz_rate(qx, Years + 1, "rate")) %>%  # Convert qx to hazard
  filter(Years >= 0) # Median age = 63, treat as mean                                                CHECK IF QANT QX or MX
Gen_haz = tibble(Time = my_time, Years = floor(my_time/12)) %>% left_join(Eng_2016, by = "Years")

PSA_Doc = PSA_Doc %>%
  mutate(Res2 = pmap(list(OS=OS, PFS=PFS, Cost=Cost, Util=QALY),
                    function(OS, PFS, Cost, Util) fun_PartSA(OS, PFS, Cost, Util, Total=TRUE, Gen_pop = Gen_haz$Hazard)))
PSA_Doc_res = do.call(rbind.data.frame, PSA_Doc$Res2) %>% bind_cols(Model = PSA_Doc$Model)
tmp = PSA_Doc_res %>% group_by(Model) %>% summarise(Cost = mean(total_costs_disc),
                                                    QALY = mean(total_qalys_disc),
                                                    LY = mean(total_lys))
write.csv(PSA_Doc_res, here("Tables", "PSA_Doc_GP.csv"))
write.csv(tmp, here("Tables", "PSA_Doc_GP2.csv"))

PSA_Niv = PSA_Niv %>%
  mutate(Res2 = pmap(list(OS=OS, HR=HR, PFS=PFS, Cost=Cost, Util=QALY),
                    function(OS, HR, PFS, Cost, Util) fun_PartSA(OS*HR, PFS, Cost, Util, Total=TRUE, Gen_pop = Gen_haz$Hazard)))
PSA_Niv_res = do.call(rbind.data.frame, PSA_Niv$Res2) %>% bind_cols(Model = PSA_Niv$Model)
tmp = PSA_Niv_res %>% group_by(Model) %>% summarise(Cost = mean(total_costs_disc),
                                                    QALY = mean(total_qalys_disc),
                                                    LY = mean(total_lys))
write.csv(PSA_Niv_res, here("Tables", "PSA_Niv_GP.csv"))
write.csv(tmp, here("Tables", "PSA_Niv2_GP.csv"))

#########################################################################
# Now just for DSMs (Local trend, damped trend) and DRSMs (local trend, damped trend, local level)
PSA_OS_Doc2 = readRDS(here("Files", "TP_PSA_Doc_GPDSM.rds")) %>% left_join(PSA_OS_Doc, by="PSA_ID") %>% select(-c("RP", "GAM"))
PSA_Doc2 = PSA_OS_Doc2 %>% mutate(PFS = PSA_Doc_PFS, Cost = to_list(PSA_Doc_cost), QALY = to_list(PSA_Doc_qalys)) %>%
  gather(key = "Model", value = "OS", -c(PSA_ID, PFS, Cost, QALY)) %>%
  mutate(Res = pmap(list(OS=OS, PFS=PFS, Cost=Cost, Util=QALY),
                    function(OS, PFS, Cost, Util) fun_PartSA(OS, PFS, Cost, Util, Total=TRUE)))
PSA_Doc_res2 = do.call(rbind.data.frame, PSA_Doc2$Res) %>% bind_cols(Model = PSA_Doc2$Model)
tmp = PSA_Doc_res2 %>% group_by(Model) %>% summarise(Cost = mean(total_costs_disc),
                                              QALY = mean(total_qalys_disc),
                                              LY = mean(total_lys))
write.csv(PSA_Doc_res2, here("Tables", "PSA_DocA.csv"))
write.csv(tmp, here("Tables", "PSA_Doc2A.csv"))

PSA_Niv2 = PSA_OS_Doc2 %>% mutate(HR = PSA_Niv_HR, PFS = PSA_Niv_PFS, Cost = to_list(PSA_Niv_cost), QALY = to_list(PSA_Niv_qalys)) %>%
  gather(key = "Model", value = "OS", -c(HR, PSA_ID, PFS, Cost, QALY)) %>%
  mutate(Res = pmap(list(OS=OS, HR=HR, PFS=PFS, Cost=Cost, Util=QALY),
                    function(OS, HR, PFS, Cost, Util) fun_PartSA(OS*HR, PFS, Cost, Util, Total=TRUE)))
PSA_Niv_res2 = do.call(rbind.data.frame, PSA_Niv2$Res) %>% bind_cols(Model = PSA_Niv2$Model)
tmp = PSA_Niv_res2 %>% group_by(Model) %>% summarise(Cost = mean(total_costs_disc),
                                                    QALY = mean(total_qalys_disc),
                                                    LY = mean(total_lys))
write.csv(PSA_Niv_res2, here("Tables", "PSA_NivA.csv"))
write.csv(tmp, here("Tables", "PSA_Niv2A.csv"))

# # CEAC
PSA_Inc = tibble(Costs = PSA_Niv_res2$total_costs_disc - PSA_Doc_res2$total_costs_disc,
                 QALYs = PSA_Niv_res2$total_qalys_disc - PSA_Doc_res2$total_qalys_disc,
                 Model = fct_relevel(PSA_Doc_res2$Model, c("ll", "LT", "DT", "Level", "Trend")))
PSA_Inc = PSA_Inc %>% mutate(Model = fct_recode(PSA_Inc$Model, `Log-logistic` = "ll", `DSM: Damped trend` = "DT", `DSM: Local trend` = "LT",
                                    `DRSM: Damped trend` = "Damped", `DRSM: Local trend` = "Trend", `DRSM: Local level` = "Level"))
f1 = ggplot(data = PSA_Inc, aes(x = QALYs, y = Costs, shape = Model, colour = Model)) +
  geom_hline(yintercept = 0, colour = "grey60") + geom_vline(xintercept = 0, colour = "grey60") + geom_point(alpha = 0.2) +
  scale_shape(solid = FALSE) + stat_ellipse() + labs(x="Incremental QALYs", y="Incremental Costs") +
  scale_y_continuous(labels=function(x) format(x, big.mark = ",", scientific = FALSE)) + theme(legend.position = "bottom")
ggsave(here("Graphs", "Scatter_all2A.pdf"), f1, height = 6, width = 8)
