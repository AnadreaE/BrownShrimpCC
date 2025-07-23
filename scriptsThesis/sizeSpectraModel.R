
#SIZE SPECTRA MODEL

library(ggplot2)
library(deSolve)
library(RColorBrewer)
library(BrownShrimp)
library(viridis)
library(dplyr)
library(tidyr)
library("colorspace")
library(lubridate)
library(latex2exp)
library(profvis)

#library(devtools)

#UPLOAD TEMPERATURE DATA

temperature_dataSet <- read.csv("./data/Temperature_data_1982-2018_MarBiol.csv", header = TRUE, sep = ';')
temperature_14_18 <- temperature_func(temperature_dataSet, "01/01/2014", "31/12/2018")
t_5years = seq(0,length(temperature_14_18$temperature)-0.1, by = 0.1)


state_r = c(P = 2, E = 0.1, L= 0.0928 *1.5 ,
            BF = BF+0.1, BM = BM+0.1,
            Pred = 0.05 )


#### TEST (1 out of 4) WITHOUT ANY KIND OF PREASURE (P NEITHER F) ####

noPreasure_params = parameters_solv
noPreasure_params$general_params$Fi = 0 #No fishery
noPreasure_params$general_params$Imax_ik = 0 #No predation


start <- Sys.time()
test_noPreasure <- solver_sizeClass.v4(t = t_5years, state = state_r, parameters = noPreasure_params, temperature_dataSet = temperature_14_18)
print(Sys.time() - start)

#Plot size spectra
plot_sizeSpectra(test_noPreasure, 2015, title = "NP")


#### TEST (2 out of 4) ONLY PREDATION ####

pred_params = parameters_solv
pred_params$general_params$Fi = 0 #No fishery


start <- Sys.time()
test_predation <- solver_sizeClass.v4(t = t_5years, state = state_r, parameters = pred_params, temperature_dataSet = temperature_14_18)
print(Sys.time() - start)

#Plot size spectra
plot_sizeSpectra(test_predation, 2015, title = "Pred")

#### TEST (3 out of 4) ONLY FISHERY ####

fishery_params = parameters_solv
fishery_params$general_params$Imax_ik = 0 #No predation
#fishery_params$general_params$Fi = 0.01 #reduce fishery

start <- Sys.time()
test_fishery <- solver_sizeClass.v4(t = t_5years, state = state_r, parameters = fishery_params, temperature_dataSet = temperature_14_18)
print(Sys.time() - start)

#Plot size spectra
plot_sizeSpectra(test_fishery, 2015, title = "Fi")

#### TEST (4 out of 4) PREDATION AND FISHERY####

bothPF_params = parameters_solv
bothPF_params$general_params$Fi = 0.1 #increased fishery

start <- Sys.time()
test_bothFP <- solver_sizeClass.v4(t = t_5years, state = state_r, parameters = bothPF_params, temperature_dataSet = temperature_14_18)
print(Sys.time() - start)

#Plot size spectra
plot_sizeSpectra(test_bothFP, 2015, title = "F&P")


#### PLOTS TO COMPARE ALL ####

#(1) L_avg over the time:
size_classes = female_cols <- paste0("BL", 1:8)

#calculate weighted average for each time step:

#first join F and M in one column
#Without preasure (1)
noPreasure_noSex = test_noPreasure %>%
  mutate(L1 = BF1 + BM1, L2 = BF2 + BM2, L3 = BF3 + BM3, L4 = BF4 + BM4,  L5 = BF5 + BM5) %>%
  mutate(ttlB = L1 + L2+ L3+ L4+ L5+ BF6+ BF7 + BF8) %>%
  select(25:29, 10:12,24, ncol(noPreasure_noSex) ) %>% #Column 24 HAS TO EB CHANGED WHEN RUNNING NEW SIM TO SELECT DATETIME
  filter(as.Date(dateTime) > as.Date("31/12/2014", format= "%d/%m/%Y" )) #delete first year 'warm up' period


L_avg_noPreasure = noPreasure_noSex %>%
  mutate(L_avg = ( L1*size_mean_F[1] + L2*size_mean_F[2]+ L3*size_mean_F[3]+ L4*size_mean_F[4] + L5*size_mean_F[5]
                  + BF6*size_mean_F[6] + BF7*size_mean_F[7]+ BF8*size_mean_F[8])/ttlB )

#Only Predation (2)
predation_noSex = test_predation %>%
  mutate(L1 = BF1 + BM1, L2 = BF2 + BM2, L3 = BF3 + BM3, L4 = BF4 + BM4,  L5 = BF5 + BM5) %>%
  mutate(ttlB = L1 + L2+ L3+ L4+ L5+ BF6+ BF7 + BF8) %>%
  select(25:29, 10:12,24, ncol(predation_noSex) ) %>% #Column 24 HAS TO EB CHANGED WHEN RUNNING NEW SIM TO SELECT DATETIME
  filter(as.Date(dateTime) > as.Date("31/12/2014", format= "%d/%m/%Y" )) #delete first year 'warm up' period

Lavg_pred = predation_noSex %>%
  mutate(L_avg = ( L1*size_mean_F[1] + L2*size_mean_F[2]+ L3*size_mean_F[3]+ L4*size_mean_F[4] + L5*size_mean_F[5]
                   + BF6*size_mean_F[6] + BF7*size_mean_F[7]+ BF8*size_mean_F[8])/ttlB )

#Only Fishery (3)
fishery_noSex = test_fishery %>%
  mutate(L1 = BF1 + BM1, L2 = BF2 + BM2, L3 = BF3 + BM3, L4 = BF4 + BM4,  L5 = BF5 + BM5) %>%
  mutate(ttlB = L1 + L2+ L3+ L4+ L5+ BF6+ BF7 + BF8) %>%
  select(25:29, 10:12,24, ncol(fishery_noSex) ) %>% #Column 24 HAS TO EB CHANGED WHEN RUNNING NEW SIM TO SELECT DATETIME
  filter(as.Date(dateTime) > as.Date("31/12/2014", format= "%d/%m/%Y" )) #delete first year 'warm up' period

Lavg_fishery = fishery_noSex %>%
  mutate(L_avg = ( L1*size_mean_F[1] + L2*size_mean_F[2]+ L3*size_mean_F[3]+ L4*size_mean_F[4] + L5*size_mean_F[5]
                   + BF6*size_mean_F[6] + BF7*size_mean_F[7]+ BF8*size_mean_F[8])/ttlB )

#Both predation & fishery (4)
bothPF_noSex = test_bothFP %>%
  mutate(L1 = BF1 + BM1, L2 = BF2 + BM2, L3 = BF3 + BM3, L4 = BF4 + BM4,  L5 = BF5 + BM5) %>%
  mutate(ttlB = L1 + L2+ L3+ L4+ L5+ BF6+ BF7 + BF8) %>%
  select(25:29, 10:12,24, ncol(bothPF_noSex) ) %>% #Column 24 HAS TO EB CHANGED WHEN RUNNING NEW SIM TO SELECT DATETIME
  filter(as.Date(dateTime) > as.Date("31/12/2014", format= "%d/%m/%Y" )) #delete first year 'warm up' period

Lavg_bothPF = bothPF_noSex %>%
  mutate(L_avg = ( L1*size_mean_F[1] + L2*size_mean_F[2]+ L3*size_mean_F[3]+ L4*size_mean_F[4] + L5*size_mean_F[5]
                   + BF6*size_mean_F[6] + BF7*size_mean_F[7]+ BF8*size_mean_F[8])/ttlB )
dev.off()

plot(noPreasure_noSex$dateTime, L_avg_noPreasure$L_avg, las = 1,  col = 'grey', cex = 0.8,
     xlab = 'years', ylab = paste(TeX("$l$"), '[cm]'), ylim = c(1, 5) ,type = 'l', lty = 2, lwd = 3, #TeX(r"($\ell$)" )
     main = 'Changes in average size ')
lines(Lavg_pred$dateTime, Lavg_pred$L_avg, col = 'indianred3', lty=1, lwd = 2.5)
lines(Lavg_fishery$dateTime, Lavg_fishery$L_avg, col = 'skyblue4', lty=1, lwd = 2.5)
lines(Lavg_bothPF$dateTime, Lavg_bothPF$L_avg, col = 'hotpink', lty=1, lwd = 2.5)

legend("topright", legend = c('no preasure', 'only predation', 'only fishery', 'both F & P'),
       fill = c('grey', 'indianred3', 'skyblue4','hotpink'), border = NA, bty = "n", y.intersp = 0.65)



plot(noPreasure_noSex$dateTime, L_avg_noPreasure$ttlB, las = 1,  col = 'grey', cex = 0.8,
     xlab = 'years', ylab = 'biomass [?]', type = 'l', lty = 2, lwd = 3,
     main = 'Changes in biomass ', ylim = c(0, 85) )
lines(Lavg_pred$dateTime, Lavg_pred$ttlB, col = 'indianred3', lty=1, lwd = 2.5)
lines(Lavg_fishery$dateTime, Lavg_fishery$ttlB, col = 'skyblue4', lty=1, lwd = 2.5)
lines(Lavg_bothPF$dateTime, Lavg_bothPF$ttlB, col = 'hotpink', lty=1, lwd = 2.5)

legend("topright", legend = c('no preasure', 'only predation', 'only fishery', 'both F & P'),
       fill = c('grey', 'indianred3', 'skyblue4','hotpink'), border = NA, bty = "n", y.intersp = 0.6)









