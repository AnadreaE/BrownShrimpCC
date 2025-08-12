
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

ble_data = read.csv("./data/BLE_Inlandslandungen_SpeiseKrabbe.csv", sep = ';')

t_5years = seq(0,length(temperature_14_18$temperature)-0.1)

#BF_p = c(1.1476,  2.739,  5.604, 11.65, 21.524, 17.976, 11.613,  5.346) #see calculation in file 'sizeSpectraMoldel.R'
#BM_p = c(1.464, 3.645, 8.855, 22.622, 37.614) #see calculation in file 'sizeSpectraMoldel.R'


#IC_parameterized = c(P = 2, E = 0.759, L= 0.38  ,
 #            BF = BF_p, BM = BM_p) #Initial conditions according to estimations with fishery landings data.

BF_p = c(3.801852,  7.063263,  12.657331, 25.166351, 49.297133, 45.876772, 34.672340,  19.396426) #see calculation in file 'sizeSpectraMoldel.R'
BM_p = c(4.573332  , 8.203546 , 13.791605 , 26.328928 , 43.371121  ) #see calculation in file 'sizeSpectraMoldel.R'
IC_parameterized.v2 = c(P = 2, E = 0.759, L= 0.38  ,
                        BF = BF_p, BM = BM_p) #Initial conditions with reduced area (only up to 20 m depth).


#### TEST V5 (1 out of 4) WITHOUT ANY KIND OF PREASURE (P NEITHER F) ####


noPreasure_params = parameters_solv
noPreasure_params$general_params$Fi = 0 #No fishery
noPreasure_params$general_params$Imax_ik = 0 #No predation

start <- Sys.time()
test_noPreasure <- solver_sizeClass.v5(t = t_5years, state = IC_parameterized.v2, parameters = noPreasure_params, temperature_dataSet = temperature_14_18)
print(Sys.time() - start)

#Plot size spectra
plot_sizeSpectra(test_noPreasure, 2017, title = "NP")


#### TEST V5 (2 out of 4) ONLY PREDATION ####

pred_params = parameters_solv
pred_params$general_params$Fi = 0 #No fishery


start <- Sys.time()
test_predation <- solver_sizeClass.v5(t = t_5years, state = IC_parameterized.v2, parameters = pred_params, temperature_dataSet = temperature_14_18)
print(Sys.time() - start)

#Plot size spectra
plot_sizeSpectra(test_predation, 2015, title = "Pred")

#### TEST V5 (3 out of 4) ONLY FISHERY ####

fishery_params = parameters_solv
fishery_params$general_params$Imax_ik = 0 #No predation
fishery_params$general_params$Fi = 0.05 #reduced fishery

start <- Sys.time()
test_fishery <- solver_sizeClass.v5(t = t_5years, state = IC_parameterized.v2, parameters = fishery_params, temperature_dataSet = temperature_14_18)
print(Sys.time() - start)

#Plot size spectra
plot_sizeSpectra(test_fishery, 2015, title = "Fi")

#### TEST V5 (4 out of 4) PREDATION AND FISHERY####

bothPF_params = parameters_solv
bothPF_params$general_params$Fi = 0.05 #reduced fishery

start <- Sys.time()
test_bothFP <- solver_sizeClass.v5(t = t_5years, state = IC_parameterized.v2, parameters = bothPF_params, temperature_dataSet = temperature_14_18)
print(Sys.time() - start)

#Plot size spectra
plot_sizeSpectra(test_bothFP, 2018, title = "F&P")

#Plot fishery catch (for only one year):
#define one year to be inspected
year_toSee = 2018
#subset data and sum up to monthly values:
sol_yearToSee = test_bothFP %>%
               filter(as.numeric(format(test_bothFP$dateTime,'%Y')) == year_toSee) %>%
              mutate(month = month(dateTime)) %>%
              group_by(month) %>%
              summarise(across(where(is.numeric), \(x) sum(x, na.rm = TRUE)))


par(mfrow = c(1,1))
#plot(sol_yearToSee$month, sol_yearToSee$catch.BF1,  type = 'b', main = 'Fished Kg per Km2')
#lines(1:length(monthly_Feffort), monthly_Feffort, col = 'red')

#calculate ttl landings in Germany
monthly_landing_year = sol_yearToSee$catch.BF1 * 7587.1 / 1000 # ttl area = 7587.1; convert Kg to t /1000; ttl area = 14701,5 (incl. till 30m)
plot(sol_yearToSee$month, monthly_landing_year,  type = 'b',  ylim = c(0,1800),
     main = paste('Monthly landings Germany', year_toSee), col = '#c994c7', lwd = 2)
lines(ble_data$month[ble_data$year == year_toSee], ble_data$t[ble_data$year == year_toSee], col = 'gray42', lty = 2, lwd = 2)



#stacked area plot

#test_pred.V5_red <- test_predation[ , c(3:18)] # delete timesteps column, plancton  and predator column
test_both.V5_red <- test_bothFP[ , c(3:17, 19)] # delete timesteps column, plancton  and predator column

#cc_long_2 <- gather(cc_2, key = "Type", value = "Value", P, E, L, J, J2, J3, J4, J5, A1, A2)
test_both_long <- test_both.V5_red %>%
  pivot_longer(-dateTime, names_to = "variable", values_to = "value")

ggplot(test_both_long, aes(x = dateTime, y = value, fill = variable)) +
  geom_area() +
  scale_fill_viridis_d() +
  labs(
    title = "System incl. predation ",
    x = "Date",
    y = "Biomass",
    fill = "Size class"
  )



#Look at plankton
par(mfrow = c(1,1))
plot(test_bothFP$time, test_bothFP$P, type = 'l', xlim = c(0,365), ylim = c(0,30) )

#### PLOTS TO COMPARE ALL ####

#(1) L_avg over the time:
size_classes = female_cols <- paste0("BL", 1:8)

#calculate weighted average for each time step:

#first join F and M in one column
#Without preasure (1)
noPreasure_noSex = test_noPreasure %>%
  mutate(L1 = BF1 + BM1, L2 = BF2 + BM2, L3 = BF3 + BM3, L4 = BF4 + BM4,  L5 = BF5 + BM5) %>%
  mutate(ttlB = L1 + L2+ L3+ L4+ L5+ BF6+ BF7 + BF8) %>%
  select(20:24, 10:12,25, 19 ) %>%
  filter(as.Date(dateTime) > as.Date("31/12/2014", format= "%d/%m/%Y" )) #delete first year 'warm up' period


L_avg_noPreasure = noPreasure_noSex %>%
  mutate(L_avg = ( L1*size_mean_F[1] + L2*size_mean_F[2]+ L3*size_mean_F[3]+ L4*size_mean_F[4] + L5*size_mean_F[5]
                   + BF6*size_mean_F[6] + BF7*size_mean_F[7]+ BF8*size_mean_F[8])/ttlB )

#Only Predation (2)
predation_noSex = test_predation %>%
  mutate(L1 = BF1 + BM1, L2 = BF2 + BM2, L3 = BF3 + BM3, L4 = BF4 + BM4,  L5 = BF5 + BM5) %>%
  mutate(ttlB = L1 + L2+ L3+ L4+ L5+ BF6+ BF7 + BF8) %>%
  select(20:24, 10:12,25, 19 ) %>% #Column 24 HAS TO EB CHANGED WHEN RUNNING NEW SIM TO SELECT DATETIME
  filter(as.Date(dateTime) > as.Date("31/12/2014", format= "%d/%m/%Y" )) #delete first year 'warm up' period

Lavg_pred = predation_noSex %>%
  mutate(L_avg = ( L1*size_mean_F[1] + L2*size_mean_F[2]+ L3*size_mean_F[3]+ L4*size_mean_F[4] + L5*size_mean_F[5]
                   + BF6*size_mean_F[6] + BF7*size_mean_F[7]+ BF8*size_mean_F[8])/ttlB )

#Only Fishery (3)
fishery_noSex = test_fishery %>%
  mutate(L1 = BF1 + BM1, L2 = BF2 + BM2, L3 = BF3 + BM3, L4 = BF4 + BM4,  L5 = BF5 + BM5) %>%
  mutate(ttlB = L1 + L2+ L3+ L4+ L5+ BF6+ BF7 + BF8) %>%
  select(20:24, 10:12,25, 19 ) %>% #Column 24 HAS TO EB CHANGED WHEN RUNNING NEW SIM TO SELECT DATETIME
  filter(as.Date(dateTime) > as.Date("31/12/2014", format= "%d/%m/%Y" )) #delete first year 'warm up' period

Lavg_fishery = fishery_noSex %>%
  mutate(L_avg = ( L1*size_mean_F[1] + L2*size_mean_F[2]+ L3*size_mean_F[3]+ L4*size_mean_F[4] + L5*size_mean_F[5]
                   + BF6*size_mean_F[6] + BF7*size_mean_F[7]+ BF8*size_mean_F[8])/ttlB )

#Both predation & fishery (4)
bothPF_noSex = test_bothFP %>%
  mutate(L1 = BF1 + BM1, L2 = BF2 + BM2, L3 = BF3 + BM3, L4 = BF4 + BM4,  L5 = BF5 + BM5) %>%
  mutate(ttlB = L1 + L2+ L3+ L4+ L5+ BF6+ BF7 + BF8) %>%
  select(20:24, 10:12,25, 19 ) %>% #Column 24 HAS TO EB CHANGED WHEN RUNNING NEW SIM TO SELECT DATETIME
  filter(as.Date(dateTime) > as.Date("31/12/2014", format= "%d/%m/%Y" )) #delete first year 'warm up' period

Lavg_bothPF = bothPF_noSex %>%
  mutate(L_avg = ( L1*size_mean_F[1] + L2*size_mean_F[2]+ L3*size_mean_F[3]+ L4*size_mean_F[4] + L5*size_mean_F[5]
                   + BF6*size_mean_F[6] + BF7*size_mean_F[7]+ BF8*size_mean_F[8])/ttlB )
#dev.off()

plot(noPreasure_noSex$dateTime, L_avg_noPreasure$L_avg, las = 1,  col = 'grey', cex = 0.8,
     xlab = 'years', ylab = paste(TeX("$l$"), '[cm]'), ylim = c(2, 8) ,type = 'l', lty = 2, lwd = 3, #TeX(r"($\ell$)" )
     main = 'Changes in average size ')
lines(Lavg_pred$dateTime, Lavg_pred$L_avg, col = 'indianred3', lty=1, lwd = 2.5)
lines(Lavg_fishery$dateTime, Lavg_fishery$L_avg, col = 'skyblue4', lty=1, lwd = 2.5)
lines(Lavg_bothPF$dateTime, Lavg_bothPF$L_avg, col = 'hotpink', lty=1, lwd = 2.5)

legend("topright", legend = c('no preasure', 'only predation', 'only fishery', 'both F & P'),
       fill = c('grey', 'indianred3', 'skyblue4','hotpink'), border = NA, bty = "n", y.intersp = 0.65)


par(mar = c(5, 4, 4, 2) + 0.5)
plot(noPreasure_noSex$dateTime, L_avg_noPreasure$ttlB, las = 1,  col = 'grey', cex = 0.8,
     xlab = 'years', ylab = TeX(" biomass Kg $Km^{-2}$"), type = 'l', lty = 2, lwd = 3,
     main = 'Changes in biomass ', ylim = c(0, 900) )
lines(Lavg_pred$dateTime, Lavg_pred$ttlB, col = 'indianred3', lty=1, lwd = 2.5)
lines(Lavg_fishery$dateTime, Lavg_fishery$ttlB, col = 'skyblue4', lty=1, lwd = 2.5)
lines(Lavg_bothPF$dateTime, Lavg_bothPF$ttlB, col = 'hotpink', lty=1, lwd = 2.5)

legend("topright", legend = c('no preasure', 'only predation', 'only fishery', 'both F & P'),
       fill = c('grey', 'indianred3', 'skyblue4','hotpink'), border = NA, bty = "n", y.intersp = 0.6)


#### now plot L_avg vs Biomass ####

#Now group by the week number and calculate avg of the size

weekly_avg_noPreasure  <- L_avg_noPreasure %>%
  mutate(
    year = year(dateTime),
    week = isoweek(dateTime)  # or week(dateTime), depending on your preference
  ) %>%
  group_by(year, week) %>%
  summarise(avg_L = mean(L_avg, na.rm = TRUE), avg_B = mean(ttlB, na.rm = TRUE), .groups = "drop")

weekly_avg_fish  <- Lavg_fishery %>%
  mutate(
    year = year(dateTime),
    week = isoweek(dateTime)  # or week(dateTime), depending on your preference
  ) %>%
  group_by(year, week) %>%
  summarise(avg_L = mean(L_avg, na.rm = TRUE), avg_B = mean(ttlB, na.rm = TRUE), .groups = "drop")

weekly_avg_pred  <- Lavg_pred %>%
  mutate(
    year = year(dateTime),
    week = isoweek(dateTime)  # or week(dateTime), depending on your preference
  ) %>%
  group_by(year, week) %>%
  summarise(avg_L = mean(L_avg, na.rm = TRUE), avg_B = mean(ttlB, na.rm = TRUE), .groups = "drop")

weekly_avg_both  <- Lavg_bothPF %>%
  mutate(
    year = year(dateTime),
    week = isoweek(dateTime)  # or week(dateTime), depending on your preference
  ) %>%
  group_by(year, week) %>%
  summarise(avg_L = mean(L_avg, na.rm = TRUE), avg_B = mean(ttlB, na.rm = TRUE), .groups = "drop")

init15 = as.POSIXct("2015-01-01", format="%Y-%m-%d", tz = "UTC")
end15 = as.POSIXct("2015-12-31", format="%Y-%m-%d", tz = "UTC")


plot_year = 2017
plot(weekly_avg_noPreasure$avg_L[weekly_avg_noPreasure$year == plot_year], weekly_avg_noPreasure$avg_B[weekly_avg_noPreasure$year == plot_year],
     lwd = 3, col = 'grey', lty=1, ylim = c(20,500), xlim = c(3.25, 5.2) ,
     xlab = 'L in cm', ylab = TeX(" biomass Kg $Km^{-2}$"), main = paste('Fi baseline ', bothPF_params$general_params$Fi, '\n ', plot_year))
points(weekly_avg_fish$avg_L[weekly_avg_fish$year == plot_year], weekly_avg_fish$avg_B[weekly_avg_fish$year == plot_year], lwd = 2, col ='skyblue4' )
points(weekly_avg_fish$avg_L[weekly_avg_fish$year == plot_year][1], weekly_avg_fish$avg_B[weekly_avg_fish$year == plot_year][1], lwd = 2, col ='skyblue', cex = 1.5, pch = 8 )

points(weekly_avg_pred$avg_L[weekly_avg_pred$year == plot_year], weekly_avg_pred$avg_B[weekly_avg_pred$year == plot_year], lwd = 2, col ='indianred4' )
points(weekly_avg_pred$avg_L[weekly_avg_pred$year == plot_year][1], weekly_avg_pred$avg_B[weekly_avg_pred$year == plot_year][1], lwd = 2, col ='indianred', cex = 1.5, pch = 8 )

points(weekly_avg_both$avg_L[weekly_avg_both$year == plot_year], weekly_avg_both$avg_B[weekly_avg_pred$year == plot_year], lwd = 2, col= 'hotpink3' )
points(weekly_avg_both$avg_L[weekly_avg_both$year == plot_year][1], weekly_avg_both$avg_B[weekly_avg_pred$year == plot_year][1], lwd = 2, col= 'hotpink', cex = 1.5, pch = 8 )


legend("topleft", legend = c('no preasure', 'only predation', 'only fishery', 'both F & P'),
       fill = c('grey', 'indianred3', 'skyblue4','hotpink'), border = NA, bty = "n", y.intersp = 0.7)



