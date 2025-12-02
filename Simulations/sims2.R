#<!--
# SPDX-FileCopyrightText: 2025 Helmholtz-Zentrum hereon GmbH
# SPDX-License-Identifier: CC0-1.0
# SPDX-FileContributor Andrea Farfan <farfanqbb@gmail.de>
#  -->


#########################
## In this file, we want to simulate time spam from 2010 to 2023, considering the
## respective changes in conditions between these years.
## (1) Predation has slightly increased / Joanna K. Bluemel et al 2021
## (2) fishery net selectivity introduced in 2016, and amendded in 2019


library(deSolve)
library(RColorBrewer)
library(BrownShrimp)
library(viridis)
library(dplyr)
library(tidyr)
#library("colorspace")
#library(lubridate)
library(latex2exp)
#library(profvis)
library(stringr)

#library(devtools)

#Read water temperature from Wadden Sea:
temperature_dataSet <- read.csv("./data/Temperature_data_1982-2018_MarBiol.csv", header = TRUE, sep = ';')
temperature_10_15 <- temperature_func(temperature_dataSet, "01/01/2010", "31/12/2015")


##########

temperature_germany <- read.csv("./data/temperature_10ger.csv")

tempG.9_21 = temperature_germany %>% filter(date_time > as.Date("31/12/2009", format = "%d/%m/%Y") & date_time < as.Date("01/01/2022", format = "%d/%m/%Y"))

ble_data = read.csv("./data/BLE_Inlandslandungen_SpeiseKrabbe.csv", sep = ';')

t_years = seq(0,length(tempG.9_21$temperature)-.1)

BF = c(8.695212, 15.349638, 21.461595, 24.672213,25.871618, 17.536662, 14.090960, 11.933151)
BM = c(10.85817, 19.93554, 29.01086, 36.67041, 58.26882)

Init.v3 = c(P = 2, E = 2, L= 1 ,
            BF = BF, BM = BM)

inti_y = 2010
fin_y = 2021
years = seq(inti_y, fin_y)

#### #Sim. baseline scenario: reality ####

Imax_ik_vals = seq(0.17, 0.19, length = length(seq(inti_y,fin_y)))
L50_vals = c(rep(3.69, length(seq(inti_y,2015))), rep(4.02, length(seq(2016,2019))),  rep(4.33, length(seq(2020,fin_y))))
SR_vals = c(rep(0.75, length(seq(inti_y,2015))), rep(.82, length(seq(2016,2019))),  rep(0.9, length(seq(2020,fin_y))))

PF_vars_df = data.frame(year=years, Imax = Imax_ik_vals, L50 = L50_vals, SR= SR_vals)


bothPF_params = parameters_solv.v2
bothPF_params$general_params$Fi = 4/365 #Temming & Hufnagl 2015

tempG.9_21 = temperature_germany %>% filter(date_time > as.Date("31/12/2009", format = "%d/%m/%Y") & date_time < as.Date("01/01/2022", format = "%d/%m/%Y"))
t_years = seq(0,length(tempG.9_21$temperature)-.1)


start <- Sys.time()
FP10_21 <- solver_sizeClass.v6(t = t_years, state = Init.v3, parameters = bothPF_params,
                               temperature_dataSet = tempG.9_21, PF_dataset = PF_vars_df)
print(Sys.time() - start)


years_to_see = c(2019, 2020, 2021)

par(mfrow = c(1,length(years_to_see)))
for (i in years_to_see){
  #Solution with cod-end 22 mm
  sol_yearToSee = FP10_21 %>%
    filter(as.numeric(format(FP10_21$dateTime,'%Y')) == i) %>%
    mutate(month = month(dateTime)) %>%
    group_by(month) %>%
    summarise(across(where(is.numeric), \(x) sum(x, na.rm = TRUE)))

  #Solution with cod-end 20 mm

  #sol_2 = FP_15_18_g.20 %>%
  #  filter(as.numeric(format(FP_15_18_g.20$dateTime,'%Y')) == i) %>%
  #  mutate(month = month(dateTime)) %>%
  #  group_by(month) %>%
  #  summarise(across(where(is.numeric), \(x) sum(x, na.rm = TRUE)))

  #Calculate the total landings per month (considering 50% discard):
  plot(ble_data$month[ble_data$year == i], ble_data$t[ble_data$year == i], col = 'gray42',  lwd = 2,  type = 'b',
       main = paste(' ', i), ylim = c(0,2250), cex.main = 1.5,
       xlab = "Month", ylab = "Landing [t]", cex.lab = 1.5, cex.axis = 1.25, las = 1)
  lines(sol_yearToSee$month, (sol_yearToSee$catch_commercial.BF6 )*7587.1/1000, col = 'green' )
  #lines(sol_yearToSee$month, (sol_yearToSee$catch_undersized.BF1)*7587.1/1000, col = '#c994c7')
  #lines(sol_2$month, (sol_2$catch_undersized.BF1)*7587.1/1000, col = 'red')

}

legend("topleft", legend = c('Landing data', 'Sim result'),#, 'bycatch 2.2', 'bycatch 2'),
       border = NA, y.intersp = 0.85, cex = 1.2, fill = c('gray42', 'green'))#, '#c994c7','red') )


#### #Sim. constant mesh size 26 m from 2010 / Pred remains same ####

Imax_ik_vals.c26 = seq(0.17, 0.19, length = length(seq(inti_y,fin_y))) #same as baseline
L50_vals.c26 = c(rep(4.64, length(seq(inti_y, fin_y)) ) )
SR_vals.c26 = c(rep(0.97, length(seq(inti_y,fin_y))))

PF_vars.c26 = data.frame(year=years, Imax = Imax_ik_vals.c26, L50 = L50_vals.c26, SR= SR_vals.c26)


bothPF_params.c26 = parameters_solv.v2
bothPF_params.c26$general_params$Fi = 4/365 #Temming & Hufnagl 2015

tempG.9_21 = temperature_germany %>% filter(date_time > as.Date("31/12/2009", format = "%d/%m/%Y") & date_time < as.Date("01/01/2022", format = "%d/%m/%Y"))
t_years = seq(0,length(tempG.9_21$temperature)-.1)


start <- Sys.time()
FP10_21.c26 <- solver_sizeClass.v6(t = t_years, state = Init.v3, parameters = bothPF_params.c26,
                               temperature_dataSet = tempG.9_21, PF_dataset = PF_vars.c26)
print(Sys.time() - start)

#### #Sim. constant mesh size 22 m from 2010 / pred remains same ####

Imax_ik_vals.c22 = seq(0.17, 0.19, length = length(seq(inti_y,fin_y))) #same as baseline
L50_vals.c22 = c(rep(4.02, length(seq(inti_y, fin_y)) ) )
SR_vals.c22 = c(rep(0.82, length(seq(inti_y,fin_y))))

PF_vars.c22 = data.frame(year=years, Imax = Imax_ik_vals.c22, L50 = L50_vals.c22, SR= SR_vals.c22)


#bothPF_params.c22 = parameters_solv.v2
#bothPF_params.c22$general_params$Fi = 4/365 #Temming & Hufnagl 2015

start <- Sys.time()
FP10_21.c22 <- solver_sizeClass.v6(t = t_years, state = Init.v3, parameters = bothPF_params.c26,
                                   temperature_dataSet = tempG.9_21, PF_dataset = PF_vars.c22)
print(Sys.time() - start)


#### #Sim. constant mesh size 22 m from 2010 / predation as of levels around '90 2/3 higher ####


Imax_ik_vals.ip = seq(0.17*5/3, 0.19*5/3, length = length(seq(inti_y,fin_y))) #same as baseline
L50_vals.ip = c(rep(4.02, length(seq(inti_y, fin_y)) ) )
SR_vals.ip = c(rep(0.82, length(seq(inti_y,fin_y))))

PF_vars.ip = data.frame(year=years, Imax = Imax_ik_vals.ip, L50 = L50_vals.ip, SR= SR_vals.ip)


#bothPF_params.c22 = parameters_solv.v2
#bothPF_params.c22$general_params$Fi = 4/365 #Temming & Hufnagl 2015

#tempG.9_21 = temperature_germany %>% filter(date_time > as.Date("31/12/2009", format = "%d/%m/%Y") & date_time < as.Date("01/01/2022", format = "%d/%m/%Y"))
#t_years = seq(0,length(tempG.9_21$temperature)-.1)


start <- Sys.time()
FP10_21.ip <- solver_sizeClass.v6(t = t_years, state = Init.v3, parameters = bothPF_params.c26,
                                   temperature_dataSet = tempG.9_21, PF_dataset = PF_vars.ip)
print(Sys.time() - start)


#### Prepare all (4) rsults to print the avg size
#(1) baseline
baseline_noSex = prep_sol_Lavg(FP10_21, "31/12/2009", "01/01/2022")

L_avg_baseline = baseline_noSex %>%
  mutate(L_avg = ( L1*size_mean_F[1] + L2*size_mean_F[2]+ L3*size_mean_F[3]+ L4*size_mean_F[4] + L5*size_mean_F[5]
                   + BF6*size_mean_F[6] + BF7*size_mean_F[7]+ BF8*size_mean_F[8])/ttlB )

#(2) constant mesh size 26 mm
c26_noSex = prep_sol_Lavg(FP10_21.c26, "31/12/2009", "01/01/2022")

L_avg_c26 = c26_noSex %>%
  mutate(L_avg = ( L1*size_mean_F[1] + L2*size_mean_F[2]+ L3*size_mean_F[3]+ L4*size_mean_F[4] + L5*size_mean_F[5]
                   + BF6*size_mean_F[6] + BF7*size_mean_F[7]+ BF8*size_mean_F[8])/ttlB )


#(3) constant mesh size 22 mm
c22_noSex = prep_sol_Lavg(FP10_21.c22, "31/12/2009", "01/01/2022")

L_avg_c22 = c22_noSex %>%
  mutate(L_avg = ( L1*size_mean_F[1] + L2*size_mean_F[2]+ L3*size_mean_F[3]+ L4*size_mean_F[4] + L5*size_mean_F[5]
                   + BF6*size_mean_F[6] + BF7*size_mean_F[7]+ BF8*size_mean_F[8])/ttlB )

#(2) constant mesh size 22 mm & increased predation
ip_noSex = prep_sol_Lavg(FP10_21.ip, "31/12/2009", "01/01/2022")

L_avg_ip = ip_noSex %>%
  mutate(L_avg = ( L1*size_mean_F[1] + L2*size_mean_F[2]+ L3*size_mean_F[3]+ L4*size_mean_F[4] + L5*size_mean_F[5]
                   + BF6*size_mean_F[6] + BF7*size_mean_F[7]+ BF8*size_mean_F[8])/ttlB )


cols = brewer.pal(3,"Dark2")

layout(matrix(c(1,2), nrow = 1), widths = c(2,1))

plot(L_avg_baseline$dateTime, L_avg_baseline$L_avg, las = 1,  col = 'black', cex.axis = 1.4, cex.lab = 1.4,
     xlab = 'years', ylab = paste(TeX("$L$"), '[cm]'), type = 'l', lty = 2, lwd = 3, #TeX(r"($\ell$)" )
     main = ' ', cex.main = 1.4) #main = "main = 'Changes in average size'"
lines(L_avg_c26$dateTime, L_avg_c26$L_avg, col = cols[1], lty=1, lwd = 2.5)
lines(L_avg_c22$dateTime, L_avg_c22$L_avg, col = cols[2], lty=1, lwd = 2.5)
lines(L_avg_ip$dateTime, L_avg_ip$L_avg, col = cols[3], lty=1, lwd = 2.5)
plot.new()
legend("topright", legend = c('baseline', '2.6 cm', '2.2 cm', '2.2 & ip'), xjust = 1,
       yjust = 0, fill = c('grey', cols[1], cols[2],cols[3]), border = NA, y.intersp = 0.75,
       bty = "n", xpd = TRUE, cex = 1.35)



cols = brewer.pal(3,"Dark2")
layout(matrix(c(1,2), nrow = 1), widths = c(2.3,1))
par(mar = c(4.5, 5, 3, 1), oma = c(0, 0, 0, 0))

plot(L_avg_baseline$dateTime, L_avg_baseline$ttlB, las = 1,  col = 'black', cex.axis = 1.4, cex.lab = 1.4,
     xlab = 'years', ylab = TeX("B [$kg km^{-2}$]"),type = 'l', lty = 2, lwd = 3, #, ylim = c(3.15, 5.1)
     main = ' ', cex.main = 1.4) #main = "main = 'Changes in average size'"
lines(L_avg_c26$dateTime, L_avg_c26$ttlB, col = cols[1], lty=1, lwd = 2.5)
lines(L_avg_c22$dateTime, L_avg_c22$ttlB, col = cols[2], lty=1, lwd = 2.5)
lines(L_avg_ip$dateTime, L_avg_ip$ttlB, col = cols[3], lty=1, lwd = 2.5)
plot.new()
par(mar = c(0,0,0,0))
legend("topright", legend = c('baseline', '2.6 cm', '2.2 cm', '2.2 & ip'), xjust = 1,
       yjust = 0, fill = c('grey', cols[1], cols[2],cols[3]), border = NA, y.intersp = 0.75,
       bty = "n", xpd = TRUE, cex = 1.35)#inset = c(-0.05, -0.05)




#### Size spectrum sheldon

# Suppose you have:

body_mass = convertL_to_W(sizes)

yearToSee = 2017

prepare_data = prep_sizeSpectrum(FP10_21, 2017)

biomassQ1 = prepare_data["Q1", ]
biomassQ3 = prepare_data["Q3", ]
biomassQ4 = prepare_data["Q4", ]

# Compute biomass *per log bin width* (Sheldon normalization)
# If bins are logarithmically spaced, no correction is needed.
# If bins are linear (like your case), you should divide by Δlog(m):
delta_log_m <- diff(log10(body_mass))
# For simplicity, assign each bin's width as the next-minus-previous
# Use average for interior bins
delta_log_m <- c(delta_log_m[1], delta_log_m, delta_log_m[length(delta_log_m)])
delta_log_m <- delta_log_m[1:length(biomassQ1)]

# Normalize biomass by bin width
B_per_log_bin <- biomassQ1 / delta_log_m

# Fit Sheldon spectrum
sheldon_model <- lm(log10(B_per_log_bin) ~ log10(body_mass))

summary(sheldon_model)

# Plot
dev.off()
plot(log10(body_mass), log10(B_per_log_bin),
     pch = 16, xlab = "log10(Body mass)", ylab = "log10(Biomass per log bin)", main = 'Q1')
abline(sheldon_model, col = "red", lwd = 2)



#Q3
B_per_log_bin <- biomassQ3 / delta_log_m

# Fit Sheldon spectrum
sheldon_model <- lm(log10(B_per_log_bin) ~ log10(body_mass))

summary(sheldon_model)

# Plot

plot(log10(body_mass), log10(B_per_log_bin),
     pch = 16, xlab = "log10(Body mass)", ylab = "log10(Biomass per log bin)", main = 'Q3')
abline(sheldon_model, col = "red", lwd = 2)

#Q4
B_per_log_bin <- biomassQ4 / delta_log_m

# Fit Sheldon spectrum
sheldon_model <- lm(log10(B_per_log_bin) ~ log10(body_mass))

summary(sheldon_model)

# Plot

plot(log10(body_mass), log10(B_per_log_bin),
     pch = 16, xlab = "log10(Body mass)", ylab = "log10(Biomass per log bin)", main = 'Q4')
abline(sheldon_model, col = "red", lwd = 2)



###teste plots for thesis :

#plot mu_fi against mu_pred
L_divers = c(2,4,6)
days = seq(1,360)

fi_f = as.vector(sapply(monthly_factors.v, function(x) rep(x, 30)))

plot(days, sel_probit(L = L_divers[1], L50 = 4.02, SR =  0.82)*fi_f, type = 'l')

plot(days, ST_predation(days, 10, 100, eta = eta_J(days), beta = 12 )*ingestion_kernel(I_max= 0.175, l_pred = 15, l_prey = L_divers[1]) )


##See how different size classes develop:
cols = brewer.pal(8,"Dark2")
plot()






