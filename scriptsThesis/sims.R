#### *NO NEED TO RUN THIS SECTION* Read original file Temp German bight & create new format #####
#The above part has to be run only once ! done on 26.09.25

start <- Sys.time()
temperature_dataSet_II <- read.csv("./data/bhv_alte_weser_leuchtturm_red.txt", header = TRUE, sep = "\t",, skip = 16, dec = '.')
print(Sys.time() - start) #4 s

colnames(temperature_dataSet_II) = c("date_time", "temperature")

start <- Sys.time()
temperature_dataSet_II$DateTime <- as.POSIXct(temperature_dataSet_II$date_time, format = "%Y-%m-%d %H:%M:%S", tz = "CET")
print(Sys.time() - start) #26 s

temperature_dataSet_II$date_time <- as.Date(temperature_dataSet_II$date_time)
temperature_dataSet_II$Time <- format(temperature_dataSet_II$date_time, format = "%H:%M:%S")

temp_germ_Davg =  temperature_dataSet_II %>% #T° as day avg
  group_by(date_time) %>% summarise(temperature = mean(temperature))

#Since file is too heavy, Filter only years from 2010
temp_10G = temp_germ_Davg %>%
  filter( date_time > as.Date("31/12/2009", format = "%d/%m/%Y")) #& date_time < as.Date("01/01/2016", format = "%d/%m/%Y"))
#save it in a new csv. file:
setwd("./data")
write.csv(temp_10G, "temperature_10ger.csv")

##### Libraries and UPLOAD TEMPERATURE DATA #####

library(deSolve)
library(RColorBrewer)
library(BrownShrimp)
library(viridis)
library(dplyr)
library(tidyr)
library("colorspace")
#library(lubridate)
library(latex2exp)
library(profvis)
library(stringr)

#library(devtools)

#Read water temperature from Wadden Sea:
temperature_dataSet <- read.csv("./data/Temperature_data_1982-2018_MarBiol.csv", header = TRUE, sep = ';')
temperature_10_15 <- temperature_func(temperature_dataSet, "01/01/2010", "31/12/2015")


##########

temperature_germany <- read.csv("./data/temperature_10ger.csv")

tempG.10_15 = temperature_germany %>% filter(date_time > as.Date("31/12/2009", format = "%d/%m/%Y") & date_time < as.Date("01/01/2016", format = "%d/%m/%Y"))

ble_data = read.csv("./data/BLE_Inlandslandungen_SpeiseKrabbe.csv", sep = ';')

t_years = seq(0,length(tempG.10_15$temperature)-.1)

#initialize Initial conditions

BF_p = c(3.801852,  7.063263,  12.657331, 25.166351, 49.297133, 45.876772, 34.672340,  19.396426)#*1e-6 #see calculation in file 'sizeSpectraMoldel.R'
BM_p = c(4.573332  , 8.203546 , 13.791605 , 26.328928 , 43.371121  )#*1e-6 #see calculation in file 'sizeSpectraMoldel.R'
IC_parameterized.v2 = c(P = 2, E = 2, L= 1 ,
                        BF = BF_p, BM = BM_p) #Initial conditions with reduced area (only up to 20 m depth).


#Initialitze parameters for comparation of 4 diff scenarios: no pressure, fishery, pred and both

#No pressure
noPreasure_params = parameters_solv
noPreasure_params$general_params$Fi = 0 #No fishery
noPreasure_params$general_params$Imax_ik = 0 #No predation

#Only predation
pred_params = parameters_solv
pred_params$general_params$Fi = 0 #No fishery
pred_params$general_params$Imax_ik = 0.185


#Only Fishery

fishery_params = parameters_solv
fishery_params$general_params$Imax_ik = 0 #No predation
fishery_params$general_params$Fi = 4/365



#### (1) THESIS VALIDATION: abundance (vs Fishery) .V5 ####

##### YEARS BEFORE 2016: 20mm cod-end #####
#Simulation with temperature measured in Germany

bothPF_params = parameters_solv
bothPF_params$general_params$Fi = 4/365 #(1-exp(-5))/365
bothPF_params$general_params$L50 = 3.69 #et. al Santos 2018
bothPF_params$general_params$SR = 0.75 #et. al Santos 2018
bothPF_params$general_params$Imax_ik = 0.16 #(1-exp(-2))/365 #reduce predation in these years as per: https://www.openseas.org.uk/news/on-thin-ices-north-sea-cod/

BF = c(8.695212, 15.349638, 21.461595, 24.672213,25.871618, 17.536662, 14.090960, 11.933151)
BM = c(10.85817, 19.93554, 29.01086, 36.67041, 58.26882)

Init.v3 = c(P = 2, E = 2, L= 1 ,
            BF = BF, BM = BM)

start <- Sys.time()
FP_10_15_g <- solver_sizeClass.v5(t = t_years, state = Init.v3, parameters = bothPF_params, temperature_dataSet = tempG.10_15)
print(Sys.time() - start)

plot_sizeSpectra(sol_df_main = FP_10_15_g, year = 2014, title = "F&P")

plot(FP_10_15_g$dateTime, FP_10_15_g$P)

years_to_see = c(2012, 2013, 2014)

par(mfrow = c(1,length(years_to_see)))
for (i in years_to_see){
  sol_yearToSee = FP_10_15_g %>%
    filter(as.numeric(format(FP_10_15_g$dateTime,'%Y')) == i) %>%
    mutate(month = month(dateTime)) %>%
    group_by(month) %>%
    summarise(across(where(is.numeric), \(x) sum(x, na.rm = TRUE)))

  #Calculate the total landings per month (considering 50% discard):
  plot(ble_data$month[ble_data$year == i], ble_data$t[ble_data$year == i], col = 'gray42',  lwd = 2,  type = 'b',
       main = paste(' ', i), ylim = c(0,2350), cex.main = 1.5,
       xlab = "Month", ylab = "Landing [tones]", cex.lab = 1.5, cex.axis = 1.25)
  lines(sol_yearToSee$month, (sol_yearToSee$catch_commercial.BF5 )*7587.1/1000, col = 'green' )
  lines(sol_yearToSee$month, (sol_yearToSee$catch_undersized.BF1)*7587.1/1000, col = '#c994c7')
}


colors_years = c('#c51b8a','#7a0177', '#7bccc4','#41b6c4')
par(mfrow = c(1,1))

for (i in 1:length(years_to_see))  {
  sim_result_month = FP_10_15_g %>%
    filter(year(FP_10_15_g$dateTime) == years_to_see[i]) %>%  #extract data from the year you want to see
    mutate(month = month(dateTime)) %>%
    group_by(month) %>%
    summarise(across(where(is.numeric), \(x) sum(x, na.rm = TRUE)))

  if (i == 1 ) {
    plot(ble_data$t[ble_data$year == years_to_see[i]], ( sim_result_month$catch_commercial.BF5)*7587.1/1000, col = colors_years[i],
         xlim = c(0,2500), ylim = c(0,2500), xlab = 'Landing - data [t]', ylab = 'Landing - simulation [t] ' ,
         lwd = 2, cex.lab = 1.3, cex.axis = 1.1)
    lines(seq(1,2500, 0.1), seq(1,2500, 0.1) )
    } else {
      points(ble_data$t[ble_data$year == years_to_see[i]], ( sim_result_month$catch_commercial.BF5)*7587.1/1000,
             col =  colors_years[i], lwd = 2)
    }


}
legend("topright", legend = years_to_see, fill = colors_years, cex = 1 , bty = 'n')



#### THESIS RESULT/DISCUSSION:  .V5 ####

#CHANGES IN AVG SIZE

#### SIM NO PREASURE ####

start <- Sys.time()
noPreasure <- solver_sizeClass.v5(t = t_years, state = IC_parameterized.v2, parameters = noPreasure_params, temperature_dataSet = tempG.10_15)
print(Sys.time() - start)

#Plot size spectra
plot_sizeSpectra(noPreasure, year = 2013, title = "NP")


#### SIM ONLY PREDATION ####

start <- Sys.time()
predation <- solver_sizeClass.v5(t = t_years, state = IC_parameterized.v2, parameters = pred_params, temperature_dataSet = tempG.10_15)
print(Sys.time() - start)

#Plot size spectra
plot_sizeSpectra(predation, year = 2015, title = "Pred")

#### SIM (3 out of 4) ONLY FISHERY ####

start <- Sys.time()
fishery <- solver_sizeClass.v5(t = t_years, state = IC_parameterized.v2, parameters = fishery_params, temperature_dataSet = tempG.10_15)
print(Sys.time() - start)

#Plot size spectra
plot_sizeSpectra(fishery, year = 2015, title = "Fi")


#### THESIS VALIDATION: Size structure .V6 (idea not finished yet) ####


start <- Sys.time()
FP_10_15_testv6 <- solver_sizeClass.v6(t = t_years, state = IC_parameterized.v2, parameters = bothPF_params, temperature_dataSet = temperature_10_15)
print(Sys.time() - start)

plot_sizeSpectra_fishery(sol_df_main = FP_10_15_testv6, year = 2015, title = "F&P")



#Calculate Nr of individuals:

weights_f = convertL_to_W(c(LJ, LJ2,  LJ3, LJ4, LJ5,LA1, LA2, LA3))
weights_m = convertL_to_W(c(LJ, LJ2,  LJ3, LJ4, LJ5))

names(weights_f) <- paste0("BF", seq_along(weights_f))   # BF1, BF2, ...
names(weights_m) <- paste0("BM", seq_along(weights_m))   # BF1, BF2, ...


par(mfrow = c(1,1))

L = seq(1, 8, by = 0.1)
plot(L, sel_probit(L, L50=4.02, SR = 0.8))


#### THESIS VALIDATION: Size structure .V5 ####

#Here we want to plot the size structure of all idividuals living between sub- and inter tidal waters

#Prepare for plot:
year = 2019
past_year <- paste0("31/12/", as.character(year - 1))
next_year <- paste0("01/01/", as.character(year + 1))

#select only relevant size classes
sol_df_year <- FP_15_18_g.22 %>%#FP_10_15_g %>%
  filter(as.Date(dateTime) > as.Date(past_year, format = "%d/%m/%Y") &
           as.Date(dateTime) < as.Date(next_year, format = "%d/%m/%Y")) %>%
  select(BF1:BF8, BM1:BM5)

cols <- c("f" = "#8856a7", "m" = "#add8e6")
#label_map <- c("BL1","BL2","BL3","BL4","BL5","BL6","BL7","BL8")
label_map  = as.character(sizes[2:9]*10)

weights = convertL_to_W(sizes)

months = c("Jan", "Feb", "Mar", "Apr")

par(mfrow = c(2, 2))
out <- list()
for (q in 1:4) {
  init_idx <-  30*(q-1)+1 #91.2 * (q - 1) + 1
  final_idx <- 30*q #91.2 * q
  #quarter_data <- sol_df_year[init_idx:final_idx, ]
  monthly_data <- sol_df_year[init_idx:final_idx, ]

  f_means <- colMeans(monthly_data[paste0("BF", 1:8)] / weights[1:8])#colMeans(monthly_data[paste0("BF", 2:6)])
  m_means <- colMeans(monthly_data[paste0("BM", 1:5)] / weights[1:5])#colMeans(monthly_data[paste0("BM", 2:5)])

  f_stack <- f_means[paste0("BF", 1:5)]
  m_stack <- m_means
  f_extra <- f_means[paste0("BF", 6:8)]

  f_vals <- c(f_stack, f_extra)
  m_vals <- c(m_stack, rep(0, length(f_extra)))

  # Calculate relative frequencies
  total_f <- sum(f_vals)
  total_m <- sum(m_vals)

  f_vals_rel <- f_vals / total_f
  m_vals_rel <- m_vals / total_m

  out[[q]] <- rbind("f" = f_vals_rel, "m" = m_vals_rel)

  mat <- out[[q]]
  bar_centers = barplot(
    mat,
    col = cols[rownames(mat)],
    main = months[q],#paste( "Q", q, "–", year),
    names.arg = label_map,
    las = 1,
    cex.main = 1.5,
    cex.axis = 1.2,
    cex.names = 1.2,
    ylim = c(0, max(colSums(mat))),
    space = 0,
    beside = FALSE,
    legend.text = rownames(mat),
    xlab = "size class",
    args.legend = list(x = "topright", bty = "n", inset = 0.02, cex = 1)
  )
}




####### PLOTS PRESENTATION ######

#PREDATION
#Function mortality due to predation

days = seq(1, 365, 0.1)
B_s = 100
eta = eta_J(days)

plot(days, ST_predation(days, Te = 10, B_s, eta, beta = 5, sigma = 0.8),
     ylab = TeX("$m_{s}$"), col = 'gray41', cex.lab = 1.5, cex.axis = 1.1)
abline(h= mean(ST_predation(days, Te = 10, B_s, eta, beta = 5, sigma = 0.8)), lty = 2)


#Size selective predation
bs_size_range = seq(0.6, 8, 0.01)

plot(bs_size_range, ingestion_kernel(I_max=0.5, l_pred = 15, bs_size_range),
    col = 'gray41', cex.lab = 1.5, cex.axis = 1.1,
     xlab = "BS size range [cm]", ylab = TeX("$SSP$"))



#FISHERY

plot(bs_size_range, sel_probit(bs_size_range, L50=4.02, SR=0.82),
     ylab = TeX("$p_{fi}$"), xlab = "BS size range [cm]",  col = 'gray41', cex.lab = 1.5, cex.axis = 1.1)


#GROWTH FUNCTION

par(mfrow = c(1,1))
plot(bs_size_range, shift_next_sizeClass(bs_size_range, 10, parameters_solv$Fem_params),
     xlab = "size range BS", ylab = "growth rate", col = 'gray41')


#### (2) THESIS DISCUSSION:  #####

#### (2.1) CHANGES IN AVG SIZE Preparation for plots (years before 2016) ####
#(1) L_avg over the time:
size_classes = paste0("BL", 1:8)

#calculate weighted average for each time step:

#first join F and M in one column
#Without preasure (1)

#noPreasure_noSex = noPreasure %>%
#  mutate(L1 = BF1 + BM1, L2 = BF2 + BM2, L3 = BF3 + BM3, L4 = BF4 + BM4,  L5 = BF5 + BM5) %>%
#  mutate(ttlB = L1 + L2+ L3+ L4+ L5+ BF6+ BF7 + BF8) %>%
#  select(21:25, 10:12,26, 20 ) %>%
#  filter(as.Date(dateTime) > as.Date("31/12/2010", format= "%d/%m/%Y" ) & as.Date(dateTime) < as.Date("01/01/2015", format= "%d/%m/%Y" )) #delete first year 'warm up' period
noPreasure_noSex = prep_sol_Lavg(noPreasure, "31/12/2010", "01/01/2015")

L_avg_noPreasure = noPreasure_noSex %>%
  mutate(L_avg = ( L1*size_mean_F[1] + L2*size_mean_F[2]+ L3*size_mean_F[3]+ L4*size_mean_F[4] + L5*size_mean_F[5]
                   + BF6*size_mean_F[6] + BF7*size_mean_F[7]+ BF8*size_mean_F[8])/ttlB )

#Only Predation (2)
predation_noSex = predation %>%
  mutate(L1 = BF1 + BM1, L2 = BF2 + BM2, L3 = BF3 + BM3, L4 = BF4 + BM4,  L5 = BF5 + BM5) %>%
  mutate(ttlB = L1 + L2+ L3+ L4+ L5+ BF6+ BF7 + BF8) %>%
  select(21:25, 10:12,26, 20 ) %>% #Column 24 HAS TO EB CHANGED WHEN RUNNING NEW SIM TO SELECT DATETIME
  filter(as.Date(dateTime) > as.Date("31/12/2010", format= "%d/%m/%Y" ) & as.Date(dateTime) < as.Date("01/01/2015", format= "%d/%m/%Y" )) #delete first year 'warm up' period

Lavg_pred = predation_noSex %>%
  mutate(L_avg = ( L1*size_mean_F[1] + L2*size_mean_F[2]+ L3*size_mean_F[3]+ L4*size_mean_F[4] + L5*size_mean_F[5]
                   + BF6*size_mean_F[6] + BF7*size_mean_F[7]+ BF8*size_mean_F[8])/ttlB )
#Only Fishery (3)
fishery_noSex = fishery %>%
  mutate(L1 = BF1 + BM1, L2 = BF2 + BM2, L3 = BF3 + BM3, L4 = BF4 + BM4,  L5 = BF5 + BM5) %>%
  mutate(ttlB = L1 + L2+ L3+ L4+ L5+ BF6+ BF7 + BF8) %>%
  select(21:25, 10:12,26, 20 ) %>% #Column 24 HAS TO EB CHANGED WHEN RUNNING NEW SIM TO SELECT DATETIME
  filter(as.Date(dateTime) > as.Date("31/12/2010", format= "%d/%m/%Y" ) & as.Date(dateTime) < as.Date("01/01/2015", format= "%d/%m/%Y" )) #delete first year 'warm up' period

Lavg_fishery = fishery_noSex %>%
  mutate(L_avg = ( L1*size_mean_F[1] + L2*size_mean_F[2]+ L3*size_mean_F[3]+ L4*size_mean_F[4] + L5*size_mean_F[5]
                   + BF6*size_mean_F[6] + BF7*size_mean_F[7]+ BF8*size_mean_F[8])/ttlB )

#Both predation & fishery (4)
bothPF_noSex = FP_10_15_g %>%
  mutate(L1 = BF1 + BM1, L2 = BF2 + BM2, L3 = BF3 + BM3, L4 = BF4 + BM4,  L5 = BF5 + BM5) %>%
  mutate(ttlB = L1 + L2+ L3+ L4+ L5+ BF6+ BF7 + BF8) %>%
  select(21:25, 10:12,26, 20 ) %>% #Column 24 HAS TO EB CHANGED WHEN RUNNING NEW SIM TO SELECT DATETIME
  filter(as.Date(dateTime) > as.Date("31/12/2010", format= "%d/%m/%Y" ) & as.Date(dateTime) < as.Date("01/01/2015", format= "%d/%m/%Y" )) #delete first year 'warm up' period

Lavg_bothPF = bothPF_noSex %>%
  mutate(L_avg = ( L1*size_mean_F[1] + L2*size_mean_F[2]+ L3*size_mean_F[3]+ L4*size_mean_F[4] + L5*size_mean_F[5]
                   + BF6*size_mean_F[6] + BF7*size_mean_F[7]+ BF8*size_mean_F[8])/ttlB )


par(mfrow = c(1,1))
plot(noPreasure_noSex$dateTime, L_avg_noPreasure$L_avg, las = 1,  col = 'grey', cex = 0.8,
     xlab = 'years', ylab = paste(TeX("$L$"), '[cm]'), ylim = c(2, 5.5) ,type = 'l', lty = 2, lwd = 3, #TeX(r"($\ell$)" )
     main = 'Changes in average size ')
lines(Lavg_pred$dateTime, Lavg_pred$L_avg, col = 'indianred3', lty=1, lwd = 2.5)
lines(Lavg_fishery$dateTime, Lavg_fishery$L_avg, col = 'skyblue4', lty=1, lwd = 2.5)
lines(Lavg_bothPF$dateTime, Lavg_bothPF$L_avg, col = 'hotpink', lty=1, lwd = 2.5)

legend("topright", legend = c('no preasure', 'only predation', 'only fishery', 'both F & P'),
       fill = c('grey', 'indianred3', 'skyblue4','hotpink'), border = NA, bty = "n", y.intersp = 0.65)



#### (2.2) MEASURING EFFECT OF CHANGING MESH SIZE OF COD-END ####

##### YEARS AFTER 2016: 22mm cod-end#####
#Simulation with temperature measured in Germany
#test first with COD-END 22mm

BF = c(8.695212, 15.349638, 21.461595, 24.672213,25.871618, 17.536662, 14.090960, 11.933151)
BM = c(10.85817, 19.93554, 29.01086, 36.67041, 58.26882)

Init.v3 = c(P = 2, E = 2, L= 1 ,
            BF = BF, BM = BM)

bothPF_params16 = parameters_solv
bothPF_params16$general_params$Fi = 4/365 #Temming & Hufnagl 2015
bothPF_params16$general_params$L50 = 4.02 #et. al Santos 2018
bothPF_params16$general_params$SR = 0.82 #et. al Santos 2018
bothPF_params16$general_params$Imax_ik = 0.18#'increasing' predation in these years as per: https://www.openseas.org.uk/news/on-thin-ices-north-sea-cod/

tempG.15_18 = temperature_germany %>% filter(date_time > as.Date("31/12/2014", format = "%d/%m/%Y") & date_time < as.Date("01/01/2020", format = "%d/%m/%Y"))
t_years16 = seq(0,length(tempG.15_18$temperature)-.1)

start <- Sys.time()
FP_15_18_g.22 <- solver_sizeClass.v5(t = t_years16, state = Init.v3, parameters = bothPF_params16, temperature_dataSet = tempG.15_18)
print(Sys.time() - start)

plot_sizeSpectra_old(sol_df_main = FP_15_18_g.22, year = 2017, title = "F&P")
dev.off()
plot(FP_15_18_g.22$dateTime, FP_15_18_g.22$P )

#test now with COD-END 20mm

bothPF_params16.2 = parameters_solv
bothPF_params16.2$general_params$Fi = 4/365 #Temming & Hufnagl 2015
bothPF_params16.2$general_params$L50 = 3.69 #et. al Santos 2018
bothPF_params16.2$general_params$SR = 0.75 #et. al Santos 2018
bothPF_params16.2$general_params$Imax_ik = 0.185 #'increasing' predation in these years as per: https://www.openseas.org.uk/news/on-thin-ices-north-sea-cod/

tempG.15_18 = temperature_germany %>% filter(date_time > as.Date("31/12/2014", format = "%d/%m/%Y") & date_time < as.Date("01/01/2020", format = "%d/%m/%Y"))
t_years16 = seq(0,length(tempG.15_18$temperature)-.1)


start <- Sys.time()
FP_15_18_g.20 <- solver_sizeClass.v5(t = t_years16, state = IC_parameterized.v2, parameters = bothPF_params16.2, temperature_dataSet = tempG.15_18)
print(Sys.time() - start)



years_to_see = c(2017, 2018, 2019)

par(mfrow = c(1,length(years_to_see)))
for (i in years_to_see){
  #Solution with cod-end 22 mm
  sol_yearToSee = FP_15_18_g.22 %>%
    filter(as.numeric(format(FP_15_18_g.22$dateTime,'%Y')) == i) %>%
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


#Plot validation

colors_years =  brewer.pal(4,"Dark2")#c('#c51b8a','#7a0177', '#7bccc4','#41b6c4')
par(mfrow = c(1,1))

for (i in 1:length(years_to_see))  {
  sim_result_month = FP_15_18_g.22 %>%
    filter(year(FP_15_18_g.22$dateTime) == years_to_see[i]) %>%  #extract data from the year you want to see
    mutate(month = month(dateTime)) %>%
    group_by(month) %>%
    summarise(across(where(is.numeric), \(x) sum(x, na.rm = TRUE)))

  if (i == 1 ) {
    plot(ble_data$t[ble_data$year == years_to_see[i]], ( sim_result_month$catch_commercial.BF6)*7587.1/1000, col = colors_years[i],
         xlim = c(0,2100), ylim = c(0,2100), xlab = 'Landing - data [t]', ylab = 'Landing - simulation [t] ' ,
         lwd = 2, cex.lab = 1.2, cex.axis = 1.2, las = 1)
    lines(seq(1,2500, 0.1), seq(1,2500, 0.1) )
  } else {
    points(ble_data$t[ble_data$year == years_to_see[i]], ( sim_result_month$catch_commercial.BF6)*7587.1/1000,
           col =  colors_years[i], lwd = 2)
  }


}
legend("topleft", legend = years_to_see, fill = colors_years, cex = 1.2 )#, bty = 'n'




#Simulate different scenarios (no fishery, only F, only P for 2016-2018)

noPreasure16 <- solver_sizeClass.v5(t = t_years16, state = IC_parameterized.v2,
                                    parameters = noPreasure_params, temperature_dataSet = tempG.15_18)

#plot_sizeSpectra_old(sol_df_main = noPreasure16, year = 2017, title = "no F&P")

predation16 <- solver_sizeClass.v5(t = t_years16, state = IC_parameterized.v2,
                                 parameters = pred_params, temperature_dataSet = tempG.15_18)

plot_sizeSpectra_old(sol_df_main = predation16, year = 2017, title = "P")

#layout(matrix(c(1,2), nrow = 1), widths = c(2,1))

plot_sizeSpectra(predation16 ,  sol_df_overlay= FP_15_18_g.20, year = 2017, title = " ", #"F&P vs P"
                 legend = c("") )
plot.new()
legend("topright", legend = c("P", "F & P", "intersect."), xjust = 1,
       yjust = 0, fill = c(adjustcolor("#43a2ca", alpha.f = 0.8), adjustcolor("gold1", alpha.f = 0.65), "#B9C86B"), border = NA, y.intersp = 0.5,
       bty = "n", xpd = TRUE, inset = c(-0.05, -0.05), cex = 1.5)

plot_sizeSpectra_freq(predation16 ,  sol_df_overlay= FP_15_18_g.20, year = 2017, title = " ", #"F&P vs P"
                 legend = c("") )

#parameters for net as of 2017
fishery_params$general_params$L50 = 4.02
fishery_params$general_params$SR = 0.82

fishery16 <- solver_sizeClass.v5(t = t_years16, state = IC_parameterized.v2,
                                 parameters = fishery_params, temperature_dataSet = tempG.15_18)

plot_sizeSpectra_old(sol_df_main = fishery16, year = 2017, title = "F")


plot_sizeSpectra(FP_15_18_g.22 , sol_df_overlay=  FP_15_18_g.20, year = 2017, title = "cod-end 2 vs 2.2 cm"
                 , legend = c("2 cm", "2.2 cm"))

plot_sizeSpectra(fishery16, sol_df_overlay = predation16, year = 2017, title = "F vs P")

plot_sizeSpectra_old(sol_df_main = FP_15_18_g.22, year = 2017, title = "2.2 cm, mesh size")



#plot_sizeSpectra(FP_15_18_g.22, sol_df_overlay= predation16 , year = 2017, title = " ", legend = c("P & F", "P"))

#plot_sizeSpectra(FP_15_18_g.22, sol_df_overlay= noPreasure16 , year = 2017, title = " ", legend = c("P & F", "noPress"))


#test now with COD-END 24mm

bothPF_params16.24 = parameters_solv
bothPF_params16.24$general_params$Fi = 4/365 #Temming & Hufnagl 2015
bothPF_params16.24$general_params$L50 = 4.33 #et. al Santos 2018
bothPF_params16.24$general_params$SR = 0.9 #et. al Santos 2018
bothPF_params16.24$general_params$Imax_ik = 0.15 #'increasing' predation in these years as per: https://www.openseas.org.uk/news/on-thin-ices-north-sea-cod/

tempG.15_18 = temperature_germany %>% filter(date_time > as.Date("31/12/2014", format = "%d/%m/%Y") & date_time < as.Date("01/01/2020", format = "%d/%m/%Y"))
t_years16 = seq(0,length(tempG.15_18$temperature)-.1)


start <- Sys.time()
FP_15_18_g.24 <- solver_sizeClass.v5(t = t_years16, state = IC_parameterized.v2, parameters = bothPF_params16.24, temperature_dataSet = tempG.15_18)
print(Sys.time() - start)

plot_sizeSpectra(FP_15_18_g.24 , sol_df_overlay=  FP_15_18_g.20, year = 2017, title = "cod-end 2 vs 2.2 cm"
                 , legend = c("2 cm", "2.4 cm"))



fishery_params$general_params$L50 = 1.5


fishery16 <- solver_sizeClass.v5(t = t_years16, state = IC_parameterized.v2,
                                 parameters = fishery_params, temperature_dataSet = tempG.15_18)

plot_sizeSpectra(fishery16, sol_df_overlay= noPreasure16 , year = 2017, title = " ",
                 legend = c("F", "noPress"))


#Prepare data for the plot
#No preasure
noPreasure_noSex16 = prep_sol_Lavg(noPreasure16, "31/12/2015", "01/01/2019")

L_avg_noPreasure16 = noPreasure_noSex16 %>%
  mutate(L_avg = ( L1*size_mean_F[1] + L2*size_mean_F[2]+ L3*size_mean_F[3]+ L4*size_mean_F[4] + L5*size_mean_F[5]
                   + BF6*size_mean_F[6] + BF7*size_mean_F[7]+ BF8*size_mean_F[8])/ttlB )

#Only Predation
pred_nosex16 = prep_sol_Lavg(predation16, "31/12/2015", "01/01/2019")

L_avg_pred_nosex16 = pred_nosex16 %>%
  mutate(L_avg = ( L1*size_mean_F[1] + L2*size_mean_F[2]+ L3*size_mean_F[3]+ L4*size_mean_F[4] + L5*size_mean_F[5]
                   + BF6*size_mean_F[6] + BF7*size_mean_F[7]+ BF8*size_mean_F[8])/ttlB )

#Only fishery
fish_nosex16 = prep_sol_Lavg(fishery16, "31/12/2015", "01/01/2019")

L_avg_fish_nosex16 = fish_nosex16 %>%
  mutate(L_avg = ( L1*size_mean_F[1] + L2*size_mean_F[2]+ L3*size_mean_F[3]+ L4*size_mean_F[4] + L5*size_mean_F[5]
                   + BF6*size_mean_F[6] + BF7*size_mean_F[7]+ BF8*size_mean_F[8])/ttlB )

#Both F&P
both_nosex16 = prep_sol_Lavg(FP_15_18_g.22, "31/12/2015", "01/01/2019")

L_avg_both_nosex16 = both_nosex16 %>%
  mutate(L_avg = ( L1*size_mean_F[1] + L2*size_mean_F[2]+ L3*size_mean_F[3]+ L4*size_mean_F[4] + L5*size_mean_F[5]
                   + BF6*size_mean_F[6] + BF7*size_mean_F[7]+ BF8*size_mean_F[8])/ttlB )

#PLOT CHANGES IN AVG SIZE
#par(mfrow = c(2,1))
cols = brewer.pal(4,"Dark2")

layout(matrix(c(1,2), nrow = 1), widths = c(2,1))

plot(L_avg_noPreasure16$dateTime, L_avg_noPreasure16$L_avg, las = 1,  col = 'grey', cex.axis = 1.4, cex.lab = 1.4,
     xlab = 'years', ylab = paste(TeX("$L$"), '[cm]'), ylim = c(3.15, 5.1) ,type = 'l', lty = 2, lwd = 3, #TeX(r"($\ell$)" )
     main = ' ', cex.main = 1.4) #main = "main = 'Changes in average size'"
lines(L_avg_pred_nosex16$dateTime, L_avg_pred_nosex16$L_avg, col = cols[1], lty=1, lwd = 2.5)
lines(L_avg_fish_nosex16$dateTime, L_avg_fish_nosex16$L_avg, col = cols[2], lty=1, lwd = 2.5)
lines(L_avg_both_nosex16$dateTime, L_avg_both_nosex16$L_avg, col = cols[3], lty=1, lwd = 2.5)
plot.new()
legend("topright", legend = c('no pressure', 'only predation', 'only fishery', 'both F & P'), xjust = 1,
       yjust = 0, fill = c('grey', cols[1], cols[2],cols[3]), border = NA, y.intersp = 0.75,
       bty = "n", xpd = TRUE, inset = c(-0.05, -0.05), cex = 1.35)


#PLOT CHANGES IN BIOMASS

#plot_sizeSpectra(FP_15_18_g.20, sol_df_overlay= FP_15_18_g.22 ,year = 2017, title = "cod-end 2 vs 2.2 cm")


#par(mfrow = c(1,2) )
layout(matrix(c(1,2), nrow = 1), widths = c(2.3,1))
par(mar = c(4.5, 5, 3, 1), oma = c(0, 0, 0, 0))
plot(noPreasure_noSex16$dateTime, log(noPreasure_noSex16$ttlB), las = 1,  col = 'grey', cex.axis = 1.2, cex.lab = 1.2,
     xlab = 'years', ylab = TeX(" biomass kg $km^{-2}$"), type = 'l', lty = 2, lwd = 3,
     main = ' ', ylim = c(5.75, 10) ) # main = 'Changes in biomass'
lines(L_avg_pred_nosex16$dateTime, log(pred_nosex16$ttlB) , col = cols[1], lty=1, lwd = 2.5)
lines(L_avg_fish_nosex16$dateTime, log(fish_nosex16$ttlB), col = cols[2], lty=1, lwd = 2.5)
lines(L_avg_both_nosex16$dateTime, log(both_nosex16$ttlB), col = cols[3], lty=1, lwd = 2.5)
par(mar = c(0, 0, 0, 0))
plot.new()
legend("center", legend = c('no pressure', 'only predation', 'only fishery', 'both F & P'),
       fill = c('grey', cols[1], cols[2],cols[3]), border = NA, bty = "n", y.intersp = 0.65, cex = 1.1)#,
       #inset = c(0.0, 0.05)) #,xpd = TRUE








####to be deleted:

# example arrays you already have
monthly_factors.v  = monthly_factors.v    # length 12
rel_landing_by_year = ttl_yearly_landing$rel_landing   # vector length 12 or function returning relat. landing per month
#sel_probit_vec         # vector of selectivity per size for 12 months or constant

# build daily fm for a single size over 365 days (simple approach: repeat monthly factors)
days = 1:365
month_of_day = rep(1:12, each = 365/12) # or better: map actual dates to months

Fi_day = parameters_solv$general_params$Fi / 365

fm_day = sapply(days, function(d) {
  m = month_of_day[d]
  Fi_day * rel_landing_by_year[m] * monthly_factors.v[m] * sel_probit(5)
})

F_year_effective = sum(fm_day)          # instantaneous rate per year
annual_survival = exp(-F_year_effective)




days = seq(1,365)

plot(days, ST_predation(days, 10, 20000 ,eta = eta_J(days) ))


