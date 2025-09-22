#library(ggplot2)
library(deSolve)
library(RColorBrewer)
library(BrownShrimp)
library(viridis)
library(dplyr)
library(tidyr)
library("colorspace")
#library(lubridate)
library(latex2exp)
#library(profvis)


#read Ble data
ble_data = read.csv("./data/BLE_Inlandslandungen_SpeiseKrabbe.csv", sep = ';')

### Determine yearly factors#####
#This yearly factor, will change the intensity of fishery according landing data.
#the factor will multiply the baseline fishery intensity estimated to be 0.05

years_landings = ble_data$year
ttl_yearly_landing = ble_data %>%
  group_by(year) %>%
  summarise(ttl_year = sum(t)) %>%
  filter(year > 2009 & year < 2025) %>% #this line can be deleted when data from this previous years are updated
  mutate(rel_landing = ttl_year / mean(ttl_year))

plot(ttl_yearly_landing$year, ttl_yearly_landing$rel_landing, type = 'b')


#Idea is to compare before and after imposition of mesh size, howwever upto today 03.08.25
#data beore 2013 is too poor to give good impressions

years = c(2011,2013,2014, 2015, 2016, 2018,2019, 2020)
#Colors: pink gradient before 2016 and blu after
cols =  c('#8856a7','#f768a1','#c51b8a','#7a0177', '#7bccc4','#41b6c4','#2c7fb8','#253494')#brewer.pal(length(years),"Dark2")
#dev.off()
par(mfrow = c(1,1))
plot(ble_data$month[ble_data$year == years[1] ], ble_data$t[ble_data$year == years[1] ], type = 'b',
     xlab = 'Months', ylab = 'landings in tones' , col = cols[1], # main = paste('Landings')
     ylim = c(0,2500), xlim = c(0,12), lwd = 2 )

#par(mfrow = c(4,1))
for (i in 2:length(years)){
  lines(ble_data$month[ble_data$year == years[i] ], ble_data$t[ble_data$year == years[i] ], type = 'b',
       main = paste('Landings', i), col = cols[i] , lwd = 2)
}

legend('topleft', legend = as.character(years), fill = cols)

landings_january_avg = mean(ble_data$t[ble_data$month ==1])


#Look at selective fishery curve (probit)


plot(size_mean_F, sel_probit(size_mean_F), type = 'b')

#### Calibrate fishery by checking simulation results with both Fi and Pred

temperature_dataSet <- read.csv("./data/Temperature_data_1982-2018_MarBiol.csv", header = TRUE, sep = ';')
temperature_14_18 <- temperature_func(temperature_dataSet, "01/01/2014", "31/12/2018")
t_5years = seq(0,length(temperature_14_18$temperature)-0.1)

#BF_p = c(1.1476,  2.739,  5.604, 11.65, 21.524, 17.976, 11.613,  5.346) #see calculation in file 'sizeSpectraMoldel.R'
#BM_p = c(1.464, 3.645, 8.855, 22.622, 37.614) #see calculation in file 'sizeSpectraMoldel.R'


#IC_parameterized = c(P = 2, E = 0.759, L= 0.38  ,
#                     BF = BF_p, BM = BM_p) #Initial conditions according to estimations with fishery landings data.

BF_p = c(3.801852,  7.063263,  12.657331, 25.166351, 49.297133, 45.876772, 34.672340,  19.396426) #see calculation in file 'sizeSpectraMoldel.R'
BM_p = c(4.573332  , 8.203546 , 13.791605 , 26.328928 , 43.371121  ) #see calculation in file 'sizeSpectraMoldel.R'
IC_parameterized.v2 = c(P = 2, E = 0.759, L= 0.38  ,
                        BF = BF_p, BM = BM_p) #Initial conditions with reduced area (only up to 20 m depth).


bothPF_params = parameters_solv
bothPF_params$general_params$Fi = 0.05 #reduced fishery

start <- Sys.time()
test_bothFP <- solver_sizeClass.v5(t = t_5years, state = IC_parameterized.v2, parameters = bothPF_params, temperature_dataSet = temperature_14_18)
print(Sys.time() - start)

#Plot size spectra
plot_sizeSpectra(test_bothFP, 2015, title = "F&P")

#Plot fishery catch for following years:

years_to_see = c(2015, 2016, 2018)

par(mfrow = c(1,length(years_to_see)))
for (i in years_to_see){
  sol_yearToSee = test_bothFP %>%
    filter(as.numeric(format(test_bothFP$dateTime,'%Y')) == i) %>%
    mutate(month = month(dateTime)) %>%
    group_by(month) %>%
    summarise(across(where(is.numeric), \(x) sum(x, na.rm = TRUE)))

  #Calculate the total landings per month (considering 50% discard):

  monthly_landing_2015 = sol_yearToSee$catch.BF1 * 7587.1 / 1000 /2 #ttl area = 7587.1 ; convert Kg to t /1000; /2 50%bycatch
  plot(sol_yearToSee$month, monthly_landing_2015,  type = 'b',
       main = paste('Monthly landings Germany', i), ylim = c(0,2250), col = '#c994c7', lwd = 2)
  lines(ble_data$month[ble_data$year == i], ble_data$t[ble_data$year == i], col = 'gray42',  lwd = 2)

}

#Look at plankton
par(mfrow = c(1,1))
plot(test_bothFP$time, test_bothFP$P, type = 'l', xlim = c(0,365), ylim = c(0,100))



