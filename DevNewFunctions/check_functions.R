#<!--
# SPDX-FileCopyrightText: 2025 Helmholtz-Zentrum hereon GmbH
# SPDX-License-Identifier: Apache-2.0
# SPDX-FileContributor Andrea Farfan <farfanqbb@gmail.de>
#  -->


library(dplyr)
library(lubridate)
library(BrownShrimp)
library(deSolve)

#Here I will check plausibility of individual functions that are used in the model

#### Some params ####
temp_range = seq(1,30, 0.1)

the.method = 'rk4'
h = 0.2 # half saturatuion constant original aprox .5
L_inf = 8.5
epsilon = 0.22
m = 3
const_c = 0.01

#### CHECK K FUNC ####

#SEE K-FUNC GAUSS VS TPC
plot(temp_range, sapply(temp_range, K_func), type = "b", col = "blue",
     xlab = "Temperature (°C)", ylab = "K(T)",
     main = "K func Gauss vs Briere")
lines(temp_range, sapply(temp_range, K_func_briere, sex_params = parameters_solv$Fem_params), type = "b", col = 'red' , cex = 0.5)
legend("topright", legend = c("Gauss", "Briere"), col = c("blue", "red"), lty = 1, pch = 1)


#### CHECK NATURAL MORTALITY ####

natural_mort_br <- sapply(temp_range, respiration_rate_b, L=5, sex_params = parameters_solv$Fem_params)
#dev.off()
plot(temp_range,  natural_mort_br, type = "b", col = "blue",
     xlab = "Temperature (°C)", ylab = " ",
     main = "natural mortality Gauss vs Briere")
#lines(temp_range,respiration_rate(temp_range, 5), type = "b", col = 'red' , cex = 0.5)
legend("topright", legend = c("Gauss", "Briere"), col = c("blue", "red"), lty = 1, pch = 1)


#Now Check how it increases with time (so test constant temperatures)

Temp = c(5, 10, 15, 20, 25)
length_range = seq(0.6, 8.5, 0.05)
colors = c('lightpink', 'navy', 'olivedrab1', 'peru', 'violet')

plot( length_range, sapply(length_range, respiration_rate_b, temperature = Temp[1], sex_params = parameters_solv$Fem_params),
      ylim = c(0,0.015), col = colors[1])
for (i in 2:5){
  lines(length_range, sapply(length_range, respiration_rate_b, temperature = Temp[i], sex_params = parameters_solv$Fem_params),
        col = colors[i], lwd = 2)
}
legend("topleft", legend = paste( 'T= ', as.character(Temp)),
       col = colors, lty = 1, pch = c(1,2,2,2,2))


#### CHECK FUNCTIONAL RESPONSE ####

func_response = function(P){
  alpha_ir = 0.1721763#8.238
  return(P/(P+alpha_ir))
}

ingestion_rate_b_test = function(temperature, L, P, sex_p){
  c_div = const_c / convertL_to_W(L)
  fr = func_response(P)
  #params = sex_parameters_func(sex)
  #L_infty = params$L_inf
  L_infty = sex_p$L_inf
  #result = m*K_func_briere(temperature, sex)*convertL_to_W(L)*L_infty*(c_div)^(1/m)*( P/(P+alpha_ir) )
  w = convertL_to_W(L)
  result = m*K_func_briere(temperature, sex_p)*(L_infty/L)*fr
  #m*convertL_to_W(L)*K_func(T)*L_inf/L*P/(P+h) #revised formula Andrea 04.04.25
  return (result)
}

B_population = 50 #biomass
Temp = 10

const_food = c(1,2, 3, 5, 15, 30)
time_range = seq(1, 365)

food_3years = new_food(time_range, Te = Temp)

#consumption = B_population*ingestion_rate_b_test(Temp, 4, const_food, parameters_solv$Fem_params)

#plot(time_range, consumption, main = paste("T =", Temp, "B =", B_population, "constant food = ", const_food) )

par(mfrow = c(2,3))

for (i in const_food){
  const_food_loop = rep(i, times = 365 )
  consumption = B_population*ingestion_rate_b_test(Temp, 4, const_food_loop, parameters_solv$Fem_params)
  plot(time_range, consumption, main = paste("T =", Temp, "B =", B_population, "constant food = ", i), ylim = c(1,2) )

}

dev.off()
food_increasing = seq(1,500)

plot(food_increasing, B_population*ingestion_rate_b_test(Temp, 4, food_increasing, parameters_solv$Fem_params))

#now check new updated instion rate function
par(mfrow = c(2,3))

for (i in const_food){
  const_food_loop = rep(i, times = length(time_range) )
  consumption = B_population*ingestion_rate_b(Temp, 4, const_food_loop, parameters_solv$general_params, parameters_solv$Fem_params)
  plot(time_range, consumption, main = paste("T =", Temp, "B =", B_population, "constant food = ", i), ylim = c(1,2) )

}




#### TEST SELECTIVE FISHERY ####

fishery = function(l){
  L50 = 4.49 #[cm]
  SR = 1.56
  nominator = exp( (1.349/SR) * (l - L50))
  denominator = 1 + nominator
  return(nominator / denominator)
}

l_range = seq(0.6, 8.5, 0.1)


plot(l_range, fishery(l_range))


#### CHECK SPAWNING ####
# check and compare share of avigerous females against Hünerlage data###
temperature_dataSet <- read.csv("./data/Temperature_data_1982-2018_MarBiol.csv", header = TRUE, sep = ';')
T_13_16 = temperature_func(temperature_dataSet, "01/01/2013", "31/12/2016")
L_avg = (5.0+8.5)/2
share_OF = T_13_16 %>%
          mutate(shareOF = sapply(T_13_16$temperature, molting_fraction, L = L_avg))

#Now simulate this time period:

state_rev1 = c(P = 2, E = 0.1, L= 0.0928 *1.5 ,
               J_f = 2, J2_f = 2, J3_f = 2 , J4_f = 1.8, J5_f = 1.8,
               J_m = 2 , J2_m = 2, J3_m = 2, J4_m = 1.8, J5_m = 1.8,
               A1_f = 5.443*0.5, A1_m = 5.443*0.5,
               A2 = 5.443*0.375, A3=5.443*0.125) #males doesn't reach this size classes

t_3years = seq(0,length(T_13_16$temperature)-0.1, by = 0.1)

parameters = parameters_solv


start <- Sys.time()
check_OF <- solver_sizeClass_sex.v2(t = t_3years, state = state_rev1, parameters = parameters, temperature_dataSet = T_13_16)
print(Sys.time() - start)

start_date <- as.POSIXct("2013-01-01", format="%Y-%m-%d", tz = "UTC")

check_OF <- mutate(check_OF, dateTime = start_date + check_OF$time * 86400) #86400 seconds in one day
check_OF <- check_OF[ , -1]

#Filter only spring data share_OF:
share_OF$date_time <- as.POSIXct(share_OF$date_time)

# Extract month
share_OF_spring = share_OF %>%
  filter(month(date_time) %in% 3:5)

#Filter only spring data simulation:
check_OF$dateTime <- as.POSIXct(check_OF$dateTime)
check_OF_spring = check_OF %>%
  filter(month(dateTime) %in% 3:5) %>%
  mutate(date = as.Date(dateTime)) %>%           # Extract date only
  group_by(date) %>%                             # Group by just the date
  summarise(across(where(is.numeric), \(x) mean(x, na.rm = TRUE)))


density_OF = check_OF_spring %>%
  mutate (densityOF = share_OF_spring$shareOF*check_OF_spring$A1_f + share_OF_spring$shareOF*check_OF_spring$A2 + share_OF_spring$shareOF*check_OF_spring$A3)

den_sum2013 = density_OF %>%
  filter(year(date) == 2013) %>%
  summarise(total_density = sum(densityOF, na.rm = TRUE))

den_sum2014 = density_OF %>%
  filter(year(date) == 2014) %>%
  summarise(total_density = sum(densityOF, na.rm = TRUE))

den_sum2015 = density_OF %>%
  filter(year(date) == 2015) %>%
  summarise(total_density = sum(densityOF, na.rm = TRUE))

den_sum2016 = density_OF %>%
  filter(year(date) == 2016) %>%
  summarise(total_density = sum(densityOF, na.rm = TRUE))


barplot(c(den_sum2013$total_density, den_sum2014$total_density, den_sum2015$total_density, den_sum2016$total_density), arg.names = c('2013', '2014', '2015', '2016'))


#Check monthly values 2013

share_OF2013 = share_OF %>%
  filter(year(date_time) == 2013)

check_OF2013 =  check_OF %>%
      filter(year(dateTime) == 2013) %>%
      mutate(date = as.Date(dateTime)) %>%           # Extract date only
      group_by(date) %>%                             # Group by just the date
      summarise(across(where(is.numeric), \(x) mean(x, na.rm = TRUE)))

density_OF2013 = share_OF2013 %>%
  mutate(density = share_OF2013$shareOF*check_OF2013$A1_f +  share_OF2013$shareOF*check_OF2013$A2 + share_OF2013$shareOF*check_OF2013$A3)

monthly2013 = density_OF2013 %>%
 # mutate(date = as.Date(date_time)) %>%
  group_by(month(date_time))%>%                             # Group by just the date
  summarise(across(where(is.numeric), \(x) sum(x, na.rm = TRUE)))

months = c('Jan', 'Feb', 'Mar', 'Apr', 'Mai', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec')

barplot(monthly2013$density, names.arg = months)

#Check monthly values 2014

share_OF2014 = share_OF %>%
  filter(year(date_time) == 2014)

check_OF2014 =  check_OF %>%
  filter(year(dateTime) == 2014) %>%
  mutate(date = as.Date(dateTime)) %>%           # Extract date only
  group_by(date) %>%                             # Group by just the date
  summarise(across(where(is.numeric), \(x) mean(x, na.rm = TRUE)))

density_OF2014 = share_OF2014 %>%
  mutate(density = share_OF2014$shareOF*check_OF2014$A1_f +  share_OF2013$shareOF*check_OF2014$A2 + share_OF2014$shareOF*check_OF2014$A3)

monthly2013 = density_OF2014 %>%
  # mutate(date = as.Date(date_time)) %>%
  group_by(month(date_time))%>%                             # Group by just the date
  summarise(across(where(is.numeric), \(x) sum(x, na.rm = TRUE)))

months = c('Jan', 'Feb', 'Mar', 'Apr', 'Mai', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec')

barplot(monthly2013$density, names.arg = months, main = '2014')

dev.off()

#Check monthly values 2015

share_OF2015 = share_OF %>%
  filter(year(date_time) == 2015)

check_OF2015 =  check_OF %>%
  filter(year(dateTime) == 2015) %>%
  mutate(date = as.Date(dateTime)) %>%           # Extract date only
  group_by(date) %>%                             # Group by just the date
  summarise(across(where(is.numeric), \(x) mean(x, na.rm = TRUE)))

density_OF2015 = share_OF2015 %>%
  mutate(density = share_OF2015$shareOF*check_OF2015$A1_f +  share_OF2013$shareOF*check_OF2015$A2 + share_OF2015$shareOF*check_OF2015$A3)

monthly2013 = density_OF2015 %>%
  # mutate(date = as.Date(date_time)) %>%
  group_by(month(date_time))%>%                             # Group by just the date
  summarise(across(where(is.numeric), \(x) sum(x, na.rm = TRUE)))

months = c('Jan', 'Feb', 'Mar', 'Apr', 'Mai', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec')

barplot(monthly2013$density, names.arg = months, main = '2015')


#Check monthly values 2016

share_OF2016 = share_OF %>%
  filter(year(date_time) == 2016)

check_OF2016 =  check_OF %>%
  filter(year(dateTime) == 2016) %>%
  mutate(date = as.Date(dateTime)) %>%           # Extract date only
  group_by(date) %>%                             # Group by just the date
  summarise(across(where(is.numeric), \(x) mean(x, na.rm = TRUE)))

density_OF2016 = share_OF2016 %>%
  mutate(density = share_OF2016$shareOF*check_OF2016$A1_f +  share_OF2013$shareOF*check_OF2016$A2 + share_OF2016$shareOF*check_OF2016$A3)

monthly2013 = density_OF2016 %>%
  # mutate(date = as.Date(date_time)) %>%
  group_by(month(date_time))%>%                             # Group by just the date
  summarise(across(where(is.numeric), \(x) sum(x, na.rm = TRUE)))

months = c('Jan', 'Feb', 'Mar', 'Apr', 'Mai', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec')

barplot(monthly2013$density, names.arg = months, main = '2016')



### test a fishery function with reltive small peaks in september

fishery = function(t) {
  toReturn = 0.2 * (1.05 + cos( (2*pi/ 365)*(t-105) ) )
  return(toReturn)
}

days = seq(1,365*5)

plot(days, fishery(days))


month_Feffort = c(0.19, 0.2, 0.86, 1.6, 1.39, 1.26, 1.19, 1.25, 1.27, 1.26, 1.09, 0.45)
x_ax = seq(1, length(month_Feffort))

plot(x_ax, month_Feffort)



#### AFTER DISCUSSION IN MEETING 03.06. ####
#Check / try to understand why diff. of biomass between Larvae and Juv. 1 is so abrupt
#First: Hatching and growth from Larvae to Juv I are only temperature dependent
#whereas growth function is size AND temperature dependent. So lets check growth ratio:

par(mfrow = c(2,2))

for (i in Temp[1:4]){
  plot(length_range,
       sapply(length_range, shift_next_sizeClass, temperatur = i, sex_params = parameters_solv$Fem_params),
       main = paste('Growth rates at T =', i ), ylim = c(0, 0.18),
       xlab = 'L [cm]', ylab = 'rate', las = 1)
  points(0.5, hatch_eggs(i), col = 'khaki4', bg = 'khaki3', pch = 21, cex = 2)
  points(0.8, shiftTo_juvenile(i), col = 'slateblue4', bg = 'slateblue', pch = 21, cex = 2 )
  legend('topright', legend= c('Eggs', 'Larvae', 'Juvenile'), col = c('khaki4', 'slateblue', 'black'),
         pch = c(1,2,2,2,2), lty = 1, cex = 1.05, lwd = 2)
}


developmentTime_larvae = 1/shiftTo_juvenile(12)
developmentTime_juv3 = 1/shift_next_sizeClass(5,12, parameters_solv$Fem_params)













