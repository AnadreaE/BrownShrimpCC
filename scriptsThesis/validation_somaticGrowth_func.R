#<!--
# SPDX-FileCopyrightText: 2025 Helmholtz-Zentrum hereon GmbH
# SPDX-License-Identifier: CC0-1.0
# SPDX-FileContributor Andrea Farfan <farfanqbb@gmail.de>
#  -->


############################################
## The growth function developed in this  ##
## thesis will be validated and tested    ##
## in this script                         ##
############################################

library(BrownShrimp)
library(deSolve)

#for this simulation, we will consider: only juvenile stage;simulation over two years (monthly values for age);
#constant temperature and following parameters:
temperature_range <- seq(0, 30, 0.1)
days_2years <- seq(1, 2*365, 1)


#FEMALES
L_as_f <- 8.5 #Asymptotic length female ~ max. length [cm]
L_as_m = 5.5

#STEP(1) TEMMINGTS SIMULATION WITH CONSTANT TEMPERATURE

#FEMALES
# create an empty vector for each temperature degree.
temming_f <- data.frame(row = days_2years)

ind <- 1
for (j in temperature_range){
  development <- c()
  initial_l <- 0.6 #cm
  for (i in days_2years) {
    new_legth <- som_growth(initial_l, j, L_as_f, 'F')
    development <- append(development, new_legth) #append the new length
    initial_l <- new_legth
  }
  temming_f[as.character(j)] <- development
}

#STEP (2) SIMULATION WITH VB l(t) with K(T) with briere:

growth_thesis_f <- data.frame(row = days_2years)

for (i in temperature_range){
  development <- c()
  initial_l <- 0.6 #cm
  for (j in days_2years){
    growth = som_growth_thesis(initial_l, j ,i , parameters_solv$Fem_params)
    development = append(development, growth)
    }
  growth_thesis_f[as.character(i)] <- development
  }


#### plot Temming vs Thesis (constant T )####

temp_to_print <- seq(6,14, 1)

par(mfrow = c(3, 3), oma = c(5, 5, 2, 1))  # outer margins: bottom, left, top, right
par(mar = c(2, 2, 2, 1))  # inner margins for subplots (smaller so things fit)

for (i in temp_to_print ){
  RSS <- sum((temming_f[[as.character(i)]] - growth_thesis_f[[as.character(i)]] )^2)
  # Calculate total sum of squares (TSS)
  TSS <- sum((temming_f[[as.character(i)]] - mean(temming_f[[as.character(i)]]))^2)
  R_squared <- 1 - (RSS / TSS)
  plot( days_2years, temming_f[[as.character(i)]], col = 'gray41', main = paste('T = ', i, ', RSQ = ', round(R_squared, 2)),
        xlab = '', ylab = '', las = 1, lwd = 2.5, cex.axis = 1.4, cex.main = 1.3) #
  lines(days_2years, growth_thesis_f[[as.character(i)]], col =  'lightseagreen', lwd = 2)
}

mtext("days", side = 1, outer = TRUE, line = 3, cex = 1.2)
mtext("L [cm]", side = 2, outer = TRUE, line = 3, cex = 1.2)

par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", axes = FALSE, xlab = "", ylab = "")
legend("bottomright", legend = c("empirical eq.", "Thesis growth func."),
       col = c("gray41", "lightseagreen"), lwd = 3, bty = "n", horiz = TRUE, cex = 1.4,
       x.intersp = 0.3)


#### plot Temming vs Thesis (with real T data )####

temperature_dataSet <- read.csv("./data/Temperature_data_1982-2018_MarBiol.csv", header = TRUE, sep = ';')

temp_15_16 = temperature_func(temperature_dataSet, "01/01/2015", "31/12/2016")

#Simulation Temming

#Fem
temming_15_16_F = c()
initial_l <- 0.6 #cm
ind <- 1

for (j in temp_15_16$temperature){
  new_legth <- som_growth(initial_l, j, L_as_f, 'F')
  temming_15_16_F <- append(temming_15_16_F, new_legth) #append the new length
  initial_l <- new_legth
}

#Male
temming_15_16_M = c()
initial_l <- 0.6 #cm
ind <- 1

for (j in temp_15_16$temperature){
  new_legth <- som_growth(initial_l, j, L_as_m, 'M')
  temming_15_16_M <- append(temming_15_16_M, new_legth) #append the new length
  initial_l <- new_legth
}


#Simulation Thesis model
#This time I compute simulation with ode because the funciton used above works only for contat T.
#with real data, 'shrinking' considerable mm is possible and therefore not good for this purpose.

system.equations <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    L_as = parameters$L_inf
    K <- K_func_briere(temp_15_16$temperature[floor(t)], parameters)  # Access temperature for the current time
    dl.dt <- K * (L_as - l)  # Differential equation

    return(list(dl.dt))  # Return the rate of change
  })
}


# Time sequence (0 to 364 represents 365 days)
t <- seq(1, length(temp_15_16$date_time))

# Initial state for l
state <- c(l = .6)  #cm


##Simulation Fem
parameters <- parameters_solv$Fem_params
cc_f <- ode(y = state, times = t, func = system.equations, parms = parameters)# Solve the ODE
cc_f <- as.data.frame(cc_f)
head(cc_f)

##Simulation Male
parameters <- parameters_solv$M_params
cc_m <- ode(y = state, times = t, func = system.equations, parms = parameters)# Solve the ODE
cc_m <- as.data.frame(cc_m)
head(cc_m)

#convert time column to dateTime so that it can be ploted:
start_time <- min(temp_15_16$date_time)  # Start time from temp_15_16
#cc_f$time <- as.POSIXct(cc_f$time, origin = start_time)
#cc_m$time <- as.POSIXct(cc_f$time, origin = start_time)
cc_f$date_time <- start_time + cc_f$time*86400
cc_m$date_time <- start_time + cc_m$time*86400
head(cc_f)

#dev.off()
par(mfrow = c(1,2))

plot(temp_15_16$date_time, temming_15_16_F, xlab = "time", ylab = "L [cm]",
     col = 'gray41', las = 1, main = 'Female', ylim = c(0.5, 8.5))
lines(cc_f$date_time, cc_f$l, col = 'maroon', lwd = 2.5 )

plot(temp_15_16$date_time, temming_15_16_M,  xlab = "time", ylab = "L [cm]",
     col = 'gray41', las = 1, main = 'Male', ylim = c(0.5, 8.5))
lines(cc_m$date_time, cc_m$l, col = 'dodgerblue4', lwd = 2.5 )

legend("topright", legend = c("empirical eq.", "Thesis growth func. F", "Thesis growth func. M"),
       col = c("gray41", 'maroon', 'dodgerblue4'), lwd = 2.5, bty = "n", cex = 1.2)
