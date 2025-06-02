library(dplyr)
library(ggplot2)
library(cowplot) #for combined plots
library(RColorBrewer)
library(deSolve)
library(patchwork)
library(minpack.lm)
library(BrownShrimp)

### UPLOAD AND FORMAT TEMPERATURE DATA ###

temperature_ds <- read.csv("C:/Users/andre/Desktop/Hereon/Data/Wadden Sea, Penning/Temperature_data_1982-2018_MarBiol.csv", header = TRUE, sep = ';')
#example to read from package <- read.csv("./data/Temperature_data_1982-2018_MarBiol.csv", header = TRUE, sep = ';')

#temperature_00_01 <- temperature_func(temperature_ds, "01/01/2000", "31/12/2001")
#temperature_05_06 <- temperature_func(temperature_ds, "01/01/2005", "31/12/2006")
#temperature_10_11 <- temperature_func(temperature_ds, "01/01/2010", "31/12/2011")
temperature_15_16 <- temperature_func(temperature_ds, "01/01/2015", "31/12/2016")
#temperature_17_18 <- temperature_func(temperature_ds, "01/01/2017", "31/12/2018")

#For the time being I will only work with data from 2010 and 2011.

#### EGG DEVELOPMENT #####
egg_development_array <- c()
hatch_status_array <- c()
initial_status_egg <- 0

#first_day_egg <- egg_dev(initial_status_egg, temperature_15_16$temperature[1])
#initial_status_egg <- first_day_egg[[1]]

for (i in 1:length(temperature_15_16$temperature) ){
  egg_development <- egg_dev(initial_status_egg, temperature_15_16$temperature[i])
  initial_status_egg <- egg_development[[1]]
  hatch_status_array <- append(hatch_status_array, egg_development[[2]] )
  egg_development_array <- append(egg_development_array, initial_status_egg)
}


#### SOMATIC GROWTH (FROM 6 MM) ####
temperature_range <- seq(0, 30, 0.1)
days <- seq(1, length(temperature_15_16$temperature))
L_as_f <- 85 #mm
L_as_m <- 55 #mm

#FEMALES
#STEP(1) Temmings simulations with constant temperatures
# create an empty DF for each temperature degree.
growth_fem <- data.frame(days=seq(1, length(temperature_15_16$temperature)))

for (j in temperature_range){
  initial_l <- 6 #mm
  development <- c()
  for (i in 1:length(temperature_15_16$temperature)) {
    new_legth <- som_growth(initial_l, j, L_as_f, 'F')
    development <- append(development, new_legth)#c(growth_fem[[ind]], new_legth) #append the new length
    initial_l <- new_legth
  }
  growth_fem[as.character(j)] <- development
}


#STEP(2) Fit Bertalanffy eq. to female simulation
K_estimations_f <- c()
all_fits_fem <- list()

for (i in 1:length(temperature_range)){
  #teperature.char <-
  set_toFit <- pull(growth_fem, as.character(temperature_range[i]))
  fit <- nls( set_toFit ~ ( L_as_f * (1 - exp(-k * days)) ),
              start = list(k = 0.5))
  all_fits_fem <- c(all_fits_fem, list(fit) )
  K_estimations_f <- append(K_estimations_f, coef(fit)["k"])
}

head(K_estimations_f)
length(K_estimations_f)

#STEP(3) plot to see the shape of the curve:
plot(temperature_range, K_estimations_f, main='K estiations female')
T_opt <- which.max(K_estimations_f)*0.1 #this is the T_opt



#MALES
#STEP(1) Temmings simulations with constant temperatures
growth_masc <- data.frame(days=seq(1, length(temperature_15_16$temperature)))

for (j in temperature_range){
  initial_l <- 6 #mm
  development <- c()
  for (i in 1:length(temperature_15_16$temperature)) {
    new_legth <- som_growth(initial_l, j, L_as_m, 'M')
    development <- append(development, new_legth)#c(growth_fem[[ind]], new_legth) #append the new length
    initial_l <- new_legth
  }
  growth_masc[as.character(j)] <- development
}


#STEP(2) Fit Bertalanffy eq. to female simulation
K_estimations_m <- c()
all_fits_masc <- list()

for (i in 1:length(temperature_range)){
  #teperature.char <-
  set_toFit <- pull(growth_masc, as.character(temperature_range[i]))
  fit <- nls( set_toFit ~ ( L_as_m * (1 - exp(-k * days)) ),
              start = list(k = 0.5))
  all_fits_masc <- c(all_fits_masc, list(fit) )
  K_estimations_m <- append(K_estimations_m, coef(fit)["k"])
}

head(K_estimations_m)
length(K_estimations_m)

#STEP(3) plot to see the shape of the curve:
plot(temperature_range, K_estimations_m, main='K estiations male')


#### STEP (4) Now apply regression to K-values 'curve' with gaussian distribution (F and M): ####
#Fem
fit_gaus_f <- nls(K_estimations_f ~  A * exp(-((temperature_range - mu)^2) / (2 * sigma^2)),
                  start = list(A = 0.0072, mu = 17, sigma = 1), trace = TRUE)

summary(fit_gaus_f)

#Detemine the goodness of the fit:
# Calculate residual sum of squares (RSS)
RSS_f <- sum((K_estimations_f - predict(fit_gaus_f))^2)
# Calculate total sum of squares (TSS)
TSS_f <- sum((K_estimations_f - mean(K_estimations_f))^2)
# Calculate R-squared
RSQ_gauss_f <- 1 - (RSS_f / TSS_f)
cat("R-squared: ", RSQ_gauss_f, "\n")


#Masc
fit_gaus_m <- nls(K_estimations_m ~  A * exp(-((temperature_range - mu)^2) / (2 * sigma^2)),
                  start = list(A = 0.016, mu = 21, sigma = 0.8), trace = TRUE)

summary(fit_gaus_m)
#plot(temperature_range, K_estimations_f)
#lines(temperature_range, predict(fit_gaus_f))

#Detemine the goodness of the fit:
# Calculate residual sum of squares (RSS)
RSS_m <- sum((K_estimations_m - predict(fit_gaus_m))^2)
# Calculate total sum of squares (TSS)
TSS_m <- sum((K_estimations_m - mean(K_estimations_m))^2)
# Calculate R-squared
RSQ_gauss_m <- 1 - (RSS_m / TSS_m)
cat("R-squared: ", RSQ_gauss_m, "\n")


dev.off()
plot(temperature_range, K_estimations_m, main = 'K - values M and F for each degree °C \n vs Fit Gauss', col = 'steelblue4',
     xlab = 'Temperature (?C)', ylab = 'K', cex.lab = 1.7, cex.main = 1.6, las = 1, lwd = 2)
points(temperature_range,K_estimations_f, col = 'orangered4', lwd = 2)
lines(temperature_range, predict(fit_gaus_f), col = 'lightsalmon2', lwd = 2)
lines(temperature_range, predict(fit_gaus_m), col = 'lightseagreen', lwd = 2)
legend("topleft",
       legend = c("Male", "Female", 'Fit M', 'Fit F'),
       col = c("steelblue4", "orangered4", 'lightseagreen', 'lightsalmon2'),
       lty = c(NA, NA, 1, 1),
       pch = c(16,16,NA, NA),
       cex = 1.25)


#### STEP (4) Again with TPC Briere instead of Gauss (F and M): ####
#TPC: Termal Performance Curve
T_min = 0.5
T_max = 30
T_range = seq(T_min, T_max, by= 0.1)

#Fem
K_estimations_f_df<- data.frame(
  temperature = temperature_range,
  performance = K_estimations_f
)

K_estimations_f_df_filtered <- subset(K_estimations_f_df, temperature > T_min & temperature < T_max) #the fit can be applied only to data within Temperature range!

fit_briere_f <- nlsLM(performance ~  c*temperature*(temperature - T_min_f)*((T_max_f - temperature)^(1/m) ),
                    data = K_estimations_f_df_filtered,
                    start = list(c = 0.000000001, T_min_f = 1, T_max_f = 35, m = 0.4),
                    lower = c(c = 0, T_min_f = 0, T_max_f = 15, m = 0.05),
                    upper = c(c = 0.1, T_min_f = 1, T_max_f = 40, m = 2),
                    control = nls.lm.control(maxiter = 500) )
summary(fit_briere_f)

RSS_f_b <- sum((K_estimations_f_df_filtered$performance - predict(fit_briere_f))^2)
# Calculate total sum of squares (TSS)
TSS_f_b <- sum((K_estimations_f_df_filtered$performance - mean(K_estimations_f_df_filtered$performance))^2)
RSQ_f_briere <- 1 - (RSS_f_b / TSS_f_b)
cat("R-squared: ", RSQ_f_briere, "\n")


#Masc
K_estimations_m_df<- data.frame(
  temperature = temperature_range,
  performance = K_estimations_m
)

K_estimations_m_df_filtered <- subset(K_estimations_m_df, temperature > T_min & temperature < T_max)  #the fit can be applied only to data within Temperature range!

fit_briere_m <- nlsLM(performance ~  c*temperature*(temperature - T_min_f)*((T_max_f - temperature)^(1/m) ),
                    data = K_estimations_m_df_filtered,
                    start = list(c = 0.000000001, T_min_f = 1, T_max_f = 35, m = 0.4),
                    lower = c(c = 0, T_min_f = 0, T_max_f = 15, m = 0.05),
                    upper = c(c = 0.1, T_min_f = 1, T_max_f = 40, m = 2),
                    control = nls.lm.control(maxiter = 500) )
summary(fit_briere_m)

RSS_m_b <- sum((K_estimations_m_df_filtered$performance - predict(fit_briere_m))^2)
# Calculate total sum of squares (TSS)
TSS_m_b <- sum((K_estimations_m_df_filtered$performance - mean(K_estimations_m_df_filtered$performance))^2)
RSQ_m_briere <- 1 - (RSS_m_b / TSS_m_b)
cat("R-squared: ", RSQ_m_briere, "\n")

plot(temperature_range, K_estimations_m, main = 'K - values M and F for each degree °C \n vs Fit Briere', col = 'steelblue4',
     xlab = 'Temperature (°C)', ylab = 'K', cex.lab = 1.7, cex.main = 1.6, las = 1, lwd = 2)
points(temperature_range,K_estimations_f, col = 'orangered4', lwd = 2)
lines(K_estimations_f_df_filtered$temperature, predict(fit_briere_f), col = 'lightsalmon2', lwd = 2)
lines(K_estimations_m_df_filtered$temperature, predict(fit_briere_m), col = 'lightseagreen', lwd = 2)
legend("topleft",
       legend = c("Male", "Female", 'Fit M', 'Fit F'),
       col = c("steelblue4", "orangered4", 'lightseagreen', 'lightsalmon2'),
       lty = c(NA, NA, 1, 1),
       pch = c(16,16,NA, NA),
       cex = 1.25)


#### STEP (4) Again with flexTPC Briere (F and M): ####
# FEMALE
r_max_f = max(K_estimations_f_df_filtered$performance)

fit_flexBriere_f <- nlsLM( performance ~  r_max_f*( (((temperature - T_min)/alpha)^alpha) * (((T_max - temperature)/ (1-alpha) )^(1-alpha)) * (1 / (T_max- T_min))  )^(alpha*(1-alpha)/beta^2) ,
                           data = K_estimations_f_df_filtered,
                           start = list(alpha = 0.5, beta = 0.25 ),
                           lower = c(alpha = 0 , beta = 0.05),
                           upper = c(alpha = 1, beta = 1),
                           control = nls.lm.control(maxiter = 500) )
summary(fit_flexBriere_f)

RSS_f_fb <- sum((K_estimations_f_df_filtered$performance - predict(fit_flexBriere_f))^2)
# Calculate total sum of squares (TSS)
TSS_f_fb <- sum((K_estimations_f_df_filtered$performance - mean(K_estimations_f_df_filtered$performance))^2)
RSQ_f_flexBriere <- 1 - (RSS_f_fb / TSS_f_fb)
cat("R-squared: ", RSQ_f_flexBriere, "\n")



# MALE
r_max_m = max(K_estimations_m_df_filtered$performance)

fit_flexBriere_m <- nlsLM( performance ~  r_max_m*( (((temperature - T_min)/alpha)^alpha) * (((T_max - temperature)/ (1-alpha) )^(1-alpha)) * (1 / (T_max- T_min))  )^(alpha*(1-alpha)/beta^2) ,
                      data = K_estimations_m_df_filtered,
                      start = list(alpha = 0.5, beta = 0.25 ),
                      lower = c(alpha = 0 , beta = 0.05),
                      upper = c(alpha = 1, beta = 1),
                      control = nls.lm.control(maxiter = 500) )
summary(fit_flexBriere_m)

RSS_m_fb <- sum((K_estimations_m_df_filtered$performance - predict(fit_flexBriere_m))^2)
# Calculate total sum of squares (TSS)
TSS_m_fb <- sum((K_estimations_m_df_filtered$performance - mean(K_estimations_m_df_filtered$performance))^2)
RSQ_m_felBriere <- 1 - (RSS_m_fb / TSS_m_fb)
cat("R-squared: ", RSQ_m_felBriere, "\n")

plot(temperature_range, K_estimations_m, main = 'K - values M and F for each degree °C \n vs Fit felxBriere', col = 'steelblue4',
     xlab = 'Temperature (°C)', ylab = 'K', cex.lab = 1.7, cex.main = 1.6, las = 1, lwd = 2)
points(temperature_range,K_estimations_f, col = 'orangered4', lwd = 2)
lines(K_estimations_f_df_filtered$temperature, predict(fit_flexBriere_f), col = 'lightsalmon2', lwd = 2)
lines(K_estimations_m_df_filtered$temperature, predict(fit_flexBriere_m), col = 'lightseagreen', lwd = 2)
legend("topleft",
       legend = c("Male", "Female", 'Fit M', 'Fit F'),
       col = c("steelblue4", "orangered4", 'lightseagreen', 'lightsalmon2'),
       lty = c(NA, NA, 1, 1),
       pch = c(16,16,NA, NA),
       cex = 1.25)


############ VALIDATION OF K(T) *Andersen with Gauss vs Temming* ###############
#pending after is decided which type of

temperature_15 <- temperature_func(temperature_ds, "01/01/2015", "31/12/2015")

#Simulation Temming size growth T data from 2015
temming_growth_f <- c()
temming_growth_m <- c()
L_as_f <- 85
L_as_m <- 55

#simulation female
initial_size <- 6 #mm
for (i in 1:length(temperature_15$date_time) ){
  development_Temming <- som_growth(initial_size, temperature_15$temperature[i], L_as_f ,"F")
  initial_size <- development_Temming
  temming_growth_f <- append(temming_growth_f, development_Temming)
}

#simulation male
initial_size <- 6 #mm
for (i in 1:length(temperature_15$date_time) ){
  development_Temming <- som_growth(initial_size, temperature_15$temperature[i], L_as_m ,'M')
  initial_size <- development_Temming
  temming_growth_m <- append(temming_growth_m, development_Temming)
}


#SIMULATE ANDERSEN
# Define the function K based on temperature
func_K <- function(T, sex) {
  if (sex == 'F'){
    alpha <- 0.007614
    T_bar <- 16.24
    T_width <- 7.098
  } else if (sex=='M'){
    alpha <- 0.01604
    T_bar <- 18.29
    T_width <- 7.529
  }

  alpha * exp(-((T - T_bar)^2 / (2 * T_width^2)))
}

# Define the system of equations (solver for dl/dt)
system.equations <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    if (sex=='F'){
      L_as = 85
    } else if (sex=='M'){
      L_as = 55
    } else {
      stop("Invalid sex value: Please ensure 'sex' is either 'F' or 'M'.")
    }
    K <- func_K(temperature_15$temperature[floor(t)], sex)  # Access temperature for the current time
    dl.dt <- K * (L_as - l)  # Differential equation

    return(list(dl.dt))  # Return the rate of change
  })
}


# Time sequence (0 to 364 represents 365 days)
t <- seq(1, 364, 0.01)

# Initial state for l
state <- c(l = 6)  #


##Simulation Fem
parameters <- list(sex='F')
cc_f <- ode(y = state, times = t, func = system.equations, parms = parameters)# Solve the ODE
cc_f <- as.data.frame(cc_f)
head(cc_f)

##Simulation Male
parameters <- list(sex='M')
cc_m <- ode(y = state, times = t, func = system.equations, parms = parameters)# Solve the ODE
cc_m <- as.data.frame(cc_m)
head(cc_m)

# Plot the results: Temming vs Andersen

days_1year <- seq(1,365)
female_plot <- ggplot() +
  geom_point(aes(x = days_1year, y = temming_growth_f, color = "Empirical function"),
             size = 2, shape = 21, fill = "gray") +  # Use points for temperature growth data
  geom_line(aes(x = t, y = cc_f$l, color = "Analytical model"),
            size = 0.75) +  # Line for cc_f
  labs(x = "Day of Year", y = "Length (l)", title = "Female Growth") +
  scale_color_manual(name = "Legend",
                     values = c("Empirical function" = "gray" , "Analytical model" = "maroon3") )+
  #labels = c("Empiric", "Analytical")) +
  theme_minimal() +
  theme(
    legend.text = element_text(size = 14),   # Increase legend text size
    legend.title = element_text(size = 14),  # Increase legend title size
    axis.title.x = element_text(size = 14),  # Increase x-axis title size
    axis.title.y = element_text(size = 14),  # Increase y-axis title size
    plot.title = element_text(size = 16, hjust = 0.5) # Increase title size and center it
  )

print(female_plot)


male_plot <- ggplot() +
  geom_point(aes(x = days_1year, y = temming_growth_m, color = "Empirical function"),
             size = 2, shape = 21, fill = "gray") +  # Use points for temperature growth data
  geom_line(aes(x = cc_m$time, y = cc_m$l, color = 'Analytical model'),
            size = 0.75) +  # Line for cc_f
  labs(x = "Day of Year", y = "Length (l)", title = "Male Growth") +
  scale_color_manual(name = "Legend",
                     values = c("Empirical function" = "gray", "Analytical model" = "darkblue") )+
  theme_minimal()  +
  theme(
    legend.text = element_text(size = 14),   # Increase legend text size
    legend.title = element_text(size = 14),  # Increase legend title size
    axis.title.x = element_text(size = 14),  # Increase x-axis title size
    axis.title.y = element_text(size = 14),  # Increase y-axis title size
    plot.title = element_text(size = 16, hjust = 0.5) # Increase title size and center it
  )

# Print the female plot
print(male_plot)

############ VALIDATION OF K(T) *Andersen with flexTCP vs Temming* ###############
func_K_briere <- function(temperature) { #refers for the time being only females
  r_max_f = 0.007638409
  T_min = 0.4
  T_max = 30
  alpha = 0.541524
  beta = 0.276430
  diff_min = temperature - T_min
  diff_max = T_max - temperature
  diff = T_max - T_min
  alpha_invert = 1 - alpha
  return(r_max_f*( ((diff_min/alpha)^alpha) * ((diff_max/ alpha_invert )^alpha_invert) * (1 / diff)  )^(alpha*alpha_invert/(beta^2) ))
}

# Define the system of equations (solver for dl/dt)
temperature_15_17 <- temperature_func(temperature_ds, "01/01/2015", "31/12/2017")


system.equations <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    if (sex=='F'){
      L_as = 85
    } else if (sex=='M'){
      L_as = 55
    } else {
      stop("Invalid sex value: Please ensure 'sex' is either 'F' or 'M'.")
    }
    K <- func_K_briere(temperature_15_17$temperature[floor(t)])  # Access temperature for the current time
    dl.dt <- K * (L_as - l)  # Differential equation

    return(list(dl.dt))  # Return the rate of change
  })
}

t2 = seq(1, 730, 0.01)
t3 = seq(1, 1095, 0.01)

##Simulation Fem
parameters <- list(sex='F')
cc_f_b <- ode(y = state, times = t3, func = system.equations, parms = parameters)# Solve the ODE
cc_f_b <- as.data.frame(cc_f_b)
head(cc_f_b)

#Simulate Temming again 2 years
temming_growth_f <- c()
temming_growth_m <- c()
L_as_f <- 85
L_as_m <- 55

#simulation female
initial_size <- 6 #mm
for (i in 1:length(temperature_15_17$date_time) ){
  development_Temming <- som_growth(initial_size, temperature_15_17$temperature[i], L_as_f ,"F")
  initial_size <- development_Temming
  temming_growth_f <- append(temming_growth_f, development_Temming)
}


##Plot Fem
female_plot_b <- ggplot() +
  geom_point(aes(x = 1: length(temperature_15_17$temperature), y = temming_growth_f, color = "Empirical function"),
             size = 2, shape = 21, fill = "gray") +  # Use points for temperature growth data
  geom_line(aes(x = cc_f_b$time, y = cc_f_b$l, color = "Analytical model"),
            size = 0.75) +  # Line for cc_f
  labs(x = "Day of Year", y = "Length (mm)", title = "Female Growth") +
  scale_color_manual(name = "Legend",
                     values = c("Empirical function" = "gray" , "Analytical model" = "maroon3") )+
  #labels = c("Empiric", "Analytical")) +
  theme_minimal() +
  theme(
    legend.text = element_text(size = 14),   # Increase legend text size
    legend.title = element_text(size = 14),  # Increase legend title size
    axis.title.x = element_text(size = 14),  # Increase x-axis title size
    axis.title.y = element_text(size = 14),  # Increase y-axis title size
    plot.title = element_text(size = 16, hjust = 0.5) # Increase title size and center it
  )

print(female_plot_b)


