#HERE SIMULATIONS WITH NEW K-FUNC tpc FLEXBRIERE
library(dplyr)
library(ggplot2)
library(cowplot) #for combined plots
library(RColorBrewer)
library("colorspace")
library(deSolve)
library(tidyr)
library(RColorBrewer)
library(BrownShrimp)
library(viridis)

library(devtools)

#UPLOAD TEMPERATURE DATA

temperature_dataSet <- read.csv("C:/Users/andre/Desktop/Hereon/Data/Wadden Sea, Penning/Temperature_data_1982-2018_MarBiol.csv", header = TRUE, sep = ';')

temperature_10_11 <- temperature_func(temperature_dataSet, "01/01/2010", "31/12/2011")
temperature_10_14 <- temperature_func(temperature_dataSet, "01/01/2010", "31/12/2014")


######### START ###########
#Def. Size Classes ranges [cm]

L_as_F <- 8.5

juvI_min <- .6
juvI_max <- 1.0

juvII_min  <- 1.0
juvII_max <- 2.0

juvIII_min  <- 2.0
juvIII_max <- 3.0

juvIV_min  <- 3.0
juvIV_max <- 4.0

juvV_min  <- 4.0
juvV_max <- 5.0

adultI_min <- 5.0
adultI_max <- 7.0

adultII_min <- 7.0
adultII_max <- L_as_F

LL = 0.3
LJ = (juvI_min+juvI_max)/2 #cm
LJ2 = (juvII_min+juvII_max)/2   # 3cm
LJ3 = (juvIII_min+juvIII_max)/2
LJ4 = (juvIV_min+juvIV_max)/2
LJ5 = (juvV_min+juvV_max)/2
LA1 = (adultI_min + adultI_max)/2
LA2 = (adultII_min + adultII_max)/2  #cm

time_range <- seq(1,1095) #3 years
temperature_range <- seq(0,30, 0.1) #(10,15)

development_df_juvI <- data.frame(row = seq(1, 1095))


#simulate development with constant temperature for each 0.1°C:
for (i in temperature_range){
  development <- c()
  initial_l <- 6 #mm
  for (j in time_range){
    new_legth <- som_growth(initial_l, i, 85, 'F')
    development <- append(development, new_legth)
    initial_l <- new_legth
  }
  development_df_juvI[as.character(i)] <- development
}

dev_time_juvI <- sapply(development_df_juvI, function(col){
  BrownShrimp::count_devDays(col, juvI_min*10, juvI_max*10) } )

dev_time_juvI <- dev_time_juvI[-c(1)] #delete first element, which doesn't make sence !


#### Following lines only for exploration purposes ####
#now simulate again, with own growth function (with K_func as TPC)

developmentAF_df_juvI <- data.frame(row = seq(1, 1095))


#simulate development with constant temperature for each 0.1°C:
for (i in temperature_range){
  development <- c()
  initial_l <- 0.6 #cm
  for (j in time_range){
    new_legth <- som_growthAF(initial_l, i, "F")
    development <- append(development, new_legth)
    initial_l <- new_legth
  }
  developmentAF_df_juvI[as.character(i)] <- development
}

dev_timeAF_juvI <- sapply(developmentAF_df_juvI, function(col){
  BrownShrimp::count_devDays(col, juvI_min, juvI_max) } )

dev_timeAF_juvI <- dev_timeAF_juvI[-c(1)]


plot(temperature_range, dev_time_juvI, type = "b", col = "blue",
     xlab = "Temperature (°C)", ylab = "Development time",
     main = "Development Time vs Temperature")
lines(temperature_range, dev_timeAF_juvI, type = "b", col = 'red' , cex = 0.5)
legend("topright", legend = c("Original", "AF model"), col = c("blue", "red"), lty = 1, pch = 1)
#Looks ok in my opinion (beside the fact that from T above ~26° due to TPC curve start to rise again)


U_function <- function(a, h, x , k){
  return (a*((x-h)^22) + k )
}

a = 1.e-22
h = 15
k = 3

plot(temperature_range, U_function(a, h, temperature_range, k))

#### continue wiht shifting funciton ####
#the parameter values where calculated in file determine_growthFunction.R

shift_next_sizeClass = function(Te, sizeClass){
  ratio = 0
  if (Te <= 1) {
    ratio = 0
  } else {
    if (sizeClass == 'juvI'){
      a = 11.73753
      b = 2.205524
      c = 1747.899
      ratio = 1 / (a + c*exp(-b*Te))
    }
    if (sizeClass == 'juvII'){
      a = 26.25511 #coef(fit_juvII)[1]
      b = 1.155226  #coef(fit_juvII)[3]
      c = 2309.817 #coef(fit_juvII)[2]
      ratio = 1 / (a + c*exp(-b*Te))
    }
    if (sizeClass == 'juvIII'){
      a = 29.84754 #coef(fit_juvIII)[1]
      b = 1.164127 #coef(fit_juvIII)[3]
      c = 4889.531 #coef(fit_juvIII)[2]
      ratio = 1 / (a + c*exp(-b*Te))
    }
    if (sizeClass == 'juvIV'){
      a = 34.40965  #coef(fit_juvIV)[1]
      b = 1.099193  #coef(fit_juvIV)[3]
      c = 9141.399  #coef(fit_juvIV)[2]
      ratio = 1 / (a + c*exp(-b*Te))
    }
    if (sizeClass == 'juvV'){
      a = 43.12777 #coef(fit_juvV)[1]
      b = 1.09148  #coef(fit_juvV)[3]
      c = 21729.47 #coef(fit_juvV)[2]
      ratio = 1 / (a + c*exp(-b*Te))
    }
    if (sizeClass == 'adultI'){
      a = 142.2711 #coef(fitII_adultI)[1]
      b = 0.9313728 #coef(fitII_adultI)[2]
      c = 122405.3  #coef(fitII_adultI)[3]
      ratio = 1 / (a + c*exp(-b*Te))
    }

    if (sizeClass == 'egg'){
      ratio = 1/ (1031.34*Te^-1.345)
    }
    if (sizeClass == 'larv' ){
      ratio = 1 / ((5.5/0.00584)*Te^-1.347)
    }

    return (ratio)
  }

}


the.method = 'rk4'
h = 0.2 # half saturatuion constant original aprox .5
L_inf = 8.5
epsilon = 0.22
m = 3
const_c = 0.01 # reference density g/cm3

state2 <- c(P = 2, E = 0.1, L= 0.4 , J = .175 , J2 = 0.125, J3 = 0.125, J4 = 0.1, J5 = 0.1, A1 = 0.35, A2 = 0.1)

t = seq(0,length(temperature_10_11$temperature)-0.1, by = 0.1)
t_5years = seq(0,length(temperature_10_14$temperature)-0.1, by = 0.1)

parameters = c()

start <- Sys.time()
test_briere <- solver_sizeClass_extended_b(t = t_5years, state = state2, parameters = parameters, temperature_10_14)
print(Sys.time() - start)

start_date <- as.POSIXct("2010-01-01", format="%Y-%m-%d", tz = "UTC")

test_briere <- mutate(test_briere, dateTime = start_date + test_briere$time * 86400) #86400 seconds in one day
test_briere <- test_briere[ , -1] # delete timesteps column

#cc_long_2 <- gather(cc_2, key = "Type", value = "Value", P, E, L, J, J2, J3, J4, J5, A1, A2)
test_briere_long <- test_briere %>%
  pivot_longer(-dateTime, names_to = "variable", values_to = "value")

#dev.off()
ggplot(test_briere_long, aes(x = dateTime, y = value, color = variable)) +
  geom_line(size = 1) +  # Plot lines
  labs(title = " extended size classes Briere aprox. each 10 mm \n 2010-2014", x = "Days", y = "Biomass") +
  scale_color_manual(
    values = c( "P" = "#d9f0a3", "E"= "#1c9099" ,"L" = "gray"  ,"J" = "#bfd3e6", "J2" = "#9ebcda", "J3" = "#8c96c6", "J4" = "#8856a7", "J5" = "#810f7c" ,
                "A1" = "#fd8d3c" ,"A2" = "#e6550d" )) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), # Center title
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        legend.position="bottom") + ylim(0,50) #+ xlim( as.POSIXct("2010-04-01", format="%Y-%m-%d", tz = "UTC") , as.POSIXct("2012-07-10", format="%Y-%m-%d", tz = "UTC") )



##### Stacked area chart with K as TPC ####
#preparation:

df_long_testBriere <- test_briere[ , -1] %>% #delete  plancton col
  pivot_longer(-dateTime, names_to = "variable", values_to = "value")

#ensure that vertical order is same as life stages
stack_order_Br <- rev(c(colnames(test_briere)[2:10]) )

# Ensure 'variable' is a factor with levels in the same order
df_long_testBriere$variable <- factor(df_long_testBriere$variable, levels = stack_order_Br)


ggplot(df_long_testBriere, aes(x = dateTime, y = value, fill = variable)) +
  geom_area() +
  scale_fill_viridis_d() +
  labs(
    title = "K func as TPC",
    x = "Date",
    y = "Biomass",
    fill = "Size class"
  )


## Second test Briere:
df_long2_testBriere <- test2Briere[ , -1] %>% #delete  plancton col
  pivot_longer(-dateTime, names_to = "variable", values_to = "value")

#ensure that vertical order is same as life stages
stack_order_Br2 <- rev(c(colnames(test2Briere)[2:10]) )

# Ensure 'variable' is a factor with levels in the same order
df_long_testBriere$variable <- factor(df_long2_testBriere$variable, levels = stack_order_Br2)

#dev.off()
ggplot(df_long2_testBriere, aes(x = dateTime, y = value, fill = variable)) +
  geom_area() +
  scale_fill_viridis_d() +
  labs(
    title = "K func as TPC with more food",
    x = "Date",
    y = "Biomass",
    fill = "Size class"
  )




