library(dplyr)
library(ggplot2)
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

temperature_95_96 <- temperature_func(temperature_dataSet, "01/01/1995", "31/12/1996")
temperature_10_11 <- temperature_func(temperature_dataSet, "01/01/2010", "31/12/2011")
temperature_10_14 <- temperature_func(temperature_dataSet, "01/01/2010", "31/12/2014")

######### START ###########

#Def. Size Classes ranges

L_as_F <- 85

juvI_min <- 6
juvI_max <- 10

juvII_min  <- 10
juvII_max <- 20

juvIII_min  <- 20
juvIII_max <- 30

juvIV_min  <- 30
juvIV_max <- 40

juvV_min  <- 40
juvV_max <- 50

adultI_min <- 50
adultI_max <- 70

adultII_min <- 70
adultII_max <- L_as_F


LL = 0.3
LJ = (juvI_min+juvI_max)/2 #cm
LJ2 = (juvII_min+juvII_max)/2 /10  # 3cm
LJ3 = (juvIII_min+juvIII_max)/2 /10
LJ4 = (juvIV_min+juvIV_max)/2 /10
LJ5 = (juvV_min+juvV_max)/2 /10
LA1 = (adultI_min + adultI_max)/2 /10
LA2 = (adultII_min + adultII_max)/2 /10 #cm

time_range <- seq(1,1095)
temperature_range <- seq(0,25, 0.1) #(10,15)

## Calculate development time for each size class

#### JUVI ####
development_df_juvI <- data.frame(row = seq(1, 1095))


#simulate development with constant temperature for each 0.1?C:
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
  BrownShrimp::count_devDays(col, juvI_min, juvI_max) } )

dev_time_juvI <- dev_time_juvI[-c(1)]

#plot( (1:251)*0.1,  dev_time_juvI, main = "Development days Juvenile \n with min_size = 6 mm and max_size = 20 mm" )

#Fitting curve for JUVI
#first: delete all values = 1095:
indexes_1095 <- which(dev_time_juvI == 1095)
dev_time_JUVI_reduced <- dev_time_juvI[dev_time_juvI != 1095]
temperature_range_red_juvI <- temperature_range[-indexes_1095]

#plot( (1:length(dev_time_JUVI_reduced))*0.1,  dev_time_JUVI_reduced, main = "Development days Juvenile \n reduced data" )

fit_juvI <- nls( dev_time_JUVI_reduced ~ ( a + c*exp(-b*temperature_range_red_juvI) ),
                 start = list(a = 40, c = 3000, b = 1))

a_juvI <- coef(fit_juvI)[1]
b_juvI <- coef(fit_juvI)[3]
c_juvI <- coef(fit_juvI)[2]

#(optional) plot results form fit:
#plot(temperature_range, dev_time_juvI, pch = 1, col = 'gray28', las = 1, xlab = 'T [?C]', ylab = 'developent time [days]', cex.lab = 1.5, main = 'Class: Juv I')
#lines(temperature_range_red_juvI, predict(fit_juvI), col = 'lightseagreen', lwd = 2.5)

#### JUVII ####
development_df_juvII <- data.frame(row = seq(1, 1095))

#simulate development with constant temperature for each 0.1?C:
for (i in temperature_range){
  development <- c()
  initial_l <- juvII_min #mm
  for (j in time_range){
    new_legth <- som_growth(initial_l, i, 85, 'F')
    development <- append(development, new_legth)
    initial_l <- new_legth
  }
  development_df_juvII[as.character(i)] <- development
}

dev_time_juvII <- sapply(development_df_juvII, function(col){
  count_devDays(col, juvII_min, juvII_max) } )

dev_time_juvII <- dev_time_juvII[-c(1)] #first row doesn't make sence

#plot( (1:251)*0.1,  dev_time_juvII, main = paste("Development days Juvenile \n with min_size =", juvII_min , "mm and max_size = ", juvII_max ,"mm") )

#Fitting curve for JUVI
#first: delete all values = 1095:
indexes_1095 <- which(dev_time_juvII == 1095)
dev_time_JUVII_reduced <- dev_time_juvII[dev_time_juvII != 1095]
temperature_range_red_juvII <- temperature_range[-indexes_1095]
#length(temperature_range_red_juvI)
#length(dev_time_JUVI_reduced)
#plot( (1:length(dev_time_JUVII_reduced))*0.1,  dev_time_JUVII_reduced, main = "Development days Juvenile II \n reduced data", xlab = "T ?C", las = 1)

fit_juvII <- nls( dev_time_JUVII_reduced ~ ( a + c*exp(-b*temperature_range_red_juvII) ),
                  start = list(a = 20, c = 1500, b = 1))

a_juvII <- coef(fit_juvII)[1]
b_juvII <- coef(fit_juvII)[3]
c_juvII <- coef(fit_juvII)[2]

#(optional) plot results form fit:
#plot(temperature_range, dev_time_juvII, pch = 1, col = 'gray28', las = 1, xlab = 'T [?C]', ylab = 'developent time [days]', cex.lab = 1.5, main = 'Class: Juv II')
#lines(temperature_range_red_juvII, predict(fit_juvII), col = 'lightseagreen', lwd = 2.5)

#### JUVIII ####

development_df_juvIII <- data.frame(row = seq(1, 1095))

#simulate development with constant temperature for each 0.1?C:
for (i in temperature_range){
  development <- c()
  initial_l <- juvIII_min #mm
  for (j in time_range){
    new_legth <- som_growth(initial_l, i, 85, 'F')
    development <- append(development, new_legth)
    initial_l <- new_legth
  }
  development_df_juvIII[as.character(i)] <- development
}

dev_time_juvIII <- sapply(development_df_juvIII, function(col){
  count_devDays(col, juvIII_min, juvIII_max) } )

dev_time_juvIII <- dev_time_juvIII[-c(1)] #first row doesn't make sence

#plot( (1:251)*0.1,  dev_time_juvIII, main = paste("Development days Juvenile \n with min_size =", juvIII_min , "mm and max_size = ", juvIII_max ,"mm") )

#Fitting curve for JUVIII
#first: delete all values = 1095:
indexes_1095 <- which(dev_time_juvIII == 1095)
dev_time_juvIII_reduced <- dev_time_juvIII[dev_time_juvIII != 1095]
temperature_range_red_juvIII <- temperature_range[-indexes_1095]

#plot( (1:length(dev_time_juvIII_reduced))*0.1,  dev_time_juvIII_reduced, main = "Development days Juvenile II \n reduced data", xlab = "T ?C", las = 1)

fit_juvIII <- nls( dev_time_juvIII_reduced ~ ( a + c*exp(-b*temperature_range_red_juvIII) ),
                   start = list(a = 20, c = 1500, b = 1))

a_juvIII <- coef(fit_juvIII)[1]
b_juvIII <- coef(fit_juvIII)[3]
c_juvIII <- coef(fit_juvIII)[2]

#(optional) plot results form fit:
#plot(temperature_range, dev_time_juvIII, pch = 1, col = 'gray28', las = 1, xlab = 'T [°C]', ylab = 'developent time [days]', cex.lab = 1.5, main = 'Class: Juv III')
#lines(temperature_range_red_juvIII, predict(fit_juvIII), col = 'lightseagreen', lwd = 2.5)

#### JUVIV ####

development_df_juvIV <- data.frame(row = seq(1, 1095))

#simulate development with constant temperature for each 0.1?C:
for (i in temperature_range){
  development <- c()
  initial_l <- juvIV_min #mm
  for (j in time_range){
    new_legth <- som_growth(initial_l, i, 85, 'F')
    development <- append(development, new_legth)
    initial_l <- new_legth
  }
  development_df_juvIV[as.character(i)] <- development
}

dev_time_juvIV <- sapply(development_df_juvIV, function(col){
  count_devDays(col, juvIV_min, juvIV_max) } )

dev_time_juvIV <- dev_time_juvIV[-c(1)] #first row doesn't make sence

#plot( (1:251)*0.1,  dev_time_juvIV, main = paste("Development days Juvenile \n with min_size =", juvIV_min , "mm and max_size = ", juvIV_max ,"mm") )

#Fitting curve for juvIV
#first: delete all values = 1095:
indexes_1095 <- which(dev_time_juvIV == 1095)
dev_time_juvIV_reduced <- dev_time_juvIV[dev_time_juvIV != 1095]
temperature_range_red_juvIV <- temperature_range[-indexes_1095]

#plot( (1:length(dev_time_juvIV_reduced))*0.1,  dev_time_juvIV_reduced, main = "Development days Juvenile II \n reduced data", xlab = "T ?C", las = 1)

fit_juvIV <- nls( dev_time_juvIV_reduced ~ ( a + c*exp(-b*temperature_range_red_juvIV) ),
                  start = list(a = 20, c = 1500, b = 1))

a_juvIV <- coef(fit_juvIV)[1]
b_juvIV <- coef(fit_juvIV)[3]
c_juvIV <- coef(fit_juvIV)[2]

#(optional) plot results form fit:
#plot(temperature_range, dev_time_juvIV, pch = 1, col = 'gray28', las = 1, xlab = 'T [?C]', ylab = 'developent time [days]', cex.lab = 1.5, main = 'Class: Juv IV')
#lines(temperature_range_red_juvIV, predict(fit_juvIV), col = 'lightseagreen', lwd = 2.5)

#### juvV ####

development_df_juvV <- data.frame(row = seq(1, 1095))

#simulate development with constant temperature for each 0.1?C:
for (i in temperature_range){
  development <- c()
  initial_l <- juvV_min #mm
  for (j in time_range){
    new_legth <- som_growth(initial_l, i, 85, 'F')
    development <- append(development, new_legth)
    initial_l <- new_legth
  }
  development_df_juvV[as.character(i)] <- development
}

dev_time_juvV <- sapply(development_df_juvV, function(col){
  count_devDays(col, juvV_min, juvV_max) } )

dev_time_juvV <- dev_time_juvV[-c(1)] #first row doesn't make sence

#plot( (1:251)*0.1,  dev_time_juvV, main = paste("Development days Juvenile \n with min_size =", juvV_min , "mm and max_size = ", juvV_max ,"mm") )

#Fitting curve for juvV
#first: delete all values = 1095:
indexes_1095 <- which(dev_time_juvV == 1095)
dev_time_juvV_reduced <- dev_time_juvV[dev_time_juvV != 1095]
temperature_range_red_juvV <- temperature_range[-indexes_1095]

#plot( (1:length(dev_time_juvV_reduced))*0.1,  dev_time_juvV_reduced, main = "Development days Juvenile II \n reduced data", xlab = "T ?C", las = 1)

fit_juvV <- nls( dev_time_juvV_reduced ~ ( a + c*exp(-b*temperature_range_red_juvV) ),
                 start = list(a = 20, c = 1500, b = 1))

a_juvV <- coef(fit_juvV)[1]
b_juvV <- coef(fit_juvV)[3]
c_juvV <- coef(fit_juvV)[2]

#(optional) plot results form fit:
#plot(temperature_range, dev_time_juvV, pch = 1, col = 'gray28', las = 1, xlab = 'T [?C]', ylab = 'developent time [days]', cex.lab = 1.5, main = 'Class: Juv V')
#lines(temperature_range_red_juvV, predict(fit_juvV), col = 'lightseagreen', lwd = 2.5)


#### adultI ####

development_df_adultI <- data.frame(row = seq(1, 1095))

#simulate development with constant temperature for each 0.1?C:
for (i in temperature_range){
  development <- c()
  initial_l <- adultI_min #mm
  for (j in time_range){
    new_legth <- som_growth(initial_l, i, 85, 'F')
    development <- append(development, new_legth)
    initial_l <- new_legth
  }
  development_df_adultI[as.character(i)] <- development
}

dev_time_adultI <- sapply(development_df_adultI, function(col){
  count_devDays(col, adultI_min, adultI_max) } )

dev_time_adultI <- dev_time_adultI[-c(1)] #first row doesn't make sence

#plot( (1:251)*0.1,  dev_time_adultI, main = paste("Development days Juvenile \n with min_size =", adultI_min , "mm and max_size = ", adultI_max ,"mm") )

#Fitting curve for adultI
#first: delete all values = 1095:
indexes_1095 <- which(dev_time_adultI == 1095)
dev_time_adultI_reduced <- dev_time_adultI[dev_time_adultI != 1095]
temperature_range_red_adultI <- temperature_range[-indexes_1095]

#plot( temperature_range_red_adultI,  dev_time_adultI_reduced, main = "Development days Juvenile II \n reduced data", xlab = "T ?C", las = 1)
#view <- development_df_adultI['20.8']

#Failed to make a good fit for the U shape, temporary solution:
#ignore all T values from < 20 (because from 20?C dev. time rises again)

temperature_range_redII_adultI <- temperature_range_red_adultI[temperature_range_red_adultI <= 20]
dev_time_adultI_reducedII <- dev_time_adultI_reduced[1: length(temperature_range_redII_adultI)] #as I know that I want to unselect the last values

fitII_adultI <- nls( dev_time_adultI_reducedII ~ ( a + c*exp(-b*temperature_range_redII_adultI) ),
                     start = list(a = 100, c = 3000, b = 1))

a_adultI <- coef(fitII_adultI)[1]
b_adultI <- coef(fitII_adultI)[3]
c_adultI <- coef(fitII_adultI)[2]

#plot(temperature_range, dev_time_adultI, pch = 1, col = 'gray28', las = 1, xlab = 'T [°C]', ylab = 'developent time [days]', cex.lab = 1.5, main = 'Class: Adult I')
#lines(temperature_range_redII_adultI, predict(fitII_adultI), col = 'lightseagreen', lwd = 2.5)

#plot(temperature_range, 1/dev_time_adultI, pch = 1, col = 'gray28', las = 1, xlab = 'T [°C]', ylab = 'developent time [days]', cex.lab = 1.5, main = 'Class: Adult I')



#Shifts of Biomass function

shift_next_sizeClass = function(Te, sizeClass){
  ratio = 0
  if (Te <= 0) {
    return (0)
  } else {
    if (sizeClass == 'juvI'){
      a = 11.73753
      b = 2.205524
      c = 1747.899
      ratio = 1 / (a + c*exp(-b*Te))
    }
    if (sizeClass == 'juvII'){
      a = coef(fit_juvII)[1]
      b = coef(fit_juvII)[3]
      c = coef(fit_juvII)[2]
      ratio = 1 / (a + c*exp(-b*Te))
    }
    if (sizeClass == 'juvIII'){
      a = coef(fit_juvIII)[1]
      b = coef(fit_juvIII)[3]
      c = coef(fit_juvIII)[2]
      ratio = 1 / (a + c*exp(-b*Te))
    }
    if (sizeClass == 'juvIV'){
      a = coef(fit_juvIV)[1]
      b = coef(fit_juvIV)[3]
      c = coef(fit_juvIV)[2]
      ratio = 1 / (a + c*exp(-b*Te))
    }
    if (sizeClass == 'juvV'){
      a = coef(fit_juvV)[1]
      b = coef(fit_juvV)[3]
      c = coef(fit_juvV)[2]
      ratio = 1 / (a + c*exp(-b*Te))
    }
    if (sizeClass == 'adultI'){
      a = coef(fitII_adultI)[1]
      b = coef(fitII_adultI)[2]
      c = coef(fitII_adultI)[3]
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


##### SYSTEM ########

#the.method = 'rk4'
#h = 0.2 # half saturatuion constant original aprox .5
#L_inf = 8.5
#epsilon = 0.22
#m = 3
#const_c = 0.01 # reference density g/cm3

t = seq(0,length(temperature_10_11$temperature)-0.1, by = 0.1)
t_5years = seq(0,length(temperature_10_14$temperature)-0.1, by = 0.1)

state = c(P = 2, E = 0.1, L= 0.4 , J = .4 , J2 = 0.3, J3 = 0.3, J4 = 0.2, J5 = 0.2, A1 = 0.35, A2 = 0.1)
#state = c(P = 100, E = 2, L= 40 , J = 70 , J2 = 50, A2 = 40)
parameters = c()

### now test new solver from library BrownShrimp AND reduce juveniles on initial contiditions

state2 <- c(P = 2, E = 0.1, L= 0.4 , J = .175 , J2 = 0.125, J3 = 0.125, J4 = 0.1, J5 = 0.1, A1 = 0.35, A2 = 0.1)

start <- Sys.time()
cc_2 <- solver_sizeClass_extended(t = t_5years, state = state2, parameters = parameters, temperature_10_14)
print(Sys.time() - start)

#add dateTime column:
start_date <- as.POSIXct("2010-01-01", format="%Y-%m-%d", tz = "UTC")
cc_2 <- mutate(cc_2, dateTime = start_date + cc_2$time * 86400) #86400 seconds in one day
cc_2 <- cc_2[ , -1] # delete timesteps column

#cc_long_2 <- gather(cc_2, key = "Type", value = "Value", P, E, L, J, J2, J3, J4, J5, A1, A2)
cc_2_long <- cc_2 %>%
  pivot_longer(-dateTime, names_to = "variable", values_to = "value")

ggplot(cc_2_long, aes(x = dateTime, y = value, color = variable)) +
  geom_line(size = 1) +  # Plot lines
  labs(title = " extended size classes aprox. each 10 mm \n 2010-2014 ", x = "Days", y = "Biomass") +
  scale_color_manual(
    values = c( "P" = "#d9f0a3", "E"= "#1c9099" ,"L" = "gray"  ,"J" = "#bfd3e6", "J2" = "#9ebcda", "J3" = "#8c96c6", "J4" = "#8856a7", "J5" = "#810f7c" ,
                "A1" = "#fd8d3c" ,"A2" = "#e6550d" )) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), # Center title
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        legend.position="bottom") + ylim(0,20)


###### Barplots #######


cols <- sequential_hcl(4, palette = "BluYl")
months <- c('Jan', 'Feb', 'Mar', 'Apr', 'Mai', 'Jun', 'Jul', 'Aug', 'Sep', 'Okt', 'Nov', 'Dez')
#legend_labels <- c("Egg", "Larv", "Juv I", "Juv II", "Adu")
x_labels <- c("0", "< 6 ", "6 - 50", "50 - L_inf")
width_vector_compact <- c((2/10)*0.2, (5/10)*0.2, ((juvV_max - juvI_min)/10)*0.2, ((adultII_max - adultI_min)/10)*0.2 )

all_juv_classes <- rowSums(cc_2[,c(4:8)])
all_adult_classes <- rowSums(cc_2[,c(9:10)])

colnames_compact <- c("time", "Egg", "Larv", "Juv", "Adu")
cc_compact_2 <- data.frame(cc_2$dateTime, cc_2$E, cc_2$L, all_juv_classes, all_adult_classes)
names(cc_compact_2) <- colnames_compact


#only life stages as size classes / Mean monthly values
#dev.off()
par(mfrow = c(3,4))
j <- 1

for (i in seq(1, 12)) {
  init_day <- 304*(j-1)
  final_day <- 304*j
  month_data <- cc_compact_2[init_day:final_day, 2:5]
  toPlot <- colMeans(month_data)
  barplot(toPlot, col = cols, main = paste(months[j],' 2010'), names.arg = x_labels, #, las = 2
          cex.main = 1.6, cex.axis = 1.4, cex.names = 1.4 ,width = width_vector_compact, ylim = c(0,11))
  j <- j + 1
}


#only juveniles / Mean monthly values
width_vector_juvs <- c(((juvI_max - juvI_min)/10)*0.2,  ((juvII_max - juvII_min)/10)*0.2, ((juvIII_max - juvIII_min)/10)*0.2,
                       ((juvIV_max - juvIV_min)/10)*0.2, ((juvV_max - juvV_min)/10)*0.2)
x_labels <- c("6-10", "10-20", "20-30", "30-40", "40-50")
#dev.off()
par(mfrow = c(3,4))
j <- 1

for (i in seq(1, 12)) {
  init_day <- 304*(j-1)
  final_day <- 304*j
  month_data <- cc_2[init_day:final_day, 5:9 ]
  toPlot <- colMeans(month_data)
  barplot(toPlot, col = cols, main = paste(months[j],' 2010'), las = 1, cex.main = 1.5, names.arg = x_labels ,
          cex.axis = 1.4, cex.names = 1.4, width = width_vector_juvs, ylim = c(0,5))
  j <- j + 1
}


#all / quarterly mean values
cols <- sequential_hcl(9, palette = "viridis")

width_vector_all <- c((2/10)*0.1, (6/10)*0.1, ((juvI_max - juvI_min)/10)*0.1,  ((juvII_max - juvII_min)/10)*0.1, ((juvIII_max - juvIII_min)/10)*0.1,
                       ((juvIV_max - juvIV_min)/10)*0.1, ((juvV_max - juvV_min)/10)*0.1, ((adultI_max - adultI_min)/10)*0.1, ((adultII_max - adultII_min)/10)*0.1)

x_labels_vec <- c( "0", "<6", "6-10", "10-20", "20-30", "30-40", "40-50", "50-70", "70-L_inf")
#dev.off()
par(mfrow = c(2,2))
j <- 1

for (i in seq(1, 4)) {
  init_day <- 912*(j-1)
  final_day <- 912*j
  month_data <- cc_1[init_day:final_day, 3:11]
  toPlot <- colMeans(month_data)
  barplot(toPlot, col = cols, main = paste('Q ', j ,' 2010'), las = 1, cex.main = 2, names.arg = x_labels_vec  ,
          width = width_vector_all, cex.axis = 2, cex.names = 1.8 , ylim = c(0,3.5), space = 0)
  j <- j + 1
}


par(mfrow = c(2,2), oma = c(3.2,3.2,1,1))
j <- 1

x_labels_vec_11 <- c( "0", "3", "8", "15", "25", "35", "45", "60", "77.5")

for (i in seq(1, 4)) {
  init_day <- 3600 + 912*(j-1)
  final_day <- 3600 + 912*j
  month_data <- cc_1[init_day:final_day, 3:11]
  toPlot <- colMeans(month_data)
  toPlot <- toPlot/sum(toPlot)
  barplot(toPlot, col = cols, main = paste('Q ', j ,' 2011'), las = 1, cex.main = 2, names.arg = x_labels_vec_11  ,
          width = width_vector_all, cex.axis = 1.4, cex.names = 1.3 , ylim = c(0,0.35), space = 0)
  j <- j + 1
}

mtext("mean of size range", side = 1, outer = TRUE, line = 2, cex = 1.5)
mtext("density", side = 2, outer = TRUE, line = 2, cex = 1.5)





##### Stacked area chart ####

df_long2 <- cc_2[ ,-1] %>%
  pivot_longer(-dateTime, names_to = "variable", values_to = "value")


# Calculate total at each time point and normalize to proportions
#df_long$Date <- as.Date(df$time, origin = "01/01/2010")  # If time starts at 0 "31/12/2009"

#ensure that vertical order is same as life stages
stack_order <- rev(c(colnames(cc_2)[2:10]) ) #rev(c(names(df_long_prev)[2:10], 'sum_all') )

# Ensure 'variable' is a factor with levels in the same order
df_long2$variable <- factor(df_long2$variable, levels = stack_order)

#dev.off()
ggplot(df_long2, aes(x = dateTime, y = value, fill = variable)) +
  geom_area() +
  scale_fill_viridis_d() +
  labs(
    title = "",
    x = "Date",
    y = "Biomass",
    fill = "Size class"
  ) #+ ylim(0, 850)


#### Stacked area chart compact version ####
all_juv_classes_2 <- rowSums(cc_2[,c(3:7)])
all_adult_classes_2 <- rowSums(cc_2[,c(8:9)])

colnames_2_compact <- c('dateTime', 'egg', 'larvae', 'juvenile', 'adult')
cc_2_compact <- data.frame(cc_2$dateTime, cc_2$E, cc_2$L, all_juv_classes_2, all_adult_classes_2)
names(cc_2_compact) <- colnames_2_compact

cc_2_compact_long <- cc_2_compact %>%
  pivot_longer(-dateTime, names_to = "variable", values_to = "value")

df_prop_2 <- cc_2_compact_long %>%
  group_by(dateTime) %>%
  mutate(prop = value / sum(value))

#Ensure same vertical order as life stages
stack_order <- rev(c('egg', 'larvae', 'juvenile', 'adult') )
df_prop_2$variable <- factor(df_prop_2$variable, levels = stack_order)


ggplot(df_prop_2, aes(x = dateTime, y = prop, fill = variable)) +
  geom_area() +
  scale_y_continuous(labels = scales::percent) +
  labs(
    title = "Proportional Stacked Area Chart \n with Temperature data 2015 and 2016",
    x = "Days",
    y = "Proportion",
    fill = "Category"
  ) +
  theme_minimal()


##### Proportional stacked area chart (not convinient)####
df_long <- cc_1 %>%
  pivot_longer(-time, names_to = "variable", values_to = "value")

# Calculate total at each time point and normalize to proportions
df_prop <- df_long %>%
  group_by(time) %>%
  mutate(prop = value / sum(value))


ggplot(df_prop, aes(x = time, y = prop, fill = variable)) +
  geom_area() +
  scale_y_continuous(labels = scales::percent) +
  labs(
    title = "Proportional Stacked Area Chart",
    x = "Time",
    y = "Proportion",
    fill = "Category"
  ) +
  theme_minimal()



####################################################
###                                             ####
###      DETERMINE G(W,T)  CONTINUOUS           ####
###                                             ####
####################################################

#Here I get one step closer to a continuous size spectrum model by determining a growth function
#that goes over a 'continuous size classes' and describe the shifts of biomass between them
temp_range <- seq(1,25, 0.1)
size_classes <- c('egg', 'larv', 'juvI', 'juvII', 'juvIII', 'juvIV', 'juvV', 'adultI')
sizes <- c( 0,0, (juvI_min+juvI_max)/2 , (juvII_min+juvII_max)/2 /10, (juvIII_min+juvIII_max)/2 /10, (juvIV_min+juvIV_max)/2 /10,
            (juvV_min+juvV_max)/2 /10, (adultI_min + adultI_max)/2 /10, (adultII_min + adultII_max)/2 /10 )

#temp_range <- seq(0,25, 0.1) #already defined above
growth <- data.frame(nrwo = c(1:length(temp_range)))

for (i in size_classes){
  g <- sapply(i, function(sc) shift_next_sizeClass( temp_range, i) )
  growth[as.character(i)] <- g
}

head(growth)

par(mfrow = c(3,3))
for (i in size_classes){
  plot(growth$nrwo*0.1, growth[[i]], main = i, las = 1, xlab = 'T', ylab = 'G(class, T)' , xlim = c(0,15), ylim = c(0, 0.085))
}

par(mfrow = c(2,2))
for (i in seq(50,201, 50) ) {
  plot(1:length(size_classes), growth[i, 2:(length(size_classes)+1)], main = paste(i/10, '°C'), xlab = "", ylab="G", xaxt = "n" )
  axis(1, at=1:8, labels = size_classes, las = 2)
}


plot(1:length(size_classes), growth[50, 2:(length(size_classes)+1)], main = paste(5, '°C'))
