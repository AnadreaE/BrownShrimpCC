######### START ###########

#Def. Size Classes ranges [cm]

Linf_F <- 8.5
Linf_M <- 5.5

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
adultI_max <- 6.0

adultII_min <- 6.0
adultII_max <- 7.0

adultIII_min <- 7.0
adultIII_max <- Linf_F

#Define mean sizes for each range [cm]
LL = 0.3
LJ = (juvI_min+juvI_max)/2
LJ2 = (juvII_min+juvII_max)/2
LJ3 = (juvIII_min+juvIII_max)/2
LJ4 = (juvIV_min+juvIV_max)/2
LJ5 = (juvV_min+juvV_max)/2
LA1 = (adultI_min + adultI_max)/2
LA2 = (adultII_min + adultII_max)/2
LA3 = (adultIII_min + adultIII_max)/2

sizeClass_names = c("LJ", "LJ2", "LJ3", "LJ4", "LJ5", "LA1", "LA2", "LA3" )
sizeClass_means = c(LJ, LJ2, LJ3, LJ4, LJ5, LA1, LA2, LA3 )

#### Def. growth function (as function of time) ####

lengthGrowth_funcTime = function (L0, time, temp, sex){ #L[cm]
  #This funciton works only in this context i.e. to simulate growth with constant temperatures in order to
  #compare and check plausibilty of the results.
  L0 = L0*10
  if (sex == "F") Linf = 85
  if (sex == "M") Linf = 55
  if (L0 > Linf ) growth = 0
  else growth = Linf- (Linf-L0)*exp(-K_func_briere(temp)*3*time)
  return(growth/10 )
}

#### Simulations ####
time_range <- seq(1,1095)
temperature_range <- seq(0,30, 0.1) #(10,15)

k_vals = sapply(temperature_range, K_func_briere)

#### JUVI ####

development.df_juvI <- data.frame(row = seq(1, 1095))

#simulate development with constant temperature for each 0.1°C & store results columns wise in the df:
for (i in temperature_range){
  development <- c()
  initial_l <- 0.6 #cm
  for (j in time_range){
    growth = lengthGrowth_funcTime(initial_l, j ,i , "F")
    #initial_l = growth #+  L_init
    development = append(development, growth)
  }
  development.df_juvI[as.character(i)] <- development
}

development.df_juvI = development.df_juvI[ , !(names(development.df_juvI) == "row")] #delete first column with timesteps, same as index

#Count how many time steps there are, where individual's size is between min and max size of its size class range:
dev_time_juvI <- sapply(development.df_juvI, function(col){
  BrownShrimp::count_devDays(col, juvI_min, juvI_max) } )

#dev_time_juvI <- dev_time_juvI[-c(1)] #delete first element cuz doesn't make sense
#dev.off()


#With following plot, we know the shape of the curve (U shape. But inverted is again TPC)
plot( (1:length(dev_time_juvI))*0.1, dev_time_juvI, main = "Development days Juvenile I \n min_size = 0.6 cm and max_size = 1 cm" )
plot( (1:length(dev_time_juvI))*0.1, 1/dev_time_juvI, main = "1/Development days Juvenile I \n min_size = 0.6 cm and max_size = 1 cm" )
lines(temperature_range, k_vals, col = 'red', ldw = 2)
legend("topright", legend = c('1/devTime', 'K_func_briere'), fill = c("black", "red"))

#Deleet values where there is no development i.e. when devTime = 1095
indexes_1095 <- which(dev_time_juvI == 1095)
dev_time_JUVI_reduced <- dev_time_juvI[dev_time_juvI != 1095]
temperature_range_red_juvI <- temperature_range[-indexes_1095]
k_vals_reduced_JUVI = k_vals[-indexes_1095]

devTime_juvI_inv_reduced = 1/dev_time_JUVI_reduced
k_vals_reduced_JUVI_inv = 1/k_vals_reduced_JUVI

fit_JuvI = lm(dev_time_JUVI_reduced ~ k_vals_reduced_JUVI_inv )
summary(fit_JuvI)

plot( (1:length(dev_time_juvI))*0.1, dev_time_juvI, main = "Development days Juvenile vs fit I \n min_size = 6 mm and max_size = 20 mm" )
lines(temperature_range_red_juvI, predict(fit_JuvI), col = 'red4', lwd = 2)
legend("topright", legend = c('1/devTime', 'K_func_briere'), fill = c("black", "red4"))

#### JUV II ####

development.df_juvII <- data.frame(row = seq(1, 1095))

#simulate development with constant temperature for each 0.1°C & store results columns wise in the df:
for (i in temperature_range){
  development <- c()
  initial_l <- 1 #cm
  for (j in time_range){
    growth = lengthGrowth_funcTime(initial_l, j ,i , "F")
    #initial_l = growth #+  L_init
    development = append(development, growth)
  }
  development.df_juvII[as.character(i)] <- development
}

development.df_juvII = development.df_juvII[ , !(names(development.df_juvII) == "row")] #delete first column with timesteps, same as index

#Count how many time steps there are, where individual's size is between min and max size of its size class range:
dev_time_juvII <- sapply(development.df_juvII, function(col){
  BrownShrimp::count_devDays(col, juvII_min, juvII_max) } )


#With following plot, we know the shape of the curve (U shape. But inverted is again TPC)
#plot( (1:length(dev_time_juvII))*0.1, dev_time_juvII, main = "Development days Juvenile I \n min_size = 0.6 cm and max_size = 1 cm" )
plot( (1:length(dev_time_juvII))*0.1, 1/dev_time_juvII, main = "1/Development days Juvenile II \n min_size = 1 cm and max_size = 2 cm" )
lines(temperature_range, k_vals, col = 'red', lwd = 2)
legend("topright", legend = c('1/devTime', 'K_func_briere'), fill = c("black", "red"))

#Deleet values where there is no development i.e. when devTime = 1095
indexes_1095 <- which(dev_time_juvII == 1095)
dev_time_juvII_reduced <- dev_time_juvII[dev_time_juvII != 1095]
temperature_range_red_juvII <- temperature_range[-indexes_1095]
k_vals_reduced_juvII = k_vals[-indexes_1095]

devTime_juvII_inv_reduced = 1/dev_time_juvII_reduced
k_vals_reduced_juvII_inv = 1/k_vals_reduced_juvII

fit_juvII = lm(dev_time_juvII_reduced ~ k_vals_reduced_juvII_inv )
summary(fit_juvII)

plot( (1:length(dev_time_juvII))*0.1, dev_time_juvII, main = "Development days Juvenile II vs fit \n min_size = 1 cm and max_size = 2 cm" )
lines(temperature_range_red_juvII, predict(fit_juvII), col = 'red4', lwd = 2)
legend("topright", legend = c('1/devTime', 'K_func_briere'), fill = c("black", "red4"))


#### JUV III ####

development.df_juvIII <- data.frame(row = seq(1, 1095))

#simulate development with constant temperature for each 0.1°C & store results columns wise in the df:
for (i in temperature_range){
  development <- c()
  initial_l <- 2 #cm
  for (j in time_range){
    growth = lengthGrowth_funcTime(initial_l, j ,i , "F")
    #initial_l = growth #+  L_init
    development = append(development, growth)
  }
  development.df_juvIII[as.character(i)] <- development
}

development.df_juvIII = development.df_juvIII[ , !(names(development.df_juvIII) == "row")] #delete first column with timesteps, same as index

#Count how many time steps there are, where individual's size is between min and max size of its size class range:
dev_time_juvIII <- sapply(development.df_juvIII, function(col){
  BrownShrimp::count_devDays(col, juvIII_min, juvIII_max) } )


#With following plot, we know the shape of the curve (U shape. But inverted is again TPC)
plot( (1:length(dev_time_juvIII))*0.1, dev_time_juvIII, main = "Development days Juvenile III \n min_size = 0.6 cm and max_size = 1 cm" )
plot( (1:length(dev_time_juvIII))*0.1, 1/dev_time_juvIII, main = "1/Development days Juvenile III \n min_size = 2 cm and max_size = 3 cm" )
lines(temperature_range, k_vals, col = 'red', lwd = 2)
legend("topright", legend = c('1/devTime', 'K_func_briere'), fill = c("black", "red"))

#Deleet values where there is no development i.e. when devTime = 1095
indexes_1095 <- which(dev_time_juvIII == 1095)
dev_time_juvIII_reduced <- dev_time_juvIII[dev_time_juvIII != 1095]
temperature_range_red_juvIII <- temperature_range[-indexes_1095]
k_vals_reduced_juvIII = k_vals[-indexes_1095]

devTime_juvIII_inv_reduced = 1/dev_time_juvIII_reduced
k_vals_reduced_juvIII_inv = 1/k_vals_reduced_juvIII

fit_juvIII = lm(dev_time_juvIII_reduced ~ k_vals_reduced_juvIII_inv )
summary(fit_juvIII)

plot( (1:length(dev_time_juvIII))*0.1, dev_time_juvIII, main = "Development days Juvenile III vs fit \n min_size = 2 cm and max_size = 3 cm" )
lines(temperature_range_red_juvIII, predict(fit_juvIII), col = 'red4', lwd = 2)
legend("topright", legend = c('1/devTime', 'K_func_briere'), fill = c("black", "red4"))



#### JUV IV ####

development.df_juvIV <- data.frame(row = seq(1, 1095))

#simulate development with constant temperature for each 0.1°C & store results columns wise in the df:
for (i in temperature_range){
  development <- c()
  initial_l <- 3 #cm
  for (j in time_range){
    growth = lengthGrowth_funcTime(initial_l, j ,i , "F")
    #initial_l = growth #+  L_init
    development = append(development, growth)
  }
  development.df_juvIV[as.character(i)] <- development
}

development.df_juvIV = development.df_juvIV[ , !(names(development.df_juvIV) == "row")] #delete first column with timesteps, same as index

#Count how many time steps there are, where individual's size is between min and max size of its size class range:
dev_time_juvIV <- sapply(development.df_juvIV, function(col){
  BrownShrimp::count_devDays(col, juvIV_min, juvIV_max) } )


#With following plot, we know the shape of the curve (U shape. But inverted is again TPC)
plot( (1:length(dev_time_juvIV))*0.1, dev_time_juvIV, main = "Development days Juvenile IV \n min_size = 3 cm and max_size = 4 cm" )
plot( (1:length(dev_time_juvIV))*0.1, 1/dev_time_juvIV, main = "1/Development days Juvenile IV \n min_size = 3 cm and max_size = 4 cm" )
lines(temperature_range, k_vals, col = 'red', lwd = 2)
legend("topright", legend = c('1/devTime', 'K_func_briere'), fill = c("black", "red"))

#Deleet values where there is no development i.e. when devTime = 1095
indexes_1095 <- which(dev_time_juvIV == 1095)
dev_time_juvIV_reduced <- dev_time_juvIV[dev_time_juvIV != 1095]
temperature_range_red_juvIV <- temperature_range[-indexes_1095]
k_vals_reduced_juvIV = k_vals[-indexes_1095]

devTime_juvIV_inv_reduced = 1/dev_time_juvIV_reduced
k_vals_reduced_juvIV_inv = 1/k_vals_reduced_juvIV

fit_juvIV = lm(dev_time_juvIV_reduced ~ k_vals_reduced_juvIV_inv )
summary(fit_juvIV)

plot( (1:length(dev_time_juvIV))*0.1, dev_time_juvIV, main = "Development days Juvenile IV vs fit \n min_size = 3 cm and max_size = 4 cm" )
lines(temperature_range_red_juvIV, predict(fit_juvIV), col = 'red4', lwd = 2)
legend("topright", legend = c('1/devTime', 'K_func_briere'), fill = c("black", "red4"))

#### JUV V ####

development.df_juvV <- data.frame(row = seq(1, 1095))

#simulate development with constant temperature for each 0.1°C & store results columns wise in the df:
for (i in temperature_range){
  development <- c()
  initial_l <- 4 #cm
  for (j in time_range){
    growth = lengthGrowth_funcTime(initial_l, j ,i , "F")
    #initial_l = growth #+  L_init
    development = append(development, growth)
  }
  development.df_juvV[as.character(i)] <- development
}

development.df_juvV = development.df_juvV[ , !(names(development.df_juvV) == "row")] #delete first column with timesteps, same as index

#Count how many time steps there are, where individual's size is between min and max size of its size class range:
dev_time_juvV <- sapply(development.df_juvV, function(col){
  BrownShrimp::count_devDays(col, juvV_min, juvV_max) } )


#With following plot, we know the shape of the curve (U shape. But inverted is again TPC)
plot( (1:length(dev_time_juvV))*0.1, dev_time_juvV, main = "Development days Juvenile V \n min_size = 4 cm and max_size = 5 cm" )
plot( (1:length(dev_time_juvV))*0.1, 1/dev_time_juvV, main = "1/Development days Juvenile V \n min_size = 4 cm and max_size = 5 cm" )
lines(temperature_range, k_vals, col = 'red', lwd = 2)
legend("topright", legend = c('1/devTime', 'K_func_briere'), fill = c("black", "red"))

#Deleet values where there is no development i.e. when devTime = 1095
indexes_1095 <- which(dev_time_juvV == 1095)
dev_time_juvV_reduced <- dev_time_juvV[dev_time_juvV != 1095]
temperature_range_red_juvV <- temperature_range[-indexes_1095]
k_vals_reduced_juvV = k_vals[-indexes_1095]

devTime_juvV_inv_reduced = 1/dev_time_juvV_reduced
k_vals_reduced_juvV_inv = 1/k_vals_reduced_juvV

fit_juvV = lm(dev_time_juvV_reduced ~ k_vals_reduced_juvV_inv )
summary(fit_juvV)

plot( (1:length(dev_time_juvV))*0.1, dev_time_juvV, main = "Development days Juvenile V vs fit \n min_size = 4 cm and max_size = 5 cm" )
lines(temperature_range_red_juvV, predict(fit_juvV), col = 'red4', lwd = 2)
legend("topright", legend = c('1/devTime', 'K_func_briere'), fill = c("black", "red4"))


#### ADULT I####

development.df_aduI <- data.frame(row = seq(1, 1095))

#simulate development with constant temperature for each 0.1°C & store results columns wise in the df:
for (i in temperature_range){
  development <- c()
  initial_l <- 5 #cm
  for (j in time_range){
    growth = lengthGrowth_funcTime(initial_l, j ,i , "F")
    #initial_l = growth #+  L_init
    development = append(development, growth)
  }
  development.df_aduI[as.character(i)] <- development
}

development.df_aduI = development.df_aduI[ , !(names(development.df_aduI) == "row")] #delete first column with timesteps, same as index

#Count how many time steps there are, where individual's size is between min and max size of its size class range:
dev_time_aduI <- sapply(development.df_aduI, function(col){
  BrownShrimp::count_devDays(col, adultI_min, adultI_max) } )


#With following plot, we know the shape of the curve (U shape. But inverted is again TPC)
plot( (1:length(dev_time_aduI))*0.1, dev_time_aduI, main = "Development days Adult I \n min_size = 5 cm and max_size = 6 cm" )
plot( (1:length(dev_time_aduI))*0.1, 1/dev_time_aduI, main = "1/Development days Adult I \n min_size = 5 cm and max_size = 6 cm" )
lines(temperature_range, k_vals, col = 'red', lwd = 2)
legend("topright", legend = c('1/devTime', 'K_func_briere'), fill = c("black", "red"))

#Deleet values where there is no development i.e. when devTime = 1095
indexes_1095 <- which(dev_time_aduI == 1095)
dev_time_aduI_reduced <- dev_time_aduI[dev_time_aduI != 1095]
temperature_range_red_aduI <- temperature_range[-indexes_1095]
k_vals_reduced_aduI = k_vals[-indexes_1095]

devTime_aduI_inv_reduced = 1/dev_time_aduI_reduced
k_vals_reduced_aduI_inv = 1/k_vals_reduced_aduI

fit_aduI = lm(dev_time_aduI_reduced ~ k_vals_reduced_aduI_inv )
summary(fit_aduI)

plot( (1:length(dev_time_aduI))*0.1, dev_time_aduI, main = "Development days Adult I vs fit \n min_size = 5 cm and max_size = 6 cm" )
lines(temperature_range_red_aduI, predict(fit_aduI), col = 'red4', lwd = 2)
legend("topright", legend = c('1/devTime', 'K_func_briere'), fill = c("black", "red4"))


#### ADULT II####

development.df_aduII <- data.frame(row = seq(1, 1095))

#simulate development with constant temperature for each 0.1°C & store results columns wise in the df:
for (i in temperature_range){
  development <- c()
  initial_l <- 6 #cm
  for (j in time_range){
    growth = lengthGrowth_funcTime(initial_l, j ,i , "F")
    #initial_l = growth #+  L_init
    development = append(development, growth)
  }
  development.df_aduII[as.character(i)] <- development
}

development.df_aduII = development.df_aduII[ , !(names(development.df_aduII) == "row")] #delete first column with timesteps, same as index

#Count how many time steps there are, where individual's size is between min and max size of its size class range:
dev_time_aduII <- sapply(development.df_aduII, function(col){
  BrownShrimp::count_devDays(col, adultII_min, adultII_max) } )


#With following plot, we know the shape of the curve (U shape. But inverted is again TPC)
plot( (1:length(dev_time_aduII))*0.1, dev_time_aduII, main = "Development days Adult II \n min_size = 6 cm and max_size = 7 cm" )
plot( (1:length(dev_time_aduII))*0.1, 1/dev_time_aduII, main = "1/Development days Adult II \n min_size = 6 cm and max_size = 7 cm" )
lines(temperature_range, k_vals, col = 'red', lwd = 2)
legend("topright", legend = c('1/devTime', 'K_func_briere'), fill = c("black", "red"))

#Deleet values where there is no development i.e. when devTime = 1095
indexes_1095 <- which(dev_time_aduII == 1095)
dev_time_aduII_reduced <- dev_time_aduII[dev_time_aduII != 1095]
temperature_range_red_aduII <- temperature_range[-indexes_1095]
k_vals_reduced_aduII = k_vals[-indexes_1095]

devTime_aduII_inv_reduced = 1/dev_time_aduII_reduced
k_vals_reduced_aduII_inv = 1/k_vals_reduced_aduII

fit_aduII = lm(dev_time_aduII_reduced ~ k_vals_reduced_aduII_inv )
summary(fit_aduII)

plot( (1:length(dev_time_aduII))*0.1, dev_time_aduII, main = "Development days Adult II vs fit \n min_size = 6 cm and max_size = 7 cm" )
lines(temperature_range_red_aduII, predict(fit_aduII), col = 'red4', lwd = 2)
legend("topright", legend = c('1/devTime', 'K_func_briere'), fill = c("black", "red4"))


#### FITS SUMMARY  ####
intercepts = c(coef(fit_JuvI)[1], coef(fit_juvII)[1], coef(fit_juvIII)[1], coef(fit_juvIV)[1], coef(fit_juvV)[1], coef(fit_aduI)[1], coef(fit_aduII)[1])
factors = c(coef(fit_JuvI)[2], coef(fit_juvII)[2], coef(fit_juvIII)[2], coef(fit_juvIV)[2], coef(fit_juvV)[2], coef(fit_aduI)[2], coef(fit_aduII)[2])

sizeClass_means_reduced = sizeClass_means[1:length(sizeClass_means)-1]

fits_df = data.frame( row.names = sizeClass_names[1:length(sizeClass_names)-1], intercept = intercepts, factor = factors)

#check weather INTERCEPTS show any relation with the mean size of the class:

#The intercept doesn't seem to follow any relytion with the mean of the size class....
#we may work only with the mean of all itnercepts as the range is not that big:
mean_intercept = mean(fits_df$intercept)

#OR linear regression
fit_intercepts = lm(fits_df$intercept ~ sizeClass_means_reduced)
summary(fit_intercepts)

plot(1:(length(sizeClass_names)-1), fits_df$intercept)
abline(h= mean(fits_df$intercept), lty=2 )
lines(1:(length(sizeClass_names)-1), predict(fit_intercepts), col = 'orange3')


#check weather FACTORS show any relation with the mean size of the class:
plot(1:length(sizeClass_means_reduced), fits_df$factor, col = 'red', ylim = c(-0.05, 0.8))
lines( 1:length(sizeClass_means_reduced), sizeClass_means_reduced/10)

plot(1:length(sizeClass_means_reduced), fits_df$factor)
lines(1:length(sizeClass_means_reduced), predict(fit_factorts), col = 'lightblue4', lwd = 2)

#fit_factorts = nls(fits_df$factor ~ a * exp(b * sizeClass_means_reduced), start = list(a = 0.05, b = 0.2))
logx = log(sizeClass_means_reduced)
logy = log(fits_df$factor)


fit_factortsLog = lm(logy ~ logx)

summary(fit_factortsLog)

log_a <- coef(fit_factortsLog)[1]
b_est <- coef(fit_factortsLog)[2]
a_est <- exp(log_a)

plot(sizeClass_means_reduced, fits_df$factor, log = 'xy')
lines(sizeClass_means_reduced, a_est * sizeClass_means_reduced^b_est, col = 'red')

lines(logx, log_a  +  logx*b_est, col = 'red')

plot(sizeClass_means_reduced, fits_df$factor)
lines(sizeClass_means_reduced, a_est * sizeClass_means_reduced^b_est, col = 'red')

#### GROWTH FUNCTION OPTION 1 ####
#here we use the intercept and factor as funciton of L mean

#check 2 Juv classes and 1 adult
growth_optionA_t = function(L_mean, temperature, sex){
  intercept = -0.4968367 #mean_intercept
  factor = a_est * L_mean^b_est #
  k = K_func_briere(temperature)
  return(1/ (intercept + factor*(1/k)) )
}



par(mfrow = c(3, 3), mar = c(4, 4, 2, 1))
ind = 1
for (i in sizeClass_names){
  plot(temperature_range,  1 / (fits_df[i, 'intercept'] + fits_df[i, 'factor']*(1/k_vals) ),
       main = paste('growth', i ,'vs function'), ylim = c(0, 0.6), ylab = "growth rate", las = 1, cex.main = 1.4, cex.lab = 1.3, cex.axis = 1.5 )
  lines(temperature_range, sapply(temperature_range, growth_optionA_t, L_mean = sizeClass_means[ind], sex='F' ) , col = 'red3' , lwd = 1.8)
  ind = ind + 1
}

plot.new()  # empty plot
legend("center", legend = c('1/devTime', 'growth function'), fill = c("black", "red3"))

