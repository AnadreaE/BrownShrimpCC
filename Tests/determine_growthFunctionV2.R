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

time_range <- seq(1,1095)
temperature_range <- seq(0,30, 0.1) #(10,15)

#### Simulations FEMALES ####

k_vals = sapply(temperature_range, K_func_briere, sex = 'F')

#### JUVI ####

development.df_juvI <- data.frame(row = seq(1, 1095))

#simulate development with constant temperature for each 0.1°C & store results columns wise in the df:
for (i in temperature_range){
  development <- c()
  initial_l <- 0.6 #cm
  for (j in time_range){
    growth = som_growth_thesis(initial_l, j ,i , "F")
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
lines(temperature_range, k_vals, col = 'red', lwd = 2)
legend("topright", legend = c('1/devTime', 'K_func_briere'), fill = c("black", "red"))

#Deleet values where there is no development i.e. when devTime = 1095
indexes_1095 <- which(dev_time_juvI == 1095 | dev_time_juvI == 0)
dev_time_JUVI_reduced <- dev_time_juvI[dev_time_juvI != 1095 & dev_time_juvI != 0]
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
    growth = som_growth_thesis(initial_l, j ,i , "F")
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
plot( (1:length(dev_time_juvII))*0.1, dev_time_juvII, main = "Development days Juvenile I \n min_size = 0.6 cm and max_size = 1 cm" )
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
    growth = som_growth_thesis(initial_l, j ,i , "F")
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
    growth = som_growth_thesis(initial_l, j ,i , "F")
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
    growth = som_growth_thesis(initial_l, j ,i , "F")
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
    growth = som_growth_thesis(initial_l, j ,i , "F")
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
    growth = som_growth_thesis(initial_l, j ,i , "F")
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

plot(sizeClass_means_reduced, fits_df$factor)
lines(sizeClass_means_reduced, a_est * sizeClass_means_reduced^b_est, col = 'red')

#### GROWTH FUNCTION OPTION 1 ####
#here we use the factor as funciton of L mean and only the mean of the intercept:

growth_optionA_t = function(L_mean, temperature, sex){
  intercept = -0.4968367 #mean_intercept
  factor = a_est * L_mean^b_est #
  k = K_func_briere(temperature, sex)
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


#### GROWTH FUNCTION OPTION 2 ####
#here we use both factor and intercept as funciton of L mean:

growth_optionB_t = function(L_mean, temperature, sex){
  intercept = -0.530907 + L_mean*0.007645 #mean_intercept
  factor = a_est * L_mean^b_est #
  k = K_func_briere(temperature, sex)
  return(1/ (intercept + factor*(1/k)) )
}



par(mfrow = c(3, 3), mar = c(4, 4, 2, 1))
ind = 1
for (i in sizeClass_names){
  plot(temperature_range,  1 / (fits_df[i, 'intercept'] + fits_df[i, 'factor']*(1/k_vals) ),
       main = paste('growth', i ,'vs function'), ylim = c(0, 0.6), ylab = "growth rate", las = 1, cex.main = 1.4, cex.lab = 1.3, cex.axis = 1.5 )
  lines(temperature_range, sapply(temperature_range, growth_optionB_t, L_mean = sizeClass_means[ind], sex='F' ) , col = 'red3' , lwd = 1.8)
  ind = ind + 1
}

plot.new()  # empty plot
legend("center", legend = c('1/devTime', 'growth function'), fill = c("black", "red3"))

#CONCLUSION COMPARATION A and B:
#No major difference, therfore so simplicity and saving computing effort I choose option A



#### Simulations MALES ####

k_vals_M = sapply(temperature_range, K_func_briere, sex = 'M')

#### M JUVI  ####

Mdevelopment.df_juvI <- data.frame(row = seq(1, 1095))

#simulate development with constant temperature for each 0.1°C & store results columns wise in the df:
for (i in temperature_range){
  development <- c()
  initial_l <- 0.6 #cm
  for (j in time_range){
    growth = som_growth_thesis(initial_l, j ,i , "M")
    #initial_l = growth #+  L_init
    development = append(development, growth)
  }
  Mdevelopment.df_juvI[as.character(i)] <- development
}

Mdevelopment.df_juvI = Mdevelopment.df_juvI[ , !(names(Mdevelopment.df_juvI) == "row")] #delete first column with timesteps, same as index

#Count how many time steps there are, where individual's size is between min and max size of its size class range:
dev_time_MjuvI <- sapply(Mdevelopment.df_juvI, function(col){
  BrownShrimp::count_devDays(col, juvI_min, juvI_max) } )

#dev_time_juvI <- dev_time_juvI[-c(1)] #delete first element cuz doesn't make sense
#dev.off()


#With following plot, we know the shape of the curve (U shape. But inverted is again TPC)
plot( (1:length(dev_time_MjuvI))*0.1, dev_time_MjuvI, main = "Development days Male Juvenile I \n min_size = 0.6 cm and max_size = 1 cm" )
plot( (1:length(dev_time_MjuvI))*0.1, 1/dev_time_MjuvI, main = "1/Development days Juvenile I \n min_size = 0.6 cm and max_size = 1 cm" )
lines(temperature_range, k_vals_M, col = '#1c9099', lwd = 2)
legend("topright", legend = c('1/devTime', 'K_func_briere'), fill = c("grey48", "#1c9099"))

#Deleet values where there is no development i.e. when devTime = 1095
indexes_1095 <- which(dev_time_MjuvI == 1095 | dev_time_MjuvI == 0)
dev_time_MJUVI_reduced <- dev_time_MjuvI[dev_time_MjuvI != 1095 & dev_time_MjuvI != 0]
temperature_range_red_MjuvI <- temperature_range[-indexes_1095]
k_vals_reduced_MJUVI = k_vals_M[-indexes_1095]

devTime_MjuvI_inv_reduced = 1/dev_time_MJUVI_reduced
k_vals_reduced_MJUVI_inv = 1/k_vals_reduced_MJUVI

fit_M.JuvI = lm(dev_time_MJUVI_reduced ~ k_vals_reduced_MJUVI_inv )
summary(fit_M.JuvI)

plot( (1:length(dev_time_MjuvI))*0.1, dev_time_MjuvI, main = "Development days Juvenile M vs fit I \n min_size = 0.6 cm and max_size = 1 cm" )
lines(temperature_range_red_MjuvI, predict(fit_M.JuvI), col = '#1c9099', lwd = 2)
legend("topright", legend = c('1/devTime', 'K_func_briere'), fill = c("grey48", "#1c9099"))

#### M JUV II ####

development.df_MjuvII <- data.frame(row = seq(1, 1095))

#simulate development with constant temperature for each 0.1°C & store results columns wise in the df:
for (i in temperature_range){
  development <- c()
  initial_l <- 1 #cm
  for (j in time_range){
    growth = som_growth_thesis(initial_l, j ,i , "M")
    #initial_l = growth #+  L_init
    development = append(development, growth)
  }
  development.df_MjuvII[as.character(i)] <- development
}

development.df_MjuvII = development.df_MjuvII[ , !(names(development.df_MjuvII) == "row")] #delete first column with timesteps, same as index

#Count how many time steps there are, where individual's size is between min and max size of its size class range:
dev_time_MjuvII <- sapply(development.df_MjuvII, function(col){
  BrownShrimp::count_devDays(col, juvII_min, juvII_max) } )


#With following plot, we know the shape of the curve (U shape. But inverted is again TPC)
plot( (1:length(dev_time_MjuvII))*0.1, dev_time_MjuvII, main = "Development days Juvenile II \n min_size = 1 cm and max_size = 2 cm" )
plot( (1:length(dev_time_MjuvII))*0.1, 1/dev_time_MjuvII, main = "1/Development days Juvenile II \n min_size = 1 cm and max_size = 2 cm" )
lines(temperature_range, k_vals_M, col = '#1c9099', lwd = 2)
legend("topright", legend = c('1/devTime', 'K_func_briere'), fill = c("grey48", "#1c9099"))

#Deleet values where there is no development i.e. when devTime = 1095
indexes_1095 <- which(dev_time_MjuvII == 1095)
dev_time_MjuvII_reduced <- dev_time_MjuvII[dev_time_MjuvII != 1095]
temperature_range_red_MjuvII <- temperature_range[-indexes_1095]
k_vals_reduced_MjuvII = k_vals_M[-indexes_1095]

devTime_MjuvII_inv_reduced = 1/dev_time_MjuvII_reduced
k_vals_reduced_MjuvII_inv = 1/k_vals_reduced_MjuvII

fit_MjuvII = lm(dev_time_MjuvII_reduced ~ k_vals_reduced_MjuvII_inv )
summary(fit_MjuvII)

plot( (1:length(dev_time_MjuvII))*0.1, dev_time_MjuvII, main = "Development days Juvenile II vs fit \n min_size = 1 cm and max_size = 2 cm" )
lines(temperature_range_red_MjuvII, predict(fit_MjuvII), col = '#1c9099', lwd = 2)
legend("topright", legend = c('1/devTime', 'K_func_briere'), fill = c("black", "#1c9099"))


#### M JUV III ####

development.df_MjuvIII <- data.frame(row = seq(1, 1095))

#simulate development with constant temperature for each 0.1°C & store results columns wise in the df:
for (i in temperature_range){
  development <- c()
  initial_l <- 2 #cm
  for (j in time_range){
    growth = som_growth_thesis(initial_l, j ,i , "M")
    #initial_l = growth #+  L_init
    development = append(development, growth)
  }
  development.df_MjuvIII[as.character(i)] <- development
}

development.df_MjuvIII = development.df_MjuvIII[ , !(names(development.df_MjuvIII) == "row")] #delete first column with timesteps, same as index

#Count how many time steps there are, where individual's size is between min and max size of its size class range:
dev_time_MjuvIII <- sapply(development.df_MjuvIII, function(col){
  BrownShrimp::count_devDays(col, juvIII_min, juvIII_max) } )


#With following plot, we know the shape of the curve (U shape. But inverted is again TPC)
#plot( (1:length(dev_time_MjuvIII))*0.1, dev_time_MjuvIII, main = "Development days Juvenile III \n min_size = 0.6 cm and max_size = 1 cm" )
plot( (1:length(dev_time_MjuvIII))*0.1, 1/dev_time_MjuvIII, main = "1/Development days Juvenile III \n min_size = 2 cm and max_size = 3 cm" )
lines(temperature_range, k_vals_M, col = '#1c9099', lwd = 2)
legend("topright", legend = c('1/devTime', 'K_func_briere'), fill = c("grey48", "#1c9099"))

#Deleet values where there is no development i.e. when devTime = 1095
indexes_1095 <- which(dev_time_MjuvIII == 1095)
dev_time_MjuvIII_reduced <- dev_time_MjuvIII[dev_time_MjuvIII != 1095]
temperature_range_red_MjuvIII <- temperature_range[-indexes_1095]
k_vals_reduced_MjuvIII = k_vals_M[-indexes_1095]

devTime_MjuvIII_inv_reduced = 1/dev_time_MjuvIII_reduced
k_vals_reduced_MjuvIII_inv = 1/k_vals_reduced_MjuvIII

fit_MjuvIII = lm(dev_time_MjuvIII_reduced ~ k_vals_reduced_MjuvIII_inv )
summary(fit_MjuvIII)

plot( (1:length(dev_time_MjuvIII))*0.1, dev_time_MjuvIII, main = "Development days Juvenile III vs fit \n min_size = 2 cm and max_size = 3 cm" )
lines(temperature_range_red_MjuvIII, predict(fit_MjuvIII), col = '#1c9099', lwd = 2)
legend("topright", legend = c('1/devTime', 'K_func_briere'), fill = c("black", "#1c9099"))



#### M JUV IV ####

development.df_MjuvIV <- data.frame(row = seq(1, 1095))

#simulate development with constant temperature for each 0.1°C & store results columns wise in the df:
for (i in temperature_range){
  development <- c()
  initial_l <- 3 #cm
  for (j in time_range){
    growth = som_growth_thesis(initial_l, j ,i , "M")
    #initial_l = growth #+  L_init
    development = append(development, growth)
  }
  development.df_MjuvIV[as.character(i)] <- development
}

development.df_MjuvIV = development.df_MjuvIV[ , !(names(development.df_MjuvIV) == "row")] #delete first column with timesteps, same as index

#Count how many time steps there are, where individual's size is between min and max size of its size class range:
dev_time_MjuvIV <- sapply(development.df_MjuvIV, function(col){
  BrownShrimp::count_devDays(col, juvIV_min, juvIV_max) } )


#With following plot, we know the shape of the curve (U shape. But inverted is again TPC)
#plot( (1:length(dev_time_MjuvIV))*0.1, dev_time_MjuvIV, main = "Development days Juvenile IV \n min_size = 3 cm and max_size = 4 cm" )
plot( (1:length(dev_time_MjuvIV))*0.1, 1/dev_time_MjuvIV, main = "1/Development days Juvenile IV \n min_size = 3 cm and max_size = 4 cm" )
lines(temperature_range, k_vals_M, col = '#1c9099', lwd = 2)
legend("topright", legend = c('1/devTime', 'K_func_briere'), fill = c("grey48", "#1c9099"))

#Deleet values where there is no development i.e. when devTime = 1095
indexes_1095 <- which(dev_time_MjuvIV == 1095)
dev_time_MjuvIV_reduced <- dev_time_MjuvIV[dev_time_MjuvIV != 1095]
temperature_range_red_MjuvIV <- temperature_range[-indexes_1095]
k_vals_reduced_MjuvIV = k_vals_M[-indexes_1095]

devTime_MjuvIV_inv_reduced = 1/dev_time_MjuvIV_reduced
k_vals_reduced_MjuvIV_inv = 1/k_vals_reduced_MjuvIV

fit_MjuvIV = lm(dev_time_MjuvIV_reduced ~ k_vals_reduced_MjuvIV_inv )
summary(fit_MjuvIV)

plot( (1:length(dev_time_MjuvIV))*0.1, dev_time_MjuvIV, main = "Development days Juvenile IV vs fit \n min_size = 3 cm and max_size = 4 cm" )
lines(temperature_range_red_MjuvIV, predict(fit_MjuvIV), col = '#1c9099', lwd = 2)
legend("topright", legend = c('1/devTime', 'K_func_briere'), fill = c("black", "#1c9099"))

#### M JUV V ####

development.df_MjuvV <- data.frame(row = seq(1, 1095))

#simulate development with constant temperature for each 0.1°C & store results columns wise in the df:
for (i in temperature_range){
  development <- c()
  initial_l <- 4 #cm
  for (j in time_range){
    growth = som_growth_thesis(initial_l, j ,i , "M")
    #initial_l = growth #+  L_init
    development = append(development, growth)
  }
  development.df_MjuvV[as.character(i)] <- development
}

development.df_MjuvV = development.df_MjuvV[ , !(names(development.df_MjuvV) == "row")] #delete first column with timesteps, same as index

#Count how many time steps there are, where individual's size is between min and max size of its size class range:
dev_time_MjuvV <- sapply(development.df_MjuvV, function(col){
  BrownShrimp::count_devDays(col, juvV_min, juvV_max) } )


#With following plot, we know the shape of the curve (U shape. But inverted is again TPC)
#plot( (1:length(dev_time_MjuvV))*0.1, dev_time_MjuvV, main = "Development days Juvenile V \n min_size = 4 cm and max_size = 5 cm" )
plot( (1:length(dev_time_MjuvV))*0.1, 1/dev_time_MjuvV, main = "1/Development days Juvenile V \n min_size = 4 cm and max_size = 5 cm" )
lines(temperature_range, k_vals_M, col = '#1c9099', lwd = 2)
legend("topright", legend = c('1/devTime', 'K_func_briere'), fill = c("grey48", "#1c9099"))

#Deleet values where there is no development i.e. when devTime = 1095
indexes_1095 <- which(dev_time_MjuvV == 1095)
dev_time_MjuvV_reduced <- dev_time_MjuvV[dev_time_MjuvV != 1095]
temperature_range_red_MjuvV <- temperature_range[-indexes_1095]
k_vals_reduced_MjuvV = k_vals_M[-indexes_1095]

devTime_MjuvV_inv_reduced = 1/dev_time_MjuvV_reduced
k_vals_reduced_MjuvV_inv = 1/k_vals_reduced_MjuvV

fit_MjuvV = lm(dev_time_MjuvV_reduced ~ k_vals_reduced_MjuvV_inv )
summary(fit_MjuvV)

plot( (1:length(dev_time_MjuvV))*0.1, dev_time_MjuvV, main = "Development days Juvenile V vs fit \n min_size = 4 cm and max_size = 5 cm" )
lines(temperature_range_red_MjuvV, predict(fit_MjuvV), col = '#1c9099', lwd = 2)
legend("topright", legend = c('1/devTime', 'K_func_briere'), fill = c("black", "#1c9099"))

#### M ADULT I#### not needed

development.df_MaduI <- data.frame(row = seq(1, 1095))

#simulate development with constant temperature for each 0.1°C & store results columns wise in the df:
for (i in temperature_range){
  development <- c()
  initial_l <- 5 #cm
  for (j in time_range){
    growth = som_growth_thesis(initial_l, j ,i , "M")
    #initial_l = growth #+  L_init
    development = append(development, growth)
  }
  development.df_MaduI[as.character(i)] <- development
}

development.df_MaduI = development.df_MaduI[ , !(names(development.df_MaduI) == "row")] #delete first column with timesteps, same as index

#Count how many time steps there are, where individual's size is between min and max size of its size class range:
dev_time_MaduI <- sapply(development.df_MaduI, function(col){
  BrownShrimp::count_devDays(col, adultI_min, Linf_M) } )


#With following plot, we know the shape of the curve (U shape. But inverted is again TPC)
plot( (1:length(dev_time_MaduI))*0.1, dev_time_MaduI, main = "Development days Adult I \n min_size = 5 cm and max_size = 6 cm" )
plot( (1:length(dev_time_MaduI))*0.1, 1/dev_time_MaduI, main = "1/Development days Adult I \n min_size = 5 cm and max_size = 6 cm" )
lines(temperature_range, k_vals_M, col = '#1c9099', lwd = 2)
legend("topright", legend = c('1/devTime', 'K_func_briere'), fill = c("gray48", "#1c9099"))

#Deleet values where there is no development i.e. when devTime = 1095
indexes_1095 <- which(dev_time_MaduI == 1095)
dev_time_MaduI_reduced <- dev_time_MaduI[dev_time_MaduI != 1095]
temperature_range_red_MaduI <- temperature_range[-indexes_1095]
k_vals_reduced_MaduI = k_vals_M[-indexes_1095]

devTime_MaduI_inv_reduced = 1/dev_time_MaduI_reduced
k_vals_reduced_MaduI_inv = 1/k_vals_reduced_MaduI

fit_MaduI = lm(dev_time_MaduI_reduced ~ k_vals_reduced_MaduI_inv )
summary(fit_MaduI)

plot( (1:length(dev_time_MaduI))*0.1, dev_time_MaduI, main = "Development days Adult I vs fit \n min_size = 5 cm and max_size = 6 cm" )
lines(temperature_range_red_MaduI, predict(fit_MaduI), col = '#1c9099', lwd = 2)
legend("topright", legend = c('1/devTime', 'K_func_briere'), fill = c("gray48", "#1c9099"))


#### Fits ####

intercepts_M = c(coef(fit_M.JuvI)[1], coef(fit_MjuvII)[1], coef(fit_MjuvIII)[1], coef(fit_MjuvIV)[1], coef(fit_MjuvV)[1] )
factors_M = c(coef(fit_M.JuvI)[2], coef(fit_MjuvII)[2], coef(fit_MjuvIII)[2], coef(fit_MjuvIV)[2], coef(fit_MjuvV)[2])

sizeClass_means_reduced_M = sizeClass_means[1:5]

fitsM_df = data.frame( row.names = sizeClass_names[1:5], intercept = intercepts_M, factor = factors_M)

#check weather INTERCEPTS show any relation with the mean size of the class:
fit_intercepts_M = lm(fitsM_df$intercept ~ sizeClass_means_reduced_M)
summary(fit_intercepts_M)

plot(1:5, fitsM_df$intercept)
abline(h= mean(fitsM_df$intercept), lty=2 )
lines(1:5, predict(fit_intercepts_M), col = 'orange3')



#check weather FACTORS show any relation with the mean size of the class:
plot(1:5, fits_df$factor)

logx_m = log(sizeClass_means_reduced_M)
logy_m = log(fitsM_df$factor)


fit_factortsLog_m = lm(logy_m ~ logx_m)

summary(fit_factortsLog_m)

log_a_m <- coef(fit_factortsLog_m)[1]
b_est_m <- coef(fit_factortsLog_m)[2]
a_est_m <- exp(log_a_m)

plot(sizeClass_means_reduced_M, fitsM_df$factor, log = 'xy')
lines(sizeClass_means_reduced_M, a_est_m * sizeClass_means_reduced_M^b_est_m, col = '#1c9099')
