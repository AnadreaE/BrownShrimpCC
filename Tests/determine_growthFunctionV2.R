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


#### Def. growth function (as function of time) ####

lengthGrowth_funcTime = function (L0, time, temp, sex){ #L[cm]
  #This funciton works only in this context i.e. to simulate growth with constant temperatures in order to
  #compare and check plausibilty of the results.
  L0 = L0*10
  if (sex == "F") Linf = 85
  if (sex == "M") Linf = 55
  if (L0 > Linf ) growth = 0
  else growth = Linf- (Linf-L0)*exp(-K_func_briere(temp)*time)
  return(growth/10 )
}

#### Simulations ####
time_range <- seq(1,1095)
temperature_range <- seq(0,30, 0.1) #(10,15)

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

dev_time_juvI <- sapply(development.df_juvI, function(col){
  BrownShrimp::count_devDays(col, juvI_min, juvI_max) } )

#dev_time_juvI <- dev_time_juvI[-c(1)] #delete first element cuz doesn't make sense

k_vals = sapply(temperature_range, K_func_briere)

devTime_juvI_inv = 1/dev_time_juvI

#With following plot, we know the shape of the curve (U shape. But inverted is again TPC)
plot( (1:length(dev_time_juvI))*0.1, devTime_juvI_inv, main = "Development days Juvenile I \n min_size = 6 mm and max_size = 20 mm" )
lines(temperature_range, k_vals, col = 'red')

#Deleet values where there is no development i.e. when devTime = 1095
indexes_1095 <- which(dev_time_juvI == 1095)
dev_time_JUVI_reduced <- dev_time_juvI[dev_time_juvI != 1095]
temperature_range_red_juvI <- temperature_range[-indexes_1095]
k_vals_red = k_vals[-indexes_1095]

devTime_juvI_inv_reduced = 1/dev_time_JUVI_reduced
test = 1/k_vals_red

fit_JuvI = lm(dev_time_JUVI_reduced ~ test )
summary(fit_JuvI)

plot( (1:length(dev_time_juvI))*0.1, 1/devTime_juvI_inv, main = "Development days Juvenile I \n min_size = 6 mm and max_size = 20 mm" )
lines(temperature_range_red_juvI, predict(fit_JuvI), col = 'red4')

