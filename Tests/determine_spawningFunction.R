library(BrownShrimp)

#APPROACH: to compare and if applicable, to adapt own spawning function to Temings result.
#For all simulations, I will consider a constant number of females (100) with an specific size (6 cm)

#### some constants needed ####
the.method = 'rk4'

m = 3 # ww-lenght scaling exponent
n = 3/4
const_c = 0.01 # reference density g/cm3
#L_inf = 8.5 #cm
#a = 0.008 #
epsilon = 0.22

avg_weight_egg = 18.e-6 #~avg 18 microgramm over the year (from table from seasonal changes on eggs Urzua)
temp_range = seq(0,30, 0.1)

fem_abundance <- 100 #number of females (we will consider first also a constant size L=6 cm)
L_test1 <- 6 #cm

#TEMING

# SI = C_l * E_l * EF_l * MF_l
# c_l: number of individuals per a size class 'l'
# E_l: Share of egg-bearing Fem in size class 'l'
# EF_l: mean number of eggs per size class
# MF_l: molting frequency

number_eggs <- function(L){
  0.01805*L^3.539 #incongruency on Temming's paper weather 0.01805 or 0.001805. The first seem more plausible for me
}


#### SIMULATIONS ####

eggs_production_temming <- c()
eggs_production_AF <- c()

spawning_rate_old <- function(L, temperature){
  s = 0
  if(L>5) s = convertL_to_W(L)*K_func(temperature)*3*epsilon #molting_fraction(L*10, T) * convertL_to_W(L)*K_func(T)*3*epsilon #here L for Temming in mm
  return(s)
}

for (i in temp_range){
  prod_temming <- fem_abundance*molting_fraction(L_test1, i)*number_eggs(L_test1*10)
  B_fe = convertL_to_W(6)
  prod_AF <- B_fe*fem_abundance*molting_fraction(L_test1, i)*spawning_rate_old(L_test1, i)/epsilon
  eggs_production_temming <- append(eggs_production_temming, prod_temming) #number of eggs
  eggs_production_AF <- append(eggs_production_AF, prod_AF) #biomass
}

plot(temp_range, eggs_production_temming, main = "Temming egg production of 100Fems \n with constant T and for L=6cm",
     xlab = 'T', ylab = 'Number Eggs')

#Now plot to compare biomass

Teming_biomass = eggs_production_temming*avg_weight_egg

plot(temp_range, eggs_production_temming*avg_weight_egg, main = "Temming egg production of 100Fems \n with constant T and for L=6cm",
     xlab = 'T', ylab = 'Biomass')
lines(temp_range, eggs_production_AF, col = 'orangered3', lwd = 2)


#We observe that my function's output is ways smaller than the one from Temming. Thus I will try to adapt my
#results to Temming's one by fitting them BUT with an important remark: Spawning from Temming doesn't consider
#decay of curve after passing optimal temperature. As I will consder the effect of thermal performance curve,
#the fit will be carried out only for all results bellow T_op = 16.2°C.


fit <- lm(Teming_biomass[1:163] ~ eggs_production_AF[1:163])
summary(fit)
intercept = 3.54616
factor = 4.40332

plot(temp_range, eggs_production_temming*avg_weight_egg, main = "Temming vs AF fitted \n with constant T and for L=6cm",
     xlab = 'T', ylab = 'Biomass')
lines(temp_range, intercept + eggs_production_AF*factor, col = 'orangered3', lwd = 2)

#END OF REVISED VERSION 30.05.25

####

temperature_dataSet <- read.csv("./data/Temperature_data_1982-2018_MarBiol.csv", header = TRUE, sep = ';')
temperature_90 <- temperature_func(temperature_dataSet, "01/01/1990", "31/12/1990")
temperature_85 <- temperature_func(temperature_dataSet, "01/01/1985", "31/12/1985")

plot(temp_range, sp_g, main = 'Spawning rate with temperature changes')

spG_90 = spawning_rate(5.5, temperature_90$temperature)
plot(temperature_90$date_time, spG_90, main = 'Spawning rate with temperature data from 1990',
     type = 'l', col = '#33638DFF')


spG_85 = spawning_rate(5.5, temperature_85$temperature)
plot(temperature_85$date_time, spG_85, main = 'Spawning rate with temperature data from 1985',
     type = 'l', col = '#33638DFF')



