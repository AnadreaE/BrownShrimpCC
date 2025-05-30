

#APPROACH: to compare and if applicable, to adapt own spawning function to Temings result.
#therefore first run a simulation with Temming's formula.

temperature_dataSet <- read.csv("C:/Users/andre/Desktop/Hereon/Data/Wadden Sea, Penning/Temperature_data_1982-2018_MarBiol.csv", header = TRUE, sep = ';')
temperature_90 <- temperature_func(temperature_dataSet, "01/01/1990", "31/12/1990")
temperature_85 <- temperature_func(temperature_dataSet, "01/01/1985", "31/12/1985")


temp_range = seq(0,30, 0.1)

sizes = c(5.5, 6.5, 7.5, 8.5)

sp_g = spawning_rate(5.5, temp_range)

plot(temp_range, sp_g, main = 'Spawning rate with temperature changes')

spG_90 = spawning_rate(5.5, temperature_90$temperature)
plot(temperature_90$date_time, spG_90, main = 'Spawning rate with temperature data from 1990',
     type = 'l', col = '#33638DFF')


spG_85 = spawning_rate(5.5, temperature_85$temperature)
plot(temperature_85$date_time, spG_85, main = 'Spawning rate with temperature data from 1985',
     type = 'l', col = '#33638DFF')


#TEMING
# For a simulation with Temming's formula, I will consider a constant number of females over the time

# SI = C_l * E_l * EF_l * MF_l
# c_l: number of individuals per a size class 'l'
# E_l: Share of egg-bearing Fem in size class 'l'
# EF_l: mean number of eggs per size class
# MF_l: molting frequency

fem_abundance <- 100 #number of females (we will consider first also a constant size L=6 cm)
L_test1 <- 6

number_eggs <- function(L){
  0.01805*L^3.539 #incongruency on Temming's paper weather 0.01805 or 0.001805. The first seem more plausible for me
}

eggs_production_temming <- c()
eggs_production_AF <- c()

for (i in temp_range){
  prod_temming <- fem_abundance*molting_fraction(L_test1*10, i)*number_eggs(L_test1)
  prod_AF <- fem_abundance*spawning_rate(L_test1, i)
  eggs_production_temming <- append(eggs_production_temming, prod_temming)
  eggs_production_AF <- append(eggs_production_AF, prod_AF)
}

plot(temp_range, eggs_production_temming, main = "Temming egg production of 100Fems \n with constant T and for L=6cm",
     xlab = 'T', ylab = 'Number Eggs')
#lines(temp_range, eggs_production_AF, col = 'orangered3',type = 'l' )
plot(temp_range, eggs_production_AF, col = 'orangered3', lwd = 2, main = 'Spawning Thesis Model',
     ylab = 'Biomass eggs produced / day')






