library(BrownShrimp)

#APPROACH: to compare and if applicable, to adapt own spawning function to Temings result.
#For all simulations, I will consider a constant number of females (100) with an specific size (6 cm)

#### MOLTING FRACTION ####
#Check behavior of molting fraction
molting_fraction_old <- function(L, temperature){
  mf = 0
  if (L >= 5){
    mf = (1 / (5.7066 * L^0.7364 * exp(temperature*-0.09363) ) ) #revision 28.05. Temming's implementation (above was Temming's formula on paper)
  }
  return ( mf )
}

temp_range_wide = seq(1,40, 0.1)

plot(temp_range_wide, molting_fraction_old(6, temp_range_wide))
abline(h = 1, lty= 2)

#We identify here that curve keeps growing with temperatures that are not realistic and that are not congruent
#with TPC theory; furthermore, the mf goes above 1, which is also completely unrealistic.
#Solution: (already implemented in bib) mf for T above T_opt is always same as mf at T_opt


#### SPAWNIGN RATE ####

avg_weight_egg = 17.725*10^-6 #[gr] ~avg 17.725 microgramm over the year (from table from seasonal changes on eggs Urzua)
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
  if(L>5) s = convertL_to_W(L)*K_func(temperature)*3
  return(s)
}

for (i in temp_range){
  prod_temming <- fem_abundance*molting_fraction(L_test1, i)*number_eggs(L_test1*10) #L*10 because Temming's fomula is in mm
  B_fe = convertL_to_W(6)
  prod_AF <- B_fe*fem_abundance*molting_fraction(L_test1, i)*spawning_rate_old(L_test1, i)#/epsilon #Biommass
  eggs_production_temming <- append(eggs_production_temming, prod_temming) #number of eggs
  eggs_production_AF <- append(eggs_production_AF, prod_AF) #biomass
}


#Now plot to compare biomass

Teming_biomass = eggs_production_temming*avg_weight_egg #[gr]

plot(temp_range, Teming_biomass, main = "Egg production Temming vs Bertalanffy \n  Nr. Fems. = 100, L=6cm",
     xlab = 'T [°C]', ylab = 'Biomass [gr]', col = 'gray41', las = 1, lwd = 2,
     cex.main = 1.5, cex.lab = 1.75, cex.axis = 1.75)
lines(temp_range, eggs_production_AF, col = 'orangered3', lwd = 3)
abline(v = temp_range[which.max(Teming_biomass)], lty = 2)
legend("topleft",
       legend = c("Temming", "v. Bertalanffy k*w", "T_opt"),
       col = c("gray41", "orangered3", "black"),
       pch = c(16, NA, NA),
       lty = c(1,1,2))


#We observe that my function's output is ways smaller than the one from Temming. Thus I will try to adapt my
#results to Temming's one by fitting them BUT with an important remark: Spawning from Temming doesn't consider
#decay of curve after passing optimal temperature. As I will consder the effect of thermal performance curve,
#the fit will be carried out only for all results bellow T_op = 16.2°C.


fit <- lm(Teming_biomass[1:163] ~ eggs_production_AF[1:163])
summary(fit)
intercept = fit$coefficients[1] #3.49198
factor = fit$coefficients[2]#0.95393

plot(temp_range, Teming_biomass, main = "Temming vs AF fitted \n with constant T and for L=6cm",
     xlab = 'T', ylab = 'Biomass')
lines(temp_range, intercept + eggs_production_AF*factor, col = 'orangered3', lwd = 2)


#END OF REVISED VERSION 30.05.25

####

temperature_dataSet <- read.csv("./data/Temperature_data_1982-2018_MarBiol.csv", header = TRUE, sep = ';')
temperature_90 <- temperature_func(temperature_dataSet, "01/01/1990", "31/12/1990")
temperature_85 <- temperature_func(temperature_dataSet, "01/01/1985", "31/12/1985")

spG_90 = sapply(temperature_90$temperature, spawning_rate_b, L = 5.5, sex_params= parameters_solv$Fem_params)
plot(temperature_90$date_time, spG_90, main = 'Spawning rate with temperature data from 1990',
     type = 'l', col = '#33638DFF')


spG_85 = sapply(temperature_85$temperature, spawning_rate_b, L = 5.5, sex_params= parameters_solv$Fem_params)
plot(temperature_85$date_time, spG_85, main = 'Spawning rate with temperature data from 1985',
     type = 'l', col = '#33638DFF')



