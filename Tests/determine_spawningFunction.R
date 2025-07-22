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
mat_size_range = seq(5.5, 8.5, by = 0.1)

fem_abundance <- 1 #number of females (we will consider first also a constant size L=6 cm)
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

spawning_rate_vB <- function(L, temperature){
  s = 0
  if(L>5.5) s = convertL_to_W(L)*K_func(temperature)*3
  return(s)
}


spawning_s = data.frame(row = seq(1, length(temp_range)))
egg_prod_T = data.frame(row = seq(1, length(temp_range)))

for (s in mat_size_range){
  eggs_production_temming = c()
  eggs_production_AF = c()
  for (i in temp_range){
    prod_temming <- fem_abundance*number_eggs(s*10) #L*10 because Temming's fomula is in mm
    B_fe = convertL_to_W(s)
    #prod_AF <- B_fe*fem_abundance*spawning_rate_vB(s, i)#, parameters_solv$Fem_params) #Biommass
    prod_AF <-fem_abundance*spawning_rate_vB(s, i)#, parameters_solv$Fem_params) #Biommass
    eggs_production_temming <- append(eggs_production_temming, prod_temming) #number of eggs
    eggs_production_AF <- append(eggs_production_AF, prod_AF) #biomass
  }
  spawning_s[as.character(s)] <- eggs_production_AF
  egg_prod_T[as.character(s)] <- eggs_production_temming
}



#Now plot to compare biomass

Teming_biomass = egg_prod_T[, 2:ncol(egg_prod_T)]*17.725*10^-6#eggs_production_temming*avg_weight_egg #[gr]

factor = c()
for (i in 1:ncol(Teming_biomass)){
  B_i = convertL_to_W(mat_size_range[i])
  factor = append(factor, Teming_biomass[[i]][1]/B_i)
}

B_sc =  convertL_to_W(mat_size_range)
logx = log(B_sc)
logy = log(factor)


fit_factortsLog = lm(logy ~ logx)

summary(fit_factortsLog)

log_a <- coef(fit_factortsLog)[1]
b_est <- coef(fit_factortsLog)[2]
a_est <- exp(log_a)

plot(B_sc, factor, log = 'xy')
lines(B_sc, a_est * B_sc^b_est, col = 'red')

plot(B_sc, factor, ylim = c(0.2, 0.4))
lines(B_sc, a_est * B_sc^b_est, col = 'red')





#FUNCTION TO CALCULATE BIOMASS OF EGGS PRODUCED by 12 individual (only in dependency of L)

egg_biomass = function(L){
  return(L*(coef(fit_fact + )))
}





plot(temp_range, Teming_biomass[['6']], main = "Egg production Temming vs Bertalanffy \n  Nr. Fems. = 100, L=6cm",
     xlab = 'T [°C]', ylab = 'Biomass [gr]', col = 'gray41', las = 1, lwd = 2,
     cex.main = 1.35, cex.lab = 1.5, cex.axis = 1.25)#, ylim = c(0,8))
lines(temp_range, spawning_s[['6']], col = 'maroon', lwd = 3)
#abline(v = temp_range[which.max(Teming_biomass)], lty = 2)
legend("topleft",
       legend = c("Benchark", "v. Bertalanffy k*w", "T_opt"),
       col = c("gray41", 'maroon', "black"),
       pch = c(16, NA, NA),
       lty = c(1,1,2), lwd = c(2,2,2))


#We observe that my function's output is ways smaller than the one from Temming. Thus I will try to adapt my
#results to Temming's one by fitting them BUT with an important remark: Spawning from Temming doesn't consider
#decay of curve after passing optimal temperature. As I will consder the effect of thermal performance curve,
#the fit will be carried out only for all results bellow T_op = 16.2°C.

intercepts = c()
factors = c()

for (i in colnames(Teming_biomass)){
  fit <- lm(Teming_biomass[[i]][1:163] ~ spawning_s[[i]][1:163])
  #summary(fit)
  intercepts = append(intercepts, fit$coefficients[1] ) #3.49198
  factors = append(factors, fit$coefficients[2])#0.95393

}

#Check relation between intercepts and factors with length

plot(mat_size_range, intercepts)

plot(mat_size_range, factors)

#fit linearly:

fit_intercepts = lm(intercepts ~ mat_size_range)
summary(fit_intercepts)
plot(mat_size_range, intercepts, ylim = c(0, 0.02))
lines(mat_size_range, predict(fit_intercepts))

fit_factors = lm(factors ~ mat_size_range)
summary(fit_factors)
plot(mat_size_range, factors)
lines(mat_size_range, predict(fit_factors))

##### conclusion ####
#Spanwning rate is dependet on temperature and lenght:

spawning_rate = function(L, Te){
  interc = coef(fit_intercepts)[1] + L*coef(fit_intercepts)[2]
  fact = coef(fit_factors)[1] + L*coef(fit_factors)[2]
  return(B_fe*fem_abundance*molting_fraction(s, i)*spawning_rate_vB(s, i))
}





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



