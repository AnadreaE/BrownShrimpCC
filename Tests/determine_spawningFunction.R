library(BrownShrimp)

#APPROACH: to compare and if applicable, to adapt own spawning function to Temings result.
#For all simulations, I will consider a constant number of females (100) with an specific size (6 cm)

#### MOLTING FRACTION ####
#Check behavior of molting fraction
molting_fraction_Temming <- function(L, temperature){
  mf = 0
  if (L >= 5.5){
    L=10*L
    mf = (1 / (5.7066 * L^0.7364 * exp(temperature*-0.09363) ) ) #revision 28.05. Temming's implementation (above was Temming's formula on paper)
  }
  return ( mf )
}

vB_term = function(L, fem_params, Te){
  w = convertL_to_W(L)
  return(3*K_func_briere(Te, parameters_solv$Fem_params)*w)
}

number_eggs <- function(L){
  L=L*10
  0.01805*L^3.539 #incongruency on Temming's paper weather 0.01805 or 0.001805. The first seem more plausible for me
}


temp_range_wide = seq(1,40, 0.1)
par(mfrow=c(1,1))
plot(temp_range_wide, molting_fraction_Temming(6, temp_range_wide), ylim = c(0,0.4), main = "Molting fraction Temming vs Bertalanffy term \n L=6cm",
     xlab = 'T [°C]', ylab = 'molting fraction')
abline(h = 1, lty= 2) #Surpassing one is very unrealistic !
lines(temp_range_wide, sapply(temp_range_wide, vB_term, L=6, fem_params=parameters_solv$Fem_params), col = 'red')


#We identify here that curve keeps growing with temperatures that are not realistic and that are not congruent
#with TPC theory; furthermore, the mf goes above 1, which is also completely unrealistic.
#Solution: piece-wise (already implemented in bib) mf for T above T_opt is always same as mf at T_opt

#### BELLOW LINES WAS NOS USED IN THE MODEL ######

K_vals = sapply(temp_range, K_func_briere, parameters_solv$Fem_params)
T_opt = which.max(K_vals)/10
#temp_less_opt = temp_range_wide[temp_range_wide <= T_opt]

#Simulate benchmark
avg_weight_egg = 17.725*10^-6 #[gr] ~avg 17.725 microgramm over the year (from table from seasonal changes on eggs Urzua)
temp_range = seq(0,30, 0.1)
mat_size_range = seq(5.5, 8.5, by = 0.1)

benchmark_egg_biomass = data.frame(row = seq(1, length(temp_range)))

fem_abundance <- 1 #number of females (we will consider first also a constant size L=6 cm)
for (s in mat_size_range){
  eggs_production_temming = c()
  #eggs_production_AF = c()
  for (i in temp_range){
    prod_temming <- fem_abundance*number_eggs(s) #L*10 because Temming's fomula is in mm
    #B_fe = convertL_to_W(s)
    #prod_AF <- B_fe*fem_abundance*spawning_rate_vB(s, i)#, parameters_solv$Fem_params) #Biommass
    #prod_AF <-fem_abundance*spawning_rate_vB(s, i)#, parameters_solv$Fem_params) #Biommass
    eggs_production_temming <- append(eggs_production_temming, prod_temming) #number of eggs
    #eggs_production_AF <- append(eggs_production_AF, prod_AF) #biomass
  }
  #spawning_s[as.character(s)] <- eggs_production_AF
  benchmark_egg_biomass[as.character(s)] <- eggs_production_temming
}

benchmark_egg_biomass = benchmark_egg_biomass[, 2:ncol(benchmark_egg_biomass)]*avg_weight_egg#eggs_production by one individual

#Now, lets try to fit molting fraction to vB term:

benchmark_moltingFraction = data.frame(row = seq(1, length(temp_range)))

for (s in mat_size_range){
  mol_fraction_temming = c()
  #eggs_production_AF = c()
  for (i in temp_range){
    mf <- molting_fraction_Temming(s, i) #L*10 because Temming's fomula is in mm
    mol_fraction_temming <- append(mol_fraction_temming, mf) #number of eggs
  }
  benchmark_moltingFraction[as.character(s)] <- mol_fraction_temming
}
#delete first column 'row'
benchmark_moltingFraction =  benchmark_moltingFraction[,2:ncol(benchmark_moltingFraction)]

plot(temp_range, benchmark_moltingFraction$`6`)
t_optstep = T_opt*10
vB_vals = sapply(temp_range, vB_term, L=6, fem_params=parameters_solv$Fem_params)
fit_test = lm(benchmark_moltingFraction$`6`[1:t_optstep] ~ vB_vals[1:t_optstep])
summary(fit_test)

#Check with different sizes
plot(temp_range, benchmark_moltingFraction$`6`, main = 'test for L = 6cm')
lines(temp_range, coef(fit_test)[1] + vB_vals * coef(fit_test)[2], col = 'red' )

vB_vals = sapply(temp_range, vB_term, L=7, fem_params=parameters_solv$Fem_params)
plot(temp_range, benchmark_moltingFraction$`7`, main = 'test for L = 7cm')
lines(temp_range, coef(fit_test)[1] + vB_vals * coef(fit_test)[2], col = 'red' )
#A fit for each size has to be carried out. Then the intercepts and factors can be fitted to the adult size clases
#to get a function with L as parameter:

intercept_fits = c()
factor_fits = c()

for (i in colnames(benchmark_moltingFraction)){
  L = as.numeric(i)
  vB_vals = sapply(temp_range, vB_term, L=L, fem_params=parameters_solv$Fem_params)
  fit = lm(benchmark_moltingFraction[[i]][1:t_optstep] ~ vB_vals[1:t_optstep])
  intercept_fits = append(intercept_fits, coef(fit)[1])
  factor_fits = append(factor_fits, coef(fit)[2])
}

plot(mat_size_range, intercept_fits)
abline(h=mean(intercept_fits), lty= 2)

fit_intercept_exp = nls(intercept_fits ~ a*exp(mat_size_range*-b), start = list(a=1 ,b=0.05) )
summary(fit_intercept_exp)
plot(mat_size_range, intercept_fits)
lines(mat_size_range,predict(fit_intercept_exp))
#since values are very small, we may just consider the average value



plot(mat_size_range, factor_fits)
#here we see the decaying exponential value and thus proceed with fit to convert it as funcion of length
fit_facto_exp = nls(factor_fits ~ a*exp(mat_size_range*-b), start = list(a=1 ,b=0.5) )
summary(fit_facto_exp)
plot(mat_size_range, factor_fits)
lines(mat_size_range, predict(fit_facto_exp), col = 'red')


spawning_energy = function(L, Te, Fem_params){
  w = convertL_to_W(L)
  vB = 3*K_func_briere(Te, parameters_solv$Fem_params)*w
  intercept = coef(fit_intercept_exp)[1]*exp(L*-coef(fit_intercept_exp)[2])
  factor = coef(fit_facto_exp)[1]*exp(L*-coef(fit_facto_exp)[2])
  return(mean(intercept_fits) + vB*factor)
}

#vB_vals = sapply(temp_range, vB_term, L=7, fem_params=parameters_solv$Fem_params)
plot(temp_range, benchmark_moltingFraction$`6`, main = 'test for L = 7cm')
lines(temp_range, sapply(temp_range, spawning_energy, L=6, Fem_params = parameters_solv$Fem_params), col = 'red' )


#### NOW TEST RESULTS ####
#Manual test for different L
#L=6
Nr_Afem = 100

molFraction_Temming = molting_fraction_Temming(6, temp_range)
ovigerous_fem_Temming = Nr_Afem*molFraction_Temming
biomass_temming = ovigerous_fem_Temming*number_eggs(6)*avg_weight_egg

molFraction_Thesis = sapply(temp_range, spawning_energy, L=6,  Fem_params = parameters_solv$Fem_params)
ovigerous_fem_Thesis = Nr_Afem*molFraction_Thesis
biomass_thesis = ovigerous_fem_Thesis*spawning_rate_b(6)

#dev.off()

plot(temp_range, molFraction_Temming)
lines(temp_range, molFraction_Thesis, col = 'red')


#L=7
Nr_Afem = 100

molFraction_Temming = molting_fraction_Temming(7, temp_range)
ovigerous_fem_Temming = Nr_Afem*molFraction_Temming
biomass_temming = ovigerous_fem_Temming*number_eggs(7)*avg_weight_egg

molFraction_Thesis = sapply(temp_range, molting_fraction, L=7,  Fem_params = parameters_solv$Fem_params)
ovigerous_fem_Thesis = Nr_Afem*molFraction_Thesis
biomass_thesis = ovigerous_fem_Thesis*spawning_rate_b(7)

#dev.off()

plot(temp_range, molFraction_Temming)
lines(temp_range, molFraction_Thesis, col = 'red')




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

temp_range = seq(1,30, 0.1)

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
    #B_fe = convertL_to_W(s)
    #prod_AF <- B_fe*fem_abundance*spawning_rate_vB(s, i)#, parameters_solv$Fem_params) #Biommass
    prod_AF <-fem_abundance*spawning_rate_vB(s, i)#, parameters_solv$Fem_params) #Biommass
    eggs_production_temming <- append(eggs_production_temming, prod_temming) #number of eggs
    eggs_production_AF <- append(eggs_production_AF, prod_AF) #biomass
  }
  spawning_s[as.character(s)] <- eggs_production_AF
  egg_prod_T[as.character(s)] <- eggs_production_temming
}



#Now plot to compare biomass

Teming_biomass = egg_prod_T[, 2:ncol(egg_prod_T)]*avg_weight_egg#eggs_production_temming*avg_weight_egg #[gr]

plot(temp_range, Teming_biomass$'6', ylim = c(0,0.7))
lines(temp_range, spawning_s$'6', col = 'red')

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
