#<!--
# SPDX-FileCopyrightText: 2025 Helmholtz-Zentrum hereon GmbH
# SPDX-License-Identifier: Apache-2.0
# SPDX-FileContributor Andrea Farfan <farfanqbb@gmail.de>
#  -->


################################################
## Following functions are meant to be        ##
## used for solving differential eq. systems  ##
## Source: Andrea F.                          ##
################################################

library(dplyr)
library(lubridate)

##### constant parameter values: #####

the.method = 'rk4'

L_as_F = 8.6
L_as_M = 5.6

juvI_min = .6
juvI_max = 1.6#1.0

juvII_min = 1.6#1.0
juvII_max = 2.6#2.0

juvIII_min = 2.6#2.0
juvIII_max = 3.6#3.0

juvIV_min = 3.6#3.0
juvIV_max = 4.6#4.0

juvV_min = 4.6#4.0
juvV_max = 5.6#5.0

adultI_min = 5.6#5.0
adultI_max = 6.6#6.0

adultII_min = 6.6#6.0
adultII_max = 7.6#7.0

adultIII_min = 7.6#7.0
adultIII_max = L_as_F

#Mean of the Size class ranges
LL = 0.3
LJ = (juvI_min+juvI_max)/2 #cm
LJ2 = (juvII_min+juvII_max)/2   # 3cm
LJ3 = (juvIII_min+juvIII_max)/2
LJ4 = (juvIV_min+juvIV_max)/2
LJ5 = (juvV_min+juvV_max)/2
LA1 = (adultI_min + adultI_max)/2
LA2 = (adultII_min + adultII_max)/2 #cm
LA3 = (adultIII_min + adultIII_max)/2  #cm

sizes = c(LL, LJ, LJ2, LJ3, LJ4, LJ5, LA1, LA2, LA3)

parameters_solv = list(
  Fem_params = list(
    L_inf = 8.5,
    #K-func:
    r_max = 0.007638409,
    alpha = 0.539550,
    beta = 0.27774,
    #Biomass shifts
    a_est_factor = 0.0947, #0.07254394,
    b_est_factor = 0.811, #0.9319859,
    intercept =  -0.5243291 #-0.4968367 #-0.522060 + L_mean*0.006774
  ),
  M_params = list(
    L_inf = 5.5,
    #K-func:
    r_max = 0.01585687,
    alpha = 0.634518,
    beta = 0.29947,
    #Biomass shifts
    intercept = -0.5046928,#-0.5066616,
    a_est_factor = 0.188, #0.1192449,
    b_est_factor = 0.842 #1.328571
    ),
  general_params = list(
    #func response
    #alpha_ir = 8.238, #1/(attac_rate*handling_time)for an are of 209 cm2 ~8.238 for 1m2 ???
    m = 3,
    const_c = 0.01,
    #Fishery
    Fi = 0.05, #only as example for the time being
    L50 = 4.02,#applies to a mesh size of cod-end of 22 mm
    SR = 0.82, #applies to a mesh size of cod-end of 22 mm
    #Ingestion rate & others
    Imax_ik = 0.1,
    a_mu = 0.000913, #aging daily mortality 1/3years for fems, 2 years for males
    L_mat = 5.5 #size at maturation fem [cm]
  )
)

# Relative fishing effort factors for each month (Temming, 2011)
#monthly_Feffort = c(0.19, 0.2, 0.86, 2.2, 1.25, 1, 1.19, 1.25, 3, 2, 1.6, 0.6)
#monthly_Feffort = c(0.19, 0.2, 0.86, 1.6, 1.39, 1.26, 1.19, 1.25, 1.27, 1.26, 1.09, 0.45)
monthly_Feffort = c(0.19, 0, 0, 1.6, 1.39, 1.26, 1.19, 1.25, 1.27, 1.26, 1.09, 0.45)


#Calculate factors for yearly variable fishery intensities according to landings data.
ble_data = read.csv("./data/BLE_Inlandslandungen_SpeiseKrabbe.csv", sep = ';')
years_landings = ble_data$year
ttl_yearly_landing = ble_data %>%
            group_by(year) %>%
            summarise(ttl_year = sum(t)) %>%
            filter(year < 2025) %>% #Landing data avaiable from 2010
            mutate(rel_landing = ttl_year / mean(ttl_year))

#agregar anho a esta tabla ! e integrar eso bien en el Eq. system (?)
monthly_factors = ble_data %>%
  filter(year < 2024) %>%  #Landing data avaiable from 2010
  group_by(month) %>%
  summarise(ttl_month = sum(t)) %>%
  mutate(rel_landing = ttl_month / mean(ttl_month))

monthly_factors.v = monthly_factors$rel_landing
#monthly_factors.v[4] =1.1
#monthly_factors.v = monthly_Feffort

##### FUNCTIONS #####


#' K is the vB growth factor transfored to flexTPC to add
#' temperature dependency
#'
#' @param temperature [°C]
#'
#' @returns  K value
#' @export
#'
#' @examples K_func_briere(10, parameters_solv$fem_params)

K_func_briere = function(temperature, sex_params){
  T_min = 0.5
  T_max = 30

  r_max = sex_params$r_max
  alpha = sex_params$alpha
  beta = sex_params$beta

  if(temperature < T_min) { temperature = T_min } #if this true, then T_min is negative and invalid to set to the power of alpha

  diff_min = temperature - T_min
  diff_max = T_max - temperature
  diff = T_max - T_min
  alpha_invert = 1 - alpha
  toReturn = r_max*( ((diff_min/alpha)^alpha) * ((diff_max/ alpha_invert )^alpha_invert) * (1 / diff)  )^(alpha*alpha_invert/(beta^2) )

  return(max(0, toReturn) )
}


#' Alpha param. used in ingestion rate function.
#' NB: this function works under the assumption that we consider 10 cm water depth as
#' to work with area instead of volume.
#' constants from kiorboe2014shifts
#'
#' @param w [g]
#'
#' @returns alpha value [gr (food source) m-2] , m2 considering 10 cm water depth as
#' @export
#'
#' @examples alpha_igr(convertL_to_W(5))
#' @examples alpha_igr(2)
#'

alpha_igr = function(w){
  return(64.32*w^(0.23))#.064327 #.1429
}


#' Growth in weight for an specific size class 'L' in dependence of
#' food resource availability.
#'
#'
#' @param temperature Temperature [°C]
#' @param L size class [cm]
#' @param Food Plancton (variable state) mg C m-2 (same as IC in I(F))
#'
#' @returns growth rate in weight [g day-1]
#' @export
#'
#' @examples ingestion_rate(temperature = 15, L = 5.5, P=2, parameters_solv$general, parameters_solv$Fem_params)

ingestion_rate_b = function(temperature, L, Food, general_params, sex_params){
  w = convertL_to_W(L)
  alpha = alpha_igr(w)
  m = general_params$m
  L_infty = sex_params$L_inf
  F_safe = max(Food, 0) #avoids negative values of P (this condition is better to put here than on the dP/dt to avoid inconsistency in the solver)
  result = m*K_func_briere(temperature, sex_params)*(L_infty/L)*( F_safe/(F_safe+alpha) )
  return (result)
}



#' Egg biomass that one adult Females release in
#' dependence of its size.
#' NB: This works for all size classes, however, as it is always
#' multiplied by molting fraction, where release only for L > maturity
#' is implemented, this funciton doesn't include this condition.
#'
#' @param L size [cm]
#'
#' @returns rate [g day^-1]
#' @export
#'
#' @examples spawning_rate_b(6)
spawning_rate_b = function(L){
  #w = convertL_to_W(L)
  #rate =  0.253161 * w^0.1796667
  L = L*10
  temming_nr_egg = 0.01805*L^3.539
  avg_weight_egg = 17.725*10^-9 #[kg] ~avg 17.725 microgramm over the year (from table from seasonal changes on eggs Urzua)
  return(temming_nr_egg*avg_weight_egg)
}



#' Respiration rate: is the energy loss.
#' NB: as the model includes other types of mortality, respiration rate
#' is kept very small, therefore multiplied by 0.1
#'
#' @param L size [cm]
#' @param temperature Temperature [°C]
#'
#' @returns [g day^-1 ]
#' @export
#'
#' @examples
respiration_rate_b = function(temperature, L, sex_params){
  m = parameters_solv$general_params$m
  mu = m*convertL_to_W(L)*K_func_briere(temperature, sex_params)*0.1# this is the rigth term of vB eq.
  return(mu)
}


#' Shift of biomass between juvenile and adult Size classes or 'growth function'
#' Valid only for Juvenile and adult
#' @param L_mean size mean of the max and min of the size class range [cm]
#' @param temperature [°C]
#' @param sex_params
#' @param size_width defines how big are the size ranges [cm] e.g. 1: then each cm one size class.
#'
#' @returns [day^-1]
#' @export
#'
#' @examples shift_next_sizeClass(6.1, 10, parameters_solv$Fe_params)

shift_next_sizeClass = function(L_mean, temperature, sex_params, size_width = 1){
  k = K_func_briere(temperature, sex_params)
  a = sex_params$a_est_factor
  b = sex_params$b_est_factor
  inter = sex_params$intercept
  factor = a * L_mean^b
  toReturn = k/ (k*inter + factor)/size_width #original: 1/ (inter + factor*(1/k))
  return(toReturn)
}


#' Shift of biomass from eggs to larvae or 'growth function for eggs'.
#'
#'
#' @param Te [°C]
#'
#' @returns rate at which eggs are hatching [day^-1]
#' @export
#'
#' @examples hatch_eggs(10)
hatch_eggs = function(Te){
  Te = max(0.00001, Te)
  return(1/ (1031.34*Te^-1.345)) #return the ratio
}


#' Shift of biomass from larvae to juvenile or 'growth function for larvae'
#'
#' @param Te [cm]
#'
#' @returns rate at which, larves change life stage to juvenile [day^-1]
#' @export
#'
#' @examples shiftTo_juvenile(10)
shiftTo_juvenile = function(Te){
  Te = max(0.00001, Te)
  return(1 / (941.78*Te^-1.347)) #return the ratio. # (5.5/0.00584)=941,78
}



#' Selectivity probability (Fishery)
#'#Rename it to Sel Logit !!! it has been changed
#'
#'
#' @param L
#' @param L50
#' @param SR
#'
#' @returns the probability that a size class L will be catched by net
#' of characteristics L50 and SR
#' @export
#'
#' @examples sel_probit(6) or sel_probit(L50 = 4.49 , SR = 1.56)

sel_probit = function(L, L50 = 4.49 , SR = 1.56){
  nominator = exp( (log(9)/SR) * (L - L50))
  denominator = 1 + nominator
  return(nominator / denominator)
}


#' Predatior population. Oversimplification of predator population with
#' cosinus preiodic function presenting peaks around sprint.
#'
#' @param t [days]
#'
#' @returns biomass of predator [@tbd]
#' @export
#'
#' @examples
new_Bpredator = function(t) {
  toReturn = 0.1 * (1.05 + cos( (2*pi/ 365)*(t-105) ) )
  return(toReturn)
}



#' optimal prey length
#' source of param values: Ovidio (@tbd correct source)
#'
#' @param l_pred [cm]
#'
#' @returns the optimal length of prey for a predator of size l_pred
#' @export
#'
#' @examples lopt(2.75)
lopt = function(l_pred){
  l = log(l_pred)
  r=-1.65
  gamma = 0.011
  return(r + l - gamma*l^2)
}


#' Predation activity
#' Ingestion kernel function for predation
#'
#' @param I_max
#' @param l_pred [cm]
#' @param l_prey [cm]
#'
#' @returns probability that size class l_prey will be eaten by predator of size l_pred [day^-1]
#' @export
#'
#' @examples ingestion_kernel
ingestion_kernel = function(I_max, l_pred, l_prey){
  l_opt = lopt(l_pred)
  return(I_max*exp(-3/2*(log(l_prey)-l_opt)^2))
}


#' Eta function
#' This function is used only internally as is nested in 'ST_predation'
#'
#' @param t_day
#' @param power_of
#' @param J
#'
#' @returns
#' @export
#'
#' @examples
eta_J <- function(t_day, power_of = 2, J = 250) {
  cos_term <- cos(2 * pi * (t_day - J) / 365)
  eta <- 0.5 + 0.5*cos_term  # gives values from 0 to 2
  return(eta^power_of)         # ^2 seasonal amplification
}



#' Spatial-Temporal mortality model
#'
#' @param t_day simulation day
#' @param Te Temperature [°C]
#' @param B_s Biomass of the prey (shrimp)
#' @param eta eta function
#' @param beta
#' @param sigma
#'
#' @returns
#' @export
#'
#' @examples
ST_predation = function(t_day, Te, B_s, eta, beta = 9, sigma = 0.8){
  mu_ref = 0.025 #day^-1 #same as mesozooplankton paper
  f_t = 3**((Te-10)/10) #Q10 for Cod metabolism, source:10.1111/j.1439-0426.2007.01004.x
  #sigma = 0.7 #0.5
  #eta = eta_J(t_day, J = J,  power_of = power_of)
  #beta = 9 # mesozooplankton paper = 18
  gamma = 3.7469 *1e-6 #Km2 (Kg d)^-1 #original 0.1 m3 (molC d)^-1
  return(mu_ref*f_t*sigma*(gamma*B_s + beta*eta ))
}


#' Somatic growth developed in this Thesis as function of time
#'
#' @param time time step
#' @param Tem temperature [°C]
#' @param sex 'F' or 'M'
#' @param L0 initial or current size [cm]
#'
#' @returns new size after gorwth dependent on temperature
#' @export
#'
#' @examples som_growth_thesis(3, 1, 5,"F") or:
#' development.df_juvI <- data.frame(row = seq(1, 1095))
#' simulate development with constant temperature for each 0.1°C & store results columns wise in the df:
#' for (i in temperature_range){
#'   development <- c()
#'   initial_l <- 0.6 #cm
#'   for (j in time_range){
#'     growth = som_growth_thesis(initial_l, j ,i , "F")
#'     development = append(development, growth) }
#'development.df_juvI[as.character(i)] <- development }
#'
som_growth_thesis = function (L0, time, temp, sex_params){ #L[cm]
  #This funciton works only in this context i.e. to simulate growth with constant temperatures in order to
  #compare and check plausibilty of the results AND always with L0 = 0 .
  Linf = sex_params$L_inf
  if (L0 > Linf ) growth = 0
  else growth = Linf- (Linf-L0)*exp(-K_func_briere(temp, sex_params)*time)
  #else growth = Linf*( 1 - exp(-K_func_briere(temp, sex)*time) )
  return(max(0,growth) )
}



########################################################################
########################################################################

size_width = 1
size_limits_F = seq(0.6, L_as_F, by=size_width)
size_mean_F = size_limits_F[1:length(size_limits_F)-1] + 0.5*size_width
BF = 0*size_mean_F
dBF.dt = 0*size_mean_F
N_max_F = length(size_mean_F)

size_limits_M = seq(0.6,L_as_M,by=size_width)
size_mean_M = size_limits_M[1:length(size_limits_M)-1] + 0.5*size_width
BM = 0*size_mean_M
dBM.dt = 0*size_mean_M
N_max_M = length(size_mean_M)

#VERSION 4
#Real size class model
solver_sizeClass.v4 = function(t, state, parameters, temperature_dataSet){
  system.equations = function(t, state, parameters) {
    #following line avoid negative values in state variables

    state[state < 0] = 0 # shift negative biomasses to 0
    list2env(as.list(state), envir = environment())  # re-assign the corrected state variables

    BF = state[grep("^BF", names(state))]#extract all elements from state whose names start with "BF"
    BM = state[grep("^BM", names(state))]

    dBF.dt = BF * 0
    dBM.dt = BM * 0

    Te = temperature_funcSolver(temperature_dataSet, t) # getting temperature for day t from temperature_dataSet
    month = month(temperature_dataSet$date_time[t+1]) # t starts in 0, but indices in R start at 1

    consumed_plankton = 0

    # predator
    dPred.dt = new_Bpredator(t) - Pred*0.15 #0.08 mortality

    #LARVAE
    gE = hatch_eggs(Te)
    IL = ingestion_rate_b(Te, LJ, P, parameters$general_params, parameters$Fem_params)#note that LJ is arbitrarily chossen as LL gives unrealistc values
    mL = respiration_rate_b(Te,LL, parameters$Fem_params)
    gL = shiftTo_juvenile(Te)
    pL = Pred*ingestion_kernel(I_max= parameters$general_params$Imax_ik, l_pred = 2.75, l_prey = LL) #predation Larvae #l_pred abg of larvae cod
    dL.dt = gE*E + IL*L - mL*L - gL*L - pL*L

    consumed_plankton = consumed_plankton + IL*L #updates plancton consumed by larvae

    promoting_L = gL*L
    produced_eggs = 0

    promoting_f = 0.5*promoting_L
    promoting_sizeClass = 0
    mol_i = 0

    mortality_eggs = 0

    #Juvenile or adult shrimp F
    for(i in 1:N_max_F){
      I_i_f = ingestion_rate_b(Te, size_mean_F[i], P, parameters$general_params, parameters$Fem_params)
      #s_i = spawning_rate_b(size_mean_F[i], Te, parameters$Fem_params)
      s_i = spawning_rate_b(size_mean_F[i])
      m_i = respiration_rate_b(Te, size_mean_F[i], parameters$Fem_params)
      g_i = shift_next_sizeClass(size_mean_F[i], Te, parameters$Fem_params,size_width=size_width)
      mol_i = molting_fraction(size_mean_F[i], Te)
      fm_i = parameters$general_params$Fi*monthly_Feffort[month]*sel_probit(size_mean_F[i]) #fishing mortality
      pL_i = Pred*ingestion_kernel(I_max= parameters$general_params$Imax_ik, l_pred = 2.75, l_prey = size_mean_F[i]) #predation Larvae #l_pred abg of larvae cod
      aging_i = 0.0
      if(i == N_max_F) aging_i = parameters$general_params$a_mu

      dBF.dt[i] = promoting_f + promoting_sizeClass + BF[i]*(I_i_f- m_i - mol_i*(1/convertL_to_W(i))*s_i - g_i - fm_i - pL_i - aging_i) #old spawning: mol_i*(1/convertL_to_W(i))*s_i
      promoting_sizeClass = g_i*BF[i]
      produced_eggs = produced_eggs + s_i*mol_i*BF[i]/convertL_to_W(i)
      mortality_eggs = s_i*mol_i*( BF[i]*fm_i + BF[i]*pL_i)
      consumed_plankton = consumed_plankton + I_i_f*BF[i]
      promoting_f = 0
    }

    promoting = 0.5*promoting_L

    #Juvenile or adult shrimp M
    for(i in 1:N_max_M){
      I_i = ingestion_rate_b(Te, size_mean_M[i], P, parameters$general_params, parameters$M_params)
      m_i = respiration_rate_b(Te, size_mean_M[i], parameters$M_params)
      g_i = shift_next_sizeClass(size_mean_M[i], Te, parameters$M_params,size_width=size_width)
      fm_i = parameters$general_params$Fi*monthly_Feffort[month]*sel_probit(size_mean_M[i]) #fishing mortality
      pL_i = Pred*ingestion_kernel(I_max= parameters$general_params$Imax_ik, l_pred = 2.75, l_prey = size_mean_M[i]) #predation Larvae #l_pred abg of larvae cod
      aging_i = 0.0
      if(i == N_max_M) aging_i = parameters$general_params$a_mu

      dBM.dt[i] = promoting + BM[i]*(I_i- m_i - g_i - fm_i - pL_i - aging_i)
      promoting = g_i*BM[i]
      consumed_plankton = consumed_plankton + I_i*BM[i]

    }

    #Eggs
    dE.dt =  produced_eggs - gE*E - mortality_eggs

    #Plankton
    #dP.dt = max(0, new_food(t) - consumed_plankton)
    dP.dt = new_food(t) - consumed_plankton

    return(list(c(dP.dt, dE.dt, dL.dt,
           dBF.dt, dBM.dt,
           dPred.dt)))

  }

  sol = euler(y = state, times = t, func = system.equations, parms = parameters)#, method=the.method)
  sol = as.data.frame(sol)

  start_date <- as.POSIXct(temperature_dataSet$date_time[1], format="%Y-%m-%d", tz = "UTC")

  sol <- mutate(sol, dateTime = start_date + sol$time * 86400) #86400 seconds in one day

  return(sol)
}

#.v5 is teh version used in the thesis
solver_sizeClass.v5 = function(t, state, parameters, temperature_dataSet){
  system.equations = function(t, state, parameters) {
    #following line avoid negative values in state variables

    state[state < 0] = 0 # shift negative biomasses to 0
    list2env(as.list(state), envir = environment())  # re-assign the corrected state variables

    BF = state[grep("^BF", names(state))]#extract all elements from state whose names start with "BF"
    BM = state[grep("^BM", names(state))]

    dBF.dt = BF * 0
    dBM.dt = BM * 0

    ttl_biomass = L + sum(BF) + sum(BM)

    Te = temperature_funcSolver(temperature_dataSet, t) # getting temperature for day t from temperature_dataSet
    month = month(temperature_dataSet$date_time[t+1]) # t starts in 0, but indices in R start at 1
    year = year(temperature_dataSet$date_time[t+1]) # t starts in 0, but indices in R start at 1

    consumed_plankton = 0 #var related to dP/dt (plancton diff eq. )
    food_consumption = 0 #var related to growth of shrimp through ingestion rate

    # predator
    #dPred.dt = new_Bpredator(t) - Pred*0.15 #0.08 mortality

    #LARVAE
    gE = hatch_eggs(Te)
    IL = ingestion_rate_b(Te, LJ, P*(L/ttl_biomass), parameters$general_params, parameters$Fem_params)#note that LJ is arbitrarily chossen as LL gives unrealistc values
    mL = respiration_rate_b(Te,LL, parameters$Fem_params)
    gL = shiftTo_juvenile(Te)
    eta = max( 0.01, eta_J(t, J = 190, power_of = 2))
    Pred = ST_predation(t, Te, L, eta, beta = 15)#eta = eta_J(t, J = 90, power_of = 4))
    pL = Pred*ingestion_kernel(I_max= parameters$general_params$Imax_ik, l_pred = 2.75, l_prey = LL) #predation Larvae #l_pred abg of larvae cod

    food_consumption = min(P*(L/ttl_biomass), IL*L)

    dL.dt = gE*E + food_consumption - mL*L - gL*L - pL*L

    consumed_plankton = consumed_plankton + food_consumption #updates plankton consumed by larvae

    promoting_L = gL*L
    produced_eggs = 0

    promoting_f = 0.5*promoting_L
    #promoting_sizeClass = 0
    mol_i = 0

    mortality_eggs = 0

    bycatch = 0 #fishery of undersized classes
    fishery_catch = 0

    food_consumption_f = 0
    #Juvenile or adult shrimp F
    for(i in 1:N_max_F){
      I_i_f = ingestion_rate_b(Te, size_mean_F[i], P*(BF[i]/ttl_biomass), parameters$general_params, parameters$Fem_params)
      s_i = spawning_rate_b(size_mean_F[i])
      m_i = respiration_rate_b(Te, size_mean_F[i], parameters$Fem_params)
      g_i = shift_next_sizeClass(size_mean_F[i], Te, parameters$Fem_params,size_width=size_width)
      mol_i = molting_fraction(size_mean_F[i], Te)
      fm_i = (parameters$general_params$Fi * ttl_yearly_landing$rel_landing[ttl_yearly_landing$year == year])*monthly_factors.v[month]*sel_probit(size_mean_F[i], L50 = parameters$general_params$L50, SR = parameters$general_params$SR) #fishing mortality
      eta = eta_J(t, J = 100, power_of = 2)
      Pred = ST_predation(t, Te, BF[i], eta, beta = 12)
      pL_i = Pred*ingestion_kernel(I_max= parameters$general_params$Imax_ik, l_pred = 15, l_prey = size_mean_F[i]) #predation Larvae #l_pred abg of larvae cod
      aging_i = 0.0
      if(i == N_max_F) aging_i = parameters$general_params$a_mu

      food_consumption_f = min( P*(BF[i]/ttl_biomass) ,food_consumption_f + I_i_f*BF[i])

      dBF.dt[i] = promoting_f + food_consumption_f - BF[i]*(m_i + mol_i*(1/convertL_to_W(size_mean_F[i]))*s_i + g_i + fm_i + pL_i + aging_i) #old spawning: mol_i*(1/convertL_to_W(i))*s_i
      #promoting_sizeClass = g_i*BF[i]
      produced_eggs = produced_eggs + s_i*mol_i*BF[i]/(convertL_to_W(size_mean_F[i])/1000) #convertLtoW result is g #
      mortality_eggs = s_i*mol_i*( BF[i]*fm_i + BF[i]*pL_i)
      consumed_plankton = consumed_plankton + food_consumption_f #I_i_f*BF[i]
      promoting_f = g_i*BF[i]

      if (size_mean_F[i] < 5) bycatch = bycatch + fm_i*BF[i]

      if (size_mean_F[i] == 5) {
        bycatch = bycatch + fm_i*BF[i]*0.3
        fishery_catch = fishery_catch + fm_i*BF[i]*0.7
      }

      if (size_mean_F[i] >= 6) fishery_catch = fishery_catch + fm_i*BF[i]
    }

    #MALES
    promoting = 0.5*promoting_L
    food_consumption_m = 0
    #Juvenile or adult shrimp M
    for(i in 1:N_max_M){
      I_i = ingestion_rate_b(Te, size_mean_M[i], P*(BM[i]/ttl_biomass), parameters$general_params, parameters$M_params)
      m_i = respiration_rate_b(Te, size_mean_M[i], parameters$M_params)
      g_i = shift_next_sizeClass(size_mean_M[i], Te, parameters$M_params,size_width=size_width)
      fm_i = (parameters$general_params$Fi * ttl_yearly_landing$rel_landing[ttl_yearly_landing$year == year])*monthly_factors.v[month]*sel_probit(size_mean_M[i], L50 = parameters$general_params$L50, SR = parameters$general_params$SR) #fishing mortality
      eta = eta_J(t, J = 100, power_of = 2)
      Pred = ST_predation(t, Te, BM[i], eta, beta = 12)
      pL_i = Pred*ingestion_kernel(I_max= parameters$general_params$Imax_ik, l_pred = 15, l_prey = size_mean_M[i]) #predation Larvae #l_pred abg of larvae cod
      aging_i = 0.0
      if(i == N_max_M) aging_i = parameters$general_params$a_mu

      food_consumption_m = min(  P*(BM[i]/ttl_biomass) ,food_consumption_m + I_i*BM[i])

      dBM.dt[i] = promoting + food_consumption_m - BM[i]*(m_i + g_i + fm_i + pL_i + aging_i)
      promoting = g_i*BM[i]
      consumed_plankton = consumed_plankton + food_consumption_m#I_i*BM[i]

      if (i < 5) bycatch = bycatch + fm_i*BM[i]

      if (i == N_max_M) {
        bycatch = bycatch + fm_i*BM[i]*0.3
        fishery_catch = fishery_catch + fm_i*BM[i]*0.7
      }

     # if (size_mean_M[i] >= 5) fishery_catch = fishery_catch + fm_i*BM[i]

    }

    #Eggs
    dE.dt =  produced_eggs - gE*E - mortality_eggs

    #Plankton
    #dP.dt = max(0, new_food(t) - consumed_plankton)
    dP.dt = new_food(t, Te, scale = 2) - consumed_plankton #1 es immigration rate ()

    return(list(c(dP.dt, dE.dt, dL.dt,
                  dBF.dt, dBM.dt) ,
                catch_undersized = bycatch, catch_commercial = fishery_catch ))

  }

  sol = euler(y = state, times = t, func = system.equations, parms = parameters)#, method=the.method)
  sol = as.data.frame(sol)

  start_date <- as.POSIXct(temperature_dataSet$date_time[1], format="%Y-%m-%d", tz = "UTC")

  sol <- mutate(sol, dateTime = start_date + sol$time * 86400) #86400 seconds in one day

  return(sol)
}

#in .v6 we want to make it possible to simulate several years (from 2010 to 2023), eventhough, the
#paramers fo predation and fishery change all over these years

parameters_solv.v2 = list(
  Fem_params = list(
    L_inf = 8.5,
    #K-func:
    r_max = 0.007638409,
    alpha = 0.539550,
    beta = 0.27774,
    #Biomass shifts
    a_est_factor = 0.0947, #0.07254394,
    b_est_factor = 0.811, #0.9319859,
    intercept =  -0.5243291 #-0.4968367 #-0.522060 + L_mean*0.006774
  ),
  M_params = list(
    L_inf = 5.5,
    #K-func:
    r_max = 0.01585687,
    alpha = 0.634518,
    beta = 0.29947,
    #Biomass shifts
    intercept = -0.5046928,#-0.5066616,
    a_est_factor = 0.188, #0.1192449,
    b_est_factor = 0.842 #1.328571
  ),
  general_params = list(
    #func response
    m = 3,
    const_c = 0.01,
    #Fishery
    Fi = 0.05, #only as example for the time being
    #L50 = 4.02,#applies to a mesh size of cod-end of 22 mm
    #SR = 0.82, #applies to a mesh size of cod-end of 22 mm
    #Ingestion rate & others
    #Imax_ik = 0.1,
    a_mu = 0.000913, #aging daily mortality 1/3years for fems, 2 years for males
    L_mat = 5.5 #size at maturation fem [cm]
  )
)



#Following version .V6 is not finished and doen'st work yet
#to represent the slightly growing population of predation:

solver_sizeClass.v6 = function(t, state, parameters, temperature_dataSet, PF_dataset){
  system.equations = function(t, state, parameters) {
    #following line avoid negative values in state variables

    state[state < 0] = 0 # shift negative biomasses to 0
    list2env(as.list(state), envir = environment())  # re-assign the corrected state variables

    BF = state[grep("^BF", names(state))]#extract all elements from state whose names start with "BF"
    BM = state[grep("^BM", names(state))]

    dBF.dt = BF * 0
    dBM.dt = BM * 0

    ttl_biomass = L + sum(BF) + sum(BM)

    Te = temperature_funcSolver(temperature_dataSet, t) # getting temperature for day t from temperature_dataSet
    month = month(temperature_dataSet$date_time[t+1]) # t starts in 0, but indices in R start at 1
    year = year(temperature_dataSet$date_time[t+1]) # t starts in 0, but indices in R start at 1

    I_max.i = PF_dataset$Imax[PF_dataset$year == year]
    L50.i = PF_dataset$L50[PF_dataset$year == year]
    SR.i = PF_dataset$SR[PF_dataset$year == year]

    consumed_plankton = 0 #var related to dP/dt (plancton diff eq. )
    food_consumption = 0 #var related to growth of shrimp through ingestion rate

    # predator
    #dPred.dt = new_Bpredator(t) - Pred*0.15 #0.08 mortality

    #LARVAE
    gE = hatch_eggs(Te)
    IL = ingestion_rate_b(Te, LJ, P*(L/ttl_biomass), parameters$general_params, parameters$Fem_params)#note that LJ is arbitrarily chossen as LL gives unrealistc values
    mL = respiration_rate_b(Te,LL, parameters$Fem_params)
    gL = shiftTo_juvenile(Te)
    eta = max( 0.01, eta_J(t, J = 190, power_of = 2))
    Pred = ST_predation(t, Te, L, eta, beta = 15)#eta = eta_J(t, J = 90, power_of = 4))
    pL = Pred*ingestion_kernel(I_max= I_max.i, l_pred = 2.75, l_prey = LL) #predation Larvae #l_pred abg of larvae cod

    food_consumption = min(P*(L/ttl_biomass), IL*L)

    dL.dt = gE*E + food_consumption - mL*L - gL*L - pL*L

    consumed_plankton = consumed_plankton + food_consumption #updates plankton consumed by larvae

    promoting_L = gL*L
    produced_eggs = 0

    promoting_f = 0.5*promoting_L
    #promoting_sizeClass = 0
    mol_i = 0

    mortality_eggs = 0

    bycatch = 0 #fishery of undersized classes
    fishery_catch = 0

    food_consumption_f = 0
    #Juvenile or adult shrimp F
    for(i in 1:N_max_F){
      I_i_f = ingestion_rate_b(Te, size_mean_F[i], P*(BF[i]/ttl_biomass), parameters$general_params, parameters$Fem_params)
      s_i = spawning_rate_b(size_mean_F[i])
      m_i = respiration_rate_b(Te, size_mean_F[i], parameters$Fem_params)
      g_i = shift_next_sizeClass(size_mean_F[i], Te, parameters$Fem_params,size_width=size_width)
      mol_i = molting_fraction(size_mean_F[i], Te)
      fm_i = (parameters$general_params$Fi * ttl_yearly_landing$rel_landing[ttl_yearly_landing$year == year])*monthly_factors.v[month]*sel_probit(size_mean_F[i], L50 = L50.i, SR = SR.i) #fishing mortality
      eta = eta_J(t, J = 100, power_of = 2)
      Pred = ST_predation(t, Te, BF[i], eta, beta = 12)
      pL_i = Pred*ingestion_kernel(I_max= I_max.i, l_pred = 15, l_prey = size_mean_F[i]) #predation Larvae #l_pred abg of larvae cod
      aging_i = 0.0
      if(i == N_max_F) aging_i = parameters$general_params$a_mu

      food_consumption_f = min( P*(BF[i]/ttl_biomass) ,food_consumption_f + I_i_f*BF[i])

      dBF.dt[i] = promoting_f + food_consumption_f - BF[i]*(m_i + mol_i*(1/convertL_to_W(size_mean_F[i]))*s_i + g_i + fm_i + pL_i + aging_i) #old spawning: mol_i*(1/convertL_to_W(i))*s_i
      #promoting_sizeClass = g_i*BF[i]
      produced_eggs = produced_eggs + s_i*mol_i*BF[i]/(convertL_to_W(size_mean_F[i])/1000) #convertLtoW result is g #
      mortality_eggs = s_i*mol_i*( BF[i]*fm_i + BF[i]*pL_i)
      consumed_plankton = consumed_plankton + food_consumption_f #I_i_f*BF[i]
      promoting_f = g_i*BF[i]

      if (size_mean_F[i] < 5) bycatch = bycatch + fm_i*BF[i]

      if (size_mean_F[i] == 5) {
        bycatch = bycatch + fm_i*BF[i]*0.3
        fishery_catch = fishery_catch + fm_i*BF[i]*0.7
      }

      if (size_mean_F[i] >= 6) fishery_catch = fishery_catch + fm_i*BF[i]
    }

    #MALES
    promoting = 0.5*promoting_L
    food_consumption_m = 0
    #Juvenile or adult shrimp M
    for(i in 1:N_max_M){
      I_i = ingestion_rate_b(Te, size_mean_M[i], P*(BM[i]/ttl_biomass), parameters$general_params, parameters$M_params)
      m_i = respiration_rate_b(Te, size_mean_M[i], parameters$M_params)
      g_i = shift_next_sizeClass(size_mean_M[i], Te, parameters$M_params,size_width=size_width)
      fm_i = (parameters$general_params$Fi * ttl_yearly_landing$rel_landing[ttl_yearly_landing$year == year])*monthly_factors.v[month]*sel_probit(size_mean_M[i], L50 = L50.i, SR = SR.i) #fishing mortality
      eta = eta_J(t, J = 100, power_of = 2)
      Pred = ST_predation(t, Te, BM[i], eta, beta = 12)
      pL_i = Pred*ingestion_kernel(I_max= I_max.i, l_pred = 15, l_prey = size_mean_M[i]) #predation Larvae #l_pred abg of larvae cod
      aging_i = 0.0
      if(i == N_max_M) aging_i = parameters$general_params$a_mu

      food_consumption_m = min(  P*(BM[i]/ttl_biomass) ,food_consumption_m + I_i*BM[i])

      dBM.dt[i] = promoting + food_consumption_m - BM[i]*(m_i + g_i + fm_i + pL_i + aging_i)
      promoting = g_i*BM[i]
      consumed_plankton = consumed_plankton + food_consumption_m#I_i*BM[i]

      if (i < 5) bycatch = bycatch + fm_i*BM[i]

      if (i == N_max_M) {
        bycatch = bycatch + fm_i*BM[i]*0.3
        fishery_catch = fishery_catch + fm_i*BM[i]*0.7
      }

      # if (size_mean_M[i] >= 5) fishery_catch = fishery_catch + fm_i*BM[i]

    }

    #Eggs
    dE.dt =  produced_eggs - gE*E - mortality_eggs

    #Plankton
    #dP.dt = max(0, new_food(t) - consumed_plankton)
    dP.dt = new_food(t, Te, scale = 2) - consumed_plankton #1 es immigration rate ()

    return(list(c(dP.dt, dE.dt, dL.dt,
                  dBF.dt, dBM.dt) ,
                catch_undersized = bycatch, catch_commercial = fishery_catch ))

  }

  sol = euler(y = state, times = t, func = system.equations, parms = parameters)#, method=the.method)
  sol = as.data.frame(sol)

  start_date <- as.POSIXct(temperature_dataSet$date_time[1], format="%Y-%m-%d", tz = "UTC")

  sol <- mutate(sol, dateTime = start_date + sol$time * 86400) #86400 seconds in one day

  return(sol)
}

