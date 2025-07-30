################################################
## Following functions are meant to be        ##
## used for solving differential eq. systems  ##
## Source: Andrea F.                          ##
################################################

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
    Imax_ik = 0.1,
    a_mu = 0.000913, #aging daily mortality 1/3years
    L_mat = 5.5 #size at maturation fem [cm]
  )
)

# Relative fishing effort factors for each month (Temming, 2011)
monthly_Feffort = c(0.19, 0.2, 0.86, 1.6, 1.39, 1.26, 1.19, 1.25, 1.27, 1.26, 1.09, 0.45)



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
#' @returns alpha value [gr C (prey) m-2] , m2 considering 10 cm water depth as
#' @export
#'
#' @examples alpha_igr(2)
alpha_igr = function(w){
  return(0.064327*w^(0.23))
}


#' Growth in weight for an specific size class 'L' in dependence of
#' food resource availability.
#'
#'
#' @param temperature Temperature [°C]
#' @param L size class [cm]
#' @param Food Plancton (variable state)
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
  F_safe = max(Food, 1e-6) #avoids negative values of P (this condition is better to put here than on the dP/dt to avoid inconsistency in the solver)
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
  avg_weight_egg = 17.725*10^-6 #[gr] ~avg 17.725 microgramm over the year (from table from seasonal changes on eggs Urzua)
  return(temming_nr_egg*avg_weight_egg)
}



spawning_rate_b_old = function(L, temperature, sex_params){
  s = 0
  intercept = 0.6407258 #3.49198
  factor = 4.336045 #0.95393
  #sex = 'F'
  #if(L>5) s = convertL_to_W(L)*K_func(temperature)*3*epsilon #molting_fraction(L*10, T) * convertL_to_W(L)*K_func(T)*3*epsilon #here L for Temming in mm
  if(L>parameters_solv$general_params$L_mat) {
    s = intercept + convertL_to_W(L)*K_func_briere(temperature, sex_params)*3*factor
  }
  #print(paste("spawning_rate_b sucsessful", s) )
  return(max(1,s))
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
  nominator = exp( (1.349/SR) * (L - L50))
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


eta_J <- function(t_day) {
  cos_term <- cos(2 * pi * (t_day - 150) / 365)
  eta <- 0.5 + 0.5*cos_term  # gives values from 0 to 2
  return(eta^2)         # seasonal amplification
}



ST_predation = function(t_day, Te, B_s){
  mu_ref = 0.025 #day^-1 #same as mesozooplankton paper
  f_t = 2
  sigma = 0.5
  eta = eta_J(t_day)
  beta = 18 # mesozooplankton paper = 18
  gamma = 3.7469 #m2 (Kg d)^-1 #original 0.1 m2 (molC d)^-1
  return(mu_ref*f_t*sigma*(gamma*B_s + beta*eta ))
}



solver_sizeClass_sex = function(t, state, parameters, temperature_dataSet){
  #params_f = sex_parameters_func('F')
  #params_m = sex_parameters_func('M')

  system.equations = function(t, state, parameters) {
    with(as.list(c(state, parameters)), {

      #following 2 lines avoid negative values in state variables
      state[state < 0] <- 0
      list2env(as.list(state), envir = environment())  # re-assign the corrected state variables

      Te = temperature_funcSolver(temperature_dataSet, t)

      #LARVAE
      gE = hatch_eggs(Te)
      IL = ingestion_rate_b(Te, LL, P, 'F')
      mL = respiration_rate_b(Te, LL, 'F')
      gL = shiftTo_juvenile(Te)
      dL.dt = gE*E + IL*L - mL*L - gL*L

      #Juv I F
      IJ_f = ingestion_rate_b(Te, LJ, P, 'F')
      mJ_f = respiration_rate_b(Te, LJ, 'F')
      gJI_f = shift_next_sizeClass(LJ, Te, 'F')
      dJ_f.dt = gL*L*0.5 + IJ_f*J_f - mJ_f*J_f - gJI_f*J_f

      #Juv I M
      IJ_m = ingestion_rate_b(Te, LJ, P, 'M')
      mJ_m = respiration_rate_b(Te, LJ, 'M')
      gJI_m = shift_next_sizeClass(LJ, Te, 'M')
      dJ_m.dt = gL*L*0.5 + IJ_m*J_m - mJ_m*J_m - gJI_m*J_m

      #Juv II F
      IJ2_f = ingestion_rate_b(Te,LJ2,P, 'F')
      mJ2_f = respiration_rate_b(Te, LJ2, 'F')
      gJ2_f = shift_next_sizeClass(LJ2, Te, 'F')
      dJ2_f.dt = gJI_f*J_f + IJ2_f*J2_f - mJ2_f*J2_f - gJ2_f*J2_f

      #Juv II M
      IJ2_m = ingestion_rate_b(Te,LJ2,P, 'M')
      mJ2_m = respiration_rate_b(Te, LJ2, 'M')
      gJ2_m = shift_next_sizeClass(LJ2, Te, 'M')
      dJ2_m.dt = gJI_m*J_m + IJ2_m*J2_m - mJ2_m*J2_m - gJ2_m*J2_m

      #Juv III F
      IJ3_f = ingestion_rate_b(Te, LJ3, P, 'F')
      mJ3_f = respiration_rate_b(Te, LJ3, 'F')
      gJ3_f = shift_next_sizeClass(LJ3, Te, 'F')
      dJ3_f.dt = gJ2_f*J2_f + IJ3_f*J3_f - mJ3_f*J3_f - gJ3_f*J3_f

      #Juv III M
      IJ3_m = ingestion_rate_b(Te, LJ3, P, 'M')
      mJ3_m = respiration_rate_b(Te, LJ3, 'M')
      gJ3_m = shift_next_sizeClass(LJ3, Te, 'M')
      dJ3_m.dt = gJ2_m*J2_m + IJ3_m*J3_m - mJ3_m*J3_m - gJ3_m*J3_m


      #Juv IV F
      IJ4_f = ingestion_rate_b(Te, LJ4, P, 'F')
      mJ4_f = respiration_rate_b(Te, LJ4, 'F')
      gJ4_f = shift_next_sizeClass(LJ4, Te, 'F')
      dJ4_f.dt = gJ3_f*J3_f + IJ4_f*J4_f - mJ4_f*J4_f - gJ4_f*J4_f

      #Juv IV M
      IJ4_m = ingestion_rate_b(Te, LJ4, P, 'M')
      mJ4_m = respiration_rate_b(Te, LJ4, 'M')
      gJ4_m = shift_next_sizeClass(LJ4, Te, 'M')
      dJ4_m.dt = gJ3_m*J3_m + IJ4_m*J4_m - mJ4_m*J4_m - gJ4_m*J4_m


      #Juv V F
      IJ5_f = ingestion_rate_b(Te, LJ5, P, 'F')
      mJ5_f = respiration_rate_b(Te, LJ5, 'F')
      gJ5_f = shift_next_sizeClass(LJ5, Te, 'F')
      dJ5_f.dt = gJ4_f*J4_f + IJ5_f*J5_f - mJ5_f*J5_f - gJ5_f*J5_f

      #Juv V M
      IJ5_m = ingestion_rate_b(Te, LJ5, P, 'M')
      mJ5_m = respiration_rate_b(Te, LJ5, 'M')
      gJ5_m = shift_next_sizeClass(LJ5, Te, 'M')
      dJ5_m.dt = gJ4_m*J4_m + IJ5_m*J5_m - mJ5_m*J5_m - gJ5_m*J5_m


      #Adult I F
      IA1_f = ingestion_rate_b(Te, LA1, P, 'F')
      sA1_f = spawning_rate_b(LA1, Te)
      mA1_f = respiration_rate_b(Te, LA1, 'F')
      gA1_f = shift_next_sizeClass(LA1, Te, 'F')
      molA1_f = molting_fraction(LA1, Te)
      dA1_f.dt = gJ5_f*J5_f + IA1_f*A1_f - mA1_f*A1_f - sA1_f*A1_f*molA1_f - gA1_f*A1_f

      #Adult I M
      IA1_m = ingestion_rate_b(Te, LA1, P, 'M')
      mA1_m = respiration_rate_b(Te, LA1, 'M')
      gA1_m = shift_next_sizeClass(LA1, Te, 'M')
      dA1_m.dt = gJ5_m*J5_m + IA1_m*A1_m - mA1_m*A1_m - gA1_m*A1_m #they actiually don't growth no a next size class, but let's say this is mortality, they growth old


      #Adult II F (only F reach this size class)
      IA2 = ingestion_rate_b(Te, LA2, P, 'F')
      sA2 = spawning_rate_b(LA2, Te)
      mA2 = respiration_rate_b(Te, LA2, 'F')
      gA2 = shift_next_sizeClass(LA2, Te, 'F')
      molA2 = molting_fraction(LA2, Te)
      dA2.dt = gA1_f*A1_f + IA2*A2 - mA2*A2 - sA2*A2*molA2 - gA2*A2



      #Adult III (only F reach this size class)
      IA3 = ingestion_rate_b(Te,LA3, P, 'F')
      sA3= spawning_rate_b(LA3, Te)
      mA3 = respiration_rate_b(Te, LA3, 'F')
      molA3 = molting_fraction(LA3, Te)
      dA3.dt = gA2*A2 + IA3*A3 - mA3*A3 - sA3*A3*molA3

      #Eggs
      dE.dt =  sA1_f*A1_f*molA1_f + sA2*A2*molA2 + sA3*A3*molA3 - gE*E # mA1*E: for adults, m equals cero because this is transfered to the spawning.therefore ake only sense to add mu of adults related to fishery (?)

      #Plancton
      dP.dt = new_food(t) - IL*L
            -  IJ_f*J_f - IJ_m*J_m  - IJ2_f*J2_f - IJ2_m*J2_m - IJ3_f*J3_f- IJ3_m*J3_m
            - IJ4_f*J4_f - IJ4_m*J4_m - IJ5_f*J5_f - IJ5_m*J5_m
            - IA1_f*A1_f - IA1_m*A1_m
            - IA2*A2 - IA3*A3

      list(c( dP.dt, dE.dt, dL.dt,
              dJ_f.dt, dJ2_f.dt, dJ3_f.dt, dJ4_f.dt, dJ5_f.dt,
              dJ_m.dt, dJ2_m.dt, dJ3_m.dt, dJ4_m.dt, dJ5_m.dt,
              dA1_f.dt, dA1_m.dt,
              dA2.dt, dA3.dt))
    })
  }


  sol = ode(y = state, times = t, func = system.equations, parms = parameters, method=the.method)
  sol = as.data.frame(sol)

  return(sol)
}


solver_sizeClass_sex.v2 = function(t, state, parameters, temperature_dataSet){
  system.equations = function(t, state, parameters) {
    #following line avoid negative values in state variables

    state[state < 0] <- 0
    list2env(as.list(state), envir = environment())  # re-assign the corrected state variables

    Te = temperature_funcSolver(temperature_dataSet, t)

    #LARVAE
    gE = hatch_eggs(Te)
    IL = ingestion_rate_b(Te, LL, P, parameters$general_params, parameters$Fem_params)
    mL = respiration_rate_b(Te,LL, parameters$Fem_params)
    gL = shiftTo_juvenile(Te)
    dL.dt = gE*E + IL*L - mL*L - gL*L

    #Juv I F
    IJ_f = ingestion_rate_b(Te, LJ, P, parameters$general_params, parameters$Fem_params)
    mJ_f = respiration_rate_b(Te, LJ, parameters$Fem_params)
    gJI_f = shift_next_sizeClass(LJ, Te, parameters$Fem_params)
    dJ_f.dt = gL*L*0.5 + IJ_f*J_f - mJ_f*J_f - gJI_f*J_f

    #Juv I M
    IJ_m = ingestion_rate_b(Te, LJ, P, parameters$general_params, parameters$M_params)
    mJ_m = respiration_rate_b(Te, LJ, parameters$M_params)
    gJI_m = shift_next_sizeClass(LJ, Te, parameters$M_params)
    dJ_m.dt = gL*L*0.5 + IJ_m*J_m - mJ_m*J_m - gJI_m*J_m

    #Juv II F
    IJ2_f = ingestion_rate_b(Te, LJ2, P, parameters$general_params, parameters$Fem_params)
    mJ2_f = respiration_rate_b(Te, LJ2, parameters$Fem_params)
    gJ2_f = shift_next_sizeClass(LJ2, Te, parameters$Fem_params)
    dJ2_f.dt = gJI_f*J_f + IJ2_f*J2_f - mJ2_f*J2_f - gJ2_f*J2_f


    #Juv II M
    IJ2_m = ingestion_rate_b(Te, LJ2, P, parameters$general_params, parameters$M_params)
    mJ2_m = respiration_rate_b(Te, LJ2, parameters$M_params)
    gJ2_m = shift_next_sizeClass(LJ2, Te, parameters$M_params)
    dJ2_m.dt = gJI_m*J_m + IJ2_m*J2_m - mJ2_m*J2_m - gJ2_m*J2_m

    #Juv III F
    IJ3_f = ingestion_rate_b(Te, LJ3, P, parameters$general_params, parameters$Fem_params)
    mJ3_f = respiration_rate_b(Te, LJ3, parameters$Fem_params)
    gJ3_f = shift_next_sizeClass(LJ3, Te, parameters$Fem_params)
    fmJ3_f = parameters$general_params$Fi*sel_probit(LJ3) #fishing mortality
    dJ3_f.dt = gJ2_f*J2_f + IJ3_f*J3_f - mJ3_f*J3_f - gJ3_f*J3_f - fmJ3_f*J3_f

    #Juv III M
    IJ3_m = ingestion_rate_b(Te, LJ3, P, parameters$general_params, parameters$M_params)
    mJ3_m = respiration_rate_b(Te, LJ3, parameters$M_params)
    gJ3_m = shift_next_sizeClass(LJ3, Te, parameters$M_params)
    fmJ3_m = parameters$general_params$Fi*sel_probit(LJ3) #fishing mortality
    dJ3_m.dt = gJ2_m*J2_m + IJ3_m*J3_m - mJ3_m*J3_m - gJ3_m*J3_m - fmJ3_m*J3_m


    #Juv IV F
    IJ4_f = ingestion_rate_b(Te, LJ4, P, parameters$general_params, parameters$Fem_params)
    mJ4_f = respiration_rate_b(Te, LJ4, parameters$Fem_params)
    gJ4_f = shift_next_sizeClass(LJ4, Te, parameters$Fem_params)
    fmJ4_f = parameters$general_params$Fi*sel_probit(LJ4) #fishing mortality
    dJ4_f.dt = gJ3_f*J3_f + IJ4_f*J4_f - mJ4_f*J4_f - gJ4_f*J4_f - fmJ4_f*J4_f

    #Juv IV M
    IJ4_m = ingestion_rate_b(Te, LJ4, P, parameters$general_params, parameters$M_params)
    mJ4_m = respiration_rate_b(Te, LJ4, parameters$M_params)
    gJ4_m = shift_next_sizeClass(LJ4, Te, parameters$M_params)
    #HERE FISHERY IS MISSING !!!
    dJ4_m.dt = gJ3_m*J3_m + IJ4_m*J4_m - mJ4_m*J4_m - gJ4_m*J4_m


    #Juv V F
    IJ5_f = ingestion_rate_b(Te, LJ5, P, parameters$general_params, parameters$Fem_params)
    mJ5_f = respiration_rate_b(Te, LJ5, parameters$Fem_params)
    gJ5_f = shift_next_sizeClass(LJ5, Te, parameters$Fem_params)
    fmJ5_f = parameters$general_params$Fi*sel_probit(LJ5) #fishing mortality
    dJ5_f.dt = gJ4_f*J4_f + IJ5_f*J5_f - mJ5_f*J5_f - gJ5_f*J5_f - fmJ5_f*J5_f

    #Juv V M
    IJ5_m = ingestion_rate_b(Te, LJ5, P, parameters$general_params, parameters$M_params)
    mJ5_m = respiration_rate_b(Te, LJ5, parameters$M_params)
    gJ5_m = shift_next_sizeClass(LJ5, Te, parameters$M_params)
    dJ5_m.dt = gJ4_m*J4_m + IJ5_m*J5_m - mJ5_m*J5_m - gJ5_m*J5_m


    #Adult I F
    IA1_f = ingestion_rate_b(Te, LA1, P, parameters$general_params, parameters$Fem_params)
    sA1_f = spawning_rate_b(LA1, Te, parameters$Fem_params)
    mA1_f = respiration_rate_b(Te, LA1, parameters$Fem_params)
    gA1_f = shift_next_sizeClass(LA1, Te, parameters$Fem_params)
    molA1_f = molting_fraction(LA1, Te)
    fmA1_f = parameters$general_params$Fi*sel_probit(A1_f) #fishing mortality
    dA1_f.dt = gJ5_f*J5_f + IA1_f*A1_f - mA1_f*A1_f - sA1_f*A1_f*molA1_f - gA1_f*A1_f - fmA1_f*A1_f

    #Adult I M
    IA1_m = ingestion_rate_b(Te, LA1, P, parameters$general_params, parameters$M_params)
    mA1_m = respiration_rate_b(Te, LA1,  parameters$M_params)
    gA1_m = shift_next_sizeClass(LA1, Te, parameters$M_params)
    dA1_m.dt = gJ5_m*J5_m + IA1_m*A1_m - mA1_m*A1_m - gA1_m*A1_m #they actiually don't growth no a next size class, but let's say this is mortality, they growth old

    #Adult II F (only F reach this size class)
    IA2 = ingestion_rate_b(Te, LA2, P, parameters$general_params, parameters$Fem_params)
    sA2 = spawning_rate_b(LA2, Te, parameters$Fem_params)
    mA2 = respiration_rate_b(Te, LA2, parameters$Fem_params)
    gA2 = shift_next_sizeClass(LA2, Te, parameters$Fem_params)
    molA2 = molting_fraction(LA2, Te)
    fmA2 = parameters$general_params$Fi*sel_probit(A2) #fishing mortality
    dA2.dt = gA1_f*A1_f + IA2*A2 - mA2*A2 - sA2*A2*molA2 - gA2*A2 - fmA2*A2
   # cat("dA2.dt:", dA2.dt, "\n")


    #Adult III (only F reach this size class)
    IA3 = ingestion_rate_b(Te, LA3, P, parameters$general_params, parameters$Fem_params)
    sA3= spawning_rate_b(LA3, Te, parameters$Fem_params)
    mA3 = respiration_rate_b(Te, LA3, parameters$Fem_params)
    molA3 = molting_fraction(LA3, Te)
    fmA3 = parameters$general_params$Fi*sel_probit(A3) #fishing mortality
    dA3.dt = gA2*A2 + IA3*A3 - mA3*A3 - sA3*A3*molA3 - fmA3*A3
   # cat("dA3.dt:", dA3.dt, "\n")
    #Eggs
    dE.dt =  sA1_f*A1_f*molA1_f + sA2*A2*molA2 + sA3*A3*molA3 - gE*E # mA1*E: for adults, m equals cero because this is transfered to the spawning.therefore ake only sense to add mu of adults related to fishery (?)
   # print(paste0("dE.dt: ", dE.dt))


    #Plancton
    dP.dt = ( new_food(t) - IL*L - IJ_f*J_f - IJ_m*J_m  - IJ2_f*J2_f - IJ2_m*J2_m - IJ3_f*J3_f - IJ3_m*J3_m -
            IJ4_f*J4_f - IJ4_m*J4_m - IJ5_f*J5_f - IJ5_m*J5_m -
            IA1_f*A1_f - IA1_m*A1_m -
            IA2*A2 - IA3*A3 )

   # print(paste0("dP.dt: ", dP.dt))

    list(c(dP.dt, dE.dt, dL.dt,
           dJ_f.dt, dJ2_f.dt, dJ3_f.dt, dJ4_f.dt, dJ5_f.dt,
           dJ_m.dt, dJ2_m.dt, dJ3_m.dt, dJ4_m.dt, dJ5_m.dt,
           dA1_f.dt, dA1_m.dt,
           dA2.dt, dA3.dt) )

  }


  sol = ode(y = state, times = t, func = system.equations, parms = parameters, method=the.method)
  sol = as.data.frame(sol)

  return(sol)
}


#VERSION 3 WITH PREDATION
solver_sizeClass_sex.v3 = function(t, state, parameters, temperature_dataSet){
  system.equations = function(t, state, parameters) {
    #following line avoid negative values in state variables

    state[state < 0] = 0 # shift negative biomasses to 0
    list2env(as.list(state), envir = environment())  # re-assign the corrected state variables

    Te = temperature_funcSolver(temperature_dataSet, t) # getting temperature for day t from temperature_dataSet
    month = month(temperature_dataSet$date_time[t+1]) # t starts in 0, but indices in R start at 1

    dPred.dt = new_Bpredator(t) - Pred*0.15 #0.08 mortality

    #LARVAE
    gE = hatch_eggs(Te)
    IL = ingestion_rate_b(Te, LL, P, parameters$general_params, parameters$Fem_params)
    mL = respiration_rate_b(Te,LL, parameters$Fem_params)
    gL = shiftTo_juvenile(Te)
    pL = Pred*ingestion_kernel(I_max= parameters$general_params$Imax_ik, l_pred = 2.75, l_prey = LL) #predation Larvae #l_pred abg of larvae cod
    dL.dt = gE*E + IL*L - mL*L - gL*L - pL*L

    #Juv I F
    IJ_f = ingestion_rate_b(Te, LJ, P, parameters$general_params, parameters$Fem_params)
    mJ_f = respiration_rate_b(Te, LJ, parameters$Fem_params)
    gJI_f = shift_next_sizeClass(LJ, Te, parameters$Fem_params, size_width = 0.4)
    dJ_f.dt = gL*L*0.5 + IJ_f*J_f - mJ_f*J_f - gJI_f*J_f

    #Juv I M
    IJ_m = ingestion_rate_b(Te, LJ, P, parameters$general_params, parameters$M_params)
    mJ_m = respiration_rate_b(Te, LJ, parameters$M_params)
    gJI_m = shift_next_sizeClass(LJ, Te, parameters$M_params, size_width = 0.4)
    dJ_m.dt = gL*L*0.5 + IJ_m*J_m - mJ_m*J_m - gJI_m*J_m

    #Juv II F
    IJ2_f = ingestion_rate_b(Te, LJ2, P, parameters$general_params, parameters$Fem_params)
    mJ2_f = respiration_rate_b(Te, LJ2, parameters$Fem_params)
    gJ2_f = shift_next_sizeClass(LJ2, Te, parameters$Fem_params)
    fmJ2_f = parameters$general_params$Fi*monthly_Feffort[month]*sel_probit(LJ2) #fishing mortality
    dJ2_f.dt = gJI_f*J_f + IJ2_f*J2_f - mJ2_f*J2_f - gJ2_f*J2_f - fmJ2_f*J2_f


    #Juv II M
    IJ2_m = ingestion_rate_b(Te, LJ2, P, parameters$general_params, parameters$M_params)
    mJ2_m = respiration_rate_b(Te, LJ2, parameters$M_params)
    gJ2_m = shift_next_sizeClass(LJ2, Te, parameters$M_params)
    fmJ2_m = parameters$general_params$Fi*monthly_Feffort[month]*sel_probit(LJ2) #fishing mortality
    dJ2_m.dt = gJI_m*J_m + IJ2_m*J2_m - mJ2_m*J2_m - gJ2_m*J2_m - fmJ2_m*J2_m

    #Juv III F
    IJ3_f = ingestion_rate_b(Te, LJ3, P, parameters$general_params, parameters$Fem_params)
    mJ3_f = respiration_rate_b(Te, LJ3, parameters$Fem_params)
    gJ3_f = shift_next_sizeClass(LJ3, Te, parameters$Fem_params)
    fmJ3_f = parameters$general_params$Fi*monthly_Feffort[month]*sel_probit(LJ3) #fishing mortality
    dJ3_f.dt = gJ2_f*J2_f + IJ3_f*J3_f - mJ3_f*J3_f - gJ3_f*J3_f - fmJ3_f*J3_f

    #Juv III M
    IJ3_m = ingestion_rate_b(Te, LJ3, P, parameters$general_params, parameters$M_params)
    mJ3_m = respiration_rate_b(Te, LJ3, parameters$M_params)
    gJ3_m = shift_next_sizeClass(LJ3, Te, parameters$M_params)
    fmJ3_m = parameters$general_params$Fi*monthly_Feffort[month]*sel_probit(LJ3) #fishing mortality
    dJ3_m.dt = gJ2_m*J2_m + IJ3_m*J3_m - mJ3_m*J3_m - gJ3_m*J3_m - fmJ3_m*J3_m


    #Juv IV F
    IJ4_f = ingestion_rate_b(Te, LJ4, P, parameters$general_params, parameters$Fem_params)
    mJ4_f = respiration_rate_b(Te, LJ4, parameters$Fem_params)
    gJ4_f = shift_next_sizeClass(LJ4, Te, parameters$Fem_params)
    fmJ4_f = parameters$general_params$Fi*monthly_Feffort[month]*sel_probit(LJ4) #fishing mortality
    dJ4_f.dt = gJ3_f*J3_f + IJ4_f*J4_f - mJ4_f*J4_f - gJ4_f*J4_f - fmJ4_f*J4_f

    #Juv IV M
    IJ4_m = ingestion_rate_b(Te, LJ4, P, parameters$general_params, parameters$M_params)
    mJ4_m = respiration_rate_b(Te, LJ4, parameters$M_params)
    gJ4_m = shift_next_sizeClass(LJ4, Te, parameters$M_params)
    fmJ4_m = parameters$general_params$Fi*monthly_Feffort[month]*sel_probit(LJ4) #fishing mortality
    dJ4_m.dt = gJ3_m*J3_m + IJ4_m*J4_m - mJ4_m*J4_m - gJ4_m*J4_m - fmJ4_m*J4_m


    #Juv V F
    IJ5_f = ingestion_rate_b(Te, LJ5, P, parameters$general_params, parameters$Fem_params)
    mJ5_f = respiration_rate_b(Te, LJ5, parameters$Fem_params)
    gJ5_f = shift_next_sizeClass(LJ5, Te, parameters$Fem_params)
    fmJ5_f = parameters$general_params$Fi*monthly_Feffort[month]*sel_probit(LJ5) #fishing mortality
    dJ5_f.dt = gJ4_f*J4_f + IJ5_f*J5_f - mJ5_f*J5_f - gJ5_f*J5_f - fmJ5_f*J5_f

    #Juv V M
    IJ5_m = ingestion_rate_b(Te, LJ5, P, parameters$general_params, parameters$M_params)
    mJ5_m = respiration_rate_b(Te, LJ5, parameters$M_params)
    gJ5_m = shift_next_sizeClass(LJ5, Te, parameters$M_params)
    fmJ5_m = parameters$general_params$Fi*monthly_Feffort[month]*sel_probit(LJ5) #fishing mortality
    dJ5_m.dt = gJ4_m*J4_m + IJ5_m*J5_m - mJ5_m*J5_m - gJ5_m*J5_m - fmJ5_m*J5_m


    #Adult I F
    IA1_f = ingestion_rate_b(Te, LA1, P, parameters$general_params, parameters$Fem_params)
    #sA1_f = spawning_rate_b(LA1, Te, parameters$Fem_params)
    sA1_f = spawning_rate_b(LA1)
    mA1_f = respiration_rate_b(Te, LA1, parameters$Fem_params)
    gA1_f = shift_next_sizeClass(LA1, Te, parameters$Fem_params)
    molA1_f = molting_fraction(LA1, Te)
    fmA1_f = parameters$general_params$Fi*monthly_Feffort[month]*sel_probit(LA1) #fishing mortality
    dA1_f.dt = gJ5_f*J5_f + IA1_f*A1_f - mA1_f*A1_f - sA1_f*A1_f*molA1_f - gA1_f*A1_f - fmA1_f*A1_f

    #Adult I M
    IA1_m = ingestion_rate_b(Te, 5.25, P, parameters$general_params, parameters$M_params)
    mA1_m = respiration_rate_b(Te, 5.25,  parameters$M_params)
    gA1_m = shift_next_sizeClass(5.25, Te, parameters$M_params)
    fmA1_m = parameters$general_params$Fi*monthly_Feffort[month]*sel_probit(LA1) #fishing mortality
    aging_mu = parameters$general_params$a_mu
    dA1_m.dt = gJ5_m*J5_m + IA1_m*A1_m - mA1_m*A1_m - aging_mu*A1_m  - fmA1_m*A1_m #they actiually don't growth no a next size class, but let's say this is mortality, they growth old

    #Adult II F (only F reach this size class)
    IA2 = ingestion_rate_b(Te, LA2, P, parameters$general_params, parameters$Fem_params)
    #sA2 = spawning_rate_b(LA2, Te, parameters$Fem_params)
    sA2 = spawning_rate_b(LA2)
    mA2 = respiration_rate_b(Te, LA2, parameters$Fem_params)
    gA2 = shift_next_sizeClass(LA2, Te, parameters$Fem_params)
    molA2 = molting_fraction(LA2, Te)
    fmA2 = parameters$general_params$Fi*monthly_Feffort[month]*sel_probit(LA2) #fishing mortality
    dA2.dt = gA1_f*A1_f + IA2*A2 - mA2*A2 - sA2*A2*molA2 - gA2*A2 - fmA2*A2


    #Adult III (only F reach this size class)
    IA3 = ingestion_rate_b(Te, LA3, P, parameters$general_params, parameters$Fem_params)
    #sA3= spawning_rate_b(LA3, Te, parameters$Fem_params)
    sA3= spawning_rate_b(LA3)
    mA3 = respiration_rate_b(Te, LA3, parameters$Fem_params)
    molA3 = molting_fraction(LA3, Te)
    fmA3 = parameters$general_params$Fi*monthly_Feffort[month]*sel_probit(LA3) #fishing mortality
    aging_mu = parameters$general_params$a_mu
    dA3.dt = gA2*A2 + IA3*A3 - mA3*A3 - sA3*A3*molA3 - fmA3*A3 - aging_mu*A3

    #Eggs
    produced_eggs =  sA1_f*A1_f*molA1_f + sA2*A2*molA2 + sA3*A3*molA3
    dE.dt = produced_eggs - gE*E # mA1*E: for adults, m equals cero because this is transfered to the spawning.therefore make only sense to add mu of adults related to fishery (?)

    dP.dt = ( new_food(t) - IL*L - IJ_f*J_f - IJ_m*J_m  - IJ2_f*J2_f - IJ2_m*J2_m - IJ3_f*J3_f - IJ3_m*J3_m -
                IJ4_f*J4_f - IJ4_m*J4_m - IJ5_f*J5_f - IJ5_m*J5_m -
                IA1_f*A1_f - IA1_m*A1_m -
                IA2*A2 - IA3*A3 )

    list(c(dP.dt, dE.dt, dL.dt,
           dJ_f.dt, dJ2_f.dt, dJ3_f.dt, dJ4_f.dt, dJ5_f.dt,
           dJ_m.dt, dJ2_m.dt, dJ3_m.dt, dJ4_m.dt, dJ5_m.dt,
           dA1_f.dt, dA1_m.dt,
           dA2.dt, dA3.dt,
           dPred.dt), produced_eggs = produced_eggs, ing_rate = IA3 , ing_rateL = IL, new_eggs=gE  )

  }


  #sol = ode(y = state, times = t, func = system.equations, parms = parameters, method=the.method)
  sol = euler(y = state, times = t, func = system.equations, parms = parameters)
  sol = as.data.frame(sol)
  start_date <- as.POSIXct(temperature_dataSet$date_time[1], format="%Y-%m-%d", tz = "UTC")
  sol <- mutate(sol, dateTime = start_date + sol$time * 86400) #86400 seconds in one day
  return(sol)
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




solver_sizeClass.v5 = function(t, state, parameters, temperature_dataSet){
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
    #dPred.dt = new_Bpredator(t) - Pred*0.15 #0.08 mortality

    #LARVAE
    gE = hatch_eggs(Te)
    IL = ingestion_rate_b(Te, LJ, P, parameters$general_params, parameters$Fem_params)#note that LJ is arbitrarily chossen as LL gives unrealistc values
    mL = respiration_rate_b(Te,LL, parameters$Fem_params)
    gL = shiftTo_juvenile(Te)
    Pred = ST_predation(t, Te, L)
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
      Pred = ST_predation(t, Te, BF[i])
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
      Pred = ST_predation(t, Te, BM[i])
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
                  dBF.dt, dBM.dt)))

  }

  sol = euler(y = state, times = t, func = system.equations, parms = parameters)#, method=the.method)
  sol = as.data.frame(sol)

  start_date <- as.POSIXct(temperature_dataSet$date_time[1], format="%Y-%m-%d", tz = "UTC")

  sol <- mutate(sol, dateTime = start_date + sol$time * 86400) #86400 seconds in one day

  return(sol)
}
