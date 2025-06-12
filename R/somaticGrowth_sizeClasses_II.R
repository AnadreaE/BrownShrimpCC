################################################
## Following functions are meant to be        ##
## used for solving differential eq. systems  ##
## Source: Andrea F.                          ##
################################################


##### constant parameter values: #####

attac_rate = 0.25
handling_time = 22
alpha_ir = 1/(attac_rate*handling_time)
#L_inf_f = 8.5
#L_inf_f = 5.5
epsilon = 0.22
m = 3
const_c = 0.01
the.method = 'rk4'

##### FUNCTIONS #####
sex_parameters_func = function(sex){
  if (sex == 'F'){
    #general:
    L_inf = 8.5
    #K-func:
    r_max = 0.007638409
    alpha = 0.539550
    beta = 0.277740
  } else if (sex == 'M'){
    #general:
    L_inf = 5.5
    #K-func:
    r_max = 0.01585687
    alpha = 0.634518
    beta = 0.299470
  }
  df = data.frame(L_inf = L_inf, r_max = r_max, alpha = alpha, beta = beta  )
  return(df)
}

#' K func based of flexTPC NB: for the time being all params only for fems
#'
#' @param temperature
#'
#' @returns  K value as function of temperature
#' @export
#'
#' @examples
K_func_briere = function(temperature, sex){
  T_min = 0.5
  T_max = 30
  if (sex == 'F'){
    r_max = 0.007638409
    alpha = 0.539550
    beta = 0.277740
  } else if (sex == 'M'){
    r_max = 0.01585687
    alpha = 0.634518
    beta = 0.299470
  }

  if(temperature < T_min) { temperature = T_min } #if this true, then T_min is negative and invalid to set to the power of alpha

  diff_min = temperature - T_min
  diff_max = T_max - temperature
  diff = T_max - T_min
  alpha_invert = 1 - alpha
  toReturn = r_max*( ((diff_min/alpha)^alpha) * ((diff_max/ alpha_invert )^alpha_invert) * (1 / diff)  )^(alpha*alpha_invert/(beta^2) )

  return(max(0, toReturn) ) #3.205227e-06 is the equivalent to K_func_briere(T_min)
}



#'Growth in weight for an specific size class 'L' in dependence of resouce availabilty. Applies for sizes between (6 and L_inf)
#'
#' @param temperature Temperature [°C]
#' @param L size class [cm]
#' @param P Plancton (variable state)
#'
#' @returns growth rate in weight [gr / time ]
#' @export
#'
#' @examples ingestion_rate(temperature = 15, L = 55, P=2, sex = 'F')

ingestion_rate_b = function(temperature, L, P, sex){
  c_div = const_c / convertL_to_W(L)
  #params = sex_parameters_func(sex)
  #L_infty = params$L_inf
  if (sex == 'F') L_infty = 8.5
  if (sex == 'M') L_infty = 5.5
  #result = m*K_func_briere(temperature, sex)*convertL_to_W(L)*L_infty*(c_div)^(1/m)*( P/(P+alpha_ir) )
  w = convertL_to_W(L)
  result = m*K_func_briere(temperature, sex)*(L_infty/L)*w^(1/4)*( P/(P+alpha_ir) )
  #m*convertL_to_W(L)*K_func(T)*L_inf/L*P/(P+h) #revised formula Andrea 04.04.25
  #1.14*K_func(T)*(L_inf/ L)^n *P/(P+h)
  return (result)
}


#' Rate of Biomass from Females that flows into egg biomass
#'
#' @param L size [cm]
#' @param T Temperature [°C]
#'
#' @returns
#' @export
#'
#' @examples
spawning_rate_b = function(L, temperature){
  s = 0
  intercept = 3.49198
  factor = 0.95393
  sex = 'F'
  #if(L>5) s = convertL_to_W(L)*K_func(temperature)*3*epsilon #molting_fraction(L*10, T) * convertL_to_W(L)*K_func(T)*3*epsilon #here L for Temming in mm
  if(L>5) {
    s = intercept + convertL_to_W(L)*K_func_briere(temperature, sex)*3*factor
  }
  return(s)
}



spawning_rate_old = function(L, temperature){
  s = 0
  if(L>5) s = convertL_to_W(L)*K_func_briere(temperature)*3*epsilon #molting_fraction(L*10, T) * convertL_to_W(L)*K_func(T)*3*epsilon #here L for Temming in mm
  return(s)
}




#' Respiration rate
#'
#' @param L size [cm]
#' @param temperature Temperature [°C]
#'
#' @returns [gr / time ]
#' @export
#'
#' @examples
respiration_rate_b = function(temperature, L, sex){

  mu = m*convertL_to_W(L)*K_func_briere(temperature, sex)*0.1# this is the rigth term of vB eq.
  #if (L>5) {
  #  mu =  convertL_to_W(L) *K_func_briere(temperature) *( 1 - molting_fraction(L, temperature)) #tbc (1-molting_fraction) all non molting fems still have a natural mortality
  #}
  return(mu)
}


respiration_rate_b_OLD = function(temperature, L){
  mu = 0
  if (temperature < 0.5){ #this funtion returns NaN when T (0.0001, 0.4)
    mu = 2.596234e-06 # equivalent to natural_mortality_b(0.5)
  } else {
    if (L>5) {} # nothing happens, mu = 0
    else{
      mu = m*convertL_to_W(L)*K_func_briere(temperature)# this is the rigth term of vB eq.
    }
  }
  return(mu)
}



#growth growth_optionA
#HERE PARAMS VALUES FOR MASCULINE ARE MISSING
shift_next_sizeClass = function(L_mean, temperature, sex){

  #factor = 0.022447 * exp(0.306237 * L_mean) #0.06734 * exp(0.30622*L_mean)
  if (sex == 'F'){
    a_est_factor = 0.07254394
    b_est_factor = 0.9319859
    intercept = -0.4968367 #-0.522060 + L_mean*0.006774
    k = K_func_briere(temperature, 'F')
  } else if (sex == 'M'){
    #calculations for a and b pending
    k = K_func_briere(temperature, 'M')
  }
  factor = a_est_factor * L_mean^b_est_factor
  return(1/ (intercept + factor*(1/k)) )
}

#Shift to next size class for Eggs
hatch_eggs = function(Te){
  Te = max(0.00001, Te)
  return(1/ (1031.34*Te^-1.345)) #return the ratio
}

#Shift to next size class for Larvae
shiftTo_juvenile = function(Te){
  Te = max(0.00001, Te)
  return(1 / ((5.5/0.00584)*Te^-1.347)) #return the ratio
}


solver_sizeClass_extended_b = function(t, state, parameters, temperature_dataSet){
  the.method = 'rk4'
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

      #Juv I
      IJ = ingestion_rate_b(Te, LJ, P, 'F')
      mJ = respiration_rate_b(Te, LJ, 'F')
      gJI = shift_next_sizeClass(LJ, Te, 'F')
      dJ.dt = gL*L + IJ*J - mJ*J - gJI*J

      #Juv II
      IJ2 = ingestion_rate_b(Te,LJ2,P, 'F')
      mJ2 = respiration_rate_b(Te, LJ2, 'F')
      gJ2 = shift_next_sizeClass(LJ2, Te, 'F')
      dJ2.dt = gJI*J + IJ2*J2 - mJ2*J2 - gJ2*J2 # - sA1*A1

      #Juv III
      IJ3 = ingestion_rate_b(Te,LJ3,P, 'F')
      mJ3 = respiration_rate_b(Te, LJ3, 'F')
      gJ3 = shift_next_sizeClass(LJ3, Te, 'F')
      dJ3.dt = gJ2*J2 + IJ3*J3 - mJ3*J3 - gJ3*J3 # - sA1*A1

      #Juv IV
      IJ4 = ingestion_rate_b(Te,LJ4,P, 'F')
      mJ4 = respiration_rate_b(Te, LJ4, 'F')
      gJ4 = shift_next_sizeClass(LJ4, Te, 'F')
      dJ4.dt = gJ3*J3 + IJ4*J4 - mJ4*J4 - gJ4*J4

      #Juv V
      IJ5 = ingestion_rate_b(Te,LJ5,P, 'F')
      mJ5 = respiration_rate_b(Te, LJ5, 'F')
      gJ5 = shift_next_sizeClass(LJ5, Te, 'F')
      dJ5.dt = gJ4*J4 + IJ5*J5 - mJ5*J5 - gJ5*J5


      #Adult I
      IA1 = ingestion_rate_b(Te,LA1,P, 'F')
      sA1= spawning_rate_b(LA1, Te)
      mA1 = respiration_rate_b(Te, LA1, 'F')
      gA1 = shift_next_sizeClass(LA1, Te, 'F')
      molA1 = molting_fraction(LA1, Te)
      dA1.dt = gJ5*J5 + IA1*A1 - mA1*A1 - sA1*A1*molA1 - gA1*A1 #+ (1/dev_tA1(T,L))*A1


      #Adult II
      IA2 = ingestion_rate_b(Te,LA2,P, 'F')
      sA2= spawning_rate_b(LA2, Te)
      mA2 = respiration_rate_b(Te, LA2, 'F')
      gA2 = shift_next_sizeClass(LA2, Te, 'F')
      molA2 = molting_fraction(LA2, Te)
      dA2.dt = gA1*A1 + IA2*A2 - mA2*A2 - sA2*A2*molA2 - gA2*A2 #+ (1/dev_tA1(T,L))*A1

      #Adult III only
      IA3 = ingestion_rate_b(Te,LA3, P, 'F')
      sA3= spawning_rate_b(LA3, Te)
      mA3 = respiration_rate_b(Te, LA3, 'F')
      molA3 = molting_fraction(LA3, Te)
      dA3.dt = gA2*A2 + IA3*A3 - mA3*A3 - sA3*A3*molA3 #+ (1/dev_tA1(T,L))*A1

      #Adult II
      dE.dt =  sA1*A1*molA1 + sA2*A2*molA2 + sA3*A3*molA3 - gE*E # mA1*E: for adults, m equals cero because this is transfered to the spawning.therefore ake only sense to add mu of adults related to fishery (?)

      #Plancton
      dP.dt = new_food(t) - IL*L - IJ*J - IJ2*J2 - IJ3*J3- IJ4*J4- IJ5*J5  - IA1*A1 - IA2*A2 - IA3*A3

      list(c( dP.dt, dE.dt, dL.dt, dJ.dt, dJ2.dt, dJ3.dt, dJ4.dt, dJ5.dt, dA1.dt, dA2.dt, dA3.dt))
    })
  }


  sol = ode(y = state, times = t, func = system.equations, parms = parameters,method=the.method)
  sol = as.data.frame(sol)

  return(sol)
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
      gJI_m = shift_next_sizeClass(LJ, Te, 'F')
      dJ_m.dt = gL*L*0.5 + IJ_m*J_m - mJ_m*J_m - gJI_m*J_m

      #Juv II F
      IJ2_f = ingestion_rate_b(Te,LJ2,P, 'F')
      mJ2_f = respiration_rate_b(Te, LJ2, 'F')
      gJ2_f = shift_next_sizeClass(LJ2, Te, 'F')
      dJ2_f.dt = gJI_f*J_f + IJ2_f*J2_f - mJ2_f*J2_f - gJ2_f*J2_f

      #Juv II M
      IJ2_m = ingestion_rate_b(Te,LJ2,P, 'M')
      mJ2_m = respiration_rate_b(Te, LJ2, 'M')
      gJ2_m = shift_next_sizeClass(LJ2, Te, 'F')
      dJ2_m.dt = gJI_m*J_m + IJ2_m*J2_m - mJ2_m*J2_m - gJ2_m*J2_m

      #Juv III F
      IJ3_f = ingestion_rate_b(Te, LJ3, P, 'F')
      mJ3_f = respiration_rate_b(Te, LJ3, 'F')
      gJ3_f = shift_next_sizeClass(LJ3, Te, 'F')
      dJ3_f.dt = gJ2_f*J2_f + IJ3_f*J3_f - mJ3_f*J3_f - gJ3_f*J3_f

      #Juv III M
      IJ3_m = ingestion_rate_b(Te, LJ3, P, 'M')
      mJ3_m = respiration_rate_b(Te, LJ3, 'M')
      gJ3_m = shift_next_sizeClass(LJ3, Te, 'F')
      dJ3_m.dt = gJ2_m*J2_m + IJ3_m*J3_m - mJ3_m*J3_m - gJ3_m*J3_m


      #Juv IV F
      IJ4_f = ingestion_rate_b(Te, LJ4, P, 'F')
      mJ4_f = respiration_rate_b(Te, LJ4, 'F')
      gJ4_f = shift_next_sizeClass(LJ4, Te, 'F')
      dJ4_f.dt = gJ3_f*J3_f + IJ4_f*J4_f - mJ4_f*J4_f - gJ4_f*J4_f

      #Juv IV M
      IJ4_m = ingestion_rate_b(Te, LJ4, P, 'M')
      mJ4_m = respiration_rate_b(Te, LJ4, 'M')
      gJ4_m = shift_next_sizeClass(LJ4, Te, 'F')
      dJ4_m.dt = gJ3_m*J3_m + IJ4_m*J4_m - mJ4_m*J4_m - gJ4_m*J4_m


      #Juv V F
      IJ5_f = ingestion_rate_b(Te, LJ5, P, 'F')
      mJ5_f = respiration_rate_b(Te, LJ5, 'F')
      gJ5_f = shift_next_sizeClass(LJ5, Te, 'F')
      dJ5_f.dt = gJ4_f*J4_f + IJ5_f*J5_f - mJ5_f*J5_f - gJ5_f*J5_f

      #Juv V M
      IJ5_m = ingestion_rate_b(Te, LJ5, P, 'M')
      mJ5_m = respiration_rate_b(Te, LJ5, 'M')
      gJ5_m = shift_next_sizeClass(LJ5, Te, 'F')
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
      gA1_m = shift_next_sizeClass(LA1, Te, 'F')
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
              dA1_f.dt,  dA1_m.dt,
              dA2.dt, dA3.dt))
    })
  }


  sol = ode(y = state, times = t, func = system.equations, parms = parameters, method=the.method)
  sol = as.data.frame(sol)

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
som_growth_thesis = function (L0, time, temp, sex){ #L[cm]
  #This funciton works only in this context i.e. to simulate growth with constant temperatures in order to
  #compare and check plausibilty of the results AND always with L0 = 0 .
  if (sex == 'F') Linf = 8.5
  if (sex == 'M') Linf = 5.5
  if (L0 > Linf ) growth = 0
  else growth = Linf- (Linf-L0)*exp(-K_func_briere(temp, sex)*time)
  #else growth = Linf*( 1 - exp(-K_func_briere(temp, sex)*time) )
  return(max(0,growth) )
}


