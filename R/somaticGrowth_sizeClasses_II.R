################################################
## Following functions are meant to be        ##
## used for solving differential eq. systems  ##
## Source: Andrea F.                          ##
################################################



##### constant parameter values: #####

attac_rate = 0.25
handling_time = 22
alpha_ir = 1/(attac_rate*handling_time)
L_inf = 8.5
epsilon = 0.22
m = 3
const_c = 0.01
the.method = 'rk4'

##### FUNCTIONS #####


#' K func based of flexTPC NB: for the time being all params only for fems
#'
#' @param temperature
#'
#' @returns  K value as function of temperature
#' @export
#'
#' @examples
K_func_briere = function(temperature){
  r_max_f = 0.007638409
  T_min = 0.5
  T_max = 30
  alpha = 0.541524
  beta = 0.276430
  if(temperature < T_min) { temperature = T_min } #if this true, then T_min is negative and invalid to set to the power of alpha

  diff_min = temperature - T_min
  diff_max = T_max - temperature
  diff = T_max - T_min
  alpha_invert = 1 - alpha
  toReturn = r_max_f*( ((diff_min/alpha)^alpha) * ((diff_max/ alpha_invert )^alpha_invert) * (1 / diff)  )^(alpha*alpha_invert/(beta^2) )

  return(max(3.205227e-06, toReturn) ) #3.205227e-06 is the equivalent to K_func_briere(T_min)
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
#' @examples ingestion_rate(temperature = 15, L = 55, P)

ingestion_rate_b = function(temperature, L, P){
  c_div = const_c / convertL_to_W(L)

  result = m*K_func_briere(temperature)*convertL_to_W(L)*L_inf*(c_div)^(1/m)*( P/(P+alpha_ir) )

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
  intercept = 3.54616
  factor = 4.40332
  #if(L>5) s = convertL_to_W(L)*K_func(temperature)*3*epsilon #molting_fraction(L*10, T) * convertL_to_W(L)*K_func(T)*3*epsilon #here L for Temming in mm
  if(L>5) {
    s = intercept + convertL_to_W(L)*K_func_briere(temperature)*3*factor
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
respiration_rate_b = function(temperature,L){
  mu = m*convertL_to_W(L)*K_func(temperature)*0.1# this is the rigth term of vB eq.
  #if (L>5) {
  #  mu =  convertL_to_W(L) *K_func_briere(temperature) *( 1 - molting_fraction(L, temperature)) #tbc (1-molting_fraction) all non molting fems still have a natural mortality
  #}
  return(mu)
}


respiration_rate_b_OLD = function(temperature,L){
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
shift_next_sizeClass = function(L_mean, temperature, sex){
  intercept = -0.4968367 #-0.522060 + L_mean*0.006774
  #factor = 0.022447 * exp(0.306237 * L_mean) #0.06734 * exp(0.30622*L_mean)
  factor = a_est * L_mean^b_est
  k = K_func_briere(temperature)
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

      Te = temperature_funcSolver(temperature_dataSet, t)

      #LARVAE
      gE = hatch_eggs(Te)
      IL = ingestion_rate_b(Te, LL, P)
      mL = respiration_rate_b(Te, LL)
      gL = shiftTo_juvenile(Te)
      dL.dt = gE*E + IL*L - mL*L - gL*L

      #Juv I
      IJ = ingestion_rate_b(Te, LJ, P)
      mJ = respiration_rate_b(Te, LJ)
      gJI = shift_next_sizeClass(LJ, Te, 'F')
      dJ.dt = gL*L + IJ*J - mJ*J - gJI*J

      #Juv II
      IJ2 = ingestion_rate_b(Te,LJ2,P)
      mJ2 = respiration_rate_b(Te, LJ2)
      gJ2 = shift_next_sizeClass(LJ2, Te, 'F')
      dJ2.dt = gJI*J + IJ2*J2 - mJ2*J2 - gJ2*J2 # - sA1*A1

      #Juv III
      IJ3 = ingestion_rate_b(Te,LJ3,P)
      mJ3 = respiration_rate_b(Te, LJ3)
      gJ3 = shift_next_sizeClass(LJ3, Te, 'F')
      dJ3.dt = gJ2*J2 + IJ3*J3 - mJ3*J3 - gJ3*J3 # - sA1*A1

      #Juv IV
      IJ4 = ingestion_rate_b(Te,LJ4,P)
      mJ4 = respiration_rate_b(Te, LJ4)
      gJ4 = shift_next_sizeClass(LJ4, Te, 'F')
      dJ4.dt = gJ3*J3 + IJ4*J4 - mJ4*J4 - gJ4*J4

      #Juv V
      IJ5 = ingestion_rate_b(Te,LJ5,P)
      mJ5 = respiration_rate_b(Te, LJ5)
      gJ5 = shift_next_sizeClass(LJ5, Te, 'F')
      dJ5.dt = gJ4*J4 + IJ5*J5 - mJ5*J5 - gJ5*J5


      #Adult I
      IA1 = ingestion_rate_b(Te,LA1,P)
      sA1= spawning_rate_b(LA1, Te)
      mA1 = respiration_rate_b(Te, LA1)
      gA1 = shift_next_sizeClass(LA1, Te, 'F')
      molA1 = molting_fraction(LA1*10, Te)
      dA1.dt = gJ5*J5 + IA1*A1 - mA1*A1 - sA1*A1*molA1 - gA1*A1 #+ (1/dev_tA1(T,L))*A1 // 0.1 Fishery


      #Adult II
      IA2 = ingestion_rate_b(Te,LA2,P)
      sA2= spawning_rate_b(LA2, Te)
      mA2 = respiration_rate_b(Te, LA2)
      molA2 = molting_fraction(LA2*10, Te)
      dA2.dt = gA1*A1 + IA2*A2 - mA2*A2 - sA2*A2*molA2 #+ (1/dev_tA1(T,L))*A1 // 0.1 Fishery

      #Adult II
      dE.dt =  sA1*A1*molA1 + sA2*A2*molA2 - gE*E # mA1*E: for adults, m equals cero because this is transfered to the spawning.therefore ake only sense to add mu of adults related to fishery (?)

      #Plancton
      dP.dt = new_food(t) - IL*L - IJ*J - IJ2*J2 - IJ3*J3- IJ4*J4- IJ5*J5  - IA1*A1 - IA2*A2

      list(c( dP.dt, dE.dt, dL.dt, dJ.dt, dJ2.dt, dJ3.dt, dJ4.dt, dJ5.dt, dA1.dt, dA2.dt))
    })
  }


  sol = ode(y = state, times = t, func = system.equations, parms = parameters,method=the.method)
  sol = as.data.frame(sol)

  return(sol)
}



#' Somatic growth in one time step of 1 day with K_func as TPC
#'
#' @param Tem temperature [°C]
#' @param sex 'F' or 'M'
#' @param L current size [cm]
#'
#' @returns new size after gorwth dependent on temperature
#' @export
#'
#' @examples som_growthAF(3, 10, "F")
som_growthAF <- function(L, Tem, sex){
  growth <- 0

  L_inf_f = 85
  L_inf_m <- 55

  if ((L > L_inf_f/10 && sex == 'F') || (L> L_inf_m/10 && sex == 'M') ){ } #nothing happens, growth = 0

  else {
    if (sex == 'F'){
      growth <- L_inf_f*(1 - exp(-K_func_briere(Tem)))
    }
    if (sex == 'M'){
      growth <- L_inf_m*(1 - exp(-K_func_briere(Tem)))
    }

  }

  return (L + growth/10)
}
