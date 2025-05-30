################################################
## Following functions are meant to be        ##
## used for solving differential eq. systems  ##
## Source: Andrea F.                          ##
################################################


#' K func based of flexTPC NB: for the time being all params only for fems
#'
#' @param temperature
#'
#' @returns  K value as function of temperature
#' @export
#'
#' @examples
K_func_briere = function(temperature){
  toReturn = 0
  if (temperature < 1.5) {#this func returns NaN when T (0.0001, 0.4) but also values extremely small bellow T=2 in comparison with original K_vals
    toReturn = 0.0001759071  } #this value is equivanlent to K_func_briere(1.4)
  else {
    r_max_f = 0.007638409
    T_min = 0.5
    T_max = 30
    alpha = 0.541524
    beta = 0.276430
    diff_min = temperature - T_min
    diff_max = T_max - temperature
    diff = T_max - T_min
    alpha_invert = 1 - alpha
    toReturn = r_max_f*( ((diff_min/alpha)^alpha) * ((diff_max/ alpha_invert )^alpha_invert) * (1 / diff)  )^(alpha*alpha_invert/(beta^2) )
  }

  return(toReturn)
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

  result = m*K_func_briere(temperature)*convertL_to_W(L)*L_inf*(c_div)^(1/m)*( P/(P+h) )

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
  if(L>5) s = convertL_to_W(L)*K_func_briere(temperature)*3*epsilon #molting_fraction(L*10, T) * convertL_to_W(L)*K_func(T)*3*epsilon #here L for Temming in mm
  return(s)
}

#' Natural mortality
#'
#' @param L size [cm]
#' @param temperature Temperature [°C]
#'
#' @returns [gr / time ]
#' @export
#'
#' @examples
natural_mortality_b = function(temperature,L){
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

natural_mortality_test_b = function(temperature,L){
  m=3
  mu = 0
  if (L>5) {} # nothing happens, mu = 0
  else{
      K_range = 0.007638361*m - 7.310261e-05*m #max_K - min_K #max(sapply(temp_range, K_func_briere))
      K_inv = (0.007638361*m - K_func_briere(temperature)*m) #/ K_range
      mu = convertL_to_W(L) * K_inv  # this is the rigth term of vB eq.
    }
  return(mu)
}



solver_sizeClass_extended_b = function(t, state, parameters, temperature_dataSet){
  the.method = 'rk4'
  system.equations = function(t, state, parameters) {
    with(as.list(c(state, parameters)), {

      Te = temperature_funcSolver(temperature_dataSet, t)

      #LARVAE
      gE = shift_next_sizeClass(Te, 'egg')
      IL = ingestion_rate_b(Te, LL, P)
      mL = natural_mortality_b(Te, LL)
      gL = shift_next_sizeClass(Te, 'larv')
      dL.dt = gE*E + IL*L - mL*L - gL*L

      #Juv I
      IJ = ingestion_rate_b(Te, LJ, P)
      mJ = natural_mortality_b(Te, LJ)
      gJI = shift_next_sizeClass(Te, 'juvI')
      dJ.dt = gL*L + IJ*J - mJ*J - gJI*J

      #Juv II
      IJ2 = ingestion_rate_b(Te,LJ2,P)
      mJ2 = natural_mortality_b(Te, LJ2)
      gJ2 = shift_next_sizeClass(Te, 'juvII')
      dJ2.dt = gJI*J + IJ2*J2 - mJ2*J2 - gJ2*J2 # - sA1*A1

      #Juv III
      IJ3 = ingestion_rate_b(Te,LJ3,P)
      mJ3 = natural_mortality_b(Te, LJ3)
      gJ3 = shift_next_sizeClass(Te, 'juvIII')
      dJ3.dt = gJ2*J2 + IJ3*J3 - mJ3*J3 - gJ3*J3 # - sA1*A1

      #Juv IV
      IJ4 = ingestion_rate_b(Te,LJ4,P)
      mJ4 = natural_mortality_b(Te, LJ4)
      gJ4 = shift_next_sizeClass(Te, 'juvIV')
      dJ4.dt = gJ3*J3 + IJ4*J4 - mJ4*J4 - gJ4*J4

      #Juv V
      IJ5 = ingestion_rate_b(Te,LJ5,P)
      mJ5 = natural_mortality_b(Te, LJ5)
      gJ5 = shift_next_sizeClass(Te, 'juvV')
      dJ5.dt = gJ4*J4 + IJ5*J5 - mJ5*J5 - gJ5*J5


      #Adult I
      IA1 = ingestion_rate_b(Te,LA1,P)
      sA1= spawning_rate_b(LA1, Te)
      mA1 = natural_mortality_b(Te, LA1)
      dA1.dt = gJ5*J5 + IA1*A1 - mA1*A1 - sA1*A1 - 0.1*A1 #+ (1/dev_tA1(T,L))*A1 // 0.1 Fishery


      #Adult II
      IA2 = ingestion_rate_b(Te,LA2,P)
      sA2= spawning_rate_b(LA2, Te)
      mA2 = natural_mortality_b(Te, LA2)
      dA2.dt = gJ2*J2 + IA2*A2 - mA2*A2 - sA2*A2 - 0.1*A2 #+ (1/dev_tA1(T,L))*A1 // 0.1 Fishery

      #Adult II
      molA1 = molting_fraction(LA1*10, Te)
      molA2 = molting_fraction(LA2*10, Te)
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
