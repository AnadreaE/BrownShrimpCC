################################################
## Following functions are meant to be        ##
## used for solving differential eq. systems  ##
## Source: Andrea F.                          ##
################################################


#' Convert size to weight W
#'
#' @param L size [cm]
#'
#' @returns equivalent weight [gr]
#' @export
#'
#' @examples convertL_to_W(L = 2.5)

convertL_to_W = function(L){
  #const_c <- 0.01 # reference density g/cm^m
  #m <- 3 # ww-lenght scaling exponent
  return (const_c*L^m)
}

#' Function to calculate K (growth parameter from vB eqs.)
#'
#' @returns K value as function of temperature
#' @export
#'
#' @examples
K_func = function(temperature){
  #following param. vals have been estimated and docum ented in Thesis Appx. A
  a = 0.007614#0.008
  mu =  16.24 # average
  sigma = 7.098
  #T_opt = 16.5 # optimal temperature
  return (a * exp(-((temperature - mu)^2) / (2 * sigma^2)))
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

ingestion_rate = function(temperature, L, P){
  c_div = const_c / convertL_to_W(L)

  result = m*K_func(temperature)*convertL_to_W(L)*L_inf*(c_div)^(1/m)*( P/(P+h) )

  #m*convertL_to_W(L)*K_func(T)*L_inf/L*P/(P+h) #revised formula Andrea 04.04.25
  #1.14*K_func(T)*(L_inf/ L)^n *P/(P+h)
  return (result)
}


#' Fraction of females that are molting and therefore releasing eggs.
#'
#' @param L size [cm]
#' @param temperature Temperature [°C]
#'
#' @returns a fraction [dimensionless]
#' @export
#'
#' @examples

molting_fraction <- function(L, temperature){
  mf = 0
  if (temperature <= 0 ) {}
  else {
    mf =  (1 / (5.7066 * L^0.7364 * temperature ^(-0.09363) ) )/2  # exp(temperature*-0.09363) ) #/2 considering
  }
  return ( mf )
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
spawning_rate = function(L, temperature){
  s = 0
  if(L>5) s = convertL_to_W(L)*K_func(temperature)*3*epsilon #molting_fraction(L*10, T) * convertL_to_W(L)*K_func(T)*3*epsilon #here L for Temming in mm
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
natural_mortality = function(temperature,L){
  mu = m*convertL_to_W(L)*K_func(temperature)# this is the rigth term of vB eq.
  if (L>5) mu = 0 #tbc
  return(mu)
}



#new_food = function(t) .2*(2+cos(t))
#' Production of resources (Plancton) - TEMPORARILY ARTIFICIAL FUNCTION WITH COSINUS; TO BE REPLACED BY A MORE REALISTIC ONE
#'
#' @param t time
#'
#' @returns
#' @export
#'
#' @examples
new_food = function(t) {
  0.2 * (2 + cos( (2*pi/ 365)*(t-172) ) )
  #( 1 + cos( (2*pi/(365/0.1) ) *(t-(100/0.1) ) )) # 2 + cos(2 * pi * t / 365
}



#' Solver for a size class model with Eggs, Larvae, 5 Juvenile classes and 2 Adult classes
#'
#' @param t
#' @param state
#' @param parameters
#' @param temperature_dataSet
#'
#' @returns a df with solution to the equation system
#' @export
#'
#' @examples
#'
solver_sizeClass_extended = function(t, state, parameters, temperature_dataSet){
  the.method = 'rk4'
  system.equations = function(t, state, parameters) {
    with(as.list(c(state, parameters)), {

      Te = temperature_funcSolver(temperature_dataSet, t)

      #LARVAE
      gE = shift_next_sizeClass(Te, 'egg')
      IL = ingestion_rate(Te, LL, P)
      mL = natural_mortality(Te, LL)
      gL = shift_next_sizeClass(Te, 'larv')
      dL.dt = gE*E + IL*L - mL*L - gL*L

      #Juv I
      IJ = ingestion_rate(Te, LJ, P)
      mJ = natural_mortality(Te, LJ)
      gJI = shift_next_sizeClass(Te, 'juvI')
      dJ.dt = gL*L + IJ*J - mJ*J - gJI*J

      #Juv II
      IJ2 = ingestion_rate(Te,LJ2,P)
      mJ2 = natural_mortality(Te, LJ2)
      gJ2 = shift_next_sizeClass(Te, 'juvII')
      dJ2.dt = gJI*J + IJ2*J2 - mJ2*J2 - gJ2*J2 # - sA1*A1

      #Juv III
      IJ3 = ingestion_rate(Te,LJ3,P)
      mJ3 = natural_mortality(Te, LJ3)
      gJ3 = shift_next_sizeClass(Te, 'juvIII')
      dJ3.dt = gJ2*J2 + IJ3*J3 - mJ3*J3 - gJ3*J3 # - sA1*A1

      #Juv IV
      IJ4 = ingestion_rate(Te,LJ4,P)
      mJ4 = natural_mortality(Te, LJ4)
      gJ4 = shift_next_sizeClass(Te, 'juvIV')
      dJ4.dt = gJ3*J3 + IJ4*J4 - mJ4*J4 - gJ4*J4

      #Juv V
      IJ5 = ingestion_rate(Te,LJ5,P)
      mJ5 = natural_mortality(Te, LJ5)
      gJ5 = shift_next_sizeClass(Te, 'juvV')
      dJ5.dt = gJ4*J4 + IJ5*J5 - mJ5*J5 - gJ5*J5


      #Adult I
      IA1 = ingestion_rate(Te,LA1,P)
      sA1= spawning_rate(LA1, Te)
      mA1 = natural_mortality(Te, LA1)
      dA1.dt = gJ5*J5 + IA1*A1 - mA1*A1 - sA1*A1 - 0.1*A1 #+ (1/dev_tA1(T,L))*A1 // 0.1 Fishery


      #Adult II
      IA2 = ingestion_rate(Te,LA2,P)
      sA2= spawning_rate(LA2, Te)
      mA2 = natural_mortality(Te, LA2)
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










