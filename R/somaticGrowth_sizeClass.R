################################################
## Following functions are meant to be        ##
## used for solving differential eq. systems  ##
## Source: Andrea F.                          ##
################################################
#source("./R/somaticGrowth_sizeClasses_II.R")

##### constant parameter values: #####
#All now in somaticGrowth_sizeClasses_II

##### FUNCTIONS #####

#' Convert size to weight W
#'
#' @param L size [cm]
#'
#' @returns equivalent weight [gr]
#' @export
#'
#' @examples convertL_to_W(L = 2.5)

convertL_to_W = function(L){
  const_c = parameters_solv$general_params$const_c #reference density g/cm^m
  m = parameters_solv$general_params$m  # ww-lenght scaling exponent
  return (const_c*L^m)
}


#' *Deprecated Function to calculate K (growth parameter from vB eqs.)
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



#' Fraction of females that are molting and therefore releasing eggs.
#'
#' @param L size [cm]
#' @param temperature Temperature [°C]
#'
#' @returns a fraction [dimensionless]
#' @export
#'
#' @examples


molting_fraction_old = function(L, Te, Fem_params){
  w = convertL_to_W(L)
  vB = 3*K_func_briere(Te, parameters_solv$Fem_params)*w
  interc = 0.01493457*exp(L*-0.10730890)
  fract = 13.96105*exp(L*-0.5684764)
  return(interc + vB*fract )
}




molting_fraction <- function(L, temperature){
  mf = 0
  L=L*10
  if (L >= parameters_solv$general_params$L_mat*10){
    if (temperature <= 0 ) {} #nothing happens, mf = 0
    else if (temperature <= 16.2 && temperature > 0){ #T_opt fems= 16.2
      #mf =  (1 / (5.7066 * L^0.7364 * temperature ^(-0.09363) ) )/2  # exp(temperature*-0.09363) ) #/2 considering
      mf = (1 / (5.7066 * L^0.7364 * exp(temperature*-0.09363) ) ) #revision 28.05. Temming's implementation (above was Temming's formula on paper)
    }
    else if (temperature > 16.2){
      mf = (1 / (5.7066 * L^0.7364 * exp(16.2*-0.09363) ) ) #this aviod that fraction increases with higher temperatures than T_opt and even with unrealistic T like 40°
    }
  }

  return ( mf )
}



#' *Deprecated Natural mortality
#'
#' @param L size [cm]
#' @param temperature Temperature [°C]
#'
#' @returns [gr / time ]
#' @export
#'
#' @examples
respiration_rate = function(temperature,L){
  mu = m*convertL_to_W(L)*K_func(temperature)*0.1 # this is the rigth term of vB eq.#0.1
  #if (L>5) {
  #  mu =  convertL_to_W(L) *K_func(temperature) *( 1 - molting_fraction(L, temperature)) #tbc (1-molting_fraction) all non molting fems still have a natural mortality
  #}
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
new_food = function(t, Te, scale = 2) {
  Q10 = 1.6**((Te-10)/10)
  toReturn = scale*Q10  * (2.2 + cos( (2*pi/ 365)*(t-105) ) )^2 #reduced to test 24.07. from 0.2 to 1.5
 # print(paste("new_food sucsessful",  0.2 * (1.2 + cos( (2*pi/ 365)*(t-172) ) ) ))
  #( 1 + cos( (2*pi/(365/0.1) ) *(t-(100/0.1) ) )) # 2 + cos(2 * pi * t / 365
  return(toReturn)
}









