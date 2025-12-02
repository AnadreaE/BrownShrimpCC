#<!--
# SPDX-FileCopyrightText: 2025 Helmholtz-Zentrum hereon GmbH
# SPDX-License-Identifier: CC0-1.0
# SPDX-FileContributor Andrea Farfan <farfanqbb@gmail.de>
#  -->

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
  #const_c = parameters_solv$general_params$const_c #reference density g/cm^m
  #m = parameters_solv$general_params$m  # ww-lenght scaling exponent
  L = L*10
  const_c = 4.625e-6
  m = 3.08
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


molting_fraction <- function(L, temperature){
  mf = 0
  L=L*10 #*10 because Teming eq. is in mm
  if (L >= parameters_solv$general_params$L_mat*10){
    if (temperature <= 0.5 ) {} #nothing happens, mf = 0
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



# new_food
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
  toReturn = scale*Q10* (7 + cos( (2*pi/ 365)*(t-215) ) )#^2 #102 spring
  return(toReturn)
}


