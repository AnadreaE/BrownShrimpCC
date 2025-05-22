################################################
## Following functions are meant to be        ##
## used in simulation for either individuals  ##
## or cohorts.                                ##
## Source: Temming                            ##
################################################


#' Development status Egg
#'
#' @param dev_status is the current dev. status of an individual or cohort. Valid values (0,1]. When status => 1, then egg hatch and change status to 'Larvae'
#' @param temperature temperature of the day [°C]
#'
#' @returns a list with 'new dev. status' a boolean wheather ready to hatch or nor in order to track history of development
#' @export
#'
#' @examples egg_dev(0.5, 15)
#' @examples for (i in 2:length(temperature2015_Q2$temperature) ){
#' egg_development <- egg_dev(current_dev_status, temperature2015_Q2$temperature[i])
#' current_dev_status <- egg_development[[1]]
#' hatch_status_array_Q2 <- append(hatch_status_array_Q2, egg_development[[2]] )
#' egg_development_array_Q2 <- append(egg_development_array_Q2,current_dev_status)}

egg_dev <- function(dev_status, temperature){
  ready_to_hatch <- FALSE
  alpha <- 1031.34
  if (dev_status < 1) {
    development_day <- 1 / (alpha*temperature^-1.354)
    dev_accumulated <- dev_status + development_day
    if (dev_accumulated >= 1) {
      ready_to_hatch <- TRUE
    }
    return (list(dev_accumulated, ready_to_hatch) ) }
  else if (dev_status >= 1)
  { ready_to_hatch <- TRUE
  development_day <- 0
  dev_accumulated <- dev_status + development_day
  return (list(dev_accumulated, ready_to_hatch))}
}


#' Development status Larvae
#'
#' @param dev_status is the current dev. status of an individual or cohort. Valid values (0,1]. When status => 1, then Larvae changes status to 'Juvenile'
#' @param temperature temperature of the day [°C]
#'
#' @returns a list with 'new dev. status' a boolean wheather ready to hatch or nor in order to track history of development
#' @export
#'
#' @examples larvae_dev(0.5, 15)

larvae_dev <- function(dev_status, temperature) {
  ready_for_juv <- FALSE
  beta <- 5 / 0.00584
  if (dev_status < 1) {
    development_day <- 1 / (beta * temperature^-1.347)
    dev_accumulated <- dev_status + development_day
    if (dev_accumulated >= 1) {
      ready_for_juv <- TRUE
      dev_accumulated <- 1  # Ensure dev_accumulated does not exceed 1
    }
  } else {
    ready_for_juv <- TRUE
    dev_accumulated <- dev_status  # No change if already at or above 1
  }
  # Return a list after both conditions
  return(list(dev_accumulated, ready_for_juv))
}


#Def function for somatic growth for Juvenile_I
#' Somatic growth for juvenile and adult
#'
#' @param b_length current length [mm]
#' @param temperature current day temperature [°C]
#' @param max_size asymptotic size [mm]
#' @param sex distinct param. values for male and females 'F' or 'M'.
#'
#' @returns new length
#' @export
#'
#' @examples sg_juv(15, 20, 85, 'F')

som_growth <- function( b_length, temperature, max_size, sex){
  #Def parameters for somatic growth
  if (sex == 'F') {
    b <- 0.04028
    c <- 0.00193
    d <- 0.0877
  }
  else {
    b <- 0.03424
    c <- 0.002
    d <- 0.0877
  }

  growth_calc <- b*temperature - exp(d*temperature)*b_length*c
  if (b_length < max_size & growth_calc >= 0) {
    growth <- growth_calc
  } else {growth <- 0}

  return(b_length + growth)
}
