
#' Create specific data set for temperature within a time range
#'
#' @param dataSet from read.csv()
#' @param startDate beginning of wished time range as Character
#' @param endDate end of wished time range as Character
#'
#' @returns temperature data set of the wished time period with two columns: 'date_time' and 'temperature'
#' @export
#'
#' @examples temperature_func(dataSet ,"01/01/2000", "31/12/2000")

temperature_func <- function(dataSet, startDate, endDate){
  #Correcting format
  start <- as.Date(startDate, format = "%d/%m/%Y")
  end <- as.Date(endDate, format = "%d/%m/%Y")
  dataSet$temperature <- as.numeric(gsub(',','.', dataSet$temperature) ) #gsub replaces ',' by '.'
  dataSet$date_time <- as.POSIXct(dataSet$date_time , format = "%d/%m/%Y %H:%M")
  #createting subset
  temperature_filteredData <- dataSet[ as.Date(dataSet$date_time) >= start & as.Date(dataSet$date_time) <= end, ]
  rownames(temperature_filteredData) <- NULL
  return (temperature_filteredData)
}


#' Count how days is needed to growth fro a customarily initial length until reach the max. L
#'
#' @param df_col column of a DF with data from simulation of somatic growth for different temperatures. The columns the temperature (as Character)
#' @param min_size min. Size to start to count
#' @param max_size max. where to stop the count
#'
#' @returns development days [days]
#' @export
#'
#' @examples sapply(development_df_juvI, function(col){ count_devDays(col, juvI_min, juvI_max) } )

count_devDays <- function(df_col, min_size, max_size){
  #counts the number of days needed to reach the max size, starting by the min size
  counter <- 0
  df_col <- as.vector(df_col)
  for (i in df_col){
    if (i < min_size) counter
    if (i >= min_size && i < max_size) counter <- counter + 1
    if (1 >= max_size) counter
  }
  return (counter)

}


#' Temperature funciton to be used within a solver
#'
#' @param df dataFrame with temperature data of a limited time frame. It must content column 'temperature'
#' @param t time step from solver
#'
#' @returns Temperature [°C]
#' @export
#'
#' @examples
temperature_funcSolver <- function(df, t){
  #df <- temperature_2015_2016
  index <- floor(t)
  return(df$temperature[1+index])
}




