
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



plot_sizeSpectra_old = function(sol_df, year, title){
  #First filter only the desired year
  past_year = paste0("31/12/", as.character(year-1))
  next_year = paste0("01/01/", as.character(year+1))
  sim_year = year - as.numeric(format(sol_df$dateTime[1], "%Y")) #e.g. first year of simulation sim_year = 0
  sol_df_year = sol_df %>%
    filter(as.Date(sol_df$dateTime) > as.Date(past_year, format= "%d/%m/%Y" ) &
             as.Date(sol_df$dateTime) < as.Date(next_year, format= "%d/%m/%Y") ) %>%
    select(BF1, BF2, BF3, BF4, BF5, BF6, BF7, BF8, BM1, BM2, BM3, BM4, BM5)


  cols <- c("f" = "#8856a7", "m" = "#add8e6")

  # Define size classes to include and their order
  class_labels <- c("BL1", "BL2", "BL3", "BL4", "BL5", "BL6", "BL7", "BL8")
  label_map <- setNames(class_labels, class_labels)

  # Set up 2x2 plot grid for 4 quarters
  par(mfrow = c(2, 2))

  for (q in 1:4) {
    # Indices for quarter (912 time steps per quarter)
    init_idx <-  91.2 * (q - 1) + 1
    final_idx <- 91.2 * q
    #prepare quarterly data:
    quarter_data <- sol_df_year[init_idx:final_idx, ]

    # 1. Extract female and male size class columns
    female_cols <- paste0("BF", 1:8)
    male_cols   <- paste0("BM", 1:5)

    # 2. Compute means
    f_means <- colMeans(quarter_data[, female_cols])
    m_means <- colMeans(quarter_data[, male_cols])

    # 3. Stackable part: BF1–BF5 with BM1–BM5
    stack_classes <- paste0("BF", 1:5)
    f_stack <- f_means[stack_classes]
    m_stack <- m_means  # BM1–BM5 assumed to match BF1–BF5 by order

    # 4. Non-stackable part: BF6–BF8 (females only)
    f_extra_classes <- paste0("BF", 6:8)
    f_extra <- f_means[f_extra_classes]

    # 5. Construct matrix for barplot: one row per sex, one column per size class
    # Columns: BF1–BF5 (stacked), BF6–BF8 (female only)
    f_vals <- c(f_stack, f_extra)
    m_vals <- c(m_stack, rep(0, length(f_extra)))  # No males for BF6–BF8

    mat <- rbind("f" = f_vals, "m" = m_vals)

    # Plot barplot V4 (background because B is higher)
    barplot(
      mat,
      col = cols[rownames(mat)],
      main = paste(title , "Q", q, "–", as.character(year)),
      names.arg = label_map,
      las = 2,
      cex.main = 1.5,
      cex.axis = 1.2,
      cex.names = 1.2,
      ylim = c(0, max(colSums(mat)) * 1.1),
      space = 0,
      beside = FALSE,
      legend.text = rownames(mat),
      args.legend = list(x = "topright", bty = "n", inset = 0.02, cex = 1)
    )
  }

}


# Helper: prepare quarterly matrices
prep_sizeSpectra <- function(sol_df, year) {
  past_year <- paste0("31/12/", as.character(year - 1))
  next_year <- paste0("01/01/", as.character(year + 1))

  sol_df_year <- sol_df %>%
    filter(as.Date(dateTime) > as.Date(past_year, format = "%d/%m/%Y") &
             as.Date(dateTime) < as.Date(next_year, format = "%d/%m/%Y")) %>%
    select(BF1:BF8, BM1:BM5)

  out <- list()

  for (q in 1:4) {
    init_idx <- 91.2 * (q - 1) + 1
    final_idx <- 91.2 * q
    quarter_data <- sol_df_year[init_idx:final_idx, ]

    f_means <- colMeans(quarter_data[paste0("BF", 1:8)])
    m_means <- colMeans(quarter_data[paste0("BM", 1:5)])

    f_stack <- f_means[paste0("BF", 1:5)]
    m_stack <- m_means
    f_extra <- f_means[paste0("BF", 6:8)]

    f_vals <- c(f_stack, f_extra)
    m_vals <- c(m_stack, rep(0, length(f_extra)))

    out[[q]] <- rbind("f" = f_vals, "m" = m_vals)
  }
  out
}


# Main plotting: can plot baseline + overlay
plot_sizeSpectra <- function(sol_df_main, sol_df_overlay = NULL, year, title) {
  cols <- c("f" = "#8856a7", "m" = "#add8e6")
  label_map <- c("BL1","BL2","BL3","BL4","BL5","BL6","BL7","BL8")

  mats_main <- prep_sizeSpectra(sol_df_main, year)
  mats_overlay <- if (!is.null(sol_df_overlay)) prep_sizeSpectra(sol_df_overlay, year) else NULL

  par(mfrow = c(2, 2))

  for (q in 1:4) {
    mat <- mats_main[[q]]
    bar_centers <- barplot(
      mat,
      col = cols[rownames(mat)],
      main = paste(title, "\n Q", q, "–", year),
      names.arg = label_map,
      las = 2,
      cex.main = 1.5,
      cex.axis = 1.2,
      cex.names = 1.2,
      ylim = if (!is.null(sol_df_overlay)) c(0, max(colSums(mat), colSums(mats_overlay[[q]]) ) * 1.1)
                else c(0, max(colSums(mat))),
      space = 0,
      beside = FALSE,
      legend.text = rownames(mat),
      args.legend = list(x = "topright", bty = "n", inset = 0.02, cex = 1)
    )

    # 2. Overlay dataset on top (gray with outline)
    if (!is.null(mats_overlay)) {
      barplot(
        mats_overlay[[q]],
        col = c("f" = adjustcolor("gray50", alpha.f = 0.1),
                "m" = adjustcolor("gray50", alpha.f = 0.1)),
        border = "darkgray",
        #lty = 2, #dashed lines
        space = 0,
        beside = FALSE,
        add = TRUE,
        axes = FALSE,
        names.arg = rep("", length(label_map))  # suppress duplicate labels
      )
    }
  }
}



