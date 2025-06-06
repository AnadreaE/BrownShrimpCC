############################################
## The growth function developed in this  ##
## thesis will be validated and tested    ##
## in this script                         ##
############################################

#for this simulation, we will consider: only juvenile stage;simulation over two years (monthly values for age);
#constant temperature and following parameters:
temperature_range <- seq(0, 30, 0.1)
days_2years <- seq(1, 2*365, 1)


#FEMALES
L_as_f <- 8.5 #Asymptotic length female ~ max. length [cm]


#STEP(1) TEMMINGTS SIMULATION WITH CONSTANT TEMPERATURE

#FEMALES
# create an empty vector for each temperature degree.
temming_f <- data.frame(row = days_2years)

ind <- 1
for (j in temperature_range){
  development <- c()
  initial_l <- 0.6 #cm
  for (i in days_2years) {
    new_legth <- som_growth(initial_l, j, L_as_f, 'F')
    development <- append(development, new_legth) #append the new length
    initial_l <- new_legth
  }
  temming_f[as.character(j)] <- development
}

#STEP (2) SIMULATION WITH VB l(t) with K(T) with briere:

growth_thesis_f <- data.frame(row = days_2years)

for (i in temperature_range){
  development <- c()
  initial_l <- 0.6 #cm
  for (j in days_2years){
    growth = som_growth_thesis(initial_l, j ,i , "F")
    development = append(development, growth)
    }
  growth_thesis_f[as.character(i)] <- development
  }


#### plot Temming vs Thesis ####

temp_to_print <- seq(6,14, 1)

par(mfrow = c(3, 3), oma = c(5, 5, 2, 1))  # outer margins: bottom, left, top, right
par(mar = c(2, 2, 2, 1))  # inner margins for subplots (smaller so things fit)

for (i in temp_to_print ){
  RSS <- sum((temming_f[[as.character(i)]] - growth_thesis_f[[as.character(i)]] )^2)
  # Calculate total sum of squares (TSS)
  TSS <- sum((temming_f[[as.character(i)]] - mean(temming_f[[as.character(i)]]))^2)
  R_squared <- 1 - (RSS / TSS)
  plot( days_2years, temming_f[[as.character(i)]], col = 'gray41', main = paste('T = ', i, ', RSQ = ', round(R_squared, 2)),
        xlab = '', ylab = '', las = 1, lwd = 2.5, cex.axis = 1.4, cex.main = 1.3) #
  lines(days_2years, growth_thesis_f[[as.character(i)]], col =  'lightseagreen', lwd = 2)
}

mtext("days", side = 1, outer = TRUE, line = 3, cex = 1.2)
mtext("L in mm", side = 2, outer = TRUE, line = 3, cex = 1.2)

par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", axes = FALSE, xlab = "", ylab = "")
legend("bottomright", legend = c("empirical", "Thesis growth func."),
       col = c("gray41", "lightseagreen"), lwd = 2, bty = "n", horiz = TRUE, cex = 1.4)

