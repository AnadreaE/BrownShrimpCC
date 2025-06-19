#HERE SIMULATIONS WITH NEW K-FUNC tpc FLEXBRIERE

library(ggplot2)
library(deSolve)
library(RColorBrewer)
library(BrownShrimp)
library(viridis)
library(dplyr)
library(tidyr)
library("colorspace")

#library(devtools)

#UPLOAD TEMPERATURE DATA

temperature_dataSet <- read.csv("./data/Temperature_data_1982-2018_MarBiol.csv", header = TRUE, sep = ';')

temperature_10_11 <- temperature_func(temperature_dataSet, "01/01/2010", "31/12/2011")
temperature_10_14 <- temperature_func(temperature_dataSet, "01/01/2010", "31/12/2014")


######### START ###########
#Def. Size Classes ranges [cm]

L_as_F <- 8.5
L_as_M = 5.5

juvI_min <- .6
juvI_max <- 1.0

juvII_min  <- 1.0
juvII_max <- 2.0

juvIII_min  <- 2.0
juvIII_max <- 3.0

juvIV_min  <- 3.0
juvIV_max <- 4.0

juvV_min  <- 4.0
juvV_max <- 5.0

adultI_min <- 5.0
adultI_max <- 6.0

adultII_min <- 6.0
adultII_max <- 7.0

adultIII_min <- 7.0
adultIII_max <- 8.5#Linf_F


LL = 0.3
LJ = (juvI_min+juvI_max)/2 #cm
LJ2 = (juvII_min+juvII_max)/2   # 3cm
LJ3 = (juvIII_min+juvIII_max)/2
LJ4 = (juvIV_min+juvIV_max)/2
LJ5 = (juvV_min+juvV_max)/2
LA1 = (adultI_min + adultI_max)/2
LA2 = (adultII_min + adultII_max)/2  #cm
LA3 = (adultIII_min + adultIII_max)/2  #cm


state <- c(P = 2, E = 0.1, L= 0.4 ,
           J_f = .175*0.5 , J2_f = 0.125*0.5, J3_f = 0.125*0.5, J4_f = 0.1*0.5, J5_f = 0.1*0.5,
           J_m = .175*0.5 , J2_m = 0.125*0.5, J3_m = 0.125*0.5, J4_m = 0.1*0.5, J5_m = 0.1*0.5,
           A1_f = 0.25*0.5, A1_m = 0.25*0.5,
           A2 = 0.1, A3=0.1) #males doesn't reach this size classes

t = seq(0,length(temperature_10_11$temperature)-0.1, by = 0.1)
t_5years = seq(0,length(temperature_10_14$temperature)-0.1, by = 0.1)

parameters = c()

start <- Sys.time()
test_sex <- solver_sizeClass_sex(t = t_5years, state = state, parameters = parameters, temperature_10_14)
print(Sys.time() - start)

start_date <- as.POSIXct("2010-01-01", format="%Y-%m-%d", tz = "UTC")

test_sex <- mutate(test_sex, dateTime = start_date + test_sex$time * 86400) #86400 seconds in one day
test_sex <- test_sex[ , -1] # delete timesteps column

#cc_long_2 <- gather(cc_2, key = "Type", value = "Value", P, E, L, J, J2, J3, J4, J5, A1, A2)
test_sex_long <- test_sex %>%
  pivot_longer(-dateTime, names_to = "variable", values_to = "value")

#dev.off()
ggplot(test_sex_long, aes(x = dateTime, y = value, color = variable)) +
  geom_line(size = 1) +  # Plot lines
  labs(title = " TPC-sex F and M \n 2010-2014", x = "Days", y = "Biomass") +
  scale_color_manual(

    values = c( "P" = "#d9f0a3", "E"= "#1c9099" ,"L" = "gray"  ,
                "J_f" = "#8856a7", "J2_f" = "#8856a7", "J3_f" = "#8856a7", "J4_f" = "#8856a7", "J5_f" = "#8856a7",
                "J_m" = "#9ebcda", "J2_m" = "#9ebcda", "J3_m" = "#9ebcda", "J4_m" = "#9ebcda", "J5_m" = "#9ebcda",
                "A1_f" = "#fd8d3c" , "A1_m" = "maroon4" ,
                "A2" = "#e6550d" , "A3" = 'red2')) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), # Center title
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        legend.position="bottom") #+ ylim(0, 100) #+ xlim( as.POSIXct("2010-04-01", format="%Y-%m-%d", tz = "UTC") , as.POSIXct("2012-07-10", format="%Y-%m-%d", tz = "UTC") )



##### Stacked area chart with K as TPC ####
#preparation:

df_long_testSex <- test_sex[ , -1] %>% #delete  plancton col
  pivot_longer(-dateTime, names_to = "variable", values_to = "value")

#ensure that vertical order is same as life stages
stack_order <- rev(c(colnames(test_sex)[2:17]) )

# Ensure 'variable' is a factor with levels in the same order
df_long_testSex$variable <- factor(df_long_testSex$variable, levels = stack_order)


ggplot(df_long_testSex, aes(x = dateTime, y = value, fill = variable)) +
  geom_area() +
  scale_fill_viridis_d() +
  labs(
    title = "TPC sex",
    x = "Date",
    y = "Biomass",
    fill = "Size class"
  )

###### Barplots #######
#See F and M separately in order to see differences and see year 2013 (third year of simulation)

#Subset desired year (2013) and by sex:
#FEM
F_2013 = test_sex %>%
    filter( as.Date(dateTime) > as.Date("31/12/2012", format= "%d/%m/%Y" ) & as.Date(dateTime) < as.Date("01/01/2014", format= "%d/%m/%Y" )  ) %>%
    select(E, L, J_f, J2_f, J3_f, J4_f, J5_f, A1_f, A2, A3)

months <- c('Jan', 'Feb', 'Mar', 'Apr', 'Mai', 'Jun', 'Jul', 'Aug', 'Sep', 'Okt', 'Nov', 'Dez')
#legend_labels <- c("Egg", "Larv", "Juv I", "Juv II", "Adu")
width_vector_compact <- c((2/10)*0.2, (5/10)*0.2, ((juvV_max - juvI_min)/10)*0.2, ((adultII_max - adultI_min)/10)*0.2 )


#all / quarterly mean values
cols <- sequential_hcl(10, palette = "viridis")

width_vector_all <- c((.2/10)*0.1, (.6/10)*0.1, ((juvI_max - juvI_min)/10)*0.1,  ((juvII_max - juvII_min)/10)*0.1, ((juvIII_max - juvIII_min)/10)*0.1,
                      ((juvIV_max - juvIV_min)/10)*0.1, ((juvV_max - juvV_min)/10)*0.1,
                      ((adultI_max - adultI_min)/10)*0.1, ((adultII_max - adultII_min)/10)*0.1, ((adultIII_max - adultIII_min)/10)*0.1)

x_labels_vec <- c( "0", "<6", "6-10", "10-20", "20-30", "30-40", "40-50", "50-60", "60-70", "70-L_inf")
dev.off()
par(mfrow = c(2,2))
j <- 1

for (i in seq(1, 4)) {#four quarters
  init_day <- 912*(j-1) #365/4 = 91,25 days per quarter, /0.1 timesteps = 912
  final_day <- 912*j
  month_data <- F_2013[init_day:final_day, ]
  toPlot <- colMeans(month_data)
  barplot(toPlot, col = cols, main = paste('F -Q ', j ,' 2013'), las = 2, cex.main = 2, names.arg = x_labels_vec  ,
          width = width_vector_all, cex.axis = 2, cex.names = 1.8 , ylim = c(0, 25), space = 0)
  j <- j + 1
}


#MALE
M_2013 = test_sex %>%
  filter( as.Date(dateTime) > as.Date("31/12/2012", format= "%d/%m/%Y" ) & as.Date(dateTime) < as.Date("01/01/2014", format= "%d/%m/%Y" )  ) %>%
  select(E, L, J_m, J2_m, J3_m, J4_m, J5_m, A1_m)

months <- c('Jan', 'Feb', 'Mar', 'Apr', 'Mai', 'Jun', 'Jul', 'Aug', 'Sep', 'Okt', 'Nov', 'Dez')
#legend_labels <- c("Egg", "Larv", "Juv I", "Juv II", "Adu")

#all / quarterly mean values
cols <- sequential_hcl(8, palette = "viridis")

width_vector_all <- c((.2/10)*0.1, (.6/10)*0.1, ((juvI_max - juvI_min)/10)*0.1,  ((juvII_max - juvII_min)/10)*0.1, ((juvIII_max - juvIII_min)/10)*0.1,
                      ((juvIV_max - juvIV_min)/10)*0.1, ((juvV_max - juvV_min)/10)*0.1,
                      ((adultI_max - adultI_min)/10)*0.1, ((adultII_max - adultII_min)/10)*0.1, ((adultIII_max - adultIII_min)/10)*0.1)

x_labels_vec <- c( "0", "<6", "6-10", "10-20", "20-30", "30-40", "40-50", "50-L_inf")
#dev.off()
par(mfrow = c(2,2))
j <- 1

for (i in seq(1, 4)) {#four quarters
  init_day <- 912*(j-1) #365/4 = 91,25 days per quarter, /0.1 timesteps = 912
  final_day <- 912*j
  month_data <- M_2013[init_day:final_day, ]
  toPlot <- colMeans(month_data)
  barplot(toPlot, col = cols, main = paste('Q ', j ,' 2013'), las = 1, cex.main = 2, names.arg = x_labels_vec,
          width = width_vector_all, cex.axis = 2, cex.names = 1.8 , ylim = c(0, 25), space = 0, las = 2)
  j <- j + 1
}
mtext("size range", side = 1, outer = TRUE, line = 2, cex = 1.5)
mtext("density", side = 2, outer = TRUE, line = 2, cex = 1.5)


#FEM for each 2 months

par(mfrow = c(2,3))
j <- 1

for (i in seq(1, 6)) {#four quarters
  init_day <- 608*(j-1) #365/4 = 91,25 days per quarter, /0.1 timesteps = 912
  final_day <- 608*j
  month_data <- F_2013[init_day:final_day, ]
  toPlot <- colMeans(month_data)
  barplot(toPlot, col = cols, main = paste('F bi-month ', j ,' 2013'), las = 2, cex.main = 2, names.arg = x_labels_vec  ,
          width = width_vector_all, cex.axis = 2, cex.names = 1.8 , ylim = c(0, 25), space = 0)
  j <- j + 1
}




