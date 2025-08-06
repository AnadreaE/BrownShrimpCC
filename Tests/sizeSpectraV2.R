
library(ggplot2)
library(deSolve)
library(RColorBrewer)
library(BrownShrimp)
library(viridis)
library(dplyr)
library(tidyr)
library("colorspace")
library(lubridate)
library(latex2exp)
library(profvis)

#library(devtools)

#UPLOAD TEMPERATURE DATA

temperature_dataSet <- read.csv("./data/Temperature_data_1982-2018_MarBiol.csv", header = TRUE, sep = ';')
temperature_14_18 <- temperature_func(temperature_dataSet, "01/01/2014", "31/12/2018")
t_5years = seq(0,length(temperature_14_18$temperature)-0.1)

BF_p = c(1.1476,  2.739,  5.604, 11.65, 21.524, 17.976, 11.613,  5.346) #see calculation in file 'sizeSpectraMoldel.R'
BM_p = c(1.464, 3.645, 8.855, 22.622, 37.614) #see calculation in file 'sizeSpectraMoldel.R'


IC_parameterized = c(P = 2, E = 0.759, L= 0.38  ,
             BF = BF_p, BM = BM_p) #Initial conditions according to estimations with fishery landings data.


#### TEST V5 (1 out of 4) WITHOUT ANY KIND OF PREASURE (P NEITHER F) ####


noPreasure_params = parameters_solv
noPreasure_params$general_params$Fi = 0 #No fishery
noPreasure_params$general_params$Imax_ik = 0 #No predation

start <- Sys.time()
test_noPreasure <- solver_sizeClass.v5(t = t_5years, state = IC_parameterized, parameters = noPreasure_params, temperature_dataSet = temperature_14_18)
print(Sys.time() - start)

#Plot size spectra
plot_sizeSpectra(test_noPreasure, 2015, title = "NP")



#### TEST V5 (2 out of 4) ONLY PREDATION ####

pred_params = parameters_solv
pred_params$general_params$Fi = 0 #No fishery


start <- Sys.time()
test_predation <- solver_sizeClass.v5(t = t_5years, state = IC_parameterized, parameters = pred_params, temperature_dataSet = temperature_14_18)
print(Sys.time() - start)

#Plot size spectra
plot_sizeSpectra(test_predation, 2015, title = "Pred")

#### TEST V5 (3 out of 4) ONLY FISHERY ####

fishery_params = parameters_solv
fishery_params$general_params$Imax_ik = 0 #No predation
fishery_params$general_params$Fi = 0.025 #reduced fishery

start <- Sys.time()
test_fishery <- solver_sizeClass.v5(t = t_5years, state = IC_parameterized, parameters = fishery_params, temperature_dataSet = temperature_14_18)
print(Sys.time() - start)

#Plot size spectra
plot_sizeSpectra(test_fishery, 2015, title = "Fi")

#### TEST V5 (4 out of 4) PREDATION AND FISHERY####

bothPF_params = parameters_solv
bothPF_params$general_params$Fi = 0.025 #reduced fishery

start <- Sys.time()
test_bothFP <- solver_sizeClass.v5(t = t_5years, state = IC_parameterized, parameters = bothPF_params, temperature_dataSet = temperature_14_18)
print(Sys.time() - start)

#Plot size spectra
plot_sizeSpectra(test_bothFP, 2015, title = "F&P")

#Plot fishery catch (for only one year):
#define one year to be inspected
year_toSee = 2015
#subset data and sum up to monthly values:
sol_yearToSee = test_bothFP %>%
               filter(as.numeric(format(test_bothFP$dateTime,'%Y')) == year_toSee) %>%
              mutate(month = month(dateTime)) %>%
              group_by(month) %>%
              summarise(across(where(is.numeric), \(x) sum(x, na.rm = TRUE)))


par(mfrow = c(1,1))
plot(sol_yearToSee$month, sol_yearToSee$catch.BF1,  type = 'b', main = 'Fished Kg per Km2')
#lines(1:length(monthly_Feffort), monthly_Feffort, col = 'red')

#stacked area plot

#test_pred.V5_red <- test_predation[ , c(3:18)] # delete timesteps column, plancton  and predator column
test_both.V5_red <- test_bothFP[ , c(3:19)] # delete timesteps column, plancton  and predator column

#cc_long_2 <- gather(cc_2, key = "Type", value = "Value", P, E, L, J, J2, J3, J4, J5, A1, A2)
test_both_long <- test_both.V5_red %>%
  pivot_longer(-dateTime, names_to = "variable", values_to = "value")

ggplot(test_both_long, aes(x = dateTime, y = value, fill = variable)) +
  geom_area() +
  scale_fill_viridis_d() +
  labs(
    title = "System incl. predation ",
    x = "Date",
    y = "Biomass",
    fill = "Size class"
  )

