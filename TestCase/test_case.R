#<!--
# SPDX-FileCopyrightText: 2025 Helmholtz-Zentrum hereon GmbH
# SPDX-License-Identifier: Apache-2.0
# SPDX-FileContributor Andrea Farfan <farfanqbb@gmail.de>
#  -->

library(deSolve)
library(RColorBrewer)
library(BrownShrimp)
#library(viridis)
library(dplyr)
library(tidyr)


################################################
## As follows a few test lines to ensure that ##
## the project was installed propperly.       ##
################################################

# (1) Test access to built-in variables:
print(parameters_solv)
#If you see the list of parameters, then test is successfully passed,
#If R doesn't find it, try to go to: Build tab/ More / Load all.

#(2) If you have passsed test (1),  now you shall have access to the functions and data too.
#As follows, we test one simulation:

#Read water temperature from Germany:
temperature_dataSet <- read.csv("./data/temperature_10ger.csv")
temperature_10_15 <- temperature_dataSet %>% filter(date_time > as.Date("31/12/2009", format = "%d/%m/%Y") & date_time < as.Date("01/01/2016", format = "%d/%m/%Y"))

#Initialize parameters
params = parameters_solv
params$general_params$Fi = 4/365 #(1-exp(-5))/365
params$general_params$L50 = 3.69 #et. al Santos 2018
params$general_params$SR = 0.75 #et. al Santos 2018
params$general_params$Imax_ik = 0.16 #(1-exp(-2))/365 #reduce predation in these years as per: https://www.openseas.org.uk/news/on-thin-ices-north-sea-cod/

#Initialize Initial conditions of the state variables:
BF = c(8.7, 15.35, 21.46, 24.67,25.87, 17.54, 14.09, 11.93)
BM = c(10.86, 19.94, 29.01, 36.67, 58.27)

Init = c(P = 2, E = 2, L= 1 ,
            BF = BF, BM = BM)

#Define time spam and resolution of time steps
t_years = seq(0,length(tempG.10_15$temperature)-.1)


#Simulate
start <- Sys.time()
test_case <- solver_sizeClass.v5(t = t_years, state = Init, parameters = params, temperature_dataSet = temperature_10_15)
print(paste('Test passed successfully, time to run simulation:', Sys.time() - start) )

#See results: quarterly average from the year you wish:
plot_sizeSpectra_old(sol_df_main = test_case, year = 2014, title = "F&P")



