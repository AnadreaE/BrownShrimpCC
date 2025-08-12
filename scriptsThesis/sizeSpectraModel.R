#SIZE SPECTRA MODEL

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


state_r = c(P = 2, E = 0.1, L= 0.0928 *1.5 ,
            BF = BF+0.1, BM = BM+0.1,
            Pred = 0.05 )


#### TEST (1 out of 4) WITHOUT ANY KIND OF PREASURE (P NEITHER F) ####

noPreasure_params = parameters_solv
noPreasure_params$general_params$Fi = 0 #No fishery
noPreasure_params$general_params$Imax_ik = 0 #No predation


start <- Sys.time()
test_noPreasure <- solver_sizeClass.v4(t = t_5years, state = state_r, parameters = noPreasure_params, temperature_dataSet = temperature_14_18)
print(Sys.time() - start)

#Plot size spectra
plot_sizeSpectra(test_noPreasure, 2015, title = "NP")


#### TEST (2 out of 4) ONLY PREDATION ####

pred_params = parameters_solv
pred_params$general_params$Fi = 0 #No fishery


start <- Sys.time()
test_predation <- solver_sizeClass.v4(t = t_5years, state = state_r, parameters = pred_params, temperature_dataSet = temperature_14_18)
print(Sys.time() - start)

#Plot size spectra
plot_sizeSpectra(test_predation, 2015, title = "Pred")

#### TEST (3 out of 4) ONLY FISHERY ####

fishery_params = parameters_solv
fishery_params$general_params$Imax_ik = 0 #No predation
fishery_params$general_params$Fi = 0.025 #reduced fishery

start <- Sys.time()
test_fishery <- solver_sizeClass.v4(t = t_5years, state = state_r, parameters = fishery_params, temperature_dataSet = temperature_14_18)
print(Sys.time() - start)

#Plot size spectra
plot_sizeSpectra(test_fishery, 2015, title = "Fi")

#### TEST (4 out of 4) PREDATION AND FISHERY####

bothPF_params = parameters_solv
bothPF_params$general_params$Fi = 0.025 #reduced fishery

start <- Sys.time()
test_bothFP <- solver_sizeClass.v4(t = t_5years, state = state_r, parameters = bothPF_params, temperature_dataSet = temperature_14_18)
print(Sys.time() - start)

#Plot size spectra
plot_sizeSpectra(test_bothFP, 2015, title = "F&P")

#stacked area plot

test_bothFP_red <- test_bothFP[ , c(3:17,19)] # delete timesteps column, plancton  and predator column

#cc_long_2 <- gather(cc_2, key = "Type", value = "Value", P, E, L, J, J2, J3, J4, J5, A1, A2)
test_pred_long <- test_bothFP_red %>%
  pivot_longer(-dateTime, names_to = "variable", values_to = "value")

ggplot(test_pred_long, aes(x = dateTime, y = value, fill = variable)) +
  geom_area() +
  scale_fill_viridis_d() +
  labs(
    title = "System incl. predation ",
    x = "Date",
    y = "Biomass",
    fill = "Size class"
  )




#### PLOTS TO COMPARE ALL ####

#(1) L_avg over the time:
size_classes = female_cols <- paste0("BL", 1:8)

#calculate weighted average for each time step:

#first join F and M in one column
#Without preasure (1)
noPreasure_noSex = test_noPreasure %>%
  mutate(L1 = BF1 + BM1, L2 = BF2 + BM2, L3 = BF3 + BM3, L4 = BF4 + BM4,  L5 = BF5 + BM5) %>%
  mutate(ttlB = L1 + L2+ L3+ L4+ L5+ BF6+ BF7 + BF8) %>%
  select(20:24, 10:12,25, 19 ) %>% #Column 24 HAS TO EB CHANGED WHEN RUNNING NEW SIM TO SELECT DATETIME
  filter(as.Date(dateTime) > as.Date("31/12/2014", format= "%d/%m/%Y" )) #delete first year 'warm up' period


L_avg_noPreasure = noPreasure_noSex %>%
  mutate(L_avg = ( L1*size_mean_F[1] + L2*size_mean_F[2]+ L3*size_mean_F[3]+ L4*size_mean_F[4] + L5*size_mean_F[5]
                  + BF6*size_mean_F[6] + BF7*size_mean_F[7]+ BF8*size_mean_F[8])/ttlB )

#Only Predation (2)
predation_noSex = test_predation %>%
  mutate(L1 = BF1 + BM1, L2 = BF2 + BM2, L3 = BF3 + BM3, L4 = BF4 + BM4,  L5 = BF5 + BM5) %>%
  mutate(ttlB = L1 + L2+ L3+ L4+ L5+ BF6+ BF7 + BF8) %>%
  select(20:24, 10:12,25, 19 ) %>% #Column 24 HAS TO EB CHANGED WHEN RUNNING NEW SIM TO SELECT DATETIME
  filter(as.Date(dateTime) > as.Date("31/12/2014", format= "%d/%m/%Y" )) #delete first year 'warm up' period

Lavg_pred = predation_noSex %>%
  mutate(L_avg = ( L1*size_mean_F[1] + L2*size_mean_F[2]+ L3*size_mean_F[3]+ L4*size_mean_F[4] + L5*size_mean_F[5]
                   + BF6*size_mean_F[6] + BF7*size_mean_F[7]+ BF8*size_mean_F[8])/ttlB )

#Only Fishery (3)
fishery_noSex = test_fishery %>%
  mutate(L1 = BF1 + BM1, L2 = BF2 + BM2, L3 = BF3 + BM3, L4 = BF4 + BM4,  L5 = BF5 + BM5) %>%
  mutate(ttlB = L1 + L2+ L3+ L4+ L5+ BF6+ BF7 + BF8) %>%
  select(20:24, 10:12,25, 19 ) %>% #Column 24 HAS TO EB CHANGED WHEN RUNNING NEW SIM TO SELECT DATETIME
  filter(as.Date(dateTime) > as.Date("31/12/2014", format= "%d/%m/%Y" )) #delete first year 'warm up' period

Lavg_fishery = fishery_noSex %>%
  mutate(L_avg = ( L1*size_mean_F[1] + L2*size_mean_F[2]+ L3*size_mean_F[3]+ L4*size_mean_F[4] + L5*size_mean_F[5]
                   + BF6*size_mean_F[6] + BF7*size_mean_F[7]+ BF8*size_mean_F[8])/ttlB )

#Both predation & fishery (4)
bothPF_noSex = test_bothFP %>%
  mutate(L1 = BF1 + BM1, L2 = BF2 + BM2, L3 = BF3 + BM3, L4 = BF4 + BM4,  L5 = BF5 + BM5) %>%
  mutate(ttlB = L1 + L2+ L3+ L4+ L5+ BF6+ BF7 + BF8) %>%
  select(20:24, 10:12,25, 19 ) %>% #Column 24 HAS TO EB CHANGED WHEN RUNNING NEW SIM TO SELECT DATETIME
  filter(as.Date(dateTime) > as.Date("31/12/2014", format= "%d/%m/%Y" )) #delete first year 'warm up' period

Lavg_bothPF = bothPF_noSex %>%
  mutate(L_avg = ( L1*size_mean_F[1] + L2*size_mean_F[2]+ L3*size_mean_F[3]+ L4*size_mean_F[4] + L5*size_mean_F[5]
                   + BF6*size_mean_F[6] + BF7*size_mean_F[7]+ BF8*size_mean_F[8])/ttlB )
dev.off()

plot(noPreasure_noSex$dateTime, L_avg_noPreasure$L_avg, las = 1,  col = 'grey', cex = 0.8,
     xlab = 'years', ylab = paste(TeX("$l$"), '[cm]'), ylim = c(1, 5) ,type = 'l', lty = 2, lwd = 3, #TeX(r"($\ell$)" )
     main = 'Changes in average size ')
lines(Lavg_pred$dateTime, Lavg_pred$L_avg, col = 'indianred3', lty=1, lwd = 2.5)
lines(Lavg_fishery$dateTime, Lavg_fishery$L_avg, col = 'skyblue4', lty=1, lwd = 2.5)
lines(Lavg_bothPF$dateTime, Lavg_bothPF$L_avg, col = 'hotpink', lty=1, lwd = 2.5)

legend("topright", legend = c('no preasure', 'only predation', 'only fishery', 'both F & P'),
       fill = c('grey', 'indianred3', 'skyblue4','hotpink'), border = NA, bty = "n", y.intersp = 0.65)



plot(noPreasure_noSex$dateTime, L_avg_noPreasure$ttlB, las = 1,  col = 'grey', cex = 0.8,
     xlab = 'years', ylab = 'biomass [?]', type = 'l', lty = 2, lwd = 3,
     main = 'Changes in biomass ', ylim = c(0, 600) )
lines(Lavg_pred$dateTime, Lavg_pred$ttlB, col = 'indianred3', lty=1, lwd = 2.5)
lines(Lavg_fishery$dateTime, Lavg_fishery$ttlB, col = 'skyblue4', lty=1, lwd = 2.5)
lines(Lavg_bothPF$dateTime, Lavg_bothPF$ttlB, col = 'hotpink', lty=1, lwd = 2.5)

legend("topright", legend = c('no preasure', 'only predation', 'only fishery', 'both F & P'),
       fill = c('grey', 'indianred3', 'skyblue4','hotpink'), border = NA, bty = "n", y.intersp = 0.6)


#### now plot L_avg vs Biomass ####

#Now group by the week number and calculate avg of the size

weekly_avg_noPreasure  <- L_avg_noPreasure %>%
  mutate(
    year = year(dateTime),
    week = isoweek(dateTime)  # or week(dateTime), depending on your preference
  ) %>%
  group_by(year, week) %>%
  summarise(avg_L = mean(L_avg, na.rm = TRUE), avg_B = mean(ttlB, na.rm = TRUE), .groups = "drop")

weekly_avg_fish  <- Lavg_fishery %>%
  mutate(
    year = year(dateTime),
    week = isoweek(dateTime)  # or week(dateTime), depending on your preference
  ) %>%
  group_by(year, week) %>%
  summarise(avg_L = mean(L_avg, na.rm = TRUE), avg_B = mean(ttlB, na.rm = TRUE), .groups = "drop")

weekly_avg_pred  <- Lavg_pred %>%
  mutate(
    year = year(dateTime),
    week = isoweek(dateTime)  # or week(dateTime), depending on your preference
  ) %>%
  group_by(year, week) %>%
  summarise(avg_L = mean(L_avg, na.rm = TRUE), avg_B = mean(ttlB, na.rm = TRUE), .groups = "drop")

weekly_avg_both  <- Lavg_bothPF %>%
  mutate(
    year = year(dateTime),
    week = isoweek(dateTime)  # or week(dateTime), depending on your preference
  ) %>%
  group_by(year, week) %>%
  summarise(avg_L = mean(L_avg, na.rm = TRUE), avg_B = mean(ttlB, na.rm = TRUE), .groups = "drop")

init15 = as.POSIXct("2015-01-01", format="%Y-%m-%d", tz = "UTC")
end15 = as.POSIXct("2015-12-31", format="%Y-%m-%d", tz = "UTC")



plot(weekly_avg_noPreasure$avg_L[weekly_avg_noPreasure$year == 2015], weekly_avg_noPreasure$avg_B[weekly_avg_noPreasure$year == 2015],
     lwd = 3, col = 'grey', lty=1, ylim = c(50,1000), xlim = c(2.5, 5) )
points(weekly_avg_fish$avg_L[weekly_avg_fish$year == 2015], weekly_avg_fish$avg_B[weekly_avg_fish$year == 2015], lwd = 2, col ='skyblue4' )
points(weekly_avg_pred$avg_L[weekly_avg_pred$year == 2015], weekly_avg_pred$avg_B[weekly_avg_pred$year == 2015], lwd = 2, col ='indianred3' )
points(weekly_avg_both$avg_L[weekly_avg_both$year == 2015], weekly_avg_both$avg_B[weekly_avg_pred$year == 2015], lwd = 2, col= 'hotpink' )

legend("topleft", legend = c('no preasure', 'only predation', 'only fishery', 'both F & P'),
       fill = c('grey', 'indianred3', 'skyblue4','hotpink'), border = NA, bty = "n", y.intersp = 0.7)


#### now plot L_avg vs Biomass ####

#Now group by the week number and calculate avg of the size

weekly_avg_noPreasure  <- L_avg_noPreasure %>%
  mutate(
    year = year(dateTime),
    week = isoweek(dateTime)  # or week(dateTime), depending on your preference
  ) %>%
  group_by(year, week) %>%
  summarise(avg_L = mean(L_avg, na.rm = TRUE), avg_B = mean(ttlB, na.rm = TRUE), .groups = "drop")

weekly_avg_fish  <- Lavg_fishery %>%
  mutate(
    year = year(dateTime),
    week = isoweek(dateTime)  # or week(dateTime), depending on your preference
  ) %>%
  group_by(year, week) %>%
  summarise(avg_L = mean(L_avg, na.rm = TRUE), avg_B = mean(ttlB, na.rm = TRUE), .groups = "drop")

weekly_avg_pred  <- Lavg_pred %>%
  mutate(
    year = year(dateTime),
    week = isoweek(dateTime)  # or week(dateTime), depending on your preference
  ) %>%
  group_by(year, week) %>%
  summarise(avg_L = mean(L_avg, na.rm = TRUE), avg_B = mean(ttlB, na.rm = TRUE), .groups = "drop")

weekly_avg_both  <- Lavg_bothPF %>%
  mutate(
    year = year(dateTime),
    week = isoweek(dateTime)  # or week(dateTime), depending on your preference
  ) %>%
  group_by(year, week) %>%
  summarise(avg_L = mean(L_avg, na.rm = TRUE), avg_B = mean(ttlB, na.rm = TRUE), .groups = "drop")

init15 = as.POSIXct("2015-01-01", format="%Y-%m-%d", tz = "UTC")
end15 = as.POSIXct("2015-12-31", format="%Y-%m-%d", tz = "UTC")


par(mfrow = c(1,1))
year_to_plot = 2015
plot(weekly_avg_noPreasure$avg_L[weekly_avg_noPreasure$year == year_to_plot], weekly_avg_noPreasure$avg_B[weekly_avg_noPreasure$year == 2015],
     lwd = 3, col = 'grey', lty=1, ylim = c(50,500), xlim = c(2.5, 5),
     xlab = "L", ylab = "Biomass", main = paste(" ", as.character(fishery_params$general_params$Fi)) )
points(weekly_avg_fish$avg_L[weekly_avg_fish$year == year_to_plot], weekly_avg_fish$avg_B[weekly_avg_fish$year == 2015], lwd = 2, col ='skyblue4' )
points(weekly_avg_pred$avg_L[weekly_avg_pred$year == year_to_plot], weekly_avg_pred$avg_B[weekly_avg_pred$year == 2015], lwd = 2, col ='indianred3' )
points(weekly_avg_both$avg_L[weekly_avg_both$year == year_to_plot], weekly_avg_both$avg_B[weekly_avg_pred$year == 2015], lwd = 2, col= 'hotpink' )

legend("topleft", legend = c('no preasure', 'only predation', 'only fishery', 'both F & P'),
       fill = c('grey', 'indianred3', 'skyblue4','hotpink'), border = NA, bty = "n", y.intersp = 0.7)



############ TEST V5 #################

state_st = c(P = 2, E = 0.1, L= 0.0928 *1.5 ,
            BF = BF+0.1, BM = BM+0.1)


#### TEST V5 (1 out of 4) WITHOUT ANY KIND OF PREASURE (P NEITHER F) ####


noPreasure_params = parameters_solv
noPreasure_params$general_params$Fi = 0 #No fishery
noPreasure_params$general_params$Imax_ik = 0 #No predation

start <- Sys.time()
test_noPreasure_stPred <- solver_sizeClass.v5(t = t_5years, state = state_st, parameters = noPreasure_params, temperature_dataSet = temperature_14_18)
print(Sys.time() - start)

#Plot size spectra
plot_sizeSpectra(test_noPreasure_stPred, 2015, title = "NP")

#### TEST V5 (2 out of 4) ONLY PREDATION ####

pred_params = parameters_solv
pred_params$general_params$Fi = 0 #No fishery


start <- Sys.time()
test_predation_stPred <- solver_sizeClass.v5(t = t_5years, state = state_st, parameters = pred_params, temperature_dataSet = temperature_14_18)
print(Sys.time() - start)

#Plot size spectra
plot_sizeSpectra(test_predation_stPred, 2018, title = "Pred")

#### TEST V5 (3 out of 4) ONLY FISHERY ####

fishery_params = parameters_solv
fishery_params$general_params$Imax_ik = 0 #No predation
fishery_params$general_params$Fi = 0.025 #reduced fishery

start <- Sys.time()
test_fishery_stPred <- solver_sizeClass.v5(t = t_5years, state = state_st, parameters = fishery_params, temperature_dataSet = temperature_14_18)
print(Sys.time() - start)

#Plot size spectra
plot_sizeSpectra(test_fishery_stPred, 2015, title = "Fi")

#### TEST V5 (4 out of 4) PREDATION AND FISHERY####

bothPF_params = parameters_solv
bothPF_params$general_params$Fi = 0.025 #reduced fishery

start <- Sys.time()
test_bothFP_stPred <- solver_sizeClass.v5(t = t_5years, state = state_st, parameters = bothPF_params, temperature_dataSet = temperature_14_18)
print(Sys.time() - start)

#Plot size spectra
plot_sizeSpectra(test_bothFP_stPred, 201, title = "F&P")

#stacked area plot

test_pred.V5_red <- test_predation_stPred[ , c(3:19)] # delete timesteps column, plancton  and predator column

#cc_long_2 <- gather(cc_2, key = "Type", value = "Value", P, E, L, J, J2, J3, J4, J5, A1, A2)
test_pred_long <- test_pred.V5_red %>%
  pivot_longer(-dateTime, names_to = "variable", values_to = "value")

ggplot(test_pred_long, aes(x = dateTime, y = value, fill = variable)) +
  geom_area() +
  scale_fill_viridis_d() +
  labs(
    title = "System incl. predation ",
    x = "Date",
    y = "Biomass",
    fill = "Size class"
  )



##### plot only January for estimation of initial conditions in Thesis#####


sol_15_jan = test_predation_stPred %>%
  filter(as.Date(test_predation_stPred$dateTime) > as.Date("31/12/2014", format= "%d/%m/%Y" ) &
           as.Date(test_predation_stPred$dateTime) < as.Date("01/02/2015", format= "%d/%m/%Y") ) %>%
  select(BF1, BF2, BF3, BF4, BF5, BF6, BF7, BF8, BM1, BM2, BM3, BM4, BM5)


cols <- c("f" = "#8856a7", "m" = "#add8e6")

# Define size classes to include and their order
class_labels <- c("BL1", "BL2", "BL3", "BL4", "BL5", "BL6", "BL7", "BL8")
label_map <- setNames(class_labels, class_labels)



# 1. Extract female and male size class columns
female_cols <- paste0("BF", 1:8)
male_cols   <- paste0("BM", 1:5)

# 2. Compute means
f_means <- colMeans(sol_15_jan[, female_cols])
m_means <- colMeans(sol_15_jan[, male_cols])

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
total_biomass <- sum(mat)
mat_prop <- mat / total_biomass

par(mfrow = c(1,1))

# Plot barplot V4 (background because B is higher)
bar_x = barplot(
  #mat,
  #col = cols[rownames(mat)],
  mat_prop,
  col = cols[rownames(mat_prop)],
  main = 'L50 = 4.49 cm',#paste("Size distribution in January"),
  names.arg =  as.character(size_mean_F), #label_map,
  ylab = 'Proportion',
  xlab = 'mean L of size class',
  las = 1,
  cex.main = 1.2,
  cex.axis = 1.2,
  cex.names = 1.2,
  cex.lab = 1.5,
  ylim = c(0, .4),#c(0,25),
  space = 0,
  beside = FALSE,
  #legend.text = rownames(mat_prop),#rownames(mat),
  #args.legend = list(x = "topleft", bty = "n", inset = 0.02, cex = 1.5)
)

#Overlap selective fishery curve (probit):
par(new = TRUE)

# Match kernel curve to barplot's x-axis midpoints
# Assume size_mean_F matches the size classes
selective_fi <- sel_probit(size_mean_F)

# Plot the kernel as a line on top
plot(x = bar_x, y = selective_fi, type = 'b', axes = FALSE, xlab = '', ylab = '', col = 'darkgreen', lwd = 2)

legend("topleft", legend = c("M", "F", "Sel. Fishery"),
       fill = c( "#add8e6", "#8856a7", "darkgreen"), cex = 1.2, bty = "n", border = NA) #, inset = 0.02 , bty = "n"
# Optional: Add area under curve with lines
polygon(c(bar_x, rev(bar_x)),
        c(rep(0, length(bar_x)), rev(selective_fi)),
        col = adjustcolor("darkgreen", alpha.f = 0.15),
        border = NA)

# Add right axis (optional)
axis(side = 4)

#### @Thesis, Methodology Biomass estimation for initial conditions ####
#Calculate distribution of 14.71  kg / km2 ... OLD:7.59 kg / km2 :

size_classes_proportions = colSums(mat_prop)

size_dist_biomass = 14.71*size_classes_proportions

est_biomass_dist = 100*size_dist_biomass/5 #considering that fishery is 5% of actual stock

F_proportions =mat["f", ] / colSums(mat)
M_proportions =mat["m", ] / colSums(mat)

est_biomass_dist_F = est_biomass_dist*F_proportions #init condition biomass Fems
est_biomass_dist_M = est_biomass_dist*M_proportions #init condition biomass Masc


#### PLOTS TO COMPARE ALL ####

#(1) L_avg over the time:
size_classes = female_cols <- paste0("BL", 1:8)

#calculate weighted average for each time step:

#first join F and M in one column
#Without preasure (1)
noPreasure_noSex_stPred = test_noPreasure_stPred %>%
  mutate(L1 = BF1 + BM1, L2 = BF2 + BM2, L3 = BF3 + BM3, L4 = BF4 + BM4,  L5 = BF5 + BM5) %>%
  mutate(ttlB = L1 + L2+ L3+ L4+ L5+ BF6+ BF7 + BF8) %>%
  select(19:23, 10:12,24, 18 ) %>% #select only the colums unisex
  filter(as.Date(dateTime) > as.Date("31/12/2014", format= "%d/%m/%Y" )) #delete first year 'warm up' period


L_avg_noPreasure_stPred = noPreasure_noSex_stPred %>%
  mutate(L_avg = ( L1*size_mean_F[1] + L2*size_mean_F[2]+ L3*size_mean_F[3]+ L4*size_mean_F[4] + L5*size_mean_F[5]
                   + BF6*size_mean_F[6] + BF7*size_mean_F[7]+ BF8*size_mean_F[8])/ttlB )

#Only Predation (2)
predation_noSex_stPred = test_predation_stPred %>%
  mutate(L1 = BF1 + BM1, L2 = BF2 + BM2, L3 = BF3 + BM3, L4 = BF4 + BM4,  L5 = BF5 + BM5) %>%
  mutate(ttlB = L1 + L2+ L3+ L4+ L5+ BF6+ BF7 + BF8) %>%
  select(19:23, 10:12,24, 18) %>% #Column 24 HAS TO EB CHANGED WHEN RUNNING NEW SIM TO SELECT DATETIME
  filter(as.Date(dateTime) > as.Date("31/12/2014", format= "%d/%m/%Y" )) #delete first year 'warm up' period

Lavg_pred_stPred = predation_noSex_stPred %>%
  mutate(L_avg = ( L1*size_mean_F[1] + L2*size_mean_F[2]+ L3*size_mean_F[3]+ L4*size_mean_F[4] + L5*size_mean_F[5]
                   + BF6*size_mean_F[6] + BF7*size_mean_F[7]+ BF8*size_mean_F[8])/ttlB )

#Only Fishery (3)
fishery_noSex_stPred = test_fishery_stPred %>%
  mutate(L1 = BF1 + BM1, L2 = BF2 + BM2, L3 = BF3 + BM3, L4 = BF4 + BM4,  L5 = BF5 + BM5) %>%
  mutate(ttlB = L1 + L2+ L3+ L4+ L5+ BF6+ BF7 + BF8) %>%
  select(19:23, 10:12,24, 18) %>% #Column 24 HAS TO EB CHANGED WHEN RUNNING NEW SIM TO SELECT DATETIME
  filter(as.Date(dateTime) > as.Date("31/12/2014", format= "%d/%m/%Y" )) #delete first year 'warm up' period

Lavg_fishery_stPred = fishery_noSex_stPred %>%
  mutate(L_avg = ( L1*size_mean_F[1] + L2*size_mean_F[2]+ L3*size_mean_F[3]+ L4*size_mean_F[4] + L5*size_mean_F[5]
                   + BF6*size_mean_F[6] + BF7*size_mean_F[7]+ BF8*size_mean_F[8])/ttlB )

#Both predation & fishery (4)
bothPF_noSex_stPred = test_bothFP_stPred %>%
  mutate(L1 = BF1 + BM1, L2 = BF2 + BM2, L3 = BF3 + BM3, L4 = BF4 + BM4,  L5 = BF5 + BM5) %>%
  mutate(ttlB = L1 + L2+ L3+ L4+ L5+ BF6+ BF7 + BF8) %>%
  select(19:23, 10:12,24, 18 ) %>% #Column 24 HAS TO EB CHANGED WHEN RUNNING NEW SIM TO SELECT DATETIME
  filter(as.Date(dateTime) > as.Date("31/12/2014", format= "%d/%m/%Y" )) #delete first year 'warm up' period

Lavg_bothPF_stPred = bothPF_noSex_stPred %>%
  mutate(L_avg = ( L1*size_mean_F[1] + L2*size_mean_F[2]+ L3*size_mean_F[3]+ L4*size_mean_F[4] + L5*size_mean_F[5]
                   + BF6*size_mean_F[6] + BF7*size_mean_F[7]+ BF8*size_mean_F[8])/ttlB )
dev.off()

plot(noPreasure_noSex_stPred$dateTime, L_avg_noPreasure_stPred$L_avg, las = 1,  col = 'grey', cex = 0.8,
     xlab = 'years', ylab = paste(TeX("$l$"), '[cm]'), ylim = c(1, 5) ,type = 'l', lty = 2, lwd = 3, #TeX(r"($\ell$)" )
     main = 'Changes in average size ')
lines(Lavg_pred_stPred$dateTime, Lavg_pred_stPred$L_avg, col = 'indianred3', lty=1, lwd = 2.5)
lines(Lavg_fishery_stPred$dateTime, Lavg_fishery_stPred$L_avg, col = 'skyblue4', lty=1, lwd = 2.5)
lines(Lavg_bothPF_stPred$dateTime, Lavg_bothPF_stPred$L_avg, col = 'hotpink', lty=1, lwd = 2.5)

legend("topright", legend = c('no preasure', 'only predation', 'only fishery', 'both F & P'),
       fill = c('grey', 'indianred3', 'skyblue4','hotpink'), border = NA, bty = "n", y.intersp = 0.65)



plot(noPreasure_noSex_stPred$dateTime, L_avg_noPreasure_stPred$ttlB, las = 1,  col = 'grey', cex = 0.8,
     xlab = 'years', ylab = 'biomass [?]', type = 'l', lty = 2, lwd = 3,
     main = 'Changes in biomass ', ylim = c(0, 600) )
lines(Lavg_pred_stPred$dateTime, Lavg_pred_stPred$ttlB, col = 'indianred3', lty=1, lwd = 2.5)
lines(Lavg_fishery_stPred$dateTime, Lavg_fishery_stPred$ttlB, col = 'skyblue4', lty=1, lwd = 2.5)
lines(Lavg_bothPF_stPred$dateTime, Lavg_bothPF_stPred$ttlB, col = 'hotpink', lty=1, lwd = 2.5)

legend("topright", legend = c('no preasure', 'only predation', 'only fishery', 'both F & P'),
       fill = c('grey', 'indianred3', 'skyblue4','hotpink'), border = NA, bty = "n", y.intersp = 0.6)


#### Plot L_avg vs Biomass ####

#Now group by the week number and calculate avg of the size

weekly_avg_noPreasure_stPred  <- L_avg_noPreasure_stPred %>%
  mutate(
    year = year(dateTime),
    week = isoweek(dateTime)  # or week(dateTime), depending on your preference
  ) %>%
  group_by(year, week) %>%
  summarise(avg_L = mean(L_avg, na.rm = TRUE), avg_B = mean(ttlB, na.rm = TRUE), .groups = "drop")

weekly_avg_fish_stPred  <- Lavg_fishery_stPred %>%
  mutate(
    year = year(dateTime),
    week = isoweek(dateTime)  # or week(dateTime), depending on your preference
  ) %>%
  group_by(year, week) %>%
  summarise(avg_L = mean(L_avg, na.rm = TRUE), avg_B = mean(ttlB, na.rm = TRUE), .groups = "drop")

weekly_avg_pred_stPred  <- Lavg_pred_stPred %>%
  mutate(
    year = year(dateTime),
    week = isoweek(dateTime)  # or week(dateTime), depending on your preference
  ) %>%
  group_by(year, week) %>%
  summarise(avg_L = mean(L_avg, na.rm = TRUE), avg_B = mean(ttlB, na.rm = TRUE), .groups = "drop")

weekly_avg_both_stPred  <- Lavg_bothPF_stPred %>%
  mutate(
    year = year(dateTime),
    week = isoweek(dateTime)  # or week(dateTime), depending on your preference
  ) %>%
  group_by(year, week) %>%
  summarise(avg_L = mean(L_avg, na.rm = TRUE), avg_B = mean(ttlB, na.rm = TRUE), .groups = "drop")

init15 = as.POSIXct("2015-01-01", format="%Y-%m-%d", tz = "UTC")
end15 = as.POSIXct("2015-12-31", format="%Y-%m-%d", tz = "UTC")



plot(weekly_avg_noPreasure_stPred$avg_L[weekly_avg_noPreasure_stPred$year == 2015], weekly_avg_noPreasure_stPred$avg_B[weekly_avg_noPreasure_stPred$year == 2015],
     lwd = 3, col = 'grey', lty=1, ylim = c(15,750), xlim = c(2.5, 5) )
points(weekly_avg_fish_stPred$avg_L[weekly_avg_fish_stPred$year == 2015], weekly_avg_fish_stPred$avg_B[weekly_avg_fish_stPred$year == 2015], lwd = 2, col ='skyblue4' )
points(weekly_avg_pred_stPred$avg_L[weekly_avg_pred_stPred$year == 2015], weekly_avg_pred_stPred$avg_B[weekly_avg_pred_stPred$year == 2015], lwd = 2, col ='indianred3' )
points(weekly_avg_both_stPred$avg_L[weekly_avg_both_stPred$year == 2015], weekly_avg_both_stPred$avg_B[weekly_avg_pred_stPred$year == 2015], lwd = 2, col= 'hotpink' )

legend("topleft", legend = c('no preasure', 'only predation', 'only fishery', 'both F & P'),
       fill = c('grey', 'indianred3', 'skyblue4','hotpink'), border = NA, bty = "n", y.intersp = 0.7)


