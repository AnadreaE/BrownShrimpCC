#HERE SIMULATIONS WITH SOLVER VERSION V3 WHERE FISHERY AND PREDATION HAS BEEN ALREADY IMPLEMENTED

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

state = c(P = 2, E = 0.1, L= 0.0928 *1.5 ,
           J_f = 2, J2_f = 2, J3_f = 2 , J4_f = 1.8, J5_f = 1.8,
           J_m = 2 , J2_m = 2, J3_m = 2, J4_m = 1.8, J5_m = 1.8,
           A1_f = 5.443*0.5, A1_m = 5.443*0.5,
           A2 = 5.443*0.375, A3=5.443*0.125, #males doesn't reach this size classes
           Pred = 0.05 )

state_rc = c(P = 2, E = 0.1, L= 0.0928 *1.5 ,
             J_f = 0.1, J2_f = 0.1, J3_f = 0.1 , J4_f = 0.1, J5_f = 0.1,
             J_m = 0.1 , J2_m = 0.1, J3_m = 0.1, J4_m = 0.1, J5_m = 0.1,
             A1_f = 0.1, A1_m = 0.1,
             A2 = 0.1, A3=0.1, #males doesn't reach this size classes
             Pred = 0.05 )


#Ref: Larv: 0.0928 kg/1000m2 but lets say that wihtout predation this is 50% bigger
#Ref Fem: 5.443 kg/1000m2 has to be distributed over the 3 Fem adult size classes: we know that A1 > A2 > A3 in biomass
#therefore lets consider following ratio: A1: 0.5, A2: 0.375, A3: 0.125
#since numbers are bigger now, increase juveniles

t_5years = seq(0,length(temperature_14_18$temperature)-0.1, by = 0.1)


#### TEST (1 out of 4) WITHOUT ANY KIND OF PREASURE (P NEITHER F) ####

noPreasure_params = parameters_solv
noPreasure_params$general_params$Fi = 0 #No fishery
noPreasure_params$general_params$Imax_ik = 0 #No predation


start <- Sys.time()
test_noPreasure <- solver_sizeClass_sex.v3(t = t_5years, state = state_rc, parameters = noPreasure_params, temperature_dataSet = temperature_14_18)
print(Sys.time() - start)


#stacked area plot

#test_noPreasure <- test_noPreasure[ , -c(1,2,19)] # delete timesteps column and predator column

#cc_long_2 <- gather(cc_2, key = "Type", value = "Value", P, E, L, J, J2, J3, J4, J5, A1, A2)
test_noPreasure_long <- test_noPreasure %>%
  pivot_longer(-dateTime, names_to = "variable", values_to = "value")

ggplot(test_noPreasure_long, aes(x = dateTime, y = value, fill = variable)) +
  geom_area() +
  scale_fill_viridis_d() +
  labs(
    title = "Without preasure of any kind, no F nor P",
    x = "Date",
    y = "Biomass",
    fill = "Size class"
  )


# Bar plot F & M together

M_F_2015 = test_noPreasure %>%
  filter( as.Date(dateTime) > as.Date("31/12/2014", format= "%d/%m/%Y" ) & as.Date(dateTime) < as.Date("01/01/2016", format= "%d/%m/%Y" )  ) %>%
  select(E, L, J_f, J2_f, J3_f, J4_f, J5_f, A1_f, J_m, J2_m, J3_m, J4_m, J5_m, A1_m, A2, A3)

# Loop through 4 quarters
par(mfrow = c(2, 2))

# Define colors for each sex: f = light blue, m = steel blue, U = gray
cols <- c("f" = "#8856a7", "m" = "#add8e6", "U" = "gray80")

# Define x-axis label order (class names)
x_labels_vec <- c("E", "L", "J", "J2", "J3", "J4", "J5", "A1", "A2", "A3")
label_map <- c(
  "E" = "E",
  "L" = "i-0.6",
  "J" = "0.6-1",
  "J2" = "1-2",
  "J3" = "2-3",
  "J4" = "3-4",
  "J5" = "4-5",
  "A1" = "5-6",
  "A2" = "6-7",
  "A3" = "7-Linf"
)
width_vector_all <- rep(1, length(x_labels_vec))  # Optional bar widths

# Loop through 4 quarters
for (q in 1:4) {
  # Indices for each quarter
  init_idx <- 912 * (q - 1) + 1
  final_idx <- 912 * q

  # Subset data for quarter
  quarter_data <- M_F_2015[init_idx:final_idx, ]

  # 1. Get sexed class columns (e.g., J_f, J2_m, A1_f, etc.)
  sexed_df <- quarter_data %>%
    select(matches("^(J|A1)\\d?_f$|^(J|A1)\\d?_m$"))

  # 2. Get unsexed class columns: E and L
  unsexed_df1 <- quarter_data %>%
    select(E, L)

  unsexed_df2 <- quarter_data %>%
    select( A2, A3)

  # 3. Compute column means
  sexed_means <- colMeans(sexed_df)
  unsexed1_means <- colMeans(unsexed_df1)
  unsexed2_means <- colMeans(unsexed_df2)

  # 4. Create long-format data frame for sexed
  sexed_df_long <- data.frame(
    class_sex = names(sexed_means),
    value = as.numeric(sexed_means)
  ) %>%
    separate(class_sex, into = c("class", "sex"), sep = "_")

  # 5. Create long-format data frame for unsexed, with sex = "U"
  unsexed1_df_long <- data.frame(
    class = names(unsexed1_means),
    value = as.numeric(unsexed1_means),
    sex = "U"
  )

  unsexed2_df_long <- data.frame(
    class = names(unsexed2_means),
    value = as.numeric(unsexed2_means),
    sex = "U"
  )

  # 6. Combine both
  means_df <- bind_rows(unsexed1_df_long , sexed_df_long, unsexed2_df_long)

  # 7. Set class as factor to preserve desired order
  means_df$class <- factor(means_df$class, levels = x_labels_vec)

  means_df$label <- label_map[as.character(means_df$class)]

  # 8. Pivot to wide format for barplot
  means_mat <- means_df %>%
    pivot_wider(names_from = sex, values_from = value, values_fill = 0 ) %>%
    arrange(class) %>%
    select(any_of(c("f", "m", "U")))  # Ensure consistent column order

  # 9. Transpose to matrix for barplot (rows = sex, columns = classes)
  mat <- t(as.matrix(means_mat))

  # 10. Draw barplot
  barplot(mat,
          col = cols[rownames(mat)],
          main = paste('NP - F & M - Q', q, '2015'),
          names.arg = label_map[levels(means_df$class)], #x_labels_vec,
          las = 2,
          cex.main = 2,
          cex.axis = 1.5,
          cex.names = 1.5,
          ylim = c(0,20), #c(0, max(colSums(mat)) * 1.1),
          width = width_vector_all,
          space = 0,
          beside = FALSE,
          legend.text = rownames(mat),
          args.legend = list(x = "topright", bty = "n", inset = 0.02, cex = 1)
  )
}



#### TEST (2 out of 4) ONLY PREDATION ####

pred_params = parameters_solv
pred_params$general_params$Fi = 0 #No fishery


start <- Sys.time()
test_pred <- solver_sizeClass_sex.v3(t = t_5years, state = state, parameters = pred_params, temperature_dataSet = temperature_14_18)
print(Sys.time() - start)

#stacked area plot

start_date <- as.POSIXct("2014-01-01", format="%Y-%m-%d", tz = "UTC")

test_pred <- mutate(test_pred, dateTime = start_date + test_pred$time * 86400) #86400 seconds in one day
test_pred <- test_pred[ , -c(1,2,19)] # delete timesteps column, plancton  and predator column

#cc_long_2 <- gather(cc_2, key = "Type", value = "Value", P, E, L, J, J2, J3, J4, J5, A1, A2)
test_pred_long <- test_pred %>%
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


# Bar plot F & M together

pred_M_F_2015 = test_pred %>%
  filter( as.Date(dateTime) > as.Date("31/12/2014", format= "%d/%m/%Y" ) & as.Date(dateTime) < as.Date("01/01/2016", format= "%d/%m/%Y" )  ) %>%
  select(E, L, J_f, J2_f, J3_f, J4_f, J5_f, A1_f, J_m, J2_m, J3_m, J4_m, J5_m, A1_m, A2, A3)

# Loop through 4 quarters
par(mfrow = c(2, 2))

# Define colors for each sex: f = light blue, m = steel blue, U = gray
cols <- c("f" = "#8856a7", "m" = "#add8e6", "U" = "gray80")

# Define x-axis label order (class names)
x_labels_vec <- c("E", "L", "J", "J2", "J3", "J4", "J5", "A1", "A2", "A3")
label_map <- c(
  "E" = "E",
  "L" = "i-0.6",
  "J" = "0.6-1",
  "J2" = "1-2",
  "J3" = "2-3",
  "J4" = "3-4",
  "J5" = "4-5",
  "A1" = "5-6",
  "A2" = "6-7",
  "A3" = "7-Linf"
)
#width_vector_all <- rep(1, length(x_labels_vec))  # Optional bar widths

# Loop through 4 quarters
for (q in 1:4) {
  # Indices for each quarter
  init_idx <- 912 * (q - 1) + 1
  final_idx <- 912 * q

  # Subset data for quarter
  quarter_data <- pred_M_F_2015[init_idx:final_idx, ]

  # 1. Get sexed class columns (e.g., J_f, J2_m, A1_f, etc.)
  sexed_df <- quarter_data %>%
    select(matches("^(J|A1)\\d?_f$|^(J|A1)\\d?_m$"))

  # 2. Get unsexed class columns: E and L
  unsexed_df1 <- quarter_data %>%
    select(E, L)

  unsexed_df2 <- quarter_data %>%
    select( A2, A3)

  # 3. Compute column means
  sexed_means <- colMeans(sexed_df)
  unsexed1_means <- colMeans(unsexed_df1)
  unsexed2_means <- colMeans(unsexed_df2)

  # 4. Create long-format data frame for sexed
  sexed_df_long <- data.frame(
    class_sex = names(sexed_means),
    value = as.numeric(sexed_means)
  ) %>%
    separate(class_sex, into = c("class", "sex"), sep = "_")

  # 5. Create long-format data frame for unsexed, with sex = "U"
  unsexed1_df_long <- data.frame(
    class = names(unsexed1_means),
    value = as.numeric(unsexed1_means),
    sex = "U"
  )

  unsexed2_df_long <- data.frame(
    class = names(unsexed2_means),
    value = as.numeric(unsexed2_means),
    sex = "U"
  )

  # 6. Combine both
  means_df <- bind_rows(unsexed1_df_long , sexed_df_long, unsexed2_df_long)

  # 7. Set class as factor to preserve desired order
  means_df$class <- factor(means_df$class, levels = x_labels_vec)

  means_df$label <- label_map[as.character(means_df$class)]

  # 8. Pivot to wide format for barplot
  means_mat <- means_df %>%
    pivot_wider(names_from = sex, values_from = value, values_fill = 0 ) %>%
    arrange(class) %>%
    select(any_of(c("f", "m", "U")))  # Ensure consistent column order

  # 9. Transpose to matrix for barplot (rows = sex, columns = classes)
  mat <- t(as.matrix(means_mat))

  # 10. Draw barplot
  barplot(mat,
          col = cols[rownames(mat)],
          main = paste('P - F & M - Q', q, '2015'),
          names.arg = label_map[levels(means_df$class)], #x_labels_vec,
          las = 2,
          cex.main = 2,
          cex.axis = 1.5,
          cex.names = 1.5,
          #ylim = c(0, max(colSums(mat)) * 1.1),
         # width = width_vector_all,
          space = 0,
          beside = FALSE,
          legend.text = rownames(mat),
         ylim = c(0,15),
          args.legend = list(x = "topright", bty = "n", inset = 0.02, cex = 1)
  )
}



#### TEST (3 out of 4) ONLY FISHERY ####

fishery_params = parameters_solv
fishery_params$general_params$Imax_ik = 0 #No predation
#fishery_params$general_params$Fi = 0.01 #reduce fishery

start <- Sys.time()
test_fishery <- solver_sizeClass_sex.v3(t = t_5years, state = state, parameters = fishery_params, temperature_dataSet = temperature_14_18)
print(Sys.time() - start)

#stacked area plot

start_date <- as.POSIXct("2014-01-01", format="%Y-%m-%d", tz = "UTC")

test_fishery <- mutate(test_fishery, dateTime = start_date + test_fishery$time * 86400) #86400 seconds in one day
test_fishery <- test_fishery[ , -c(1,2,19)] # delete timesteps column and predator column

#cc_long_2 <- gather(cc_2, key = "Type", value = "Value", P, E, L, J, J2, J3, J4, J5, A1, A2)
test_fishery_long <- test_fishery %>%
  pivot_longer(-dateTime, names_to = "variable", values_to = "value")

ggplot(test_fishery_long, aes(x = dateTime, y = value, fill = variable)) +
  geom_area() +
  scale_fill_viridis_d() +
  labs(
    title = "System incl. fishery ",
    x = "Date",
    y = "Biomass",
    fill = "Size class"
  )


# Bar plot F & M together

fishery_M_F_2015 = test_fishery %>%
  filter( as.Date(dateTime) > as.Date("31/12/2014", format= "%d/%m/%Y" ) & as.Date(dateTime) < as.Date("01/01/2016", format= "%d/%m/%Y" )  ) %>%
  select(E, L, J_f, J2_f, J3_f, J4_f, J5_f, A1_f, J_m, J2_m, J3_m, J4_m, J5_m, A1_m, A2, A3)

# Loop through 4 quarters
par(mfrow = c(2, 2))

# Define colors for each sex: f = light blue, m = steel blue, U = gray
cols <- c("f" = "#8856a7", "m" = "#add8e6", "U" = "gray80")

# Define x-axis label order (class names)
x_labels_vec <- c("E", "L", "J", "J2", "J3", "J4", "J5", "A1", "A2", "A3")
label_map <- c(
  "E" = "E",
  "L" = "i-0.6",
  "J" = "0.6-1",
  "J2" = "1-2",
  "J3" = "2-3",
  "J4" = "3-4",
  "J5" = "4-5",
  "A1" = "5-6",
  "A2" = "6-7",
  "A3" = "7-Linf"
)
#width_vector_all <- rep(1, length(x_labels_vec))  # Optional bar widths

# Loop through 4 quarters
for (q in 1:4) {
  # Indices for each quarter
  init_idx <- 912 * (q - 1) + 1
  final_idx <- 912 * q

  # Subset data for quarter
  quarter_data <- fishery_M_F_2015[init_idx:final_idx, ]

  # 1. Get sexed class columns (e.g., J_f, J2_m, A1_f, etc.)
  sexed_df <- quarter_data %>%
    select(matches("^(J|A1)\\d?_f$|^(J|A1)\\d?_m$"))

  # 2. Get unsexed class columns: E and L
  unsexed_df1 <- quarter_data %>%
    select(E, L)

  unsexed_df2 <- quarter_data %>%
    select( A2, A3)

  # 3. Compute column means
  sexed_means <- colMeans(sexed_df)
  unsexed1_means <- colMeans(unsexed_df1)
  unsexed2_means <- colMeans(unsexed_df2)

  # 4. Create long-format data frame for sexed
  sexed_df_long <- data.frame(
    class_sex = names(sexed_means),
    value = as.numeric(sexed_means)
  ) %>%
    separate(class_sex, into = c("class", "sex"), sep = "_")

  # 5. Create long-format data frame for unsexed, with sex = "U"
  unsexed1_df_long <- data.frame(
    class = names(unsexed1_means),
    value = as.numeric(unsexed1_means),
    sex = "U"
  )

  unsexed2_df_long <- data.frame(
    class = names(unsexed2_means),
    value = as.numeric(unsexed2_means),
    sex = "U"
  )

  # 6. Combine both
  means_df <- bind_rows(unsexed1_df_long , sexed_df_long, unsexed2_df_long)

  # 7. Set class as factor to preserve desired order
  means_df$class <- factor(means_df$class, levels = x_labels_vec)

  means_df$label <- label_map[as.character(means_df$class)]

  # 8. Pivot to wide format for barplot
  means_mat <- means_df %>%
    pivot_wider(names_from = sex, values_from = value, values_fill = 0 ) %>%
    arrange(class) %>%
    select(any_of(c("f", "m", "U")))  # Ensure consistent column order

  # 9. Transpose to matrix for barplot (rows = sex, columns = classes)
  mat <- t(as.matrix(means_mat))

  # 10. Draw barplot
  barplot(mat,
          col = cols[rownames(mat)],
          main = paste('Fi - F & M - Q', q, '2015'),
          names.arg = label_map[levels(means_df$class)], #x_labels_vec,
          las = 2,
          cex.main = 2,
          cex.axis = 1.5,
          cex.names = 1.5,
          ylim = c(0,8), #c(0, max(colSums(mat)) * 1.1),
          # width = width_vector_all,
          space = 0,
          beside = FALSE,
          legend.text = rownames(mat),
          args.legend = list(x = "topright", bty = "n", inset = 0.02, y.intersp = 0.5)
  )
}


#### TEST (4 out of 4) PREDATION AND FISHERY####

bothPF_params = parameters_solv
bothPF_params$general_params$Fi = 0.1 #increased fishery

state_rc = c(P = 2, E = 0.1, L= 0.0928 *1.5 ,
          J_f = 0.1, J2_f = 0.1, J3_f = 0.1 , J4_f = 0.1, J5_f = 0.1,
          J_m = 0.1 , J2_m = 0.1, J3_m = 0.1, J4_m = 0.1, J5_m = 0.1,
          A1_f = 0.1, A1_m = 0.1,
          A2 = 0.1, A3=0.1, #males doesn't reach this size classes
          Pred = 0.05 )

#profvis({
  start <- Sys.time()
  test_bothPF <- solver_sizeClass_sex.v3(t = t_5years, state = state_rc, parameters = bothPF_params, temperature_dataSet = temperature_14_18)
  print(Sys.time() - start)
#})



#stacked area plot

test_bothPF <- test_bothPF[ , -c(1,2,19)] # delete timesteps column and predator column

#cc_long_2 <- gather(cc_2, key = "Type", value = "Value", P, E, L, J, J2, J3, J4, J5, A1, A2)
test_bothPF_long <- test_bothPF %>%
  pivot_longer(-dateTime, names_to = "variable", values_to = "value")

ggplot(test_bothPF_long, aes(x = dateTime, y = value, fill = variable)) +
  geom_area() +
  scale_fill_viridis_d() +
  labs(
    title = "System incl. predation and fishery ",
    x = "Date",
    y = "Biomass",
    fill = "Size class"
  )


# Bar plot F & M together

fishery_M_F_2015 = test_bothPF %>%
  filter( as.Date(dateTime) > as.Date("31/12/2014", format= "%d/%m/%Y" ) & as.Date(dateTime) < as.Date("01/01/2016", format= "%d/%m/%Y" )  ) %>%
  select(E, L, J_f, J2_f, J3_f, J4_f, J5_f, A1_f, J_m, J2_m, J3_m, J4_m, J5_m, A1_m, A2, A3)

# Loop through 4 quarters
par(mfrow = c(2, 2))

# Define colors for each sex: f = light blue, m = steel blue, U = gray
cols <- c("f" = "#8856a7", "m" = "#add8e6", "U" = "gray80")

# Define x-axis label order (class names)
x_labels_vec.v3 <- c("E", "L", "J", "J2", "J3", "J4", "J5", "A1", "A2", "A3")
label_map.v3 <- c(
  "E" = "E",
  "L" = "i-0.6",
  "J" = "0.6-1",
  "J2" = "1-2",
  "J3" = "2-3",
  "J4" = "3-4",
  "J5" = "4-5",
  "A1" = "5-6",
  "A2" = "6-7",
  "A3" = "7-Linf"
)
#width_vector_all <- rep(1, length(x_labels_vec))  # Optional bar widths

# Loop through 4 quarters
for (q in 1:4) {
  # Indices for each quarter
  init_idx <- 912 * (q - 1) + 1
  final_idx <- 912 * q

  # Subset data for quarter
  quarter_data <- fishery_M_F_2015[init_idx:final_idx, ]

  # 1. Get sexed class columns (e.g., J_f, J2_m, A1_f, etc.)
  sexed_df <- quarter_data %>%
    select(matches("^(J|A1)\\d?_f$|^(J|A1)\\d?_m$"))

  # 2. Get unsexed class columns: E and L
  unsexed_df1 <- quarter_data %>%
    select(E, L)

  unsexed_df2 <- quarter_data %>%
    select( A2, A3)

  # 3. Compute column means
  sexed_means <- colMeans(sexed_df)
  unsexed1_means <- colMeans(unsexed_df1)
  unsexed2_means <- colMeans(unsexed_df2)

  # 4. Create long-format data frame for sexed
  sexed_df_long <- data.frame(
    class_sex = names(sexed_means),
    value = as.numeric(sexed_means)
  ) %>%
    separate(class_sex, into = c("class", "sex"), sep = "_")

  # 5. Create long-format data frame for unsexed, with sex = "U"
  unsexed1_df_long <- data.frame(
    class = names(unsexed1_means),
    value = as.numeric(unsexed1_means),
    sex = "U"
  )

  unsexed2_df_long <- data.frame(
    class = names(unsexed2_means),
    value = as.numeric(unsexed2_means),
    sex = "U"
  )

  # 6. Combine both
  means_df <- bind_rows(unsexed1_df_long , sexed_df_long, unsexed2_df_long)

  # 7. Set class as factor to preserve desired order
  means_df$class <- factor(means_df$class, levels = x_labels_vec)

  means_df$label <- label_map[as.character(means_df$class)]

  # 8. Pivot to wide format for barplot
  means_mat <- means_df %>%
    pivot_wider(names_from = sex, values_from = value, values_fill = 0 ) %>%
    arrange(class) %>%
    select(any_of(c("f", "m", "U")))  # Ensure consistent column order

  # 9. Transpose to matrix for barplot (rows = sex, columns = classes)
  mat <- t(as.matrix(means_mat))

  # 10. Draw barplot
  barplot(mat,
          col = cols[rownames(mat)],
          main = paste('P&Fi - F & M - Q', q, '2015'),
          names.arg = label_map[levels(means_df$class)], #x_labels_vec,
          las = 2,
          cex.main = 2,
          cex.axis = 1.5,
          cex.names = 1.5,
          ylim = c(0,3),#c(0, max(colSums(mat)) * 1.1),
          # width = width_vector_all,
          space = 0,
          beside = FALSE,
          legend.text = rownames(mat),
          args.legend = list(x = "topright", bty = "n", inset = 0.02, cex = 1, y.intersp = 0.5)
  )
}


### YEARS TEMPERATURE AVG ####

t_2014 = temperature_dataSet %>%
  filter(as.Date(date_time, format= "%d/%m/%Y %H:%M") > as.Date("31/12/2013", format= "%d/%m/%Y" ) &
           as.Date(date_time, format= "%d/%m/%Y %H:%M") < as.Date("01/01/2015", format= "%d/%m/%Y" ) )
T_avg_2014 = mean(as.numeric(gsub(',','.', t_2014$temperature) ))

t_2015 = temperature_dataSet %>%
  filter(as.Date(date_time, format= "%d/%m/%Y %H:%M") > as.Date("31/12/2014", format= "%d/%m/%Y" ) &
           as.Date(date_time, format= "%d/%m/%Y %H:%M") < as.Date("01/01/2016", format= "%d/%m/%Y" ) )
T_avg_2015 = mean(as.numeric(gsub(',','.', t_2015$temperature) ))

t_2016 = temperature_dataSet %>%
  filter(as.Date(date_time, format= "%d/%m/%Y %H:%M") > as.Date("31/12/2015", format= "%d/%m/%Y" ) &
           as.Date(date_time, format= "%d/%m/%Y %H:%M") < as.Date("01/01/2017", format= "%d/%m/%Y" ) )
T_avg_2016 = mean(as.numeric(gsub(',','.', t_2016$temperature) ))


t_2017 = temperature_dataSet %>%
  filter(as.Date(date_time, format= "%d/%m/%Y %H:%M") > as.Date("31/12/2016", format= "%d/%m/%Y" ) &
           as.Date(date_time, format= "%d/%m/%Y %H:%M") < as.Date("01/01/2018", format= "%d/%m/%Y" ) )
T_avg_2017 = mean(as.numeric(gsub(',','.', t_2017$temperature) ))


t_2018 = temperature_dataSet %>%
  filter(as.Date(date_time, format= "%d/%m/%Y %H:%M") > as.Date("31/12/2017", format= "%d/%m/%Y" ) &
           as.Date(date_time, format= "%d/%m/%Y %H:%M") < as.Date("01/01/2019", format= "%d/%m/%Y" ) )
T_avg_2018 = mean(as.numeric(gsub(',','.', t_2018$temperature) ))


#### PLOTS TO COMPARE ALL ####

#(1) L_avg over the time:

#calculate weighted average for each time step:
size_classes = c(LL, LJ, LJ2, LJ3, LJ4, LJ5, LA1, LA2, LA3)
#first join F and M in one column
#Without preasure (1)
noPreasure_noSex = test_noPreasure %>%
  mutate(J1 = J_f + J_m, J2 = J2_f + J2_m, J3 = J3_f + J3_m, J4 = J4_f + J4_m,  J5 = J5_f + J5_m,
         A1 = A1_f + A1_m) %>%
  mutate(ttlB = L + J1+ J2+ J3+ J4+ J5+ A1+ A2 + A3) %>%
  select(1,2, 15:24) %>%
  filter(as.Date(dateTime) > as.Date("31/12/2014", format= "%d/%m/%Y" )) #delete first year 'warm up' period

noPreasure_noSex = noPreasure_noSex[, c(1,2,6,7,8,9,10,11,3,4,12,5)] #rearrange order of columns

Lavg_noPreasure = noPreasure_noSex %>%
  mutate(L_avg = (L*LL + J1*LJ + J2*LJ2+ J3*LJ3+ J4*LJ4 + J5*LJ5
                  + A1*LA1 + A2*LA2+ A3*LA2)/ttlB )

#Only Predation (2)
predation_noSex = test_pred %>%
  mutate(J1 = J_f + J_m, J2 = J2_f + J2_m, J3 = J3_f + J3_m, J4 = J4_f + J4_m,  J5 = J5_f + J5_m,
         A1 = A1_f + A1_m) %>%
  mutate(ttlB = L + J1+ J2+ J3+ J4+ J5+ A1+ A2 + A3) %>%
  select(1,2, 15:24) %>%
  filter(as.Date(dateTime) > as.Date("31/12/2014", format= "%d/%m/%Y" )) #delete first year 'warm up' period

predation_noSex = predation_noSex[, c(1,2,6,7,8,9,10,11,3,4,12,5)] #rearrange order of columns

Lavg_pred = predation_noSex %>%
  mutate(L_avg = (L*LL + J1*LJ + J2*LJ2+ J3*LJ3+ J4*LJ4 + J5*LJ5
                  + A1*LA1 + A2*LA2+ A3*LA2)/ttlB )

#Only Predation (3)
fishery_noSex = test_fishery %>%
  mutate(J1 = J_f + J_m, J2 = J2_f + J2_m, J3 = J3_f + J3_m, J4 = J4_f + J4_m,  J5 = J5_f + J5_m,
         A1 = A1_f + A1_m) %>%
  mutate(ttlB = L + J1+ J2+ J3+ J4+ J5+ A1+ A2 + A3) %>%
  select(1,2, 15:24) %>%
  filter(as.Date(dateTime) > as.Date("31/12/2014", format= "%d/%m/%Y" )) #delete first year 'warm up' period

fishery_noSex = fishery_noSex[, c(1,2,6,7,8,9,10,11,3,4,12,5)] #rearrange order of columns

Lavg_fishery = fishery_noSex %>%
  mutate(L_avg = (L*LL + J1*LJ + J2*LJ2+ J3*LJ3+ J4*LJ4 + J5*LJ5
                  + A1*LA1 + A2*LA2+ A3*LA2)/ttlB )

#Both predation & fishery (4)
bothPF_noSex = test_bothPF %>%
  mutate(J1 = J_f + J_m, J2 = J2_f + J2_m, J3 = J3_f + J3_m, J4 = J4_f + J4_m,  J5 = J5_f + J5_m,
         A1 = A1_f + A1_m) %>%
  mutate(ttlB = L + J1+ J2+ J3+ J4+ J5+ A1+ A2 + A3) %>%
  select(1,2, 15:24) %>%
  filter(as.Date(dateTime) > as.Date("31/12/2014", format= "%d/%m/%Y" )) #delete first year 'warm up' period

bothPF_noSex = bothPF_noSex[, c(1,2,6,7,8,9,10,11,3,4,12,5)] #rearrange order of columns

Lavg_bothPF = bothPF_noSex %>%
  mutate(L_avg = (L*LL + J1*LJ + J2*LJ2+ J3*LJ3+ J4*LJ4 + J5*LJ5
                  + A1*LA1 + A2*LA2+ A3*LA2)/ttlB )
dev.off()

plot(noPreasure_noSex$dateTime, Lavg_noPreasure$L_avg, las = 1,  col = 'grey', cex = 0.8,
     xlab = 'years', ylab = paste('l', '[cm]'), ylim = c(1, 5) ,type = 'l', lty = 2, lwd = 3, #TeX(r"($\ell$)" )
     main = 'Changes in average size ')
lines(Lavg_pred$dateTime, Lavg_pred$L_avg, col = 'indianred3', lty=1, lwd = 2.5)
lines(Lavg_fishery$dateTime, Lavg_fishery$L_avg, col = 'skyblue4', lty=1, lwd = 2.5)
lines(Lavg_bothPF$dateTime, Lavg_bothPF$L_avg, col = 'hotpink', lty=1, lwd = 2.5)

legend("topright", legend = c('no preasure', 'only predation', 'only fishery', 'both F & P'),
       fill = c('grey', 'indianred3', 'skyblue4','hotpink'), border = NA, bty = "n", y.intersp = 0.65)



plot(noPreasure_noSex$dateTime, Lavg_noPreasure$ttlB, las = 1,  col = 'grey', cex = 0.8,
     xlab = 'years', ylab = 'biomass [?]', type = 'l', lty = 2, lwd = 3,
     main = 'Changes in biomass ', ylim = c(0, 110) )
lines(Lavg_pred$dateTime, Lavg_pred$ttlB, col = 'indianred3', lty=1, lwd = 2.5)
lines(Lavg_fishery$dateTime, Lavg_fishery$ttlB, col = 'skyblue4', lty=1, lwd = 2.5)
lines(Lavg_bothPF$dateTime, Lavg_bothPF$ttlB, col = 'hotpink', lty=1, lwd = 2.5)

legend("topright", legend = c('no preasure', 'only predation', 'only fishery', 'both F & P'),
       fill = c('grey', 'indianred3', 'skyblue4','hotpink'), border = NA, bty = "n", y.intersp = 0.6)



#### now plot L_avg vs Biomass ####

#Now group by the week number and calculate avg of the size

weekly_avg_noPreasure  <- Lavg_noPreasure %>%
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
     lwd = 3, col = 'grey', lty=1, ylim = c(5,100), xlim = c(1.5, 4.35) )
points(weekly_avg_fish$avg_L[weekly_avg_fish$year == 2015], weekly_avg_fish$avg_B[weekly_avg_fish$year == 2015], lwd = 2, col ='skyblue4' )
points(weekly_avg_pred$avg_L[weekly_avg_pred$year == 2015], weekly_avg_pred$avg_B[weekly_avg_pred$year == 2015], lwd = 2, col ='indianred3' )
points(weekly_avg_both$avg_L[weekly_avg_both$year == 2015], weekly_avg_both$avg_B[weekly_avg_pred$year == 2015], lwd = 2, col= 'hotpink' )

legend("topleft", legend = c('no preasure', 'only predation', 'only fishery', 'both F & P'),
       fill = c('grey', 'indianred3', 'skyblue4','hotpink'), border = NA, bty = "n", y.intersp = 0.7)




############


params_test_r = parameters_solv
params_test_r$general_params$Fi = 0.1


state_r = c(P = 2, E = 0.1, L= 0.0928 *1.5 ,
          BF = BF+0.1, BM = BM+0.1,
          Pred = 0.05 )

start <- Sys.time()
test_r = solver_sizeClass.v4(t = t_5years, state = state_r, parameters = noPreasure_params, temperature_dataSet = temperature_14_18)
print(Sys.time() - start)

plot(test_noPreasure$time, test_noPreasure[, 4], xlim = c(0,50), ylim = c(0,1000))
for (i in 5:11){
  lines(test_noPreasure$time, test_noPreasure[, i])
}

#plot(test_r$time,test_r[,2],col=8,type='n',ylim=c(-1,50))
#dev.off()
#for(i in 2:18)lines(test_r$time,test_r[,i],col=i)




### Check eeg produciton of V3 and V4:
dev.off()
plot(as.Date(test_noPreasure$dateTime, format= "%d/%m/%Y" ), test_noPreasure$produced_eggs, col = 'pink',
     lwd = 2, xlim = c(as.Date("01/04/2014", format= "%d/%m/%Y" ), as.Date("01/08/2014", format= "%d/%m/%Y" )),
    main = 'comparision egg proudction'  ,ylim = c(0, 1.00)) #, ylim = c(0,0.01))
lines(as.Date(test_r$dateTime, format= "%d/%m/%Y") , test_r$produced_eggs.BF1 )

###Check new eggs:

plot(as.Date(test_noPreasure$dateTime, format= "%d/%m/%Y" ), test_noPreasure$new_eggs, col = 'pink',
     lwd = 2, xlim = c(as.Date("01/04/2014", format= "%d/%m/%Y" ), as.Date("01/12/2014", format= "%d/%m/%Y" )),
     main = 'comparision new eggs',ylim = c(0, .600)) #, ylim = c(0,0.01))
lines(as.Date(test_r$dateTime, format= "%d/%m/%Y") , test_r$new_eggs )


#Check weather high egg production is related to overgrowth on fem adults

Bfem_2014.V3 = test_noPreasure %>%
  filter(as.Date(dateTime, format = '%d/%m/%Y') < as.Date("30/05/2014", format = '%d/%m/%Y')
         & as.Date(dateTime, format = '%d/%m/%Y') > as.Date("15/04/2014", format = '%d/%m/%Y')) %>%
  mutate(Bfem = A1_f + A2 + A3) %>%
  select(c(Bfem, dateTime) )

Bfem_2014.V4 = test_r %>%
  filter(as.Date(dateTime, format = '%d/%m/%Y') < as.Date("30/05/2014", format = '%d/%m/%Y')
         & as.Date(dateTime, format = '%d/%m/%Y') > as.Date("15/04/2014", format = '%d/%m/%Y')) %>%
  mutate(Bfem = BF6 + BF7 + BF8) %>%
  select(c(Bfem, dateTime) )


plot(Bfem_2014.V3$dateTime, Bfem_2014.V3$Bfem, lwd = 2, col = 'pink',
     main = 'comparision Biomass Fems V3 vs V4', ylim = c(0,50) )
lines(Bfem_2014.V4$dateTime, Bfem_2014.V4$Bfem)

test = spawning_rate_b(8.1)

#compare juvI :
BjuvI_2014.V3 = test_noPreasure %>%
  filter(as.Date(dateTime, format = '%d/%m/%Y') < as.Date("01/04/2014", format = '%d/%m/%Y')) %>%
  mutate(BjuvI = J_f + J_m) %>%
  select(c(BjuvI, dateTime, P) )

BjuvI_2014.V4 = test_r %>%
  filter(as.Date(dateTime, format = '%d/%m/%Y') < as.Date("01/04/2014", format = '%d/%m/%Y')) %>%
  mutate(BJuvI = BF1 + BM1) %>%
  select(c(BJuvI, dateTime, P) )


plot(BjuvI_2014.V3$dateTime, BjuvI_2014.V3$BjuvI, lwd = 2, col = 'pink',
     main = 'comparision Biomass JuvI V3 vs V4', ylim = c(0,10) )
lines(BjuvI_2014.V4$dateTime, BjuvI_2014.V4$BJuvI)

#Check Ingestion rate larvae:
plot(as.Date(test_noPreasure$dateTime, format= "%d/%m/%Y" ), test_noPreasure$ing_rateL, col = 'pink',
     lwd = 3, xlim = c(as.Date("01/01/2014", format= "%d/%m/%Y" ), as.Date("01/02/2014", format= "%d/%m/%Y" )),
     main = 'comparision ingestion rate'  ,ylim = c(0, 0.3190)) #, ylim = c(0,0.01))
points(as.Date(test_r$dateTime, format= "%d/%m/%Y") , test_r$ing_rateL, lwd=2 )


#check plankton
plot(as.Date(test_noPreasure$dateTime, format= "%d/%m/%Y" ), test_noPreasure$P, col = 'pink',
     lwd = 2, xlim = c(as.Date("01/01/2014", format= "%d/%m/%Y" ), as.Date("01/04/2014", format= "%d/%m/%Y" )),
     main = 'comparision Plankton Biomass'  )#,ylim = c(0, 0.0150)) #, ylim = c(0,0.01))
lines(as.Date(test_r$dateTime, format= "%d/%m/%Y") , test_r$P )

#check larvae
plot(as.Date(test_noPreasure$dateTime, format= "%d/%m/%Y" ), test_noPreasure$L, col = 'pink',
     lwd = 2, xlim = c(as.Date("01/01/2014", format= "%d/%m/%Y" ), as.Date("01/04/2014", format= "%d/%m/%Y" )),
     main = 'comparision Larvae'  )#,ylim = c(0, 0.0150)) #, ylim = c(0,0.01))
lines(as.Date(test_r$dateTime, format= "%d/%m/%Y") , test_r$L )


###
#test_r <- test_r[ , -c(1,2)] # delete timesteps, P and predator column

test_r_2015 = test_r %>%
  filter(as.Date(test_r$dateTime) > as.Date("31/12/2014", format= "%d/%m/%Y" ) &
           as.Date(test_r$dateTime) < as.Date("31/12/2016", format= "%d/%m/%Y") ) %>%
  select(E, L, BF1, BF2, BF3, BF4, BF5, BF6, BF7, BF8, BM1, BM2, BM3, BM4, BM5)

# Bar plot F & M together


# Define colors for each sex: f = light blue, m = steel blue, U = gray
# Define x-axis label order (class names)
x_labels_vec.v3 <- c("E", "L", "J", "J2", "J3", "J4", "J5", "A1", "A2", "A3")

cols <- c("f" = "#8856a7", "m" = "#add8e6")
bg_col <- "grey90"  # light grey for background

# Define x-axis label order (class names)
x_labels_vec <- c( "BF1", "BF2", "BF3", "BF4", "BF5", "BF6", "BF7", "BF8")
label_map <- c(
  "BF1" = "1-2",
  "BF2" = "2-3",
  "BF3" = "3-4",
  "BF4" = "4-5",
  "BF5" = "5-6",
  "BF6" = "6-7",
  "BF7" = "7-8",
  "BF8" = "8.5"
)
#width_vector_all <- rep(1, length(x_labels_vec))  # Optional bar widths



# Define size classes to include and their order
class_labels <- c("BF1", "BF2", "BF3", "BF4", "BF5", "BF6", "BF7", "BF8")
label_map <- setNames(class_labels, class_labels)

# Set up 2x2 plot grid for 4 quarters
par(mfrow = c(2, 2))

for (q in 1:4) {
  # Indices for quarter (912 time steps per quarter)
  init_idx <- 912 * (q - 1) + 1
  final_idx <- 912 * q
  #prepare test_r V4:
  quarter_data <- test_r_2015[init_idx:final_idx, ]

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

  #PREPARE test V3
  quarter_data <- M_F_2015[init_idx:final_idx, ]

  # 1. Get sexed class columns (e.g., J_f, J2_m, A1_f, etc.)
  sexed_df <- quarter_data %>%
    select(matches("^(J|A1)\\d?_f$|^(J|A1)\\d?_m$"))

  fem_df2 <- quarter_data %>%
    select( A2, A3)

  # 3. Compute column means
  sexed_means <- colMeans(sexed_df)
  fem_means <- colMeans(fem_df2)

  # 4. Create long-format data frame for sexed
  sexed_df_long <- data.frame(
    class_sex = names(sexed_means),
    value = as.numeric(sexed_means)
  ) %>%
    separate(class_sex, into = c("class", "sex"), sep = "_")

  # 5. Create long-format data frame for unsexed, with sex = "U"
  unsexed1_df_long <- data.frame(
    class = names(fem_means),
    value = as.numeric(fem_means),
    sex = "f"
  )

  # 6. Combine both
  means_df <- bind_rows(sexed_df_long, unsexed1_df_long)

  # 7. Set class as factor to preserve desired order
  means_df$class <- factor(means_df$class, levels = x_labels_vec.v3)

  means_df$label <- label_map[as.character(means_df$class)]

  # 8. Pivot to wide format for barplot
  means_mat <- means_df %>%
    pivot_wider(names_from = sex, values_from = value, values_fill = 0 ) %>%
    arrange(class) %>%
    select(any_of(c("f", "m")))  # Ensure consistent column order

  # 9. Transpose to matrix for barplot (rows = sex, columns = classes)
  matv3 <- t(as.matrix(means_mat))


  # Plot barplot V4 (background because B is higher)
  barplot(
    mat,
    col = cols[rownames(mat)],
    main = paste("Quarter", q, "– 2015"),
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

  # Draw barplot V3
  barplot(matv3,
          col = rgb(0.9, 0.9, 0.9, alpha = 0.75), #col gray a bit transparent#cols[rownames(matv3)],
          main = paste('P&Fi - F & M - Q', q, '2015'),
          names.arg = label_map[levels(means_df$class)], #x_labels_vec,
          las = 2,
          #cex.main = 2,
          #cex.axis = 1.5,
          #cex.names = 1.5,
          #ylim = c(0,3),#c(0, max(colSums(mat)) * 1.1),
          # width = width_vector_all,
          space = 0,
          add = TRUE,
          border = NA,
          axisname = FALSE,
          beside = FALSE,
          #legend.text = rownames(mat),
         # args.legend = list(x = "topright", bty = "n", inset = 0.02, cex = 1, y.intersp = 0.5)
  )
}



