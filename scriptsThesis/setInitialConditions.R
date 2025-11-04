library(BrownShrimp)
library(deSolve)
library(RColorBrewer)

#Read water temperature from Germany:

temperature_germany <- read.csv("./data/temperature_10ger.csv")

tempG.10_15 = temperature_germany %>% filter(date_time > as.Date("31/12/2009", format = "%d/%m/%Y") & date_time < as.Date("01/01/2016", format = "%d/%m/%Y"))

ble_data = read.csv("./data/BLE_Inlandslandungen_SpeiseKrabbe.csv", sep = ';')

t_years = seq(0,length(tempG.10_15$temperature)-.1)


#THE INITIAL CONDITIONS WERE DETERMINED AS FOLLOWS

# (1) Run a simulation of scenario 'Only predation'

pred_params = parameters_solv
pred_params$general_params$Fi = 0 #No fishery
pred_params$general_params$Imax_ik = 0.185


bothPF_params = parameters_solv
bothPF_params$general_params$Imax_ik = 0.185
bothPF_params$general_params$Fi = 4/365 #Temming & Hufnagl 2015

state_st = c(P = 2, E = 0.1, L= 0.0928 *1.5 ,
             BF = BF+0.1, BM = BM+0.1)

start <- Sys.time()
test_predation <- solver_sizeClass.v5(t = t_years, state = state_st, parameters = pred_params, temperature_dataSet = tempG.10_15)
print(Sys.time() - start)


##### plot only January for estimation of initial conditions in Thesis#####


sol_15_jan = test_predation %>%
  filter(as.Date(test_predation$dateTime) > as.Date("31/12/2014", format= "%d/%m/%Y" ) &
           as.Date(test_predation$dateTime) < as.Date("01/02/2015", format= "%d/%m/%Y") ) %>%
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
  mat,
  col = cols[rownames(mat)],
  main = 'L50 = 3.69 cm',#paste("Size distribution in January"),
  names.arg =  as.character(size_mean_F), #label_map,
  ylab = 'Biomass kg km^-2',
  xlab = 'mean L of size class',
  las = 1,
  cex.main = 1.2,
  cex.axis = 1.2,
  cex.names = 1.2,
  cex.lab = 1.5,
  #ylim = c(0, .4),#c(0,25),
  space = 0,
  beside = FALSE,
  #legend.text = rownames(mat_prop),#rownames(mat),
  #args.legend = list(x = "topleft", bty = "n", inset = 0.02, cex = 1.5)
)

#Overlap selective fishery curve (probit):
par(new = TRUE)

# Match kernel curve to barplot's x-axis midpoints
# Assume size_mean_F matches the size classes
selective_fi <- sel_probit(size_mean_F, L50 = 3.69, SR = 0.75)

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
axis(side = 4, cex.axis = 1.2)

#### @Thesis, Methodology Biomass estimation for initial conditions ####
#Calculate distribution of 14.71  kg / km2:

size_classes_proportions = colSums(mat_prop)

size_dist_biomass = 14.71*size_classes_proportions #kg/km-2

est_biomass_dist = size_dist_biomass*(1-(4/365))/(4/365) #considering that fishery mosrtality equals 4y-1

F_proportions =mat["f", ] / colSums(mat)
M_proportions =mat["m", ] / colSums(mat)

est_biomass_dist_F = est_biomass_dist*F_proportions #init condition biomass Fems
est_biomass_dist_M = est_biomass_dist*M_proportions #init condition biomass Masc




