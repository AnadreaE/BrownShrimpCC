#Here I will check plausibility of individual funciton that are used in the model

#### Some params ####
temp_range = seq(1,30, 0.1)

the.method = 'rk4'
h = 0.2 # half saturatuion constant original aprox .5
L_inf = 8.5
epsilon = 0.22
m = 3
const_c = 0.01

#### CHECK K FUNC ####

#SEE K-FUNC GAUSS VS TPC
plot(temp_range, sapply(temp_range, K_func), type = "b", col = "blue",
     xlab = "Temperature (°C)", ylab = "K(T)",
     main = "K func Gauss vs Briere")
lines(temp_range, sapply(temp_range, K_func_briere), type = "b", col = 'red' , cex = 0.5)
legend("topright", legend = c("Gauss", "Briere"), col = c("blue", "red"), lty = 1, pch = 1)


#### CHECK NATURAL MORTALITY ####

natural_mort_br <- sapply(temp_range, respiration_rate_b, 3)
#dev.off()
plot(temp_range, respiration_rate(temp_range, 3) , type = "b", col = "blue",
     xlab = "Temperature (°C)", ylab = " ",
     main = "natural mortality Gauss vs Briere")
lines(temp_range, natural_mort_br, type = "b", col = 'red' , cex = 0.5)
legend("topright", legend = c("Gauss", "Briere"), col = c("blue", "red"), lty = 1, pch = 1)


#Checking when funciton returns NaN or inf or similat invalid vals
respiration_rate_b_check = function(temperature,L){
  mu = m*convertL_to_W(L)*K_func(temperature)# this is the rigth term of vB eq.
  if (L>5) {
    mu =  convertL_to_W(L) *K_func_briere(temperature) *( 1 - molting_fraction(L, temperature)) #tbc (1-molting_fraction) all non molting fems still have a natural mortality
  }
  return(mu)
}


#### bellow lines where the try to invert nat. mortality curve, but this was proven wrong approach ####
natural_mortality_test = function(temperature,L){
  m=3
  K_range = 0.007613879*m - 0.0005557642*m #max_K - min_K #max(sapply(temp_range, K_func))
  K_inv = (0.007613879*m - K_func(temperature)*m) #/ K_range
  mu = convertL_to_W(L) * K_inv  # this is the rigth term of vB eq.
  if (L>5) mu = 0 #tbc
  return(mu)
}

sizes <- c(1.8, 2.8, 3.8, 4.8) #cm

par(mfrow = c(2,2))
for (i in sizes){
  plot(temp_range, natural_mortality_test(temp_range, i) , type = "b", col = "blue",
       xlab = "Temperature (°C)", ylab = " ",
       main = paste("L = ", i, " cm mortality with K_func inverese \n mu = convertL_to_W(L)* K_inv"),
       cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5)#, ylim = c(0, 0.6) )
  lines(temp_range, sapply(temp_range, natural_mortality_test_b, i) , type = "b", col = "red", cex = 0.8)
}

#######

K_func = function(temperature){
  #following param. vals have been estimated and docum ented in Thesis Appx. A
  a = 0.008
  W = 7 # Twidth
  T_opt = 16 # optimal temperature
  return (a*exp(-(temperature-T_opt)/(2*W)))
}

K_func_briere_t = function(temperature){
  r_max_f = 0.007638409
  T_min = 0.4
  T_max = 30
  alpha = 0.541524
  beta = 0.276430
  diff_min = temperature - T_min
  diff_max = T_max - temperature
  diff = T_max - T_min
  alpha_invert = 1 - alpha
  toReturn = r_max_f*( ((diff_min/alpha)^alpha) * ((diff_max/ alpha_invert )^alpha_invert) * (1 / diff)  )^(alpha*alpha_invert/(beta^2) )

  return(toReturn)
}

#### CHECK FUNCTIONAL RESPONSE ####
#**PENDING !!!




