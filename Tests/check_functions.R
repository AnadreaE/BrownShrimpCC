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
lines(temp_range, sapply(temp_range, K_func_briere, sex_params = parameters_solv$Fem_params), type = "b", col = 'red' , cex = 0.5)
legend("topright", legend = c("Gauss", "Briere"), col = c("blue", "red"), lty = 1, pch = 1)


#### CHECK NATURAL MORTALITY ####

natural_mort_br <- sapply(temp_range, respiration_rate_b, L=3, sex_params = parameters_solv$Fem_params)
#dev.off()
plot(temp_range, respiration_rate(temp_range, 3) , type = "b", col = "blue",
     xlab = "Temperature (°C)", ylab = " ",
     main = "natural mortality Gauss vs Briere")
lines(temp_range, natural_mort_br, type = "b", col = 'red' , cex = 0.5)
legend("topright", legend = c("Gauss", "Briere"), col = c("blue", "red"), lty = 1, pch = 1)



#### CHECK FUNCTIONAL RESPONSE ####

func_response = function(P){
  alpha_ir = 0.1721763#8.238
  return(P/(P+alpha_ir))
}

ingestion_rate_b_test = function(temperature, L, P, sex_p){
  c_div = const_c / convertL_to_W(L)
  fr = func_response(P)
  #params = sex_parameters_func(sex)
  #L_infty = params$L_inf
  L_infty = sex_p$L_inf
  #result = m*K_func_briere(temperature, sex)*convertL_to_W(L)*L_infty*(c_div)^(1/m)*( P/(P+alpha_ir) )
  w = convertL_to_W(L)
  result = m*K_func_briere(temperature, sex_p)*(L_infty/L)*fr
  #m*convertL_to_W(L)*K_func(T)*L_inf/L*P/(P+h) #revised formula Andrea 04.04.25
  return (result)
}


const_food = c(1,2, 3, 5, 15, 30)
time_range = seq(1, 365)

food_3years = new_food(time_range)

B_population = 50 #biomass
Temp = 10
#consumption = B_population*ingestion_rate_b_test(Temp, 4, const_food, parameters_solv$Fem_params)

#plot(time_range, consumption, main = paste("T =", Temp, "B =", B_population, "constant food = ", const_food) )

par(mfrow = c(2,3))

for (i in const_food){
  const_food_loop = rep(i, times = 365 )
  consumption = B_population*ingestion_rate_b_test(Temp, 4, const_food_loop, parameters_solv$Fem_params)
  plot(time_range, consumption, main = paste("T =", Temp, "B =", B_population, "constant food = ", i), ylim = c(1,2) )

}
dev.off()
food_increasing = seq(1,500)

plot(food_increasing, B_population*ingestion_rate_b_test(Temp, 4, food_increasing, parameters_solv$Fem_params))

#now check new updated instion rate function
par(mfrow = c(2,3))

for (i in const_food){
  const_food_loop = rep(i, times = 365 )
  consumption = B_population*ingestion_rate_b(Temp, 4, const_food_loop, parameters_solv$general_params, parameters_solv$Fem_params)
  plot(time_range, consumption, main = paste("T =", Temp, "B =", B_population, "constant food = ", i), ylim = c(1,2) )

}




#### TEST SELECTIVE FISHERY ####

fishery = function(l){
  L50 = 4.49 #[cm]
  SR = 1.56
  nominator = exp( (1.349/SR) * (l - L50))
  denominator = 1 + nominator
  return(nominator / denominator)
}

l_range = seq(0.6, 8.5, 0.1)


plot(l_range, fishery(l_range))





