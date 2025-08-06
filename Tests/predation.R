library(latex2exp)

#Here Predation function will be determined

#Optimal length of prey for a predator of size l_pred
lopt = function(l_pred){
  l = log(l_pred)
  r=-1.65
  gamma = 0.011
  return(r + l - gamma*l^2)
}


ingestion_kernel = function(I_max, l_pred, l_prey){
  l_opt = lopt(l_pred)
  return(I_max*exp(-3/2*(log(l_prey)-l_opt)^2))
}


cod_lRange = seq(0.5, 120, 0.1) #[cm]
cod_larvRange = seq(2.5, 3, 0.01) #[cm]
shrimp_lRange = seq(0.155, 8.5, 0.006981605) #(8.5-0.15)/1176 steps to meet same length as cod L range


#dev.off()
plot(cod_lRange, exp(lopt(cod_lRange)), main = 'Optimal prey lenght', xlab= 'predator L [cm]',
     ylab = 'prey (shrimp) L [cm]' , col = 'gray41', las = 1, cex.axis=1.5, cex.lab=1.5, cex.main = 1.4,
     xlim = c(0,60), ylim = c(0,10), type = "n")#, cex.main=1.3, xlim = c(1,4), ylim = c(0,2))
rect(xleft = 0, ybottom = 0, xright = 62, ytop = 1,
     col = 'skyblue3', border = NA) #larvae and juv1 #color Offshore
rect(xleft = 0, ybottom = 1, xright = 62, ytop = 4,
     col = 'lightskyblue1', border = NA) # zoea 4 and 5 (juv2,3 and 4) #color intertidal
rect(xleft = 0, ybottom = 4, xright = 62, ytop = 6,
     col = 'papayawhip', border = NA) #juv V #color bottom, sea floor
rect(xleft = 0, ybottom = 6, xright = 62, ytop = 8.7,
     col = 'lightskyblue1', border = NA) #Adult #color bottom, sea floor
abline(h=0.2, lty= 2, col = 'lightsalmon', lwd=2)
abline(h=8.5, lty=2, col='tomato3', lwd = 2)
points(cod_lRange, exp(lopt(cod_lRange)), col = 'gray41', lwd = 1.5)
legend("topleft", c(TeX("$\\L_{opt}$"), 'min L shrimp', 'max L shrimp') ,col=c('gray41','lightsalmon', 'tomato3'),
       lty = c(1,2,2), lwd=3)
legend("bottomright", c('open sea >10m', 'subtidal>1m', 'sea bottom<1m') ,fill=c('skyblue3','lightskyblue1', 'papayawhip'))



plot(log(cod_lRange), log(shrimp_lRange), main = 'Log-log size range predator Cod vs prey',
     xlab= 'log Cod L', ylab = 'log shrimp l')

fit = lm(log(shrimp_lRange) ~ log(cod_lRange))
summary(fit)
dev.off()

Imax = 1
plot(shrimp_lRange, ingestion_kernel(I_max = Imax, 8.75,shrimp_lRange), main = paste("Imax=", Imax), #l_pred =cod_lRange[20]
     xlab= 'Shrimp size range [cm]', ylab = 'ingestion kernel')#, xlim= c(0,20))


plot(cod_larvRange, ingestion_kernel(I_max = Imax, cod_larvRange), main = paste("Imax=", Imax),
     xlab= 'Cod size range [cm]', ylab = 'ingestion kernel')#, xlim= c(0,20))



#See with one single L_opt for the avg size of larvae cod:
cod_avg = (2.5+3)/2
l_opt_cod_avg = lopt(cod_avg)

#first see range of l_opt for larval cods:
range_lopt_larvCod = exp(lopt(cod_larvRange))
min(range_lopt_larvCod)
max(range_lopt_larvCod)

shrimpLarv_range = seq(min(range_lopt_larvCod), max(range_lopt_larvCod), length.out = length(cod_larvRange))

ingestion_kernel_c = function(I_max, l_prey){
  l_opt = -0.6496558
  return(I_max*exp(-3/2*(log(l_prey)-l_opt)^2))
}

Imax=1
plot(cod_larvRange, ingestion_kernel_c(I_max = Imax, shrimpLarv_range), main = paste("Imax=", Imax),
     xlab= 'Cod size range [cm]', ylab = 'ingestion kernel')



#Cod reproduciton 'peak' in spring :
Bcod = function(t) {
  toReturn = 0.2 * (1.05 + cos( (2*pi/ 365)*(t-105) ) )
  return(toReturn)
}

t = seq(1,365*5)
plot(t, Bcod(t))


#test l - l_opt against Ing. kernel
Imax = 10
plot(log(shrimp_lRange) -  l_opt_cod_avg, ingestion_kernel(I_max = Imax, shrimp_lRange), main = paste("Imax=", Imax),
     xlab= 'l - l_opt', ylab = 'clearance efficiency')


#######PREDATIO WITH: #############
#spatio-temporal distribution of biomass specific mortality

#mu_shrimp = mu_ref * f_t * sigma * (gamma*B_s + beta*n(J) )
#mu_ref: reference value = 0.025 day-1
#f_t: Q10 temperature coefficient
#sigma: function that depends on depth and salinity.-> will be a constant
#gamma*B_s: shrimp concentration
#beta*n(J): oscilation component for seasonality.

eta_J <- function(t_day) {
  cos_term <- cos(2 * pi * (t_day - 150) / 365)
  eta <- 0.5 + 0.5*cos_term  # gives values from 0 to 2
  return(eta^2)         # seasonal amplification
}

ST_predation = function(t_day, Te, B_s){
  mu_ref = 0.025 #day^-1 #same as mesozooplankton paper
  f_t = 2
  sigma = 0.5
  eta = eta_J(t_day)
  beta = 18 # mesozooplankton paper = 18
  gamma = 3.7469 #m2 (Kg d)^-1 #original 0.1 m2 (molC d)^-1
  return(mu_ref*f_t*sigma*(gamma*B_s + beta*eta ))
}

#test
days_3years = seq(1, 365*3)

B_s = 2 + cos(2 * pi * (days_3years - 250) / 365)

plot(days_3years, B_s)

predatio_preasure = ST_predation(days_3years, 10, B_s)
plot(days_3years, predatio_preasure)


#See how new_predation function look like:
plot(days_3years, new_Bpredator(days_3years), main = 'see preador funciton')


