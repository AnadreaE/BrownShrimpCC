
#Here Predation function will be determined

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


cod_lRange = seq(2.5, 120, 0.1) #[cm]
cod_larvRange = seq(2.5, 3, 0.01) #[cm]
shrimp_lRange = seq(0.155, 8.5, 0.0071) #(8.5-0.15)/1176 steps to meet same length as cod L range



plot(cod_lRange, exp(lopt(cod_lRange)), main = 'Optimal prey lenght of cod', xlab= 'Cod L [cm]',
     ylab = 'prey l [cm]' , col = 'gray41', las = 1, cex.axis=1.5, cex.lab=1.5)#, cex.main=1.3, xlim = c(1,4), ylim = c(0,2))
abline(h=0.2, lty= 2, col = 'lightsalmon', lwd=2)
abline(h=8.5, lty=2, col='tomato3', lwd = 2)
legend("topleft", c('l_opt', 'min l BS', 'max l BS') ,col=c('gray41','lightsalmon', 'tomato3'),
       lty = c(1,2,2), lwd=2)




plot(log(cod_lRange), log(shrimp_lRange), main = 'Log-log size range predator Cod vs prey',
     xlab= 'log Cod L', ylab = 'log shrimp l')

fit = lm(log(shrimp_lRange) ~ log(cod_lRange))
summary(fit)
dev.off()

Imax = 1
plot(shrimp_lRange, ingestion_kernel(I_max = Imax, cod_lRange[20],shrimp_lRange), main = paste("Imax=", Imax),
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
