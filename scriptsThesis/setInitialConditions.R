#DETERMINE INITAL CONDITIONS BASED ON REFERENCES FOUND IN LITERATURE
# Spring 2015

### 1) ovigerous fems : 40 individuals / 1000 m2
avg_l_OF = (5.0+8.5)/2
avg_w_OF = convertL_to_W(avg_l_OF)

OF_biomass = 40*avg_w_OF #[gr / 1000 m2]

#with avg T from spring 2015, we can estimate the share of OF from the total F population
#which is OF_biomass, thus we can estimate total F population

temp_data =  read.csv("./data/Temperature_data_1982-2018_MarBiol.csv", header = TRUE, sep = ';')
T_spring2015 = temperature_func(temp_data, "20/03/2015", "21/06/2015")

avg_T = mean(T_spring2015$temperature)

share = molting_fraction(avg_l_OF, avg_T) # = 0.0226

#0.0226 of all F is equivalent to OF_biomass, then total F :

fem_biomass = 100*OF_biomass/2.26
#-> fem biomass = 5443.3 gr / m2


### 2) Larve 145 individuals / 1000 m2
avg_l_larv = 0.4 # cm
avg_w_larv = convertL_to_W(0.4)

larv_biomass = 145*avg_w_larv #[gr/1000m2]




