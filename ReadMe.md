<!--
SPDX-FileCopyrightText: 2025 Helmholtz-Zentrum hereon GmbH
SPDX-License-Identifier: CC0-1.0
SPDX-FileContributor Andrea Farfan <farfanqbb@gmail.de>
-->

# BrownShrimpCC

Brown shrimp size class model  

## Description 

The BrownShrimpCC is a package that simulates the growth and development of the Brown Shrimp (Crangon crangon). The model used for brown shrimp is a size-spectrum model coupled with a stage-based component that describes the development of egg and larval stages. 
With this model, we aim to investigate how predation and fisheries affect the size structure and biomass of the brown shrimp in the North Sea.  


## Requirements 

To install and run this program, you just need to have R and RStudio installed in your computer. 
All libraries needed will be installed when downloading the project.  


## Installation 

You can either download the project from GitHub [BrownShrimpCC] (https://github.com/AnadreaE/BrownShrimpCC.git) by running following lines in the console: 

install.packages("devtools")  # install devtools if necessary
devtools::install_github("AnadreaE/BrownShrimpCC")
library(BrownShrimpCC)


Important:
Please install the package using the command above (or via RStudio → Build → Install).
If the project folder is downloaded manually and used without installation, dependencies will not install automatically.
Once installed, go to the build tab in the top right window of the RStudio GUI, click on `More`, `clean and install` to ensure that all functions and data is available. Then it shall be ready to use.   

To install [BrownShrimpCC] from the USB stick that was handed over with the Thesis you may have to install all required libraries manually:

deSolve,
    RColorBrewer,
    viridis,
    dplyr,
    tidyr,
    colorspace,
    latex2exp,
    stringr,
    lubridate

and additionally, it is recommended to install "devtools".


## License 

Apache-2.0

## Usage 

Structure: 
Folder R: contains files where all the functions that we need to run the simulations are defined. 
Folder scriptsThesis: contains the files with the steps followed to get to the functions (refer to chapter `Model development` from the Thesis work related to this project).
Folder Simulations: contains the files where we apply the functions and run simulations of different scenarios to address the research question. 
 
As an inspiration, I recommend to take a look to the file `Sims.R` under the folder `Simulations`, where some simulations of different scenarios are carried out. All simulations and plots included in the thesis are found there. 
But in case you want to start from scratch, you first need the following: 
-	Temperature data: The package includes temperature data from the Federal Waterways and Shipping Administration (original name in German: Wasserstraße und Schifffahrtsverwaltung des Bundes) from the station `LEUCHTTURM ALTE WESER' (georeference ETRS89/UTM32N). This file named `temperature_10ger.csv` can be imported from the folder `data` and contains the daily averaged temperature from 01.01.2010 to 01.11.2022.  
Another temperature data set from the Wadden Sea can be found in the same folder, which contains temperatures measured at 8:00 every day from 01.01.1982 to 31.12.2018. 
Important: the data set containing temperature must have two columns named `date_time` and `temperature`, otherwise the solver and functions nested in the solver will not work. 
-	Initial conditions: a vector with the initial values of the state variables of the system is required. The state variables represent the abundance in biomass per area, kg km-2, and must be named: `E` for egg biomass, `L` for larva biomass,  and for the size classes we differentiate males M from Females. The size class names must begin with `BM` and `BF` for males and females, respectively, and be followed by the size class number, depending on how many size classes are defined. By default, the system and solver use 8 size classes for females and 5 for males, with a size resolution or size width of 1 cm, but these can be customized. NB: In case the size resolution is changed, the parameter `size_width` in the function `shift_next_sizeClass()` must be accordingly amended.
-	List of parameters: There already exists a list of default parameters `parameters_solv`, in case any parameter has to be amended for any specific scenario, it can be easily copied, and then the changes in the copy can be used for the simulation.  
A simulation can be computed by only running the `solver_sizeClass.v5`  function: 
test <- solver_sizeClass.v5(t = time_steps, state = initi_conditions, parameters = parameters_solv, temperature_dataSet = temperature_ds)
This function includes all growth, reproduction, predation, and fishery functions developed in the Master's thesis `Changes in body size and biomass of brown shrimp (Crangon crangon) in response to predation and fishery in the Southern North Sea`. 


## Contribution Guidelines 

To contribute to BrownShrimpCC, Clone this repo locally and commit your code on a separate branch.  


## Support
Helmholz-Zentrum Hereon
 
## Authors and acknowledgment

Helmholz-Zentrum Hereon
Andrea Farfan Aragon

## Project status  
Under development

