# rAedesSim
R package for population mosquito modeling Progetto REDLAV IBIMET CNR ASL2 Lucca 


To install package 

```R
if (!require(devtools)) { install.packages("devtools")}
devtools::install_github("alfcrisci/rAedesSim")
```

An example of operative work-chain 

```R
library(rAedesSim)

###########################################################################################################################
# Load some meteorological data obtained by Weather Local Model simulations for 2012 in REDLAV coastal locations in Tuscany.

data("redlav_2012_meteo")
data("redlav_2012_monitoring")

#################################################################################################
# Load different weather based watermodels for two different mosquito habitat 
# trap indicates a standard recipient 
# tombino indicates a urban manhole
# watermodel objects are a lm/gma models  with formula = tmedwater ~ daylength + tmed (+ tmax + tmin). 
# Length of day and daily mean temperature are required.


data(trap_tosc_wmodel)
data(tombino_tosc_wmodel)

#################################################################################################
# Create biocontainer objects and define a single BS (breeding site) number ( nrecipients=1 ).
# loanding  different model to assess water temperature from weather data.

i_biocontainer_tomb=biocontainer(nrecipients=1,
                                 watermodel=tombino_tosc_wmodel, 
                                 lat=42.76090556,
                                 lon=10.88788889,
                                 elevation=5)
                                 
i_biocontainer_trap=biocontainer(nrecipients=1,
                                 watermodel=trap_tosc_wmodel,
                                 lat=42.76090556,
                                 lon=10.88788889,
                                 elevation=5)

#################################################################################################
# To perform simulation and to fit the  best parameters for a selected location.
# Castiglione_della_Pescaia

# Retrieve biodata and meteodata objects respectively obtained from local mosquito eggs monitoring 
# and meteorological simulated weather data  

C_della_Pescaia_P4_monitoring=redlav_2012_monitoring[[4]]

# monitoring egg data are weekly and dividing the ones by the number of days of monitoring interval 
# a daily egg laying activity is obtained. It is useful to find best model's fitting parameters.

C_della_Pescaia_P4_monitoring=C_della_Pescaia_P4_monitoring$ts_data/7


C_della_Pescaia_P4_meteo_2012=redlav_2012_meteo[[4]]


##################################################################################################
# How to perform a simulation. 

iniegg=50 # winter diapause initial eggs
maxjd=140 # maximal jd of diapausant eggs exclosion.
varianceday=15 # temporal variance of diapausant eggs exclosion courbe  
inistep=15 # first steps used for parameter fitting.
starting_day_simulation=61 # 1 march for italy

#######################################################################################################################
# The information to assess local environmental parameters linked to the water are contained in biocontainer object.
# Hence biometeo objects merge the information of geograpical/topographical location features and meteorological data obtained 
# from field monitoring or simulated by using weather/environmental models.


C_Pescaia_P4_bio_tombino=biometeo(C_della_Pescaia_P4_meteo_2012,
                                   startday=starting_day_simulation,
				   i_biocontainer_tomb)

C_Pescaia_P4_bio_trap=biometeo(C_della_Pescaia_P4_meteo_2012,
                                startday=starting_day_simulation,
				i_biocontainer_trap)

ini_population=biopopulation(eggs=0,larvae=0,pupae=0,adults=0,eggs_diap=iniegg)

guess_simulation_tomb=check_fit_initial(C_Pescaia_P4_bio_tombino,
                                        C_della_Pescaia_P4_monitoring,
                                        i_biocontainer_tomb,
                                        larv_alpha=0.5,
                                        adu_alpha=5,
                                        eggs_diap=iniegg,
                                        egn=95,
                                        ini_rmse = 1,
                                        end_rmse=15
                                        )


guess_simulation_tomb$resfig_all


simulation=biomodel(i_biometeo=C_Pescaia_P4_bio_tombino,
                    i_biocontainer=i_biocontainer_tomb,                            
                    i_biopopulation=ini_population,
                    i_bioparameters= bioparameters(alfa_l=0.5,alfa_a=5)
                    )


					  


##################################################################################################
# viewwhere is a function to perform a fast visualisation of simulation in its urban context.
# work of rAedesSim object 

viewwhere(simulation)


simulation_fit=biofitmodel(i_biometeo=i_biometeo,
                           i_biopopulation=i_biopopulation,
                           i_biocontainer=i_biocontainer,
                           i_monitoring=Castiglione_della_Pescaia_P4_monitoring,
                           range_alpha_a=c(0,seq(0,5,1)),
                           range_alpha_l=seq(0.6,1.6,0.2),
                           plotresults=TRUE
			   )	


simulation_fit

##################################################################################################
```
