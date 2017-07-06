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

# install the packages from CRAN if is necessary.

if (!require(dygraphs) || !require(XLConnect) ) { install.packages(c("dygraphs","XLConnect")) }

library(dygraph) # htmlwidget plots
library(XLConnect) # excel writing

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

# monitoring egg data in ts_data sub_object are weekly and dividing the ones by the number of days of monitoring interval 
# a daily egg laying activity is obtained. It is useful to find best model's fitting parameters.

C_della_Pescaia_P4_monitoring$ts_data=C_della_Pescaia_P4_monitoring$ts_data/7


C_della_Pescaia_P4_meteo_2012=redlav_2012_meteo[[4]]


##################################################################################################
# How to perform a simulation. 

iniegg=400 # winter diapause initial eggs
maxjd=140 # maximal jd of diapausant eggs exclosion.
varianceday=15 # temporal variance of diapausant eggs exclosion courbe  
inistep=15 # first steps used for parameter fitting.
mean_egg=95 # average eggs given by a female during 1 gonotrophyc cycle
starting_day_simulation=61 # 1 march for italy


ini_population=biopopulation(eggs=0,larvae=0,pupae=0,adults=0,eggs_diap=iniegg)


#######################################################################################################################
# The information to assess local environmental parameters linked to the water are contained in biocontainer object.
# Hence biometeo objects merge the information of geograpical/topographical location features and meteorological 
# data were obtained from field monitoring orsimulated by using a weather/environmental local area model.


C_Pescaia_P4_bio_tombino=biometeo(C_della_Pescaia_P4_meteo_2012,
                                  i_biocontainer_tomb,
                                  startday = starting_day_simulation,
                                  datemax=maxjd,
                                  varjd=varianceday)

C_Pescaia_P4_bio_trap=biometeo(C_della_Pescaia_P4_meteo_2012,
                               i_biocontainer_trap,
                               startday = starting_day_simulation,
                               datemax=maxjd,
                               varjd=varianceday)


#########################################################################################################
# The context of simulation considered for Castiglione della Pescaia :

# season startday: 1 march jday 61
# spring time  period (jday) : 61 to [61 +140]
# autumn time  period (jday) : [365-130] to 365 
# thermal bias 3 celsius degrees for air and water.
# first 15 observative step of monitoring
# mean egg for female:  95 eggs

guess_simulation_tomb=check_fit_initial(C_Pescaia_P4_bio_tombino,
                                        C_della_Pescaia_P4_monitoring,
                                        i_biocontainer_tomb,
                                        larv_alpha=0.8,
                                        adu_alpha=2,
                                        eggs_diap=iniegg,
                                        egn=85,
                                        Ndayscenario_prim=150,
                                        Ndayscenario_aut=130,
                                        ini_rmse = 1,
                                        end_rmse=15,
                                        deltatmed_prim=3,
                                        deltawater_prim=3,
                                        deltatmed_aut=3,
                                        deltawater_aut=3
)



guess_simulation_tomb$resfig_all   # results 


############################################################################################################
# Write results

datasim=data.frame(date=rownames(as.data.frame(guess_simulation_tomb$simulation$ts_population)),
                   as.data.frame(guess_simulation_tomb$simulation$ts_population),
                   as.data.frame(guess_simulation_tomb$simulation$ts_parameter),
                   as.data.frame(i_biometeo_tomb$timeseries)
                   )

obs=data.frame(date=rownames(as.data.frame(C_della_Pescaia_P4_monitoring$ts_data[,1])),
           observed=C_della_Pescaia_P4_monitoring$ts_data[,1],
	   row.names = NULL)

datasim=merge(datasim,obs)

name_out_excel="C_della_Pescaia_P4_simulation.xls"

if (file.exists(name_out_excel)) {file.remove(name_out_excel)}

XLConnect::writeWorksheetToFile(name_out_excel,datasim,"C_della_Pescaia_P4")

#################################################################################################
# Plot paramter courbes

guess_simulation_tomb$simulation$ts_parameter$index_day=NULL
guess_simulation_tomb$simulation$ts_parameter$d_emergency=NULL
guess_simulation_tomb$simulation$ts_parameter$d_induction=NULL
guess_simulation_tomb$simulation$ts_parameter$inib_state=NULL

guess_simulation_tomb$simulation$ts_parameter$mu=guess_simulation_tomb$simulation$ts_parameter$mu/100 
# reduce eggs mortality to 0- 1 range

#################################################################################################
#  recursive grid search of alpha adults and alpha larvae parameters

dygraph(guess_simulation_tomb$simulation$ts_parameter, xlab = "Year", ylab = "Rate")

dygraph(C_Pescaia_P4_bio_tombino$timeseries[,c("emergency","d_induction")], xlab = "Year", ylab = "Rate")

#################################################################################################
#  recursive grid search of alpha adults and alpha larvae parameters

simulation_fit=biofitmodel(i_biometeo=i_biometeo,
                           i_biopopulation=i_biopopulation,
                           i_biocontainer=i_biocontainer_tomb,
                           i_monitoring=C_della_Pescaia_P4_monitoring,
                           range_alpha_a=c(0,seq(1,2.5,0.2)),
                           range_alpha_l=seq(0.8,1.2,0.2),
                           plotresults=TRUE
			   )	


simulation_fit

##################################################################################################
# viewwhere is a function to perform a fast visualisation of simulation in its urban context.
# work of rAedesSim object 

viewwhere(i_biocontainer_tomb)



```
