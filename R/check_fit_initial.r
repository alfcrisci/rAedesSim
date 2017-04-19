#' check_fit_initial
#'
#' @description check_fit_initial is a function to identificate in a unknow location 
#' which are the best initial thermal condition partial to fit better the start of season and in spring and autmn the thermal bias due to  air and water temperature.
#' @param  i_biometeo object rAedesSim \code{ biometeo} object   
#' @param  i_moitoring object rAedesSim \code{biodata} object concerning mosquito eggs field observations.
#' @param  i_biocontainer object rAedesSim \code{biocontainer} object .  
#' @param  larv_alpha numeric alpha larvae intra-competition factor.
#' @param  adu_alpha numeric adult larvae intra-competition factor.
#' @param  eggs_diap numeric Number of diapausant eggs considered.  
#' @param  egn numeric Number of eggs pour adults.  
#' @param  Ndayscenario_prim numeric Indicate the number of days considered after the first date of meteo data series. Defalut is 120.
#' @param  Ndayscenario_aut numericIndicate the number of days considered before the last date of meteo data series. Defalut is 90.
#' @param  deltatmed_prim numeric Mean air temperature bias cosidered for spring period.
#' @param  deltawater_prim numeric Mean water temperature bias cosidered for spring period.
#' @param  deltatmed_aut numeric Mean air temperature bias cosidered for spring period.
#' @param  deltawater_aut numeric Mean water temperature bias cosidered for spring period.
#' @param  ini_rmse numeric Starting index to calculate RMSE. Defalut is 1 that means all monitoring data are considered. 
#' @param  end_rmse numeric Ending index to calculate RMSE. Defalut is NULL that means all monitoring data are considered. 
#' @param  stocastic logical performs stocastic elaboration by poisson sampling at the end of simulation for each state.
#' @param  n_sampling numeric number  of resampling if stochastic is implemented see in \code{biomodel}. Default is 10.
#' @param  saveparameter logical  Save results about parameters. Default is FALSE. 
#' @return List of results with graph, esclosion graph ,the initial rmse and the full simulation object.
#' @author Istituto di Biometeorologia Firenze Italy  Alfonso crisci \email{a.crisci@@ibimet.cnr.it} ASL 2 LUCCA Marco Selmi \email{m.selmi@@usl2.toscana.it} 
#' @keywords biofitmodel
#' 
#' @import dygraphs
#' @import lubridate
#' @import xts
#' @importFrom aTSA accurate
#' 
#' @export



check_fit_initial=function(i_biometeo,
                           i_monitoring,
                           i_biocontainer,
                           larv_alpha,
                           adu_alpha,
                           eggs_diap,
                           egn=16,
                           Ndayscenario_prim=120,
                           Ndayscenario_aut=90,
                           deltatmed_prim=0,
                           deltawater_prim=0,
                           deltatmed_aut=0,
                           deltawater_aut=0,
                           ini_rmse=1,
                           end_rmse=NULL,
                           stocastic=FALSE,
                           n_sampling=10,
                           saveparameter=TRUE) 
                           {

                            ini_population=biopopulation(eggs=0,
                                                        larvae=0,
                                                        pupae=0,
                                                        adults=0,
                                                        eggs_diap=eggs_diap
                                                        )
                             ################################################################################################
                             
                           
                             i_bioparameters=bioparameters(alfa_l=larv_alpha,alfa_a=adu_alpha,egn=egn)
  
                             l_season=length(i_biometeo$tmed_est)
                             
                             if (l_season < (Ndayscenario_aut+1)) { stop("Length of time series is shorter than limit for autumn.")}
                             if (l_season < (Ndayscenario_prim+1)) { stop("Length of time series is shorter than limit for spring.")}
   
                             ini_aut=l_season-Ndayscenario_aut
                              
                             i_biometeo$tmed_est[1:Ndayscenario_prim]=i_biometeo$tmed_est[1:Ndayscenario_prim]+deltatmed_prim
  
                             i_biometeo$twater_est[1:Ndayscenario_prim]=i_biometeo$twater_est[1:Ndayscenario_prim]+deltawater_prim
  
                             i_biometeo$tmed_est[ini_aut:l_season]=i_biometeo$tmed_est[ini_aut:l_season]+deltatmed_aut
  
                             i_biometeo$twater_est[ini_aut:l_season]=i_biometeo$twater_est[ini_aut:l_season]+deltawater_aut
  
                             simulation=biomodel(i_biometeo=i_biometeo,
                                                 i_biocontainer=i_biocontainer,                            
                                                 i_biopopulation=ini_population,
                                                 i_bioparameters= i_bioparameters,
                                                 saveparameter = saveparameter,
                                                 stocastic = stocastic,
                                                 n_sampling= n_sampling)
  
                             fig_eggs=merge(simulation$ts_population[,"eggs"],
                                       i_monitoring$ts_data[,1],
                                       join = "inner")
                             observed=i_monitoring$ts_data[,1]      
                             fig_all=merge(simulation$ts_population,
                                           observed,
                                       join = "inner")
                             fig_all$index_day=NULL
                             fig_all$diapausant_eggs=NULL
                             rmse=accurate(fig_eggs[ini_rmse:end_rmse,1],
                                           fig_eggs[ini_rmse:end_rmse,2],
                                           k=1,output=FALSE)
  
                             names(fig_eggs)=c("Simulated Eggs","Observed Eggs")
  
                             resfig=dygraph(fig_eggs,
                                            xlab = "Year", 
                                            ylab = "Eggs Counts", 
                                            main=paste(i_monitoring$location,
                                                       "RMSE first", end_rmse,":",round(rmse[4],2),
                                                       "N_D_Eggs: ",eggs_diap,
                                                       "MaxJD: ",i_biometeo$datemax,
                                                       "Dw_spring: ",deltawater_prim,
                                                       "Dw_aut: ",deltawater_aut))
                             
                             resfig_all=dygraph(fig_all,
                                            xlab = "Year", 
                                            ylab = "Counts", 
                                            main=paste(i_monitoring$location,
                                                       "RMSE first", end_rmse,":",round(rmse[4],2),
                                                       "N_D_Eggs: ",eggs_diap,
                                                       "MaxJD: ",i_biometeo$datemax,
                                                       "Dw_spring: ",deltawater_prim,
                                                       "Dw_aut: ",deltawater_aut))
         
                             esclosion_graph=dygraph(i_biometeo$timeseries[,c("emergency")]*100,
                                                     xlab = "Year", 
                                                     ylab = "Rate for emergent eggs",main=paste(i_monitoring$location," - Courbe of Eggs exclosion "))
                             res=list(resfig_eggs=resfig,
                                      resfig_all=resfig_all,
                                      esclosion_graph=esclosion_graph,
                                      rmse=rmse,
                                      simulation=simulation)
  return(res)
}
