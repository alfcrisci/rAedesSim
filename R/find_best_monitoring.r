#' find_best_monitoring
#'
#' @description 
#' \code{find_best_monitoring} Function to semplify  model tuning by using \code{biofitmodel} rAedesSim Class.
#' 
#' @param i_meteo object: rAedesSim meteo  concerning weather data.  
#' @param i_biocointaner object: rAedesSim object concerning \code{biocontainer} object used in simulation.  
#' @param i_monitoring object: rAedesSim \code{biodata} object concerning mosquito eggs field observations.
#' @param initial_eggs numeric: initial values for \code{biodata} object concerning mosquito eggs field observations.
#' @param range_alpha_a numeric: rAedesSim vector of sorted  guess' values of female adult competition. Default is obtained by c(0,seq(0,0.002,0.001)).
#' @param range_alpha_l numeric: rAedesSim vector of sorted guess' values of larval competition.  Default is obtained by seq(0.6,1.6,0.2).
#' @param range_density_l numeric: rAedesSim object guess' values of maximaum larval competion in function of waer volume in cointaner. Default is obtained by 70.
#' @param stocastic logical: if stocastic simulation by \code{biomodel} are carried out in simulation.Defalut is True.
#' @param n_sampling numeric: number  of resampling if stochastic is implemented by \code{biomodel}. Default is 10.
#' @param inibition logical: if larval density is considered. Defalut is FALSE.
#' @return rAedesSim object biofitmodel object.
#' @author  Istituto di Biometeorologia Firenze Italy  Alfonso crisci \email{a.crisci@@ibimet.cnr.it} ASL 2 LUCCA Marco Selmi \email{m.selmi@@usl2.toscana.it} 
#' @keywords   rAedesSim modeling
#'
#' 
#' 
#' @export


find_best_monitoring=function(i_meteo,
                              i_biocontainer,
                              i_monitoring,
                              initial_eggs=100,
                              range_alpha_a=c(0,seq(0,0.002,0.001)),
                              range_alpha_l=seq(0.6,1.6,0.2),
                              range_density_l=70,
                              stocastic=TRUE,
                              n_sampling=10,
                              inibition=FALSE) {
							 
                                                i_biometeo=biometeo(i_meteo,i_biocontainer)
                                                i_biopopulation=biopopulation(eggs=initial_eggs,larvae=0,pupae=0,adults=0,eggs_diap=initial_eggs)

                                                return(biofitmodel(i_biometeo=i_biometeo,
                                                                   i_biopopulation=i_biopopulation,
								   i_biocontainer=i_biocontainer,
								   i_monitoring=i_monitoring,
								   range_alpha_a=range_alpha_a,
								   range_alpha_l=range_alpha_l,
								   range_density_l=range_density_l,
							           stocastic=stocastic,
								   n_sampling=n_sampling,
								   inibition=inibition
			                                    ))				
                           }
