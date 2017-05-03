#' biofitmodel
#'
#' @description Biofitmodel is the function to find,  troughout a grid scheme, which are the best intra competition parameters Alpha Adults and Alpha Larvs 
#' by using the minimal RMSE criterion, when environmental  and population data are fixed.To build grid the function needs the ranges of alpha parameters. 
#' It is possible to esplicit the numerosity  of data used to perform the validation indicating the index of starting and ending by using observed data as reference.
#' 
#' @param i_biometeo object rAedesSim \code{biometeo} object   
#' @param i_biopopulation object  rAedesSim  \code{biopopulation} object .  
#' @param i_biocontainer object  rAedesSim \code{biocontainer} object .  
#' @param i_monitoring object  rAedesSim \code{biodata} object concerning mosquito eggs field observations.
#' @param range_alpha_a numeric rAedesSim range of female adult competition alpha. Default values are done by seq(0.005,0.01,0.001).
#' @param range_alpha_l numeric rAedesSim range of intra-larval competition alpha  .  Default values are done by seq(0.005,0.01,0.001)
#' @param range_density_l numeric rAedesSim object range of maximum larval density in liter of water volume.Default value is 100.
#' @param n_sampling numeric number  of resampling if stochastic is implemented see in \code{biomodel}. Default is 10.
#' @param inibition logical if larval density is considered in \code{biomodel}.Default is FALSE.
#' @param plotresults logical if is true a plot is done. Default is FALSE.
#' @param testRMSE logical if test the root mean square error of simualtions. Default is FALSE.
#' @param ini_rmse numeric Starting position index to calculate RMSE on the observed time series.Defalut is 1.
#' @param end_rmse numeric Ending position index to calculate RMSE on observed data.Default is NULL and means that all observed data are considered if ini_rmse = 1.
#' @return biofitmodel
#' @author  Istituto di Biometeorologia Firenze Italy  Alfonso crisci \email{a.crisci@@ibimet.cnr.it} ASL LUCCA Marco Selmi \email{marco.selmi@@uslnordovest.toscana.it }
#' @keywords biofitmodel
#' 
#' @import deSolve
#' @import lubridate
#' @import xts
#' @importFrom aTSA accurate
#' @importFrom rms.gof rms.pval
#' 
#' @export

biofitmodel  <- function(i_biometeo,
                         i_biopopulation,
                         i_biocontainer,
                         i_monitoring,
                         range_alpha_a=seq(0.005,0.01,0.001),
                         range_alpha_l=seq(0.005,0.01,0.001),
                         range_density_l=100,
                         stocastic=TRUE,
                         n_sampling=10,
                         inibition=FALSE,
                         plotresults=FALSE,
                         testRMSE=FALSE,
                         ini_rmse=1,
                         end_rmse=NULL
                        )
	             { 
			   
			   
			    #########################################################################################à
			    # Check arguments of object.
			 
			     if (class(i_biometeo) != "biometeo") { stop(" Object argument must be a  rAedesSim biometeo class." )};
			     if (class(i_biopopulation) != "biopopulation") { stop(" Object argument must be a rAedesSim bioparameters class." )};
			     if (class(i_biocontainer) != "biocontainer") { stop(" Object argument must be a rAedesSim biocontainer class." )};
  			 #########################################################################################
			    # Create matrix and list for testing parameters
				
			    
			    replies=cbind(merge(range_alpha_a,range_alpha_l),z=rep(range_density_l,nrow(merge(range_alpha_a,range_alpha_l))))
			    replies=replies[,1:3]
			    biopar_list=list()
			    for ( i in 1:nrow(replies)) { biopar_list[[i]]=bioparameters(alfa_l=replies$y[i], alfa_a=replies$x[i],l_density=replies$z[i])}
			    #########################################################################################
			    # Create list for outcomes
			    simul_ts=list();
			    simul_RMSE=numeric(nrow(replies));
			    simul_test=list();
			    simul_RMSE_pval=numeric(nrow(replies));
			    success_vector=logical(nrow(replies));
			    #######################################################################################################
				
                             for ( i in seq_along(biopar_list)){ message(paste("Working on:", i));
                                                                 success_vector[i] = FALSE;
                                                                                                 
                                                                  tryCatch({simulation = biomodel(i_biometeo=i_biometeo,
                                                                                                  i_biocontainer=i_biocontainer,
                                                                                                  i_biopopulation=i_biopopulation,
                                                                                                  i_bioparameters=biopar_list[[i]],
                                                                                                  stocastic = stocastic,
                                                                                                  n_sampling = n_sampling,
                                                                                                  inibition = inibition);
                                                                                                  message(paste("Processed case:", i,"Simulation ok!"));
                                                                                                  success_vector[i] = TRUE;
                                                                                                   },
                                                                                error=function(cond) 
										                  {
                                                                                                  success_vector[i] = FALSE
                                                                                                  simul_ts[[i]] = NA
                                                                                                  simul_RMSE[i]=NA
                                                                                                  simul_test[[i]]=NA
                                                                                                  simul_RMSE_pval[i]=NA
                                                                                                  message(paste("Processed case:", i,"Simulation aborted!"))
                                                                                                  },
                                                                                                  finally=
									                          {
											            message(paste("Processed case:", i,"Done!"))
                                                                                                  
                                                                                                  }
                                                                                                  ); #end of try-catch
						                                                 
								                                 if (success_vector[i] == TRUE) {
                                                                                                 Eggs=simulation$ts_population$full_eggs
                                                                                                 Eggs_obs=i_monitoring$ts_data
                                                                                                 merged=merge.xts(Eggs,Eggs_obs,join = "inner");
                                                                                                 simul_ts[[i]]=merged; 
												 end_rmse=ifelse(is.null(end_rmse),nrow(merged),end_rmse)
                                                                                                 if ( ini_rmse >= end_rmse )  { stop("Monitoring data are unavailable!")
															      }
                                                                                                 df_rmse=data.frame(observed=as.vector(merged$Eggs_obs[ini_rmse:end_rmse]),
                                                                                                                    predicted=as.vector(merged$full_eggs[ini_rmse:end_rmse]))
                                                                                                 df_rmse=na.omit(df_rmse) 
												 if (nrow(df_rmse)>1) 
												       {  simul_test[[i]]=accurate(df_rmse[,1],df_rmse[,2],k=1,output=FALSE)
													  simul_RMSE[i]=simul_test[[i]][4];
                                                                                                          if ( testRMSE == TRUE) 
													     {
													      simul_RMSE_pval[i]=rms.pval(df_rmse[,1],df_rmse[,2],num_sim =500)
													      }
                                                                                                       }													  
								                                 }	
		   
				                                                             }                        
                
				#########################################################################################################################################
				# Fill results
	
				res_simul_RMSE=NA;
				res_simul_test=NA;
				res_simul_RMSE_pval=NA;
				res_simul=NA;
				replies_best=NA; 
	
				best=head(which.min(simul_RMSE),1)
				names(replies)<-c("alpha_a","alpha_l","density_max_l");
	
			    	if (length(best) == 0) 
				    {best=NA;res_simul=NA;replies_best=c(NA,NA,NA);
				    } 
				else {
				     res_simul=simul_ts[[best]];
				     res_simul_RMSE=simul_RMSE[best];
				     res_simul_test=simul_test[[best]];
				     res_simul_RMSE_pval=simul_RMSE_pval[best];
				     replies_best=replies[best,];	     
				     }
				
				
				#########################################################################################################################################
				# Fill spatial objects
				
				i_biocontainer$sp_obj$alpha_a=as.numeric(replies_best[1])
				i_biocontainer$sp_obj$alpha_l=as.numeric(replies_best[2])
				i_biocontainer$sp_obj$density_max_l=as.numeric(replies_best[3])
				
				plot_ts=NULL;
				
	      if ( plotresults == TRUE)   { plot_ts=plot(simul_ts[[best]],
	                                                 main = paste("Observed (red) & Assessed (Black) - ",as.character(i_monitoring$location),"-",as.character(i_biocontainer$type),"-","Stage's competivity index: Larvae=",as.character(replies_best$alpha_l)," Adults=",as.character(replies_best$alpha_a)," Larval MaxDensity=",replies_best[3]),
	                                                 cex.axis = 1.2,
	                                                 cex.main = 2.5,
	                                                 legend.loc = "bottomright", 
	                                                 legend.pars = list(bty = "n",cex=2,horiz=TRUE),
	                                                 legend.names = c("Observed","Assessed")) 
				                             }
				#########################################################################################################################################
			   
        object  <-  list(name_model="rAedesSim",
                         location=as.character(i_monitoring$location),
                         guess_parameter=replies,
                         best_simul = res_simul,
                         best_simul_RMSE =res_simul_RMSE ,
                         best_simul_test =res_simul_test ,
                         best_simul_RMSE_pval =res_simul_RMSE_pval,
                         par_fitted_best=replies_best,
                         simul_RMSE=simul_RMSE,
                         simul_test=simul_test,
                         simul_RMSE_pval=simul_RMSE_pval,
                         n_replies=length(na.omit(simul_RMSE)),
                         stocastic=stocastic,
                         n_sampling=n_sampling,
                         inibition=inibition,
                         ID=i_biocontainer$ID,
                         sp_obj=i_biocontainer$sp_obj,
                         lat=i_biocontainer$lat,
                         lon=i_biocontainer$lon,
                         plot_ts=plot_ts,
                         ini_rmse=ini_rmse,
                         end_rmse=end_rmse
                        );
			            									  
                #########################################################################################################################################
				             	    
        attr(object,"name_model") <- "Model's name"
				attr(object,"location") <- "Location's name."
        attr(object,"guess_parameter") <- "Matrix of guess values."
        attr(object,"best_simul") <- "Timeseries object: Best simulation taking into account diapause."
        attr(object,"best_simul_RMSE") <- "Best Root mean square Error value in the grid scheme of simulations."
        attr(object,"best_simul_test") <- "Statistics associated to best RMSE in the grid scheme of simulations."
        attr(object,"best_simul_RMSE_pval") <- "P.val of RMSE in the grid scheme of simulations."
        attr(object,"par_fitted_best") <- "Parameter fitted."
        attr(object,"simul_RMSE") <- "Root mean square error value for  all simulations."
        attr(object,"simul_test") <- "Statistics of reliability for all simulations."
        attr(object,"simul_RMSE_pval") <- "Root mean square error pval for all simulations.Significant under p<.005 "
        attr(object,"n_replies") <- "Number of simulations."
        attr(object,"stocastic") <- "If stocasticity in simulation are considered."
        attr(object,"n_sampling") <- "Number of resampling."
        attr(object,"inibition") <- "Logical if larval inibition are taken into account in simulation."
        attr(object,"ID")<-"ID label of container set."
        attr(object,"sp_obj")<-"SpatialPointDataFrame of location."
        attr(object,"lat")<-"latitude coordinates of simulations."
        attr(object,"lon")<-"longitude coordinates of simulations."
        attr(object,"plot_ts")<-"Plot fitted vs observed eggs."
        attr(object,"ini_rmse")<-"Starting index to evaluate RMSE in monitoring data time series."
        attr(object,"end_rmse")<-"Ending index to evaluate RMSE in monitoring data time series."
        class(object) <- "biofitmodel"
        return(object)
}

