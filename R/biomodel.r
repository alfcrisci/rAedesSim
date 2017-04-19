#' biomodel
#'
#' @description The biomodel function is the core function that performs simulation.It consists in a S3 object where by using the others objects the simulation is performed.
#'              calculating the transitional rate of population by using daily environmental data. With these a a Runge-Kutta 
#' @param i_biometeo object rAedesSim object concerning biometeorological parameter of location where mosquito population are simulated.  
#' @param i_biocontainer object rAedesSim object concerning trap or habitat considered.  
#' @param i_biopopulation object rAedesSim object concerning mosquito population.  
#' @param i_bioparameters object rAedesSim object concerning biological parameter of mosquito population.
#' @param state numeric The state vector used in time simulation.The values are number of individuals L1 : Ovideposition rate, L3 : Eggs Mortality rate, 
#' L4: Eggs2larvae transition rate, L5 : Pupae2Larvae transition rate, L6 Larvae Mortality,L7  Pupae Mortality
#' L8 Pupae2Adult transition rate, L10 Adult mortality Default: L1=0,L3=0,L4=0,L5=0,L6=0,L7=0,L8=0,L10=0
#' @param time_interval numeric Interval times of simulation.
#' @param stocastic logical performs stocastic elaboration by poisson sampling at the end of simulation for each state.
#' @param n_sampling=10 numeric : number of resampling.
#' @param inibition logical taking into account density larvae inibition.Default is FALSE.
#' @param ID_sim character  ID of simulation.Default is FALSE.
#' @param saveparameter logical  Save results about paraemters. Default is FALSE. 
#' @return biomodel
#' @author  Istituto di Biometeorologia Firenze Italy  Alfonso crisci \email{a.crisci@@ibimet.cnr.it} ASL 2 LUCCA Marco Selmi \email{m.selmi@@usl2.toscana.it} 
#' @keywords  modeling
#'
#' @import deSolve
#' @import lubridate
#' @import xts
#' 
#' @export

biomodel  <- function(i_biometeo,
                      i_biocontainer,
                      i_biopopulation,
                      i_bioparameters,
                      state=c(L1=0,L2=0,L3=0,L4=0,L5=0,L6=0,L7=0,L8=0),
                      time_interval=24,
                      stocastic=TRUE,
                      n_sampling=10,
                      inibition=FALSE,
                      ID_sim="",
                      saveparameter=FALSE
		      )
	              { 
			   
			  
			    #########################################################################################à
			    # Check arguments of object.
			 
			   if (class(i_biopopulation)!="biopopulation") { stop("Object i_biopopulation argument must be of class biopopulation" )};
 			   if (class(i_bioparameters)!="bioparameters") { stop("Object i_bioparameters argument must be of class bioparameters" )};
                           if (class(i_biometeo)!="biometeo") { stop("Object i_biometeo argument must be of class biocontainer" )};
  			   if (class(i_biocontainer)!="biocontainer") { stop("Object i_biometeo don't have a  biocontainer object defined" )};
  			   
			    months=as.numeric(format(i_biometeo$dates,"%m"));
	                   
	                    #########################################################################################
			    # Define time vector interval and instance del loop to build outcome.
			   
			    Time <- c(0:(time_interval-1));
		            df_outcome_pop=data.frame(index_day=as.vector(1:i_biometeo$ndays),
				                      eggs_active=NA,
				                      diapausant_eggs=NA,
				                      larvae=NA,
				                      pupae=NA,
				                      adult=NA,
						      tot_eggs=NA,
						      neggs_diap=NA
				                      );
	
			    df_outcome_pop$eggs[1]=i_biopopulation$eggs;
			    df_outcome_pop$diapausant_eggs[1]=i_biopopulation$eggs_diap;
			    df_outcome_pop$larvae[1]=i_biopopulation$larvae;
			    df_outcome_pop$pupae[1]=i_biopopulation$pupae;
		            df_outcome_pop$adult[1]=i_biopopulation$adults;
			    df_outcome_pop$full_eggs[1]=i_biopopulation$eggs+i_biopopulation$eggs_diap;

				
			   df_outcome_par=data.frame(index_day=as.vector(1:i_biometeo$ndays),
                                                     f_ovo_a = NA,
			                             f_trans_u2l = NA,
			                             f_trans_l2p = NA,
			                             f_trans_p2a = NA,
			                             ma = NA,
			                             mu = NA,
			                             mp = NA,
			                             ml = NA,
						     d_emergency=NA,
						     d_induction=NA,
			                             inib_state=NA
			                             );
									
                           #########################################################################################
	                   # Main loop.
			   						  
			   for ( i_day in 2:i_biometeo$ndays) 
				
			   {
	                   #########################################################################################
			    # Define necessary i_bioparameters in the loop. Here eventual time serial modifications.
			    # Eggs - Larval - Pupae - Adult  are the stage considered in differential system
			   
			     alfa_lar = i_bioparameters$alfa_l; # Density dependent competition parameter for larvae 
			     alfa_ad = i_bioparameters$alfa_a; # Density dependent competition parameter for adults
			     par_egn = i_bioparameters$egn; # Mean eggs for individual female cycle
                             par_ef = i_bioparameters$ef;  # Pupal emergency factor
				   
			     inib_state = flag_inib(N_larvae=i_biopopulation$larvae,
			                            vol_water=i_biocontainer$pooled_volume_current,
			                            critical_density=i_bioparameters$l_density); # inibition state
			     inib_state=ifelse(inibition==FALSE,0,inib_state)	   
			   
			     f_ovo_a = trans_rates(as.numeric(i_biometeo$tmed_est[i_day]),
						   as.numeric(i_biometeo$twater_est[i_day]),
						   transition = "Ovideposition"); 
                             f_trans_u2l = trans_rates(as.numeric(i_biometeo$tmed_est[i_day]),
						       as.numeric(i_biometeo$twater_est[i_day]),
						       transition = "Eggs2larvae"); 
                             f_trans_l2p = trans_rates(as.numeric(i_biometeo$tmed_est[i_day]),
						       as.numeric(i_biometeo$twater_est[i_day]),
						       transition = "Larvae2Pupae");
                             f_trans_p2a = trans_rates(as.numeric(i_biometeo$tmed_est[i_day]),
						       as.numeric(i_biometeo$twater_est[i_day]),
						       transition = "Pupae2Adult"); 
                
			     urel=ifelse(is.na(i_biometeo$rhum[i_day]),60,i_biometeo$rhum[i_day]);
			    
			     ma = mortality_rate(as.numeric(i_biometeo$tmed_est[i_day]),
						 as.numeric(i_biometeo$twater_est[i_day]),
						 urel,
						 stage="Adult")    
                             mu = mortality_rate(as.numeric(i_biometeo$tmed_est[i_day]),
						 as.numeric(i_biometeo$twater_est[i_day]),
						 urel, 
						 stage="Eggs")  
                             mp = mortality_rate(as.numeric(i_biometeo$tmed_est[i_day]),
						 as.numeric(i_biometeo$twater_est[i_day]),
						 urel, 
						 stage="Pupae")
                             ml = mortality_rate(as.numeric(i_biometeo$tmed_est[i_day]),
						 as.numeric(i_biometeo$twater_est[i_day]),
						 urel,
						 stage="Larval")
			   
			     ##############################################################################################################################
			     # Define parameter for daily simulation
		             Parameter <- c(f_ovo_a=f_ovo_a,
                                             f_trans_u2l=f_trans_u2l,  
                                             f_trans_l2p=f_trans_l2p, 
                                             f_trans_p2a=f_trans_p2a,  
                                             ma=ma,        
                                             mu=mu,  
                                             mp=mp,  
                                             ml=ml,  
		                             alfa_l=alfa_lar, 
					     alfa_a=alfa_ad, 
		                             ef=par_ef,
					     egn=par_egn, 
			                     N_ini_a=i_biopopulation$adults,
			                     N_ini_p=i_biopopulation$pupae,
			                     N_ini_l=i_biopopulation$larvae, 
			                     N_ini_u=i_biopopulation$eggs,
			                     sexratio=i_bioparameters$sex_ratio,
					     deltatime=time_interval
					     );
							  
				#########################################################################################################################################
				# Model calling C interface library
							   
			        aedesmodelday <- suppressWarnings(ode(func = "derivs", 
						     y = state, 
						     parms = Parameter,
						     times=Time,
						     initfunc = "initmod",
						     dllname = "rAedesSim",
						     method = rkMethod("ode45"),
						     nout=1
						    ))
                
					
				##########################################################################################################################################
				# Take sum of  daily derivative 
				
				aedesmodelday=as.list(aedesmodelday[nrow(aedesmodelday),][2:9])

				if (stocastic == TRUE)
                                    { 
					for ( j in 2:(length(aedesmodelday)-1)) { 
						     aedesmodelday[[j]]=ifelse(aedesmodelday[j]<=0,0,round(mean(rpois(n_sampling,aedesmodelday[[j]]))));
					}
				    }
				
				
				
				#########################################################################################################################################
				
				if (inib_state == 1) { par_egn=ifelse(inib_state == 1,0,par_egn)}
				
				dadults=as.numeric(aedesmodelday$L7 - aedesmodelday$L8);
				dlarvae=as.numeric(aedesmodelday$L3-aedesmodelday$L5-aedesmodelday$L4);
				dpupae=as.numeric(aedesmodelday$L4-aedesmodelday$L6-aedesmodelday$L7);
			        deggs_pop=as.numeric(par_egn*aedesmodelday$L1-aedesmodelday$L2)+i_biopopulation$eggs_diap*i_biometeo$d_emergency[i_day]
				neggs_pop=ifelse(deggs_pop<0,0,deggs_pop);   
			     
				deggs_diap=as.numeric(ifelse(i_biometeo$d_induction[i_day]>0,neggs_pop*i_biometeo$d_induction[i_day],-1*i_biopopulation$eggs_diap*i_biometeo$d_emergency[i_day]));
				
			        neggs_diap_pop=ifelse(deggs_diap<0,0,deggs_diap);    
				   
				#########################################################################################################################################
				# update population
				   
				i_biopopulation$larvae=i_biopopulation$larvae+dlarvae;
				i_biopopulation$pupae=i_biopopulation$pupae+dpupae;
				i_biopopulation$adults=i_biopopulation$adults+dadults;          		
			        i_biopopulation$eggs=i_biopopulation$eggs+deggs_pop-deggs_diap;
				i_biopopulation$eggs_diap=i_biopopulation$eggs_diap+deggs_diap;
				i_biopopulation$eggs_diap=ifelse(as.numeric(format(i_biometeo$dates[i_day],"%j"))==i_biometeo$datemax,0,i_biopopulation$eggs_diap)
				#########################################################################################################################################
				# to avoid complete extiction
				   
				i_biopopulation$eggs_diap =ifelse(i_biopopulation$eggs_diap<0,0,i_biopopulation$eggs_diap)
				i_biopopulation$eggs=ifelse(i_biopopulation$eggs<0,0.1,i_biopopulation$eggs);
				i_biopopulation$eggs=ifelse(i_biopopulation$eggs<5 &as.numeric(format(i_biometeo$dates[i_day],"%j")) >274,i_biopopulation$pupae,i_biopopulation$eggs)
				i_biopopulation$larvae=ifelse(i_biopopulation$larvae<0,0,i_biopopulation$larvae);
			        i_biopopulation$pupae=ifelse(i_biopopulation$pupae<0,0,i_biopopulation$pupae);
				i_biopopulation$adults=ifelse(i_biopopulation$adults<0,0,i_biopopulation$adults);    
				   
			        #########################################################################################################################################
				# update outcomes
				
				df_outcome_pop[i_day,2]=i_biopopulation$eggs/i_biocontainer$nrecipients;
				df_outcome_pop[i_day,3]=i_biopopulation$eggs_diap/i_biocontainer$nrecipients;
				df_outcome_pop[i_day,4]=i_biopopulation$larvae/i_biocontainer$nrecipients;
				df_outcome_pop[i_day,5]=i_biopopulation$pupae/i_biocontainer$nrecipients;
			        df_outcome_pop[i_day,6]=i_biopopulation$adults/i_biocontainer$nrecipients;
				df_outcome_pop[i_day,7]=df_outcome_pop[i_day,2]+sum(neggs_pop,neggs_diap_pop)/i_biocontainer$nrecipients;
				df_outcome_pop[i_day,8]=neggs_diap_pop/i_biocontainer$nrecipients;
				   
				################################################################# 
				# update daily parameters
				if (saveparameter == TRUE) 
				{
                               
				df_outcome_par[i_day,2]=as.numeric(f_ovo_a);
				df_outcome_par[i_day,3]=as.numeric(f_trans_u2l);
				df_outcome_par[i_day,4]=as.numeric(f_trans_l2p);
				df_outcome_par[i_day,5]=as.numeric(f_trans_p2a);
				df_outcome_par[i_day,6]=as.numeric(ma);
                                
				df_outcome_par[i_day,7]=as.numeric(mu);
				df_outcome_par[i_day,8]=as.numeric(mp);
				df_outcome_par[i_day,9]=as.numeric(ml);
				df_outcome_par[i_day,10]=as.numeric(i_biometeo$d_emergency[i_day]);
				df_outcome_par[i_day,11]=as.numeric(i_biometeo$d_induction[i_day]);
				df_outcome_par[i_day,12]=as.numeric(inib_state);
				}
			    
			        }; # end loop days 	
				#########################################################################################################################################
			        # create time series objects
	
				ts_population=as.xts(zoo(df_outcome_pop,i_biometeo$dates))
				ts_parameter=as.xts(zoo(df_outcome_par,i_biometeo$dates))
				
				#########################################################################################################################################
			   
                                object  <-  list(name_model="Rbiosim",
				                 ID_sim=ID_sim,
                                                 timestep_integration = length(Time),
					         ts_population=ts_population,
					         ts_parameter=ts_parameter,
				   	         stocastic=stocastic,
					         n_sampling=n_sampling,
					         ID=i_biocontainer$ID,
					         site_name=i_biocontainer$site_name,
					         sp_obj=i_biocontainer$sp_obj,
                                                 lat=i_biocontainer$lat,
			                         lon=i_biocontainer$lon
					         );
			            									  
				               

                #########################################################################################################################################
				             	    
                attr(object,"name_model") <- "Model's name"
                attr(object,"ID_sim") <- "ID label of simulation"
                attr(object,"timestep_integration") <- "Time step of daily integration"
                attr(object,"ts_population") <- "Model outcomes of mosquito as multivariate timeseries xts object"
                attr(object,"ts_parameter") <- "Model outcomes of simulation parameters  xts object"
                attr(object,"ID")<-"ID label of container set"
                attr(object,"site_name")<-"Name of sites"
                attr(object,"sp_obj")<-"SpatialPointDataFrame of location"
                attr(object,"lat")<-"latitude coordinate of cointainer"
                attr(object,"lon")<-"longitude coordinates of cointainer"
   
                class(object) <- "biomodel"
                return(object)
}
