#' trans_rates
#'
#' @description Calculate daily transition rates for Aedes in function to trasition among the stage considered.
#' Daily rate (day^-1)  consist in a number in day. It is a unitarian parameter. It is based on Poletti et al work (2011) and Magori et al work.
#' 
#' @param tair   Numeric   Temperature in Celsius degree  for "Ovidepostion" transition.
#' @param twater  Numeric Water Temperature in Celsius degree for  any transition except "Ovidepostion".   
#' @param transition   Character Name of transition :  "Ovidepostion" - "Eggs2larvae" - "Larvae2Pupae" - "Pupae2Adult". Default is "Ovidepostion".
#' @return trans_rates
#' @references Poletti P., Messeri G, Ajelli M., Vallorani R., Rizzo C., et al. (2011) Transmission Potential of Chikungunya Virus and Control Measures: The Case of Italy. PLoS ONE 6(5): e18860. doi:10.1371/journal.pone.0018860
#'            Magori K, Legros M, Puente ME, Focks DA, Scott TW, Lloyd AL, Gould F. 2009. Skeeter Buster: A Stochastic, Spatially Explicit Modeling Tool for Studying Aedes aegypti Population Replacement and Population Suppression Strategies. PLoS Negl Trop Dis 3(9): e508. doi:10.1371/journal.pntd.0000508
#' @author  Istituto di Biometeorologia Firenze Italy  Alfonso crisci \email{a.crisci@@ibimet.cnr.it} ASL 2 LUCCA Marco Selmi \email{m.selmi@@usl2.toscana.it} 
#' @keywords  Container
#'
#' 
#' 
#' @export
 
trans_rates<-function( tair = 24, twater = 24,flag_I=FALSE, inib=0.63,transition = "Ovideposition") {

ovodepos=0
transferlp=0
transferhl=0
transferpa=0

if ( transition == "Ovideposition" ) 
   
    { ovodepos=( 1 / ((0.046 * tair^2) - (2.770 * tair) + 45.300));
                              return( ovodepos );
    }
                            
						
if ( transition == "Eggs2larvae" ) 
                            { transferlp=1 / ((0.120*twater*twater) - (6.600*twater) + 98.000);
                             return(transferlp);
                             }

							 
if ( transition == "Larvae2Pupae" ) 
                            { transferhl=1 / (6.900 - (4.000*exp(-(((twater-20)/4.100)*((twater-20)/4.100)))))
							  if(flag_I==TRUE)
                                          { transferhl=transferhl*(1-inib);
                                          }
                              return(transferhl);
                             }
							 
if ( transition == "Pupae2Adult" ) 
                            { transferpa=1 / ((0.027*twater*twater) - (1.700*twater) + 27.700);	
                             return(transferpa);
                             }	

							 
							 
}							 
