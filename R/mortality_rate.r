#' mortality_rate
#'
#' @description Calculate daily mortality rates for Aedes in function to the stage considered. 
#' 
#' @param tair   numeric:   Temperature in Celsius Degree of (I) air: for "Adult" stage - (II)  for "Diapausant Eggs" the Minimal winter seasonal temperature  .
#' @param twater  numeric: Water Temperature in Celsius degree.    
#' @param urel numeric: Relative humidity in \% ( 0-100). 
#' @param stage   char string: Name of stage ( Diapausant Eggs - Eggs - Larval - Pupae - Adult). Default is "Eggs".
#' @param air_dryness  logical: if  the saturation deficit of air  is taken into account as supplementary agent for adult mortality.
#' @return m_rate numeric: daily mortality rate for the stage. 
#' @references 
#'   Poletti et al,2011, Transmission Potential of Chikungunya Virus and Control Measures: The Case of Italy. PLoS ONE 6(5).
#'   Magori et al,2009, Skeeter Buster: A Stochastic, Spatially Explicit Modeling Tool for Studying Aedes aegypti Population Replacement and Population Suppression Strategies. PLoS Negl Trop Dis 3(9). 
#'   Brady et al,2013, Modelling adult Aedes aegypti and Aedes albopictus survival at different temperatures in laboratory 
#'   and field settings. Parasites & Vectors 2013, 6:351  
#'
#' @author  Istituto di Biometeorologia Firenze Italy  Alfonso crisci \email{a.crisci@@ibimet.cnr.it} ASL 2 LUCCA Marco Selmi \email{m.selmi@@usl2.toscana.it} 
#' @keywords  Mortality rates
#'
#' 
#' 
#' @export
 
mortality_rate<-function( tair = 24, twater = 24, urel =60, stage = "Eggs", air_dryness = FALSE) {

sature_deficit= 6.1078 * 10**(( 7.5 * tair ) / ( 237.7 + tair ) ) *( 1- urel/100 );
sature_deficit=ifelse( sature_deficit > 30 , 30 , sature_deficit);
	
if ( stage == "Eggs" ) 
   {
	me=ifelse((twater < 35.0)||(twater > 15.0),(506-(506*exp(-(((twater-25)/27.3)**6)))),1.220816);
	me=ifelse(me<0.05,0.05,me)
	return(me)
	
	}

if ( stage == "Larval" ) 
    {
	return(0.029 + (858 * exp( twater -43.4 )));
	}
	
	
if ( stage == "Pupae") 
   {
	return(0.021 + ( 37.000* exp( twater - 36.800 ) ));	
   }


if ( stage == "Adult" && air_dryness == FALSE) 
   { 

	return(0.031 + (95820*exp(tair - 50.400)));
	
   }

if ( stage == "Diapausant Eggs")
   {   
	tmin=round(tmin,0);
	if(tmin  > 2){
		return(0);
	}
	else if (tmin > -10 && tmin <= 2){
		return(-1.2346 * tmin + 8.179);
	}
	else if (tmin <= -10 && tmin > -14){
		return(-18.75 * tmin - 166.67);
	}
	else if(tmin < -14){
		return(100);
	}	
   }
   
if ( stage=="Adult" && air_dryness==TRUE) 
    {
	return(0.031 + (95820*exp(tair - 50.400))+0.04*(sd-10)/(30-10));
    }

}

