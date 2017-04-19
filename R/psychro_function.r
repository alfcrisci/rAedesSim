#' pvapsat
#'
#' Saturation saturation pressure of water (hPa) as function of air temperature (degC).
#'
#' The saturation pressure of water is the pressure at or above which water
#' would condense from from the air. It is a function of air temperature.
#'
#'
#' @param t numeric air temperature in Celsius
#' @return pvapsat 
#' @export
#' 

pvapsat<-function(tair=24)
{
   es=6.1078*10^((7.5*tair)/(237.7+tair));
  
  return(es);
}




#' delta_vapor_density
#'
#' @description  Return  vapour density at surface air-water.
#'
#'
#'
#'
#' @param tair numeric  Air temperature in Celsius
#' @param twater numeric Water temperature in Celsius
#' @param urel numeric relative humidity in percentage.
#' @return delta_vapor_density
#' @export
#' 

delta_vapor_density<-function(tair=24,twater=21,urel=60)
{  
   R=461.48; #J K-l Kg-l.
   g_standard=9.80665;
   deltad=((100*pvapsat(twater))/(R*(twater+273)))-(100*(urel/100)*pvapsat(tair))/(R*(tair+273));
   return(deltad)
}

#' delta_vapor_pressure
#'
#' @description  Difference water vapour pressure between cointaner water surface an free air. 
#'
#'
#' @param tair numeric  Air temperature in Celsius.
#' @param twater numeric Water temperature in Celsius.
#' @param urel numeric  Air Relative umidity.
#' @return delta_vapor_pressure
#' @export
#' 

delta_vap_pressure<-function(tair=24,twater=21,urel=60)
{   
   delta_p=((100*pvapsat(twater))-(100*(urel/100)*pvapsat(tair)))

   return(delta_p)
}


#' delta_air_density
#'
#'  @description  Gradient or difference in air density between cointaner water surface an free air. 
#'
#'
#' @param tair numeric: Air temperature in Celsius
#' @param twater numeric: Water temperature in Celsius
#' @param urel numeric: Air Relative umidity.
#'
#'
#' @return delta_vapor_pressure
#' @export
#' 

delta_air_density<-function(tair=24,twater=21,rh=60)
{   
    R=461.48; #J K-l Kg-l.
	  g_standard=9.80665
    delta=((100*pvapsat(twater))/(R*(twater+273)))+(100*(rh/100)*pvapsat(tair))/(R*(tair+273))/2;
    return(delta)
}



#' mass_diffusion_vap
#'
#'  @description  Diffusion coeffcient of vapor (mq/s).in function to  temperature for Air-Water Vapor Mixtures  Bolz and Tuve (1976). 
#'
#'
#' @param tair numeric: Air temperature in degC.
#'
#'
#' @return mass_diffusion_vap
#' @export
#' 

mass_diffusion_vap<-function(tair=24)
{  
  diff_t=(2.07*(10**-6))+4.479*(10**-8)*(tair+273)+1.656*(10**-10)*(tair+273)^2;
  return(diff_t);
}


#' mass_transport_coef
#'
#' @description  Mass transfer coefficient in function to air temperature in forced convection  (v_air > 0.1).
#'
#'
#' @param tair numeric Air temperature in degC.
#' @param twater numeric Water temperature in degC.
#' @param urel numeric Air Relative umidity.
#' @param v_air numeric Velocity of air movement (meter per second).
#'
#' @return mass_transport_coef
#' @export
#' 

mass_transport_coef<-function(tair=23,twater=21,urel=50,v_air=0.5)

{  ktrasp_vap=3.4*(10**-8);
   k_trasp=ktrasp_vap*delta_vap_pressure(tair,twater,urel);
   return(k_trasp);
}


#' evaporation_params
#'
#' @description The function return a list of parameters in order to assess the evaporative losses 
#' for a generic cylndric cointainer under forced convetion (v_air > 0.1). Parameters done are in order:the  Reynolds number,Raylegh number,
#' the Schmidt and Sherwood numbers,the hm coeffcient of vapor trasport, the mv evaporative loss in time in in kg/sec and mvday the evaporative loss in a day (mg/day).
#'
#'
#' @param tair numeric  Air temperature in degC.
#' @param twater numeric Water temperature in degC.
#' @param urel numeric  Air Relative umidity.
#' @param v_air numeric Velocity of air movement.
#' @param L numeric Dimension factor to caracterize the cointainer.A cyliner model is considered. Generally L is a diameter lenght measured in meter. Default is 0.1.
#' @param A numeric Area of air-water interface. Default is NULL.
#' @return evaporation_params
#' @export
#' 

evaporation_params<-function(tair=24,twater=21,urel=60,v_air=0.5,L=0.1,A=NULL)
{ 
             A = ifelse(is.null(A),3.1415*(L/2)**2,A); # Area of water container
             mv = 18.016; #kg on kmole
             d_aria_20 = 0.8216; #kg on m3
             cp_aria = 1.005; #KJ/kg K # Specific heat of air
             visc_aria = 1.583*(10**-5);# air viscosity
             g_standard = 9.80665;
             dilataz_aria = 1/(tair+273.16);
             Rzero=8314/mv; #J K-l kg-l.
             entalp_vap = 40.8; #kJ/mol
	     
	           res = list()
	
             res$reynolds = ( v_air * L )/ visc_aria;
             res$raylegh = g_standard * ((delta_air_density(tair,twater,urel) * (L**3))/(delta_vap_pressure(tair,twater,urel) * visc_aria * mass_diffusion_vap(tair)));
             res$schmidt = 0.6;
             res$sherwood = 0.230*((res$schmidt)**(1/3)) * res$raylegh**(0.321);
             res$hm = (res$sherwood * mass_diffusion_vap(tair)) / L;
             res$mev = res$hm*A*(delta_vap_pressure(tair,twater,urel)/(461 * (twater + 273 ))); # kg/sec
             res$mevday = res$mev * 3600 * 24 * 1000; # mg day integrate from a seconds
	
             return(res)
}

 


