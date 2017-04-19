#' day_length
#'
#' @description Calculate the raw lenght of day in hours related to the date
#' 
#' @param dates  Date in  YMD format   
#' @param lon Numerical Geographical Longitude of site in decimal degrees (negative == West) . Default is 11.10.
#' @param lat Numeric Geographical Latitude of site in decimal degrees. Default is 43.77.
#' @return Return sunrise time, suset time,length of day in hours decimal.
#' @references  Teets, D.A. 2003. Predicting sunrise and sunset times. The College Mathematics Journal 34(4):317-321.
#' @seealso \code{\link{weigthdry}}
#'
#' @author  Istituto di Biometeorologia Firenze Italy  Alfonso crisci \email{a.crisci@@ibimet.cnr.it} ASL 2 LUCCA Marco Selmi \email{m.selmi@@usl2.toscana.it} 
#' @keywords  dry days
#' 
#' 
#' @export

day_length<-function(dates,lon=11.10,lat=43.77){
  
  require(lubridate)
  
  d=yday(dates);
  
  ## Function to convert degrees to radians
  
  rad<-function(x)pi*x/180					   
 
  
  ##Radius of the earth (km)

  R=6378
  
  ##Radians between the xy-plane and the ecliptic plane

  epsilon=rad(23.45)

  ##Convert observer's latitude to radians

  L=rad(lat)

  
  timezone = -4*(abs(lon)%%15)*sign(lon)

  ## The earth's mean distance from the sun (km)

  r = 149598000

  theta = 2*pi/365.25*(d-80)

  z.s = r*sin(theta)*sin(epsilon)
  r.p = sqrt(r^2-z.s^2)

  t0 = 1440/(2*pi)*acos((R-z.s*sin(L))/(r.p*cos(L)))
  
  ##a kludge adjustment for the radius of the sun

  that = t0+5 

  ## Adjust "noon" for the fact that the earth's orbit is not circular:

  n = 720-10*sin(4*pi*(d-80)/365.25)+8*sin(2*pi*d/365.25)

  ## now sunrise and sunset are:

  sunrise = (n-that+timezone)/60
  sunset = (n+that+timezone)/60
  daylenght=sunset-sunrise
  
  return(list("sunrise" = sunrise,"sunset" = sunset,"daylength" = daylenght))
}
