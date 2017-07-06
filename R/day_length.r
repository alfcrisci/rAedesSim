#' day_length
#'
#' @description Calculate the raw lenght of day in hours related to the date
#' 
#' @param dates  Date in  YMD format   
#' @param lon Numerical Geographical Longitude of site in decimal degrees (negative == West) . Default is 11.10.
#' @param lat Numeric Geographical Latitude of site in decimal degrees. Default is 43.77.
#' @param timezone character  Default is "Europe/Rome".
#' @return Return sunrise time, suset time,length of day in hours decimal.
#' @references  Teets, D.A. 2003. Predicting sunrise and sunset times. The College Mathematics Journal 34(4):317-321.
#' @seealso \code{\link{weigthdry}}
#'
#' @author  Istituto di Biometeorologia Firenze Italy  Alfonso crisci \email{a.crisci@@ibimet.cnr.it} ASL LUCCA Marco Selmi \email{marco.selmi@@uslnordovest.toscana.it }
#' @keywords  dry days
#' 
#' @import lubridate
#' @export

day_length<-function(dates,lon=11.10,lat=43.77,timezone="Europe/Rome"){
                     d=yday(as.Date(dates,tz=timezone));
                     rad<-function(x)pi*x/180					   ## Internal function to convert degrees to radians
                     R=6378 ##Radius of the earth (km)
                     epsilon=rad(23.45) ##Radians between the xy-plane and the ecliptic plane
                     L=rad(lat) ##Convert observer's latitude to radians
                     timezonec = -4*(abs(lon)%%15)*sign(lon)
                     r = 149598000 ## The earth's mean distance from the sun (km)
                     theta = 2*pi/365.25*(d-80)
                     z.s = r*sin(theta)*sin(epsilon)
                     r.p = sqrt(r^2-z.s^2)
                     t0 = 1440/(2*pi)*acos((R-z.s*sin(L))/(r.p*cos(L)))
                     that = t0+5  ##a kludge adjustment for the radius of the sun
                     n = 720-10*sin(4*pi*(d-80)/365.25)+8*sin(2*pi*d/365.25) ## Adjust "noon" for the fact that the earth's orbit is not circular.
                     sunrise = (n-that+timezonec)/60
                     sunset = (n+that+timezonec)/60
                     daylenght=sunset-sunrise   ## now sunrise and sunset are:
                     return(list("sunrise" = sunrise,"sunset" = sunset,"daylength" = daylenght))
}
