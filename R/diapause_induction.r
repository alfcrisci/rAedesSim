#' diapause_induction
#'
#' @description Calculate the daily diapause incidence  in autumnal eggs of Aedes albopictus. 
#' It is based on the work published by Lacour G. at al (2015), Seasonal Synchronization of Diapause Phases in Aedes albopictus (Diptera: Culicidae),
#' in PLOS. Date limits are imposed by default. Julian dates 237 for diapause starting and 357 for diapause ending. It needs location coordinates where throughout  \code{"day_length"} function astronomical daily sunlight duration is calculated.
#' @param dates  Date in  YMD format   
#' @param lon Numeric Geographical Longitude of site in decimal degrees. Default is 7.216.
#' @param lat Numeric Geographical Latitude of site in decimal degrees. Default is 43.77.
#' @return Return Pourcentage of diapausant eggs. 
#' @references   Lacour G, Chanaud L, L’Ambert G, Hance T (2015) Seasonal Synchronization of Diapause Phases in Aedes albopictus (Diptera: Culicidae). PLoS ONE 10 12 : e0145311. doi:10.1371 journal.pone.0145311
#' @author  Istituto di Biometeorologia Firenze Italy  Alfonso crisci \email{a.crisci@@ibimet.cnr.it} ASL 2 LUCCA Marco Selmi \email{m.selmi@@usl2.toscana.it} 
#' @keywords  diapause,autumn
#' 
#' 
#' @export

diapause_induction=function(dates,lon=7.216,lat=43.658) {
  data(diapause_model)
  photo=rAedesSim::day_length(dates,lon,lat)$daylength
  jd=as.numeric(format(as.Date(dates),"%j"))
  falling=predict(diapause_model, data.frame(photo_lacour = photo))
  temp=data.frame(jd=jd,falling=falling)
  temp$falling[which(as.numeric(temp$jd) > 350)] <- 0
  temp$falling[which(as.numeric(temp$jd) < 237)] <- 0
  return(temp$falling/100)
}
