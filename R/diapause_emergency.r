#' diapause_emergency
#'
#' @description Calculate the daily diapause incidence of egg springtime hatching in probabilistic terms. 
#' 
#' @param dates  Date in  YMD format   
#' @param lon Numeric Geographical Longitude of site in decimal degrees. Default is 7.216.
#' @param lat Numeric Geographical Latitude of site in decimal degrees. Default is 43.77.
#' @param CPPmin  Numeric Number of hours of critical photoperiod. Default is 11.5
#' @param datemax Date Date of spring diapause awekaning in MM-DD. Default is 04-20.
#' @param bounddays Numeric Days of incertain for diapause emergency bounds. Default is 10.
#' @param varjd Numeric Days of variance for  diapause emergency courbes in days. Default is 10.
#' @param truncated logical Modeling with truncated gaussian respecting the limits. Default is FALSE.
#' @return Return diapause emergency in percentage.
#' @references   Lacour G, Chanaud L, L’Ambert G, Hance T (2015) Seasonal Synchronization of Diapause Phases in Aedes albopictus Diptera: Culicidae. PLoS ONE 10-12: e0145311. doi:10.1371 journal.pone.0145311
#' @seealso \code{diapause_induction}
#' @author  Istituto di Biometeorologia Firenze Italy  Alfonso crisci \email{a.crisci@@ibimet.cnr.it} ASL 2 LUCCA Marco Selmi \email{m.selmi@@usl2.toscana.it} 
#' @keywords  diapause emergency ,spring
#' 
#' @importFrom truncnorm dtruncnorm 
#' @export

diapause_emergency=function(dates,
                            lon=7.216,
                            lat=43.658,
                            CPPmin=11.5,
                            jdmax=110,
                            bounddays=10,
                            varjd=10,
                            truncated=F) {
   set.seed(0.5)
  
   jd=as.numeric(format(as.Date(dates),"%j"))
   photo=rAedesSim::day_length(dates,lon,lat)$daylength
   
  jdmax=ifelse(jdmax>max(jd),max(jd),jdmax)
  jdmax=ifelse(jdmax>min(jd),jdmax,min(jd)) 
  
  
   datemaxid=which(jd == jdmax)
   jd_ini=jd[which(photo>CPPmin)][1]
   jd_ini=ifelse(jd_ini > jdmax,jdmax,jd_ini)
   mid_jd=(jd_ini+jdmax)/2
   
   # d_dis_day=dnorm((jd_ini-bounddays):(jdmax+bounddays), mean = mid_jd, sd=varjd)
  
   d_dis_day=dnorm((jd_ini-bounddays):(jdmax), mean = mid_jd-bounddays/2, sd=varjd)
   
   scaled_dis_day=d_dis_day/sum(d_dis_day);
   
   meanpp=unlist(tapply(scaled_dis_day, as.factor((jd_ini-bounddays):(jdmax)), function(x) c(sum(x))))
   
                 
   if ( truncated == T) {d_dis_day=truncnorm::dtruncnorm(as.numeric(jd_ini:jdmax),jd_ini,jdmax,mean=mid_jd,sd=varjd);
                         scaled_dis_day=d_dis_day/sum(d_dis_day);
                         meanpp=unlist(tapply(scaled_dis_day, as.numeric(jd_ini:jdmax), function(x) c(sum(x))))
                         
   }
   
   meanpp_df=data.frame(jd=as.numeric(names(meanpp)),emergency=meanpp)
   temp=data.frame(jd=jd,photoperiod=photo)
   final_res=merge(temp,meanpp_df,all.x=T)
   final_res$emergency[which(is.na(final_res$emergency))]=0
   return(final_res)
}

