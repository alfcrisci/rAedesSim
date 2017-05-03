#' weigthdry
#'
#' @description Calculate a exponential weigth to perform the impact of dry spell in function to its length. 
#' It is used to evaluate the impact of drought on mosquito population. A rain treshshold is considered. Default value is 5 mm.
#' 
#' @param drylength numeric vector: length of dry spell as number of consecutive dry days.    
#' @param tresh nuemric: threshold in days 
#' @return weigthdry
#' @seealso \code{"drydaycons"}
#' @author  Istituto di Biometeorologia Firenze Italy  Alfonso crisci \email{a.crisci@@ibimet.cnr.it} ASL LUCCA Marco Selmi \email{marco.selmi@@uslnordovest.toscana.it }
#' @keywords  drought, dry spell
#'
#' 
#' 
#' @export
						   
weigthdry<-function(drylength,tresh=5) {
                                       res=NULL;
                                       res=round(1-exp(drylength/tresh-drylength),digits=2) 
                                       return(res)
				       }
