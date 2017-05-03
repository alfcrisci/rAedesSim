#' drydaycons
#'
#' @description Calculate the sequence of consecutive dry days in function to a certain  treshshold.
#' 
#' @param prec numeric  vector of daily cumulated rainfall in mm.    
#' @param tresh numeric Rainfall threshold in mm
#' @return numeric Numerical vector of cumulate sequence

#' @seealso \code{weigthdry}
#' @author  Istituto di Biometeorologia Firenze Italy  Alfonso crisci \email{a.crisci@@ibimet.cnr.it} ASL 2 LUCCA Marco Selmi \email{marco.selmi@@uslnordovest.toscana.it }
#' @keywords  dry days
#' 
#' 
#' @export





drydaycons<-function(prec,tresh=1) {
                                   prec[prec<tresh]<-0
                                   prec[prec>=tresh]<-1
                                   prec.rle=rle(prec);
                                   res=NULL;
                                   for (k in 1:length(prec.rle$lengths)) { res=append(res,1:prec.rle$lengths[k])}
                                   res[which(prec>=1)]=0
                                   return(res)
		                   }
