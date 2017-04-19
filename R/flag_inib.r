#' flag_inib
#'
#' @description Calculate if in the container exists a egg hatch inibition due to larval density.  
#' Needs volume water average for unitarian container.Return a logical value ( TRUE/FALSE) if inibition of hatch occurs.
#' 
#' @param N_larvae Number of larvae.
#' @param vol_water Default is 1000.  
#' @param critical_density Larval critical density to establish if hatch is inibited.
#' @return . 
#' @references Janice S. Edgerly - Michelle S. Willey - Todd P. Livdahl - Ecological Entomology - Vol. 18 - Issue 2 - 1993 - pp. 123-128  \url{http://www.clarku.edu/faculty/tlivdahl/publications/livdahl-edgerly87.pdf}
#' @author  Istituto di Biometeorologia Firenze Italy  Alfonso crisci \email{a.crisci@@ibimet.cnr.it} ASL 2 LUCCA Marco Selmi \email{m.selmi@@usl2.toscana.it} 
#' @keywords  Egg Hatch inibition
#'
#' 
#' 
#' @export
 
  
flag_inib<-function(N_larvae,vol_water=1000,critical_density=72){
	 flag_inib_l=0;
	if ( is.na(N_larvae)) {return(flag_inib_l)};
	if( critical_density <= N_larvae/(vol_water/1000))
                { flag_inib_l=1;
                }  
	return(flag_inib_l)
	}
