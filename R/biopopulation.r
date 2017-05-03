#' biopopulation
#'
#' @description S3 class to describe biomodelday object to be used in rAedesim model. 
#' 
#' @param eggs numeric  Eggs Count.  
#' @param larvae numeric Larvae Count.  
#' @param pupae numeric Pupae Count. 
#' @param adults numeric Mosquito Female Actual Count. 
#' @param eggs_diap numeric Diapausant Eggs Count.
#' @param ID character ID label Population. 
#' @return object Return a rAedesim biopopulation object.
#' @author  Istituto di Biometeorologia Firenze Italy  Alfonso crisci \email{a.crisci@@ibimet.cnr.it} ASL LUCCA Marco Selmi \email{marco.selmi@@uslnordovest.toscana.it }
#'@keywords  container
#'
#' 
#' 
#' @export

biopopulation<- function( eggs=10000,
                          larvae=0,	
                          pupae=0,
                          adults=0,
                          eggs_diap=10,
                          ID="NA"
			                    ) 
                         {object <- list(
                          eggs=eggs,
                          larvae=larvae,	
                          pupae=pupae,
                          adults=adults,
                          eggs_diap=eggs_diap,
                          ID=ID)                          
                          class(object) <- "biopopulation"
                          attr(object,"eggs") <- "Eggs Count"
                          attr(object,"larvae") <- "Larvae Count"
                          attr(object,"pupae") <- "Pupae Count"
                          attr(object,"adults") <- "Female Actual Count"
                          attr(object,"eggs_diap") <- "Diapausant Eggs Count"
                          attr(object,"ID")<-"ID population"
                          class(object) <- "biopopulation"
                          return(object)
}