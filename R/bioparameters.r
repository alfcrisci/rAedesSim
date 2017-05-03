#' bioparameters
#'
#' @description Bioparameters is a function to instantiate an S3 object  containing biological parameters used in simulation for rAedesSim .
#' @param alfa_l numeric Coefficient of competition between individuals at larval stage.
#' @param alfa_a numeric  Coefficient of competition between individuals at adult stage.
#' @param exp_lar numeric Exponential Gompertz parameter for competition modulation. Default=1.
#' @param l_density numeric Expected larval density. Default is 70.
#' @param sexratio numeric   Sex ratio of population. Default is 0.5.
#' @param egn  numeric the Mean eggs during ovideposition cycle. Default is 63.
#' @param inib numeric Inibition rate parameter ( 0-1) Default is 0.
#' @param sspp  character Name of the species. Default is "Albopictus".
#' @param genus_sspp  character  Name of the genus of the species. Default is "Aedes".
#' @param order_sspp  character  Name of the order of the species. Default is "Diptera".
#' @param geo_area   character Name of geographic area. Default is "Tuscany".
#' @param name_location  character  Name of location.
#' @param name_location numeric mean elevation  of container location. Default is missing.
#' @return S3 object Bioparameters object
#' @author  Istituto di Biometeorologia Firenze Italy  Alfonso crisci \email{a.crisci@@ibimet.cnr.it} ASL LUCCA Marco Selmi \email{marco.selmi@@uslnordovest.toscana.it }
#' @keywords  bioparameters
#'
#' 
#' 
#' @export

bioparameters <- function(alfa_l=1,
                          alfa_a=0,
                          exp_lar=1,
                          l_density=100,
                          sex_ratio=0.5,
                          egn=16,
                          ef=0.83,
                          inib=0, 					 
                          sspp="Albopictus",
                          genus_sspp="Aedes",
                          order_sspp="Diptera",
                          geo_area="Tuscany",			    
                          name_location=NA
                          ) 
							{
  object <- list(alfa_l=alfa_l,
                 alfa_a=alfa_a,
                 exp_lar=exp_lar,
                 l_density=l_density,
                 sex_ratio=sex_ratio,
                 egn=egn,
                 ef=ef,
                 inib=inib,					
                 sspp=sspp,
                 genus_sspp=genus_sspp,
                 order_sspp=order_sspp,
                 geo_area=geo_area,			    
                 name_location=name_location
                 );
				 
  attr(object,"alfa_l") <- "Environmental Parameter of larvae competition "
  attr(object,"alfa_a") <- "Environmental Parameter of adult competition"
  attr(object,"l_density") <- "Critical larval density N/cm^3"
  attr(object,"sex_ratio") <- "Ratio between sex"
  attr(object,"egn")<-"Individual mean hatch rate"
  attr(object,"ef")<-"Pupal effective Initibition factor"
  attr(object,"inib")<-"Larval Inibibition factor"
  attr(object,"sspp")<-"Name of species"
  attr(object,"genus")<-"Genus of species"
  attr(object,"order_sspp")<-"Order of species"
  attr(object,"geo_area")<-"Geographical area names"			    
  attr(object,"name_location")<-"Name of site"
  class(object) <- "bioparameters"
  return(object)
}
