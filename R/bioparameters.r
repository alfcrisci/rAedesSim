#' bioparameters
#'
#' @description Bioparameters is a function to instantiate an S3 object  containing biological parameters used in simulation for rAedesSim .
#' @param numeric alfa_l Coefficient of competition between individuals at larval stage.
#' @param numeric alfa_a  Coefficient of competition between individuals at adult stage.
#' @param numeric exp_lar Exponential Gompertz parameter for competition modulation. Default=1.
#' @param numeric l_density Expected larval density. Default is 70.
#' @param numeric sexratio   Sex ratio of population. Default is 0.5.
#' @param numeric egn  Mean eggs . Default is 63.
#' @param numeric inib Inibition rate parameter ( 0-1) Default is 0.
#' @param string sspp Name of the species. Default is "Albopictus".
#' @param string genus_sspp  Name of the genus of the species. Default is "Aedes".
#' @param string order_sspp  Name of the order of the species. Default is "Diptera".
#' @param string geo_area  Name of geographic area. Default is "Tuscany".
#' @param string name_location  Name of location.
#' @param numeric name_location mean elevation  of container location. Default is 100.
#' @return S3 object Bioparameters object
#' @author  Istituto di Biometeorologia Firenze Italy  Alfonso crisci \email{a.crisci@@ibimet.cnr.it} ASL 2 LUCCA Marco Selmi \email{m.selmi@@usl2.toscana.it} 
#' @keywords  biparameters
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
