#' biodata
#'
#' @description  biodata is a function to instantiate a rAedesSim S3 object that contains to describes data collected during a field monitoring and ready to use to perform simulation verification.
#' @param ID char character ID of monitoring. 
#' @param parameter_definition character Name or description of observative parameter collected. Default are Eggs.
#' @param units character Indicate the measure unity of the paraemter observed.  
#' @param location character Name of site of observations. 
#' @param instrument character Name of engine to provide observation
#' @param network  character If data is collected inside an observative network.
#' @param common_name character Data tipology - Observation - Sensor monitoring - Simulation
#' @param scientific_name character Data tipology - Observation - Sensor monitoring - Simulation
#' @param phenology character Data tipology - Observation - Sensor monitoring - Simulation
#' @param obs_standard logical States if data belong to a phenology Data standard class es SYNOP 
#' @param sourcedata Matrix or Data.frame or ascii file of raw data.
#' @param data_provider character Institution / Private data manager.
#' @param data_maintaner character maintainer's  name or contact or its contact . 
#' @param data_licence character Licence of data.
#' @param sourcedata data.frame Or matrix or ascii file of raw data. Field required are dates and parameter.
#' @param field_delimiter character field delimiter of file. Default is comma ",". 				                   
#' @param lat numeric latitude coordinates of  where data are collected.
#' @param lon numeric longitude coordinates of  where data are collected.
#' @param timezone character Timezone. Default is Europe/Rome.
#' @param CRS character Projection of coordinate in proj4 format.
#' @param elevation numeric Elevation of the Default is 40.
#' @param geonotes character Annotations in regard to the contest of observation.
#' @param urban logical Flag indicating if data belong to urban area. Automatic detection is done in according to UMZ Urban Morfological Zone EAA http://database.espon.eu. Default is TRUE.
#' @param nasa_radiance numeric Night radiance value. Is a proxy for urbanity.
#' @param census_level numeric Population density estimation by EEA data sources.
#' @return object  Return a rAedesim object of class biodata  
#' @author  Istituto di Biometeorologia Firenze Alfonso Crisci \email{a.crisci@@ibimet.cnr.it} ASL 2 LUCCA  Marco Selmi \email{marco.selmi@@uslnordovest.toscana.it }
#' @keywords  data, modeling 
#'
#' @import sp
#' 
#' @export

biodata <- function(ID="",
	            parameter_definition = "Eggs", 
                    units = "Integer count", 
		    location="", 
                    instrument = "Trappola REDLAV ITA",
                    network="Redlav",
                    common_name="Zanzara Tigre",
                    scientific_name="Aedes Albopictus",
                    phenology="Eggs Hatching",
                    obs_standard=TRUE,
		    sourcedata=NULL,
                    data_provider="ASL Lucca",
                    data_maintaner="Marco Selmi - m.selmi@usl2.toscana.it",
                    data_licence="",
                    field_delimiter=",",
                    lat,
                    lon,
		    timezone="Europe/Rome",
                    CRS="+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs",
                    geonotes="",
		    urban=TRUE,
                    nasa_radiance=NULL,
                    census_level=NULL)
                   {	
       		   
				  #################################################################################################################
				  # location conversion
				   
				   epgs4386="+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs";
				   splocation=data.frame(lon=lon,lat=lat)
				   coordinates(splocation) <- c('lon','lat')
				   proj4string(splocation) <- CRS
				   
				   #################################################################################################################
				   newsp <-SpatialPointsDataFrame(splocation,data.frame(Elevation=elevation))
                    
				   if (CRS != "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
				      { newsp <-spTransform(splocation,CRS(epgs4386))
					      newsp <-SpatialPointsDataFrame(newsp,data.frame(Elevation=elevation))
          	  }
				     
				    lon_geo=as.numeric(coordinates(newsp))[1]
				    lat_geo=as.numeric(coordinates(newsp))[2]
				   
				   #################################################################################################################
				   	
					
				   if ( is.null(sourcedata)) { stop( "To build meteodata valid object almost a data source is required.\n 
                                   			              Date of day in YYYY-MM-DD format (dates),parameter field. 
                                                                      Checking header data names to avoid variable error.") 
					                       } 
				   
				   #################################################################################################################
				  
                                    if (is.data.frame(sourcedata) || is.matrix(sourcedata))
				                                { filedata=as.data.frame(sourcedata)
					                        }
                                    else 
				                                {
                                           filedata=read.table(sourcedata, header=TRUE, sep=field_delimiter,na.strings="NA", dec=".", strip.white=TRUE)
				                                }
					
				    ts_zoo=NULL;
				    ts_zoo=try(as.xts(zoo(filedata$parameter,as.Date(as.character(filedata$dates,tz=timezone)))));

				    if ( is.null(ts_zoo))
				       { stop( "Timeseries creation was invalid! Check data and dates in data sources")
				       };
				  
   				  #################################################################################################################
				  						 
				   object <- list(ID = ID,
				                  parameter_definition =  parameter_definition,
				                  unity_measure = units,
				                  location = location, 
				                  instrument = instrument,
				                  network = network,
				                  common_name = common_name,
				                  scientific_name = scientific_name,
				                  phenology = phenology,
				                  sourcedata=sourcedata,
				                  obs_standard = obs_standard,
				                  data_provider = data_provider,
				                  data_maintaner = data_maintaner,
				                  data_licence = data_licence,      
				                  lat = lat_geo,
				                  lon = lon_geo,
				                  elevation = elevation,
						  timezone = "Europe/Rome",
				                  CRS = epgs4386,
				                  geonotes = geonotes,
				                  urban = urban,
				                  nasa_radiance = nasa_radiance,
				                  feature_population = census_level,
				                  ts_data = ts_zoo,
				                  sp_obj=newsp
				                  );
                 attr(object,"ID")<-"ID of monitoring/data collection"
                 attr(object,"parameter_definition")<-"Name of site of observations" 
                 attr(object,"unity_measure")<-"Measure/Observation units of parameter" 
                 attr(object,"location")<-"Name of site of observations" 
                 attr(object,"instrument")<-"Name of engine to provide observation"
                 attr(object,"network")<-"Name of observative network of data"
                 attr(object,"common_name")<-"Indicate the common name of animal  simulated"
                 attr(object,"scientific_name")<-"Indicate the scientific name of animal of simulated"
                 attr(object,"phenology")<-"Data tipology - Observation - Sensor monitoring - Simulation"
                 attr(object,"sourcedata")<-"Indicate the provenience of raw data" 
                 attr(object,"obs_standard")<-"If data belong to a phenology data standard class" 
                 attr(object,"data_provider")<-"Institution / Private data manager"
                 attr(object,"data_maintaner")<-"maintainer's  name or contact or its contact"
                 attr(object,"data_licence")<-"Licence of data"                   
                 attr(object,"lat")<-"latitude coordinates of  site where data are collected"
                 attr(object,"lon")<-"longitude coordinates of  site where data are collected"
	         attr(object,"timezone")<-"longitude coordinates of  where data are collected"
                 attr(object,"CRS")<-"Projection of coordinate in proj4 format"
                 attr(object,"elevation")<-"Elevation in meters"
	         attr(object,"timezone")<-"longitude coordinates of  where data are collected"
                 attr(object,"geonotes ")<-"Annotations in regard to the contest of observation"
                 attr(object,"urban")<-"Flag indicating if data belong to urban area"
                 attr(object,"nasa_radiance ")<-"NASA Night radiance value. Is a proxy for urbanization"
                 attr(object,"feature_population")<-"Population density estimation in the neighbour area"
                 attr(object,"ts_data")<-"Time series object of data"
                 attr(object,"sp_obj")<-"SpatialPointDataFrame"   
                 class(object) <-"biodata"
                 return(object)
}  
					  
