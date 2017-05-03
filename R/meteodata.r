#' meteodata
#'
#' @description meteodata  is the function to store raw weather information in a S3 object in order to performs the simulation.
#'
#' @param station_name character name of station or location
#' @param network  character name or project or network which data corresponds.
#' @param data_type character weather data class: forecast,observation, sensor micrometeorological data, physical parameter simulation
#' @param standard character weather data standard class. Default is SYREP. 
#' @param data_provider character institution -Private data manager.
#' @param data_maintaner character maintainer's  name or contact or its contact . 
#' @param data_licence character Licence of data.    
#' @param lat numeric  latitude coordinate.
#' @param lon numeric  longitude coordinate.
#' @param CRS character Projection of coordinates in proj4 format.
#' @param elevation numeric Elevation corresponding to data collected. Default is 40.
#' @param timeformat  numeric Time period of data aggregation. Default is "Daily".
#' @param sourcedata  matrix or data.frame or ascii file of raw data. Frame must be daily with at least 6 fields. dates ( in YYYY-MM-DD format ),tmed,tmax,tmin,rhum prec.
#' @param field_delimiter chararter field delimiter of file. Default is comma ",".
#' @param date_format character Format of the dates in raw data. Default is YYYY-MM-DD. 
#' @param timeseries  object xts R Timeseries of data
#' @return meteodata
#' @author  Istituto di Biometeorologia Firenze Italy  Alfonso crisci \email{a.crisci@@ibimet.cnr.it} ASL 2 LUCCA Marco Selmi \email{marco.selmi@@uslnordovest.toscana.it } 
#' @keywords  metorological data 
#' @import sp
#' 
#' @export


meteodata<-function(station_name="Pisa San Giusto",
                    network="Aeronautica Militare",
                    data_type="Simulation",
                    standard="SYREP",
                    data_provider="IBIMET CNR",	
                    data_maintainer="",
                    data_licence="",
                    date_format="YMD",
                    lat=43.0,	
                    lon=11.0,
                    elevation=40,
                    timeformat="daily",
                    sourcedata=NULL,
                    field_delimiter=",",
                    timeseries=NULL,
                    CRS="+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
		    )
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
				      {
					   newsp <- spTransform(splocation,CRS(epgs4386))
					   newsp <-SpatialPointsDataFrame(newsp,data.frame(Elevation=elevation))
                    
					  }
				     
				    lon_geo=as.numeric(coordinates(newsp))[1]
				    lat_geo=as.numeric(coordinates(newsp))[2]
				    
					
				    if ( is.null(sourcedata)) { 
				                               stop( "To build meteodata valid object a data source is required in one of this two ways: \n
				                                A) data.frame 
				                                B) csv ascii  formatted file.\n
                                                                Format must be daily with fields : 
                                                                [dates] Date of day better in YYYY-MM-DD format . \n 
                                                                [tmed]  Mean temperature. \n
                                                                [tmax]  Maximum temperature.\n 
                                                                [tmin]  Minimum temperature.\n 
                                                                [rhum]  Mean relative humidity.\n 
                                                                [prec]  Cumulated rainfall.\n
                                                                Header name are indicated to avoid errors. 
                                                                Mean temperature and minimum temperature fields are required."
						                            ) 
							                     	} 
				   
				   #################################################################################################################
				  
				    if (is.data.frame(sourcedata) || is.matrix(sourcedata))
				      { filemeteo=as.data.frame(sourcedata)
				      }
				    else 
				      {
				        filemeteo=read.table(sourcedata, header=TRUE, sep=field_delimiter,na.strings="NA", dec=".", strip.white=TRUE)
				      }
                   
				   #################################################################################################################
				   
				   if ( length(grep("tmed",names(filemeteo)))==0 || length( grep("dates",names(filemeteo)))== 0 || length( grep("rhum",names(filemeteo)))== 0)
                                          { stop( "Name field [dates] [tmed] [rhum] in datasource are requested. Please change the names of datafield including these ones.")
										  }				  

				   
				   if ( length(grep("anno",names(filemeteo)))>0 && length(grep("mese",names(filemeteo)))>0 && length(grep("mese",names(filemeteo))>0 ))
                                          { names(filemeteo)<-gsub("anno","year",names(filemeteo))
					    names(filemeteo)<-gsub("mese","month",names(filemeteo))
					    names(filemeteo)<-gsub("giorno","day",names(filemeteo))
					    names(filemeteo)<-gsub("urel","rhum",names(filemeteo))
					    names(filemeteo)<-gsub("data","dates",names(filemeteo))
										  
					   }	
				   
				   if ( !length(grep("year",names(filemeteo)))>0 && !length(grep("month",names(filemeteo)))>0 && length(grep("day",names(filemeteo))>0 ))
					  { filemeteo$year<-year(filemeteo$dates)
					    filemeteo$month<-month(filemeteo$dates)
					    filemeteo$day<-day(filemeteo$dates)
					  }  
				   #################################################################################################################
				   
				   filemeteo$dates=lubridate::ymd(filemeteo$dates)
                                 
	         if ( date_format == "DMY") {filemeteo$dates=lubridate::dmy(filemeteo$dates)}
	
           if ( date_format == "MDY") {filemeteo$dates=lubridate::mdy(filemeteo$dates)}
                   
				   #################################################################################################################
				   rownames(filemeteo)<-1:nrow(filemeteo)
				
				        
				   period=as.numeric(max(as.Date(filemeteo$dates))-min(as.Date(filemeteo$dates)))
				    
				   length_data_ini=nrow(filemeteo)
				   row_na=nrow(filemeteo[!complete.cases(filemeteo),])
				   perc_missing=(row_na/nrow(filemeteo))*100
				   
				   if ( period > length_data_ini-1) { 
				                          fill_dates=data.frame(dates=seq(min(as.Date(filemeteo$dates)),max(as.Date(filemeteo$dates)),1));
				                          filemeteo=merge(fill_dates,filemeteo,by=c("dates"),all = T)
				                          }
				   
				   filemeteo$daylength = day_length(as.Date(filemeteo$dates),lon_geo,lat_geo)$daylength;
				   
				   #################################################################################################################
				  
				   if (period < (length_data_ini-1)) { 
				     stop( "Data source is not correcly time indexed! Suspect date duplicates.")
				   }				  

				   
				   continuity=ifelse(period==length_data_ini-1,TRUE,FALSE)
				   
				   #################################################################################################################
				  
				   variables=names(filemeteo)[grep("tmed|tmax|tmin|rhum|prec",names(filemeteo))]
				   
				   if ( is.null(timeseries) )
				   {
				   ts_zoo=try(as.xts(zoo(filemeteo[variables],as.Date(as.character(filemeteo$dates)))))
				   
				   }
				   
				   if ( !exists("ts_zoo"))
				   {
				    warning( "Timeseries creation invalid! Check data and dates in data sources")
				   }
				  
   				  #################################################################################################################
				  						 
				   object <- list(station_name=station_name,
				                  network=network,
				                  data_type=data_type, 
				                  standard=standard,
				                  data_provider=data_provider,
				                  data_maintainer=data_maintainer,
				                  data_licence=data_licence,
				                  lat=lat_geo,	
				                  lon=lon_geo,
				                  CRS=epgs4386,
				                  elevation=elevation,
				                  timeformat=timeformat,
				                  tmed=filemeteo$tmed,
				                  tmax=filemeteo$tmax,
				                  tmin=filemeteo$tmin,
				                  urel=filemeteo$rhum,
				                  prec=filemeteo$prec,
				                  dates=filemeteo$dates,
				                  length_data_ini=length_data_ini,
				                  ndays=length(filemeteo$dates),
				                  daylenght=filemeteo$daylength,
				                  continuity=continuity,
				                  perc_missing_data=perc_missing,
				                  timeseries=ts_zoo,
				                  sp_obj=newsp)
 
                 attr(object,"station_name") <- "Name of station"
                 attr(object,"network") <- "Network of station"
                 attr(object,"data_type")<-"Station Type"
                 attr(object,"standard")<-"Data class or relative reference standard"
                 attr(object,"data_provider")<-"Istitution or Private data manager"
                 attr(object,"data_maintainer")<-"Name or contact of data maintainer"
                 attr(object,"data_licence")<-"Licence of data"
                 attr(object,"lat")<-"Latitude in decimal degrees. Datum WGS 84"
                 attr(object,"lon")<-"Longitude in decimal degrees . Datum WGS 84"	
                 attr(object,"CRS")<-"Projection used for the coordinate in proj4 string format"	
                 attr(object,"elevation")<-"Elevation (m)"	
                 attr(object,"timeformat")<-"Period of aggregation"
                 attr(object,"tmed")<-"Mean daily temperature"
                 attr(object,"tmax")<-"Maximum daily temperature"
                 attr(object,"tmin")<-"Minimum daily temperature"
                 attr(object,"urel")<-"Relative humidity daily average"
                 attr(object,"prec")<-"Rainfall cumulated in a day"
                 attr(object,"dates")<-"Dates of meteorological data matrix"
                 attr(object,"length_data_ini")<-"Initial data length of raw data"
                 attr(object,"ndays")<-"Number of days"
                 attr(object,"daylenght")<-"Day length for each date"
                 attr(object,"continuity")<-"If the temporal continuity of raw data is detected"
                 attr(object,"perc_missing_data")<-"Percentage of missing values"
                 attr(object,"timeseries")<-"Data timeseries as R xts object"
                 attr(object,"sp_obj")<-"SpatialPointDataFrame"  
                 class(object) <-"meteodata"
               return(object)
}
