#' biocontainer
#'
#' @description biocontainer is the function able to instantiate a S3 object used in rAedesSim simulation. The one is referred to the main features of the container where eggs are layed and larvae/pupae growth.  
#' 
#' @param type character  Container type.  The default trap used in REDLAV project is indicated: volume 0.75 750 cm^3 lt, heigth 13 cm, Area 50 cm^2.Describe the tipology of cointainer considered. For example (i) cointaner (ii) manhole (iii) others. 
#' @param nrecipients numeric  Represent the numerosity of container used in simulation and correspond  to BS the number of breeding sites . Default is 50.
#' @param container_shape numeric Numeric level indicates if the surface container the shape . Parallelepiped plate (1)  Cilynder plate(2) Unspecified (3). Default is 3. 
#' @param capacity numeric  The maximum capacity of water in trap/container. Default ( 1 lt 1000 cm^3). 
#' @param initVol numeric  The unitarian initial volume of water in trap/container. Default ( 0.75 lt 750 ml).
#' @param frac_initVol numeric Fraction of volume present in trap/location. Default is 0.
#' @param surface_area  numeric The water-air surface of container in squared centimeter (cmq).
#' @param underplate  logical  Presence of underplate. Default is FALSE.
#' @param frac_underplate numeric Fraction of underplate engaged (0 - 1) Default is 0.
#' @param daily_evap_loss numeric Volume of water losses daily due to evaporation. Default is 0.
#' @param frac_covering numeric Fraction of container covering  (0 - 1) Default is 0.
#' @param lat numeric Mean geographical latitude of container location. Default is 43.5. 
#' @param lon numeric Mean geographical longitude of container location. Default is 11.27.
#' @param elevation numeric  Mean elevation  of container location. Default is 100.
#' @param CRS  character Projection of coordinate in proj4 format.
#' @param timezone sring Timezone. Default is "Europe/Rome". 
#' @param sourcedata dataframe  Or ascii file consisting in   Water temperature measurements. Fields required are dates ( YYYY-MM-DD) format, Daily mean of water tmedwater.
#' @param meteodata object  rAedesim  object with the data useful for model water calibration.
#' @param watermodel logical Statistical model oject ( gam - lm - glm).Default is NULL.
#' @param watermodel_fit logical Fit a statistical model oject ( gam - lm - glm).Default is FALSE.
#' @param model_type character  Statistical model class.
#' @param date_format character Format of the dates in raw data. Default is YMD. 
#' @param field_delimiter character field delimiter of file. Default is comma ",". 
#' @param timeseries object xts R Timeseries of data 
#' @param ID  character ID of container's set.
#' @param site_name character Name of place or location.
#' @return biocontainer 
#' @author  Istituto di Biometeorologia Firenze Italy  Alfonso crisci \email{a.crisci@@ibimet.cnr.it} ASL 2 LUCCA Marco Selmi \email{marco.selmi@@uslnordovest.toscana.it }
#' @keywords  container
#'
#' @import sp
#' @import mgcv
#' @import xts
#' @import lubridate
#' @export


biocontainer<-function( type="Trap  REDLAV ITA",
			                     nrecipients=1,
			                     container_shape=3,
			                     capacity=1000,
			                     initVol=750,
			                     frac_initVol=0.75,
			                     surface_area=50,
			                     underplate=FALSE,
			                     frac_underplate=0,
			                     frac_covering=0,
			                     lat=43.5,
			                     lon=11.27,
			                     elevation=100,
			                     timezone="Europe/Rome",
			                     CRS="+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs",
			                     sourcedata=NULL,
			                     meteodata=NULL,
			                     watermodel=NULL,
			                     watermodel_fit=FALSE,
			                     model_type="gam",
			                     date_format="YMD",
			                     field_delimiter=",", 
			                     ID=c("NA") ,
			                     site_name=c("NA")			   
			                     ) {
                
				
				#################################################################################################################
				# coordinates conversion if coordinates are not geographic.
				
				epgs4386="+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs";
				splocation=data.frame(lon=lon,lat=lat)
				coordinates(splocation) <- c('lon','lat')
				proj4string(splocation) <- CRS
				   
				#################################################################################################################
				newsp <-SpatialPointsDataFrame(splocation,data.frame(Elevation=elevation))
                    
				if (CRS != "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
				      {
					   newsp <-spTransform(splocation,CRS(epgs4386))
					   newsp <-SpatialPointsDataFrame(newsp,data.frame(Elevation=elevation))
                    
					  }
				     
				lon_geo=as.numeric(coordinates(newsp))[1]
				lat_geo=as.numeric(coordinates(newsp))[2]
				
				#################################################################################################################
				
				if ( (is.null(meteodata) || is.null(sourcedata) ) &&  watermodel_fit ==TRUE)
                                          { warning( "To fit a watermodel a meteodata object and a datasource of water temperatures are needed.")
					  }	
				
			  
				#################################################################################################################
				# Instance an objects
				
				full_ts=NULL;
				
				#################################################################################################################
				# Check row index
				
				if  ( watermodel_fit == TRUE && is.null(watermodel)) 
				    {
				     if (is.data.frame(sourcedata) || is.matrix(sourcedata))
				      { filedata = as.data.frame(sourcedata)
				      }
                                 else 
				      {
				      filedata = read.table(sourcedata, header=TRUE,sep=field_delimiter,na.strings="NA", dec=".", strip.white=TRUE)
				      }
                   
				      if ( !length(grep("tmedwater",names(filedata))) > 0 || !length( grep("dates",names(filedata))) > 0 )
                                       { stop( "Field dates and tmedwater in datasource are needed.Check header file and change names.")
				       }
				      if ( date_format == "DMY") {filedata$dates=dmy(filedata$dates)};
				      if ( date_format == "MDY") {filedata$dates=mdy(filedata$dates)};
                   
				      rownames(filedata)<-1:nrow(filedata)
				
				     filedata$daylength=day_length(as.Date(as.character(filedata$dates),tz=timezone),lon_geo,lat_geo)$daylength;
				   
				     full_ts = try(as.xts(zoo(data.frame(daylength=filedata$daylength,tmedwater=filedata$tmedwater),as.Date(as.character(filedata$dates)))))
				
				     if   (!exists("full_ts"))
				         {
				          stop( "Timeseries creation invalid! Check data and dates in data sources")
			           }
				
				
				###################################################################################################################
				# Prepare 
			   
				    full_ts=merge.xts(full_ts, meteodata$timeseries)
				    full_ts_df=na.omit(as.data.frame(full_ts))
				    rownames(full_ts_df) <- NULL;
				
				    if (model_type == "lin") 
				       {
				        watermodel=try(lm(tmedwater~daylength+tmed+tmin,data=full_ts_df))
					};
				
				    if (model_type == "gam") 
				       {
				         if ( nrow(full_ts_df) <45)
						    {
				                     stop( "\n\n Gam modeling require more data to avoid overfitting. \nTry with linear models (lin) in argument." )
						    }
					    
				            watermodel = try(gam(tmedwater~s(daylength)+s(tmed)+s(tmin),data=full_ts_df))
				        }	  
										  
                                       }						  
       #################################################################################################################
				
				object <- list(type=type,
				               nrecipients = nrecipients,
				               container_shape = container_shape,
				               surface_area = surface_area,
				               underplate = underplate,
				               frac_underplate = frac_underplate,
				               frac_covering = frac_covering,
				               initVol = initVol,
				               frac_initVol = frac_initVol,
				               capacity = capacity,
				               pooled_volume_init = initVol*frac_initVol*nrecipients,
				               pooled_volume_current = initVol*frac_initVol*nrecipients,
				               lat = lat_geo,
				               lon = lon_geo,
				               CRS = epgs4386,
				               timezone=timezone,
				               elevation = elevation,
				               meteodata = meteodata,
				               watermodel = watermodel,
				               model_type = model_type,
				               date_model = Sys.Date(),
				               timeseries = full_ts,
				               sp_obj = newsp,
				               ID = ID,
				               site_name = site_name
				               );
				 
				               attr(object,"type") <- "Describe the tipology of container"
				               attr(object,"nrecipients") <- "Breeding sites numerosity"
				               attr(object,"container_shape") <- "Container's shape. Parallelepiped plate (1) Cilynder plate(2) Unspecified (3)."
				               attr(object,"surface_area") <- "Surface at water-air interface of container (cm^2) "
				               attr(object,"underplate") <- " Underplate presence"
				               attr(object,"frac_underplate") <- "Underplate fraction occupated by plate"
				               attr(object,"frac_covering") <- "Fraction of initial pooled volume(0 if open, 0.9 full covered)"
				               attr(object,"initVol") <- "Fraction of initial unitarian volume"
				               attr(object,"frac_initVol") <- "Fraction of initial pooled volume"
				               attr(object,"capacity")<-"Maximum capacity of water in trap/container (cm^3)"
				               attr(object,"pooled_volume_current")<-"Current Pooled volume of container's set (cm^3)"
				               attr(object,"pooled_volume_init")<-"Pooled volume of container's set (cm^3)"
				               attr(object,"lat")<-"latitude coordinates of  biocointainer"
				               attr(object,"lon")<-"longitude coordinates of  biocointainer"
				               attr(object,"CRS")<-"Projection used for the coordinate in proj4 string format"	
				               attr(object,"elevation")<-"Mean elevation"
				               attr(object,"metameteo")<-"Rbiomsim object to fit model"
				               attr(object,"watermodel")<-"Model fitted from raw data of water temperature measurements"
				               attr(object,"watermodel_fit")<-"If water model was fitted"
				               attr(object,"model_type")<-"Statistical model class used for data fitting"
				               attr(object,"date_model")<-"Model creation date" 
				               attr(object,"timeseries")<-"Timeseries object obtained from data"
				               attr(object,"timezone")<-"Timezone"
				               attr(object,"sp_obj")<-"SpatialPointDataFrame" 
				               attr(object,"ID")<-"ID label of container set"
				               attr(object,"site_name")<-"Name of sites"
	                     class(object) <-"biocontainer"
				               return(object)
}
