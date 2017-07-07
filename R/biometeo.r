#' biometeo
#'
#' @description Biometeo is the function to instantiate a S3 object related to environmental parameter that to be used in biomodel function.
#' It require firstly a meteodata object and biocontainer object as arguments.
#'
#' @param meteodata   object rAedesim meteodata object.
#' @param biocontainer  object rAedesim biocontainer object
#' @param startday  integer Index of date when  simulation starts. Default is 1.
#' @param deltatmax numeric Mean Error Bias considered for maximum air temperature. Default is 0.
#' @param deltatmin numeric Mean Error Bias  considered for minimum air temperature. Default is 0.
#' @param deltatmed numeric Mean Error Bias  considered for mean air temperature. Default is 0.
#' @param tresh_rain_dry numeric Rain threshold for effective precipitation. Default is 4.
#' @param weigth_k numeric Weighting Paraemter to take into account dry spell impact See also \code{weigthdry}. Default is 5.
#' @param datemax integer Julian day of maxim date of spring diapause awekaning . Default is 110 4th April.
#' @param bounddays numeric Days of incertainty for diapause emergency period bounds. Default is 10.
#' @param varjd numeric Days of variance for  diapause emergency courbes. Default is 10.
#' @param truncated logical Modeling with truncated gaussian respecting the limits. Default is FALSE.
#' @param timezone character Timezone. Default is Europe/Rome.
#' @return  Return a biometeo object.
#' @author  Istituto di Biometeorologia Firenze Italy  Alfonso crisci \email{a.crisci@@ibimet.cnr.it} ASL LUCCA Marco Selmi \email{marco.selmi@@uslnordovest.toscana.it }
#' @keywords  container
#'
#'
#'
#' @export


biometeo = function (meteo,
                     i_biocontainer,
                     startday = 1,
                     deltatmax = 0,
                     deltatmin = 0,
                     deltatmed = 0,
                     tresh_rain = 4,
                     weigth_k = 5,
                     CPPmin = 11.5,
                     datemax = 110,
                     bounddays = 10,
                     varjd = 10,
                     truncated = FALSE,
                     timezone = "Europe/Rome")
{
  if (class(meteo) != "meteodata")
  {
    stop("argument must be an meteodata object ")
  }
  
  if (class(i_biocontainer) != "biocontainer")
  {
    stop("argument must be an biocontainer object")
  }
  
  len_ts = as.numeric(length(meteo$dates))
  
  if (len_ts < 1)
  {
    stop("Time series is singular.")
  }
  
  jd = as.numeric(format(as.Date(meteo$dates, tz = timezone), "%j"))
  
  startday_jd = jd[startday]
  
  if (length(jd) <= (startday - 1)) {
    stop("Startday not correct for simulation.")
  }
  
  meteo$tmax = meteo$tmax[startday:len_ts]
  meteo$tmin = meteo$tmin[startday:len_ts]
  meteo$tmed = meteo$tmed[startday:len_ts]
  meteo$urel = meteo$urel[startday:len_ts]
  meteo$prec = meteo$prec[startday:len_ts]
  meteo$dates = meteo$dates[startday:len_ts]
  meteo$daylenght = meteo$daylenght[startday:len_ts]
  
  tmax = meteo$tmax + deltatmax
  tmin = meteo$tmin + deltatmin
  tmed = meteo$tmed + deltatmed
  
  dates = as.Date(as.character(meteo$dates), tz = timezone)
  
  day_len = day_length(dates,
                       as.numeric(i_biocontainer$lon),
                       as.numeric(i_biocontainer$lat))$daylength
  
  new_data = data.frame(daylength = day_len,
                        tmed = tmed,
                        tmin = tmin)
  w_tmed = predict(i_biocontainer$watermodel, newdata = new_data)
  
  id.na <- which(is.na(w_tmed))
  w_tmed[id.na] <- tmed[id.na]
  
  tot_evap = NA
  weight_dry = NA
  prevdrydays = NA
  
  if (!is.null(meteo$urel) ||
      length(which(is.na(as.numeric(meteo$urel)))) == length(as.numeric(meteo$urel)))
  {
    evap <-
      evaporation_params(tmed, w_tmed, meteo$urel, A = i_biocontainer$base_area)
    tot_evap = evap$mevday * i_biocontainer$nrecipients
  }
  
  diapause_emerg = diapause_emergency(
    dates,
    as.numeric(i_biocontainer$lon),
    as.numeric(i_biocontainer$lat),
    CPPmin,
    datemax,
    bounddays,
    varjd,
    truncated
  )
  
  
  diapause_induct = diapause_induction(dates,
                                       as.numeric(i_biocontainer$lon),
                                       as.numeric(i_biocontainer$lat))
  
  if (meteo$perc_missing_data < 20)
  {
    prevdrydays = drydaycons(meteo$prec, tresh_rain)
    weight_dry = weigthdry(prevdrydays, weigth_k)
  }
  
  
  df_data = data.frame(
    tmedwater = w_tmed,
    tmed = tmed,
    diapause_emerg,
    d_induction = diapause_induct,
    lengthday = day_len,
    cdrydays = prevdrydays
  )
  
  rownames(df_data) = dates
  ts_zoo = as.xts(df_data, order.by = dates)
  
  object <- list(
    tmed_est = tmed,
    twater_est = w_tmed,
    prec = meteo$prec,
    rhum = meteo$urel,
    tot_evap = tot_evap,
    tresh_rain_dry = tresh_rain,
    daylength = day_len,
    ndays = len_ts,
    d_emergency = as.numeric(df_data$emergency),
    d_induction = diapause_induct,
    tresh_rain_dry = tresh_rain,
    prevdrydays = prevdrydays,
    weight_k = weigth_k,
    weight_dry = weight_dry,
    timeformat = meteo$timeformat,
    dates = dates,
    timezone = timezone,
    timeseries = ts_zoo,
    datemax = datemax,
    bounddays = bounddays,
    varjd = varjd
  )
  
  attr(object, "tmed_est") <-
    "Daily mean air temperature estimated"
  attr(object, "twater_est") <-
    "Daily mean water temperature estimated"
  attr(object, "prec") <- "Daily precipitation"
  attr(object, "rhum") <- "Daily evaporation"
  attr(object, "tot_evap") <-
    "Mass of evaporative losses for day in mg"
  attr(object, "tresh_rain_dry") <-
    "Rain threshold to consider effective precipitation"
  attr(object, "prevdrydays") <-
    "Consecutive dry days previously accerted"
  attr(object, "weight_k") <-
    "Weighting Paraemter to take into account dryness impact"
  attr(object, "daylength") <- "Daylength"
  attr(object, "ndays") <- "Number of days"
  attr(object, "d_emergency") <-
    "Percentage of exclosion of daily diapausant eggs in springtime"
  attr(object, "d_induction") <-
    "Percentage of daily autumnal diapausant eggs"
  attr(object, "dates") <- "Dates"
  attr(object, "timezone") <- "Timezone"
  attr(object, "timeformat") <- "Period of data aggregation"
  attr(object, "timeseries") <- "Data as xts object"
  attr(object, "datemax") <-
    "Julian day of maximal date of spring diapause awekaning"
  attr(object, "bounddays") <-
    "Days of incertainty for diapause emergency period bounds"
  attr(object, "varjd") <-
    "Days of variance for diapause emergency courbes"
  class(object) <- "biometeo"
  return(object)
}
