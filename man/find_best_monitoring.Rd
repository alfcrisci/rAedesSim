% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/find_best_monitoring.r
\name{find_best_monitoring}
\alias{find_best_monitoring}
\title{find_best_monitoring}
\usage{
find_best_monitoring(i_meteo, i_biocontainer, i_monitoring,
  initial_eggs = 100, range_alpha_a = c(0, seq(0, 0.002, 0.001)),
  range_alpha_l = seq(0.6, 1.6, 0.2), range_density_l = 70,
  stocastic = TRUE, n_sampling = 10, inibition = FALSE)
}
\arguments{
\item{i_meteo}{object: rAedesSim meteo  concerning weather data.}

\item{i_monitoring}{object: rAedesSim \code{biodata} object concerning mosquito eggs field observations.}

\item{initial_eggs}{numeric: initial values for \code{biodata} object concerning mosquito eggs field observations.}

\item{range_alpha_a}{numeric: rAedesSim vector of sorted  guess' values of female adult competition. Default is obtained by c(0,seq(0,0.002,0.001)).}

\item{range_alpha_l}{numeric: rAedesSim vector of sorted guess' values of larval competition.  Default is obtained by seq(0.6,1.6,0.2).}

\item{range_density_l}{numeric: rAedesSim object guess' values of maximaum larval competion in function of waer volume in cointaner. Default is obtained by 70.}

\item{stocastic}{logical: if stocastic simulation by \code{biomodel} are carried out in simulation.Defalut is True.}

\item{n_sampling}{numeric: number  of resampling if stochastic is implemented by \code{biomodel}. Default is 10.}

\item{inibition}{logical: if larval density is considered. Defalut is FALSE.}

\item{i_biocointaner}{object: rAedesSim object concerning \code{biocontainer} object used in simulation.}
}
\value{
rAedesSim object biofitmodel object.
}
\description{
\code{find_best_monitoring} Function to semplify  model tuning by using \code{biofitmodel} rAedesSim Class.
}
\author{
Istituto di Biometeorologia Firenze Italy  Alfonso crisci \email{a.crisci@ibimet.cnr.it} ASL 2 LUCCA Marco Selmi \email{marco.selmi@uslnordovest.toscana.it}
}
\keyword{modeling}
\keyword{rAedesSim}
