\name{dataset.airplane.2017}
\alias{dataset.airplane.2017}
\docType{data}
\title{Dataset of commercial airplanes from 1965 to 2017.}
\description{Commercial airplanes from 1965 to 2017.}
\section{Columns}{
  \code{[,1] Name       }{Airplane name} \cr
  \code{[,2] EIS        }{Entry into service} \cr
  \code{[,3] Range      }{Maximum range at full payload in 1,000km} \cr
  \code{[,4] P.cap      }{Passenger capacity} \cr
  \code{[,5] PFE        }{Passenger fuel efficiency in passengers*km/L (log scale)} \cr
  \code{[,6] C.speed    }{Cruising speed in km/hr} \cr
  \code{[,7] M.speed    }{Maximum speed in km/hr}
}
\usage{data(dataset.airplane.2017)}
\source{http://www.airbus.com/aircraftfamilies \cr http://www.boeing.com/commercial}
\author{Dong-Joon Lim, PhD}
\references{
Lim, Dong-Joon, and Timothy R. Anderson. Time series benchmarking analysis for new product scheduling: who are the competitors and how fast are they moving forward?. Advances in DEA Theory and Applications: with Examples in Forecasting Models. \emph{Wiley (forthcoming)}, 2016.
}
\examples{
# Load dataset
  data(dataset.airplane.2017)
}