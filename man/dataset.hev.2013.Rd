\name{dataset.hev.2013}
\alias{dataset.hev.2013}
\docType{data}
\title{Dataset of hybrid electric vehicles from MY1997 to MY2013.}
\description{Hybrid electric vehicles from MY1997 to MY2015.}
\section{Columns}{
  \code{[,1] Name          }{Vehicle name} \cr
  \code{[,2] MY            }{Model year} \cr
  \code{[,3] MSRP_2013     }{MSRP converted to 2013 value} \cr
  \code{[,4] Acc           }{Acceleration (0-100km) in km/h/s} \cr
  \code{[,5] MPG           }{MPG in mile/gallon} \cr
  \code{[,6] MPGe          }{MPG equivalence for PHEV in mile/gallon}
}
\usage{data(dataset.hev.2013)}
\source{http://www.fueleconomy.gov}
\author{Dong-Joon Lim, PhD}
\references{
Lim, Dong-Joon, et al. "Comparing technological advancement of hybrid electric vehicles (HEV) in different market segments." \emph{Technological Forecasting and Social Change} 97 (2015): 140~153.
}
\examples{
# Load dataset
  data(dataset.hev.2013)
}