\name{dataset.engine.2015}
\alias{dataset.engine.2015}
\docType{data}
\title{Dataset of auto engines from MY2005 to MY2015.}
\description{Auto engines from MY2005 to MY2015.}
\section{Columns}{
  \code{[,1] Name          }{Vehicle name} \cr
  \code{[,2] MY            }{Model year} \cr
  \code{[,3] Cylinder      }{The number of cylinder} \cr
  \code{[,4] Displacement  }{Displacement in liter} \cr
  \code{[,5] CO2           }{CO2 emission in gram/mile} \cr
  \code{[,6] Power         }{Engine power in HP} \cr
  \code{[,7] Torque        }{Engine torque in lb.ft} \cr
  \code{[,8] Type          }{Engine system and fuel types}
}
\usage{data(dataset.engine.2015)}
\source{http://www.fueleconomy.gov \cr http://www.autoevolution.com}
\author{Dong-Joon Lim, PhD}
\references{
D.-J. Lim,  Internal combustion engine race: naturally aspirated vs turbo/super-charged, \emph{working paper} (2015).
}
\examples{
# Load dataset
  data(dataset.engine.2015)
}