\name{dm.hdf}
\alias{dm.hdf}
\title{Distance measure using HDF}
\description{Implements \emph{Fare}'s hyperbolic distance function (semi-radial & non-oriented measure).}
\usage{dm.hdf(xdata, ydata, rts="crs", 
       wd=NULL, se=FALSE, sg="ssm", date=NULL, cv="convex", o=NULL)}
\arguments{
  \item{xdata}{Input(s) vector (\emph{n} by \emph{m})}
  \item{ydata}{Output(s) vector (\emph{n} by \emph{s})}
  \item{rts}{Returns to scale assumption \cr
    \code{"crs"  }{Constant RTS (default)} \cr  
    \code{"vrs"  }{Variable RTS} \cr
    \code{"irs"  }{Increasing RTS} \cr
    \code{"drs"  }{Decreasing RTS}
  }
  \item{wd}{Weak disposability vector indicating (an) undesirable output(s) (\emph{1} by \emph{s})}
  \item{se}{Implements super-efficiency model alike \emph{Anderson & Peterson}'s model if \code{TRUE}}
  \item{sg}{Employs second-stage optimization \cr
    \code{"ssm"  }{Slack-sum maximization (default)} \cr
    \code{"max"  }{Date-sum maximization (only if \code{date} is defined)} \cr
    \code{"min"  }{Date-sum minimization (only if \code{date} is defined)}
  }
  \item{date}{Production date (\emph{n} by \emph{1})}
  \item{cv}{Convexity assumption \cr
    \code{"convex" }{Convexity holds (default)} \cr
    \code{"fdh"    }{Free disposal hull (this will override \code{rts})}
  }
  \item{o}{DMU index to calc. \code{NULL}(default) will calc for all}
}
\value{
  \item{\code{$eff}}{Efficiency score}
  \item{\code{$lambda}}{Intensity vector}
  \item{\code{$mu}}{Secondary intensity vector for weak disposability under VRS}
  \item{\code{$xslack}}{Input slack}
  \item{\code{$yslack}}{Output slack}
  \item{\code{$iteration}}{The number of iteration to obtain the hyperbolic efficiency score}
}
\references{
Fare, R., Shawna Grosskopf, and CA Knox Lovell. The Measurement of Efficiency of Production. \emph{Boston: Kulwer-Nijhoff} (1985).

Fare, Rolf, et al. "Estimating the hyperbolic distance function: A directional distance function approach." \emph{European Journal of Operational Research} 254.1 (2016): 312~319.
}
\author{Dong-Joon Lim, PhD}
\seealso{
  \code{\link{dm.ddf}} Distance measure using DDF \cr
  \code{\link{dm.dea}} Distance measure using DEA \cr
  \code{\link{dm.hdf}} Distance measure using HDF \cr
  \code{\link{dm.sbm}} Distance measure using SBM \cr
  \code{\link{dm.sf}} Distance measure using SF
}
\examples{
# Reproduce Table 2 in Fare et al.(2016)
  # ready
  x <- data.frame(x1 = c(2, 4, 9, 6.5, 10, 6, 9))
  y <- data.frame(y1 = c(3, 7, 10, 8.5, 4, 2, 8))

  # go
  sf  <- dm.sf (x, y, "vrs")$eff
  hdf <- dm.hdf(x, y, "vrs")$eff
  matrix(t(cbind(sf, hdf)), 2, 7, 
         dimnames = list(c("SF", "HDF"), 
                         paste0("DMU_", c(letters[1:4], "o", "p", "q")))) 

}