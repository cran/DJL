\name{dm.sf}
\alias{dm.sf}
\title{Distance measure using SF}
\description{Implements \emph{Luenberger}'s shortage (benefit) function (radial & non-oriented measure).}
\usage{dm.sf(xdata, ydata, rts="crs", g=NULL,
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
  \item{g}{Directional vector indicating a measurement direction (\emph{n} by \emph{(m+s)})\cr
           By default (\code{NULL}), \code{xdata} & \code{ydata} will be used}
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
  \item{\code{$w}}{Input (dual) weight}
  \item{\code{$p}}{Output (dual) weight}
  \item{\code{$u}}{Free (dual) variable}
}
\references{
Luenberger, David G. "Benefit functions and duality." \emph{Journal of mathematical economics} 21.5 (1992): 461~481.

Chambers, Robert G., Yangho Chung, and Rolf Fare. "Profit, directional distance functions, and Nerlovian efficiency." \emph{Journal of optimization theory and applications} 98.2 (1998): 351~364.
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
# Additive form shortage function
  # ready
  x <- matrix(c(5, 1, 4), ncol = 1)
  y <- matrix(c(8, 3, 5, 6, 4, 1), ncol = 2)
  g <- matrix(c(1), nrow = 3, ncol = 3) 
  w <- matrix(c(1, 0), ncol = 2)
  # go
  dm.sf(x, y, "crs", g, w)

# Multiplicative form shortage function
  # ready
  g <- cbind(x, y)
  # go
  dm.sf(x, y, "crs", g, w)
}