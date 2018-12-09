\name{roc.malmquist}
\alias{roc.malmquist}
\title{Malmquist Index: time-series productivity analysis}
\description{Employs distance measure over time to calculate the productivity changes.}
\usage{roc.malmquist(xdata, ydata, tm=NULL, dm="dea", rts="crs", orientation,
              g=NULL, wd=NULL, ncv=NULL, env=NULL, cv="convex")}
\arguments{
  \item{xdata}{Input(s) array (\emph{n} by \emph{m} by \emph{t})}
  \item{ydata}{Output(s) array (\emph{n} by \emph{s} by \emph{t})}
  \item{tm}{Tick mark of production dates (a vector length of \emph{t})}
  \item{dm}{Distance measure to calculate the productivity \cr
    \code{"dea"  }{Data Envelopment Analysis (default)} \cr  
    \code{"sbm"  }{Slack Based Model} \cr
    \code{"ddf"  }{Directional Distance Function} \cr
    \code{"hdf"  }{Hyperbolic Distance Function} \cr
    \code{"sf"   }{Shortage Function}
  }
  \item{rts}{Returns to scale assumption \cr
    \code{"crs"  }{Constant RTS (default)} \cr  
    \code{"vrs"  }{Variable RTS} \cr
    \code{"irs"  }{Increasing RTS} \cr
    \code{"drs"  }{Decreasing RTS}
  }
  \item{orientation}{Orientation of the measurement \cr
    \code{"n"  }{Non-orientation (default)} \cr
    \code{"i"  }{Input-orientation} \cr
    \code{"o"  }{Output-orientation}
  }
  \item{g}{Directional vector indicating a measurement direction (\emph{n} by \emph{(m+s)})\cr
           By default (\code{NULL}), \code{xdata} & \code{ydata} will be used}
  \item{wd}{Weak disposability vector indicating (an) undesirable output(s) (\emph{1} by \emph{s})}
  \item{ncv}{Non-controllable variable index(binary) for internal NDF (\emph{1} by \emph{(m+s)})}
  \item{env}{Environment index for external NDF (\emph{n} by \emph{1})}
  \item{cv}{Convexity assumption \cr
    \code{"convex" }{Convexity holds (default)} \cr
    \code{"fdh"    }{Free disposal hull (this will override \code{rts})}
  }
}
\value{
  \item{\code{$cu}}{Catching Up (aka technical efficiency change: TEC) index}
  \item{\code{$fs}}{Frontier Shift (FS) Index}
  \item{\code{$mi}}{Malmquist Index}
}
\references{
R. Fare, S. Grosskopf, and C. A. K. Lovell, Production Frontiers. \emph{Cambridge University Press}, 1994.
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
# Load data
  df <- array(c(4,     3,   9,  10,  7, 4,  3, 5,
                5,    12,   3,   8,  1, 4, 14, 3,
                1,     1,   1,   1,  1, 1,  1, 1,
                3.4,   2,  10,   8, 10, 4,  1, 5,
                6,    10, 3.5,   7,  2, 4, 12, 3,
                1,     1,   1,   1,  1, 1,  1, 1,
                2.8, 1.8,   8,   7, 10, 3,  1, 5,
                5.7, 8.8, 2.8,   5,  2, 5,  9, 3,
                1,     1,   1,   1,  1, 1,  1, 1,
                2.2, 1.5,   8,   5,  8, 3,  1, 5,
                6,     8, 2.3, 3.5,  2, 5,  7, 3,
                1,     1,   1,   1,  1, 1,  1, 1),
              c(8, 3, 4))

# Run
  roc.malmquist(df[,1:2,], df[,3,], dm = "sbm", orientation = "n")
}