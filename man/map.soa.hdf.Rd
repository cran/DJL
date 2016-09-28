\name{map.soa.hdf}
\alias{map.soa.hdf}
\title{SOA mapping using HDF}
\description{Employs \code{dm.hdf} over time to generate a state-of-the-art map.}
\usage{map.soa.hdf(xdata, ydata, date,
            rts="crs", wd=NULL, sg="ssm", cv="convex", mk="dmu")}
\arguments{
  \item{xdata}{Input(s) vector (\emph{n} by \emph{m})}
  \item{ydata}{Output(s) vector (\emph{n} by \emph{s})}
  \item{date}{Production date (\emph{n} by \emph{1})}
  \item{rts}{Returns to scale assumption \cr
    \code{"crs"  }{Constant RTS (default)} \cr  
    \code{"vrs"  }{Variable RTS} \cr
    \code{"irs"  }{Increasing RTS} \cr
    \code{"drs"  }{Decreasing RTS}
  }
  \item{wd}{Weak disposability vector indicating (an) undesirable output(s) (\emph{1} by \emph{s})}
  \item{sg}{Employs second-stage optimization \cr
    \code{"ssm"  }{Slack-sum maximization (default)} \cr
    \code{"max"  }{Date-sum maximization (only if \code{date} is defined)} \cr
    \code{"min"  }{Date-sum minimization (only if \code{date} is defined)}
  }
  \item{cv}{Convexity assumption \cr
    \code{"convex" }{Convexity holds (default)} \cr
    \code{"fdh"    }{Free disposal hull (this will override \code{rts})}
  }
  \item{mk}{Marker on the map \cr
    \code{"dmu"  }{DMU index (default)} \cr
    \code{"eff"  }{Efficiency score } 
  }
}
\references{
D.-J. Lim,  Internal combustion engine race: naturally aspirated vs turbo/super-charged, \emph{working paper} (2015).
}
\author{Dong-Joon Lim, PhD}
\seealso{
  \code{\link{map.soa.ddf}} SOA mapping using DDF \cr
  \code{\link{map.soa.dea}} SOA mapping using DEA \cr
  \code{\link{map.soa.hdf}} SOA mapping using HDF \cr
  \code{\link{map.soa.sbm}} SOA mapping using SBM \cr
  \code{\link{map.soa.sf}} SOA mapping using SF
}
\examples{
# Load engine dataset
  df <- dataset.engine.2015

# Subset for SC/TC 8 cylinder engines
  stc.8 <- subset(df, grepl("^.C..", df[, 8]) & df[, 3] == 8)

# Parameters
  x <- subset(stc.8, select = 4)
  y <- subset(stc.8, select = 5:7)
  d <- subset(stc.8, select = 2)

# Generate an SOA map
  map.soa.hdf(x, y, d, "vrs")
}