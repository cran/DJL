\name{map.soa.dea}
\alias{map.soa.dea}
\title{Generate SOA map using DEA}
\description{Employs \code{dm.dea} over time to generate a state-of-the-art map.}
\usage{map.soa.dea(xdata, ydata, date, rts, orientation,
            sg="ssm", ncv=NULL, env=NULL, mk="dmu")}
\arguments{
  \item{xdata}{Input(s) vector (\emph{n} by \emph{m} matrix)}
  \item{ydata}{Output(s) vector (\emph{n} by \emph{s} matrix)}
  \item{date}{Production date (\emph{n} by \emph{1} matrix)}
  \item{rts}{Returns to scale assumption \cr
    \code{"crs"  }{Constant RTS} \cr  
    \code{"vrs"  }{Variable RTS} \cr
    \code{"irs"  }{Increasing RTS} \cr
    \code{"drs"  }{Decreasing RTS}
  }
  \item{orientation}{Orientation of the measurement \cr
    \code{"i"  }{Input-orientation} \cr
    \code{"o"  }{Output-orientation}
  }
  \item{sg}{Employs second-stage optimization \cr
    \code{"ssm"  }{Slack-sum maximization (default)} \cr
    \code{"max"  }{Date-sum maximization (only if \code{date} is defined)} \cr
    \code{"min"  }{Date-sum minimization (only if \code{date} is defined)}
  }
  \item{ncv}{Non-controllable variable index for internal NDF (\emph{1} by \emph{(m+s)} binary matrix)}
  \item{env}{Environment index for external NDF (\emph{n} by \emph{1} matrix)}
  \item{mk}{Marker on the map \cr
    \code{"dmu"  }{DMU index (default)} \cr
    \code{"eff"  }{Efficiency score } 
  }
}
\author{Dong-Joon Lim, PhD}
\seealso{
  \code{\link{map.soa.ddf}} SOA map using DDF \cr
  \code{\link{map.soa.dea}} SOA map using DEA \cr
  \code{\link{map.soa.sbm}} SOA map using SBM \cr
  \code{\link{map.soa.sf}} SOA map using SF \cr
}
\examples{
# Load engine dataset
  data(dataset.engine.2015)

# Subset for forced induction systems
  fis<-subset(dataset.engine.2015,grepl("^.C..",dataset.engine.2015[,8]))

# Parameters
  d<-as.matrix(fis[,2])
  x<-as.matrix(fis[,4])
  y<-as.matrix(fis[,6:7])

# Generate an SOA map
  map.soa.dea(x,y,d,"crs","o")
}