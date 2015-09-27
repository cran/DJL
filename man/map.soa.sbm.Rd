\name{map.soa.sbm}
\alias{map.soa.sbm}
\title{Generate SOA map using SBM}
\description{Employs \code{dm.sbm} over time to generate a state-of-the-art map.}
\usage{map.soa.sbm(xdata, ydata, date, rts, orientation, sg="ssm")}
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
    \code{"n"  }{Non-orientation (default)} \cr
    \code{"i"  }{Input-orientation} \cr
    \code{"o"  }{Output-orientation}
  }
  \item{sg}{Employs second-stage optimization \cr
    \code{"ssm"  }{Slack-sum maximization (default)} \cr
    \code{"max"  }{Date-sum maximization (only if \code{date} is defined)} \cr
    \code{"min"  }{Date-sum maximization (only if \code{date} is defined)}
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

# Say, efficiency is defined as (HP & Torque)/Displacement
  x<-matrix(as.numeric(fis[,4]),ncol=1)
  y<-matrix(as.numeric(fis[,6:7]),ncol=2)
  d<-matrix(as.numeric(fis[,2]),ncol=1)

# generate an SOA map
  map.soa.sbm(x,y,d,"crs")
}