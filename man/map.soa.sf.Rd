\name{map.soa.sf}
\alias{map.soa.sf}
\title{Generate SOA map using SF}
\description{Employs \code{dm.sf} over time to generate a state-of-the-art map.}
\usage{map.soa.sf(xdata, ydata, date, rts, g, w=0, sg="ssm")}
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
  \item{g}{Directional vector indicating a measurement direction (\emph{n} by \emph{(m+s)} matrix)}
  \item{w}{Weak disposability vector indicating (an) undesirable output(s) (\emph{1} by \emph{s} matrix)}
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
  g<-matrix(c(x,y),nrow=nrow(x),ncol=3)

# generate an SOA map
  map.soa.sf(x,y,d,"crs",g)
}