\name{map.soa.sf}
\alias{map.soa.sf}
\title{SOA mapping using SF}
\description{Employs \code{dm.sf} over time to generate a state-of-the-art map.}
\usage{map.soa.sf(xdata, ydata, date, rts, g, w=NULL, sg="ssm", mk="dmu")}
\arguments{
  \item{xdata}{Input(s) vector (\emph{n} by \emph{m})}
  \item{ydata}{Output(s) vector (\emph{n} by \emph{s})}
  \item{date}{Production date (\emph{n} by \emph{1})}
  \item{rts}{Returns to scale assumption \cr
    \code{"crs"  }{Constant RTS} \cr  
    \code{"vrs"  }{Variable RTS} \cr
    \code{"irs"  }{Increasing RTS} \cr
    \code{"drs"  }{Decreasing RTS}
  }
  \item{g}{Directional vector indicating a measurement direction (\emph{n} by \emph{(m+s)})}
  \item{w}{Weak disposability vector indicating (an) undesirable output(s) (\emph{1} by \emph{s})}
  \item{sg}{Employs second-stage optimization \cr
    \code{"ssm"  }{Slack-sum maximization (default)} \cr
    \code{"max"  }{Date-sum maximization (only if \code{date} is defined)} \cr
    \code{"min"  }{Date-sum minimization (only if \code{date} is defined)}
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
  \code{\link{map.soa.sbm}} SOA mapping using SBM \cr
  \code{\link{map.soa.sf}} SOA mapping using SF \cr
}
\examples{
# Reproduce Table 2 in Lim, D-J. (2015)
  # Load engine dataset
    data(dataset.engine.2015)

  # Subset for 4 cylinder engines
    fce<-subset(dataset.engine.2015,dataset.engine.2015[,3]==4)
    
  # Parameters
    x<-subset(fce,select=4)
    y<-subset(fce,select=5:7)
    d<-subset(fce,select=2)
    g<-data.frame(0,y)
    w<-matrix(c(1,0,0),ncol=3)

  # Generate an SOA map
    map.soa.sf(x,y,d,"crs",g,w,mk="eff")
}