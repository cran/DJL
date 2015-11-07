\name{map.soa.sf}
\alias{map.soa.sf}
\title{Generate SOA map using SF}
\description{Employs \code{dm.sf} over time to generate a state-of-the-art map.}
\usage{map.soa.sf(xdata, ydata, date, rts, g, w=NULL, sg="ssm", mk="dmu")}
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
  \code{\link{map.soa.ddf}} SOA map using DDF \cr
  \code{\link{map.soa.dea}} SOA map using DEA \cr
  \code{\link{map.soa.sbm}} SOA map using SBM \cr
  \code{\link{map.soa.sf}} SOA map using SF \cr
}
\examples{
# Reproduce Table 2 in Lim, D-J. (2015)
  # Load engine dataset
    data(dataset.engine.2015)
  
  # Parameters
    c<-matrix(as.numeric(dataset.engine.2015[,3]),ncol=1)
    x<-matrix(as.numeric(dataset.engine.2015[,4]),ncol=1)
    y<-matrix(as.numeric(dataset.engine.2015[,5:7]),ncol=3)
    d<-matrix(as.numeric(dataset.engine.2015[,2]),ncol=1)
    g<-matrix(c(rep(0,nrow(x)),y),nrow=nrow(x),ncol=4)
    w<-matrix(c(1,0,0),ncol=3)

  # Subset for 4 cylinder engines
    x_s<-subset(x,c==4)
    y_s<-subset(y,c==4)
    d_s<-subset(d,c==4)
    g_s<-subset(g,c==4)

  # Generate an SOA map
    map.soa.sf(x_s,y_s,d_s,"crs",g_s,w,mk="eff")
}