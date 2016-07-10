\name{map.soa.ddf}
\alias{map.soa.ddf}
\title{SOA mapping using DDF}
\description{Employs \code{dm.ddf} over time to generate a state-of-the-art map.}
\usage{map.soa.ddf(xdata, ydata, date, rts="crs", g, 
            wd=NULL, sg="ssm", cv="convex", mk="dmu")}
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
  \item{g}{Directional vector indicating a measurement direction (\emph{n} by \emph{(m+s)})}
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
\author{Dong-Joon Lim, PhD}
\seealso{
  \code{\link{map.soa.ddf}} SOA mapping using DDF \cr
  \code{\link{map.soa.dea}} SOA mapping using DEA \cr
  \code{\link{map.soa.sbm}} SOA mapping using SBM \cr
  \code{\link{map.soa.sf}} SOA mapping using SF \cr
}
\examples{
# Load engine dataset
  df<-dataset.engine.2015

# Subset for forced induction systems
  fis<-subset(df,grepl("^.C..",df[,8]))

# Parameters
  x<-subset(fis,select=4)
  y<-subset(fis,select=6:7)
  d<-subset(fis,select=2)
  g<-matrix(c(1),nrow=nrow(x),ncol=3)

# Generate an SOA map
  map.soa.ddf(x,y,d,"crs",g)
}