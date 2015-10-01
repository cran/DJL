\name{dm.sf}
\alias{dm.sf}
\title{Distance measure using SF}
\description{Implements \emph{Luenberger}'s shortage function (radial & non-oriented measure).}
\usage{dm.sf(xdata, ydata, rts, g, w=NULL, se=0, sg="ssm", date=NULL)}
\arguments{
  \item{xdata}{Input(s) vector (\emph{n} by \emph{m} matrix)}
  \item{ydata}{Output(s) vector (\emph{n} by \emph{s} matrix)}
  \item{rts}{Returns to scale assumption \cr
    \code{"crs"  }{Constant RTS} \cr  
    \code{"vrs"  }{Variable RTS} \cr
    \code{"irs"  }{Increasing RTS} \cr
    \code{"drs"  }{Decreasing RTS}
  }
  \item{g}{Directional vector indicating a measurement direction (\emph{n} by \emph{(m+s)} matrix)}
  \item{w}{Weak disposability vector indicating (an) undesirable output(s) (\emph{1} by \emph{s} matrix)}
  \item{se}{Implements super-efficiency model alike \emph{Anderson & Peterson}'s model if \code{1}}
  \item{sg}{Employs second-stage optimization \cr
    \code{"ssm"  }{Slack-sum maximization (default)} \cr
    \code{"max"  }{Date-sum maximization (only if \code{date} is defined)} \cr
    \code{"min"  }{Date-sum maximization (only if \code{date} is defined)}
  }
  \item{date}{Production date (\emph{n} by \emph{1} matrix)}
}
\value{
  \item{\code{$eff}}{Efficiency score}
  \item{\code{$lambda}}{Intensity vector}
  \item{\code{$mu}}{Secondary intensity vector for weak disposability under VRS}
  \item{\code{$xslack}}{Input slack}
  \item{\code{$yslack}}{Output slack}
}
\references{
Luenberger, David G. "Benefit functions and duality." \emph{Journal of mathematical economics} 21.5 (1992): 461~481.

Chambers, Robert G., Yangho Chung, and Rolf Fare. "Profit, directional distance functions, and Nerlovian efficiency." \emph{Journal of optimization theory and applications} 98.2 (1998): 351~364.
}
\author{Dong-Joon Lim, PhD}
\seealso{
  \code{\link{dm.ddf}} Distance measure using DDF \cr
  \code{\link{dm.dea}} Distance measure using DEA \cr
  \code{\link{dm.sbm}} Distance measure using SBM \cr
  \code{\link{dm.sf}} Distance measure using SF
}
\examples{
# Additive form shortage function
  # ready
  x<-matrix(c(5,1,4),ncol=1)
  y<-matrix(c(8,3,5,6,4,1),ncol=2)
  g<-matrix(c(1),nrow=3,ncol=3) 
  w<-matrix(c(1,0),ncol=2)
  # go
  dm.sf(x,y,"crs",g,w)

# Ratio form shortage function
  # ready
  g<-matrix(c(x,y),nrow=3,ncol=3)
  # go
  dm.sf(x,y,"crs",g,w)
}