\name{dm.sbm}
\alias{dm.sbm}
\title{Distance measure using SBM}
\description{Implements \emph{Tone}'s slack-based model (non-radial & (non-)oriented measure).}
\usage{dm.sbm(xdata, ydata, rts, orientation="n", se=0, sg="ssm", date=NULL)}
\arguments{
  \item{xdata}{Input(s) vector (\emph{n} by \emph{m})}
  \item{ydata}{Output(s) vector (\emph{n} by \emph{s})}
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
  \item{se}{Implements super-efficiency model alike \emph{Anderson & Peterson}'s model if \code{1}}
  \item{sg}{Employs second-stage optimization \cr
    \code{"ssm"  }{Slack-sum maximization (default)} \cr
    \code{"max"  }{Date-sum maximization (only if \code{date} is defined)} \cr
    \code{"min"  }{Date-sum minimization (only if \code{date} is defined)}
  }
  \item{date}{Production date (\emph{n} by \emph{1})}
}
\value{
  \item{\code{$eff}}{Efficiency score}
  \item{\code{$lambda}}{Intensity vector}
  \item{\code{$xslack}}{Input slack}
  \item{\code{$yslack}}{Output slack}
  \item{\code{$xtarget}}{Input target}
  \item{\code{$ytarget}}{Output target}
}
\references{
Tone, Kaoru. "A slacks-based measure of efficiency in data envelopment analysis." \emph{European journal of operational research} 130.3 (2001): 498~509.

Tone, Kaoru. "A slacks-based measure of super-efficiency in data envelopment analysis." \emph{European journal of operational research} 143 (2002): 32~41.
}
\author{Dong-Joon Lim, PhD}
\seealso{
  \code{\link{dm.ddf}} Distance measure using DDF \cr
  \code{\link{dm.dea}} Distance measure using DEA \cr
  \code{\link{dm.sbm}} Distance measure using SBM \cr
  \code{\link{dm.sf}} Distance measure using SF
}
\examples{
# Reproduce Table 2 in Tone.(2001)
  # ready
  x<-matrix(c(4,6,8,8,2,3,3,1,1,4),ncol=2)
  y<-matrix(c(2,2,6,6,1,3,3,2,1,4),ncol=2)
  # go
  dm.sbm(x,y,"crs")

# Reproduce Table 1 in Tone.(2002)
  # Published input slacks are alternate optima (confirmed by Tone)
  # ready
  x<-matrix(c(4,7,8,4,2,10,12,3,3,1,2,4,1,1),ncol=2)
  y<-matrix(c(1,1,1,1,1,1,1),ncol=1)
  # go
  dm.sbm(x,y,"crs",se=1)

# Reproduce Table 4 in Tone.(2002)
  # ready
  x<-matrix(c(80,65,83,40,52,94,600,200,400,1000,600,700,54,97,72,75,20,36,8,1,4,7,3,5),
            ncol=4)
  y<-matrix(c(90,58,60,80,72,96,5,1,7,10,8,6),ncol=2)
  # go
  dm.sbm(x,y,"crs","i",se=1)
}