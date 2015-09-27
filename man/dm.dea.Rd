\name{dm.dea}
\alias{dm.dea}
\title{Distance measure using DEA}
\description{Implements \emph{data envelopment analysis} (radial & oriented measure).}
\usage{dm.dea(xdata, ydata, rts, orientation,
       se=0, sg="ssm", date=NULL, ncv=NULL, env=NULL)}
\arguments{
  \item{xdata}{Input(s) vector (\emph{n} by \emph{m} matrix)}
  \item{ydata}{Output(s) vector (\emph{n} by \emph{s} matrix)}
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
  \item{se}{Implements \emph{Anderson & Peterson}'s super-efficiency model if \code{1}}
  \item{sg}{Employs second-stage optimization \cr
    \code{"ssm"  }{Slack-sum maximization (default)} \cr
    \code{"max"  }{Date-sum maximization (only if \code{date} is defined)} \cr
    \code{"min"  }{Date-sum maximization (only if \code{date} is defined)}
  }
  \item{date}{Production date (\emph{n} by \emph{1} matrix)}
  \item{ncv}{Non-controllable variable index for internal NDF (\emph{1} by \emph{(m+s)} binary matrix)}
  \item{env}{Environment index for external NDF (\emph{n} by \emph{1} matrix)}
}
\value{
  \item{\code{$eff}}{Efficiency score}
  \item{\code{$lambda}}{Intensity vector}
  \item{\code{$xslack}}{Input slack(s)}
  \item{\code{$yslack}}{Output slack(s)}
  \item{\code{$vx}}{Input (dual) weight}
  \item{\code{$uy}}{Output (dual) weight}
  \item{\code{$w}}{Free (dual) variable}
}
\references{
Charnes, Abraham, William W. Cooper, and Edwardo Rhodes. "Measuring the efficiency of decision making units." \emph{European journal of operational research} 2.6 (1978): 429~444.

Banker, Rajiv D., Abraham Charnes, and William Wager Cooper. "Some models for estimating technical and scale inefficiencies in data envelopment analysis." \emph{Management science} 30.9 (1984): 1078~1092.

Banker, Rajiv D., and Richard C. Morey. "Efficiency analysis for exogenously fixed inputs and outputs." \emph{Operations Research} 34.4 (1986): 513~521.

Ruggiero, John. "On the measurement of technical efficiency in the public sector." \emph{European Journal of Operational Research} 90.3 (1996): 553~565.

Fried, Harold O., CA Knox Lovell, and Shelton S. Schmidt, eds. The measurement of productive efficiency and productivity growth. \emph{Oxford University Press}, 2008.
}
\author{Dong-Joon Lim, PhD}
\seealso{
  \code{\link{dm.ddf}} Distance measure using DDF \cr
  \code{\link{dm.dea}} Distance measure using DEA \cr
  \code{\link{dm.sbm}} Distance measure using SBM \cr
  \code{\link{dm.sf}} Distance measure using SF
}
\examples{
# Reproduce Table 3.9 (p.348) in Fried, H.O. et al.(2008)
  # ready
  x<-matrix(c(8,6,3,10,6,8,8,4,8,4.6,1.9,9,3.6,3.6,9,1.9),ncol=2)
  y<-matrix(c(8,5,2,9,4.5,4.5,7,2),ncol=1)
  c<-matrix(c(0,1,0),nrow=1)
  # go
  dm.dea(x,y,"crs","i")$eff
  dm.dea(x,y,"vrs","i")$eff
  dm.dea(x,y,"crs","i",ncv=c)$eff
  dm.dea(x,y,"vrs","i",ncv=c)$eff
}