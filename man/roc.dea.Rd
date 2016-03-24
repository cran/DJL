\name{roc.dea}
\alias{roc.dea}
\title{Rate of change (RoC) calculation using DEA}
\description{Employs \code{dm.dea} over time to calculate RoCs.}
\usage{roc.dea(xdata, ydata, date, t, rts, orientation,
        sg="ssm", ftype="d", ncv=NULL, env=NULL)}
\arguments{
  \item{xdata}{Input(s) vector (\emph{n} by \emph{m})}
  \item{ydata}{Output(s) vector (\emph{n} by \emph{s})}
  \item{date}{Production date (\emph{n} by \emph{1})}
  \item{t}{A vantage point from which the RoC is captured}
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
    \code{"max"  }{Date-sum maximization} \cr
    \code{"min"  }{Date-sum minimization}
  }
  \item{ftype}{Frontier type \cr
    \code{"d"  }{Dynamic frontier (default)} \cr  
    \code{"s"  }{Static frontier}
  }
  \item{ncv}{Non-controllable variable index(binary) for internal NDF (\emph{1} by \emph{(m+s)})}
  \item{env}{Environment index for external NDF (\emph{n} by \emph{1})}
}
\value{
  \item{\code{$eff_r}}{Efficiency at release (i.e., at each production date)}
  \item{\code{$eff_t}}{Efficiency at \code{t}}
  \item{\code{$lambda_t}}{Intensity vector at \code{t}}
  \item{\code{$eft_date}}{Effective date}
  \item{\code{$roc_past}}{RoC observed from the obsolete DMUs in the past}
  \item{\code{$roc_avg}}{Average RoC}
  \item{\code{$roc_local}}{Local RoC}
}
\references{
Lim, Dong-Joon, Timothy R. Anderson, and Oliver Lane Inman. "Choosing effective dates from multiple optima in Technology Forecasting using Data Envelopment Analysis (TFDEA)." \emph{Technological Forecasting and Social Change} 88 (2014): 91~97.

Lim, Dong-Joon, et al. "Comparing technological advancement of hybrid electric vehicles (HEV) in different market segments." \emph{Technological Forecasting and Social Change} 97 (2015): 140~153.

Lim, Dong-Joon, et al. Technometrics Study Using DEA on Hybrid Electric Vehicles (HEVs). Handbook of Operations Analytics Using Data Envelopment Analysis. \emph{Springer (forthcoming)}, 2016.
}
\author{Dong-Joon Lim, PhD}
\seealso{
  \code{\link{dm.dea}} Distance measure using DEA \cr
  \code{\link{roc.dea}} RoC calculation using DEA \cr
  \code{\link{target.arrival.dea}} Arrival target setting using DEA \cr
  \code{\link{target.spec.dea}} Spec target setting using DEA
}
\examples{
# Reproduce Table 3 in Lim, D-J. et al.(2014)
  # Load airplane dataset
  data(dataset.airplane.2017)
  
  # ready
  x<-data.frame(Flew=rep(1,28))
  y<-subset(dataset.airplane.2017,select=3:7)
  d<-subset(dataset.airplane.2017,select=2)

  # go
  roc.dea(x,y,d,2007,"vrs","o","min","d")$roc_past

# Reproduce Table 3 in Lim, D-J. et al.(2015)
  # Load hev dataset
  data(dataset.hev.2013)
  
  # ready
  x<-subset(dataset.hev.2013,select=3)
  y<-subset(dataset.hev.2013,select=4:6)
  d<-subset(dataset.hev.2013,select=2)
  c<-subset(dataset.hev.2013,select=7)

  # go
  results<-roc.dea(x,y,d,2013,"vrs","o","min","d",env=c)
  hev<-which(results$roc_local>0)
  data.frame(Class=c[hev,],SOA=hev,LocalRoC=results$roc_local[hev,])[order(c[hev,]),]
  # NOTE: the published results include a typo on roc_local[82,]
  #       this will be corrected in forthcoming book chapter(Lim, D-J. et al., 2016).
}