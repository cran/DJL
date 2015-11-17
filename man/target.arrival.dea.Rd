\name{target.arrival.dea}
\alias{target.arrival.dea}
\title{Arrival target setting using DEA}
\description{Employs \code{dm.dea} over time to estimate the arrival of known specifications.}
\usage{target.arrival.dea(xdata, ydata, date, t, rts, orientation, sg="ssm", ftype="d")}
\arguments{
  \item{xdata}{Input(s) vector (\emph{n} by \emph{m} matrix)}
  \item{ydata}{Output(s) vector (\emph{n} by \emph{s} matrix)}
  \item{date}{Production date (\emph{n} by \emph{1} matrix)}
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
}
\value{
  \item{\code{$eff_t}}{Efficiency at \code{t}}
  \item{\code{$lambda_t}}{Intensity vector at \code{t}}
  \item{\code{$eft_date}}{Effective date}
  \item{\code{$roc_avg}}{Average RoC}
  \item{\code{$roc_local}}{Local RoC}
  \item{\code{$roc_ind}}{Individualized RoC}
  \item{\code{$arrival_avg}}{Estimated arrival using \code{roc_avg}}
  \item{\code{$arrival_seg}}{Estimated arrival using \code{roc_ind}}
}
\references{
Lim, Dong-Joon, Timothy R. Anderson, and Oliver Lane Inman. "Choosing effective dates from multiple optima in Technology Forecasting using Data Envelopment Analysis (TFDEA)." \emph{Technological Forecasting and Social Change} 88 (2014): 91~97.

Lim, Dong-Joon, and Timothy R. Anderson. Time series benchmarking analysis for new product scheduling: who are the competitors and how fast are they moving forward?. Advances in DEA Theory and Applications: with Examples in Forecasting Models. \emph{Wiley (forthcoming)}, 2016.
}
\author{Dong-Joon Lim, PhD}
\seealso{
  \code{\link{dm.dea}} Distance measure using DEA \cr
  \code{\link{roc.dea}} RoC calculation using DEA \cr
  \code{\link{target.arrival.dea}} Arrival target setting using DEA \cr
  \code{\link{target.spec.dea}} Spec target setting using DEA
}
\examples{
# Reproduce Table 4 in Lim, D-J., and Timothy R. Anderson.(2016)
  # Load airplane dataset
  data(dataset.airplane.2017)
  
  # ready
  d<-as.matrix(dataset.airplane.2017[,2])
  x<-matrix(c(1),nrow=28,ncol=1)
  y<-as.matrix(dataset.airplane.2017[,3:7])
  
  # go
  target.arrival.dea(x,y,d,2007,"vrs","o","min","d")$arrival_seg
}