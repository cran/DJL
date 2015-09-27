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
    \code{"min"  }{Date-sum maximization}
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
  \item{\code{$roc_ind}}{RoC expected for each target}
  \item{\code{$arrival_avg}}{Estimated arrival using \code{roc_avg}}
  \item{\code{$arrival_seg}}{Estimated arrival using \code{roc_ind}}
}
\references{
Lim, Dong-Joon, Timothy R. Anderson, and Oliver Lane Inman. "Choosing effective dates from multiple optima in Technology Forecasting using Data Envelopment Analysis (TFDEA)." \emph{Technological Forecasting and Social Change} 88 (2014): 91~97.
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
  # ready
  x<-matrix(c(1),nrow=28,ncol=1)
  y<-matrix(c(9.205, 9.62, 9.8, 12.7, 10.01, 9.265, 10.2, 12.4, 12.2, 11.065, 13.45,
              12.27, 10.5, 15, 13.7, 13.408, 14.305, 11.12, 12.5, 14.6, 16.7, 14.685,
              17.37, 15.2, 15.7, 14.816, 15.75, 15, 132, 159, 366, 366, 250, 250,
              234, 412, 181, 218, 416, 293, 295, 261, 295, 293, 301, 365, 253, 380,
              313, 365, 301, 525, 242, 467, 280, 315, 2.619, 2.812, 2.973, 3.15,
              2.901, 2.824, 2.987, 3.245, 3.192, 3.28, 3.251, 3.203, 3.462, 3.229,
              3.308, 3.216, 3.225, 3.166, 3.124, 3.344, 3.192, 3.387, 3.362, 3.205,
              3.399, 3.349, 3.55, 3.549, 870, 870, 893, 893, 870, 870, 892, 902, 849,
              849, 902, 870, 870, 870, 870, 870, 892, 892, 870, 881, 881, 892, 892,
              902, 902, 908, 902, 902, 933, 965, 945, 945, 934, 934, 955, 945, 913,
              913, 977, 934, 913, 913, 913, 934, 945, 945, 913, 913, 913, 945, 945,
              945, 945, 977, 945, 945), ncol=5)
  d<-matrix(c(1965, 1966, 1969, 1971, 1972, 1973, 1979, 1983, 1984, 1988, 1989, 1990,
              1993, 1993, 1993, 1996, 1997, 1998, 1998, 2002, 2003, 2004, 2006, 2007,
              2010, 2011, 2013, 2013),ncol=1)
  # go
  target.arrival.dea(x,y,d,2007,"d","vrs","o","min")$arrival_avg
}