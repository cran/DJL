\name{target.arrival.sf}
\alias{target.arrival.sf}
\title{Arrival target setting using SF}
\description{Employs \code{dm.sf} over time to estimate the arrival of known specifications. This function is valid only when multiplicative form of directional vector is used.}
\usage{target.arrival.sf(xdata, ydata, date, t, rts="crs", g=NULL, 
                  wd=NULL, sg="ssm", ftype="d", cv="convex", anc=FALSE)}
\arguments{
  \item{xdata}{Input(s) vector (\emph{n} by \emph{m})}
  \item{ydata}{Output(s) vector (\emph{n} by \emph{s})}
  \item{date}{Production date (\emph{n} by \emph{1})}
  \item{t}{A vantage point from which the RoC is captured}
  \item{rts}{Returns to scale assumption \cr
    \code{"crs"  }{Constant RTS (default)} \cr  
    \code{"vrs"  }{Variable RTS} \cr
    \code{"irs"  }{Increasing RTS} \cr
    \code{"drs"  }{Decreasing RTS}
  }
  \item{g}{Directional vector indicating a measurement direction (\emph{n} by \emph{(m+s)})\cr
           By default (\code{NULL}), \code{xdata} & \code{ydata} will be used}
  \item{wd}{Weak disposability vector indicating (an) undesirable output(s) (\emph{1} by \emph{s})} 
  \item{sg}{Employs second-stage optimization \cr
    \code{"ssm"  }{Slack-sum maximization (default)} \cr
    \code{"max"  }{Date-sum maximization} \cr
    \code{"min"  }{Date-sum minimization}
  }
  \item{ftype}{Frontier type \cr
    \code{"d"  }{Dynamic frontier (default)} \cr  
    \code{"s"  }{Static frontier}
  } 
  \item{cv}{Convexity assumption \cr
    \code{"convex" }{Convexity holds (default)} \cr
    \code{"fdh"    }{Free disposal hull (this will override \code{rts})}
  }
  \item{anc}{Implements stepwise RoC computation if \code{TRUE}}
}
\value{
  \item{\code{$eff_t}}{Efficiency at \code{t}}
  \item{\code{$lambda_t}}{Intensity vector at \code{t}}
  \item{\code{$eft_date}}{Effective date}
  \item{\code{$roc_avg}}{Average RoC}
  \item{\code{$roc_anc}}{Local RoCs across the periods}
  \item{\code{$roc_local}}{Local RoC}
  \item{\code{$roc_ind}}{Individualized RoC}
  \item{\code{$arrival_avg}}{Estimated arrival using \code{roc_avg}}
  \item{\code{$arrival_seg}}{Estimated arrival using \code{roc_ind}}
}
\references{
Lim, Dong-Joon, et al. "Comparing technological advancement of hybrid electric vehicles (HEV) in different market segments." \emph{Technological Forecasting and Social Change} 97 (2015): 140~153.

Lim, Dong-Joon, and Timothy R. Anderson. Time series benchmarking analysis for new product scheduling: who are the competitors and how fast are they moving forward?. Advances in DEA Theory and Applications: with Examples in Forecasting Models. \emph{Wiley (forthcoming)}, 2016.
}
\author{Dong-Joon Lim, PhD}
\seealso{
  \code{\link{dm.sf}} Distance measure using SF \cr
  \code{\link{roc.sf}} RoC calculation using SF \cr
  \code{\link{map.soa.sf}} SOA mapping using SF \cr
  \code{\link{target.arrival.sf}} Arrival target setting using SF
}
\examples{
# Estimate arrivals of MY2013 hev models
  # Load hev dataset
  df <- dataset.hev.2013
  
  # ready
  x <- subset(df, select = 3)
  y <- subset(df, select = 4 : 6)
  d <- subset(df, select = 2)
  g <- data.frame(x, y)
  
  # go
  target.arrival.sf(x, y, d, 2012, "vrs", g)$arrival_seg
}