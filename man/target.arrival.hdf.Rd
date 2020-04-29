\name{target.arrival.hdf}
\alias{target.arrival.hdf}
\title{Arrival target setting using HDF}
\description{Employs \code{dm.hdf} over time to estimate the arrival of known specifications.}
\usage{target.arrival.hdf(xdata, ydata, date, t, rts="crs", 
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
  \item{anc}{Implements a stepwise RoC computation if \code{TRUE}}
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

Lim, Dong-Joon, and Timothy R. Anderson. Time series benchmarking analysis for new product scheduling: who are the competitors and how fast are they moving forward?. \emph{Advances in DEA Theory and Applications: with Examples in Forecasting Models}. (2017): 443-458.
}
\author{Dong-Joon Lim, PhD}
\seealso{
  \code{\link{dm.hdf}} Distance measure using HDF \cr
  \code{\link{roc.hdf}} RoC calculation using HDF \cr
  \code{\link{map.soa.hdf}} SOA mapping using HDF \cr
  \code{\link{target.arrival.hdf}} Arrival target setting using HDF
}
\examples{
# Estimate arrivals of MY2015 SC/TC 8 cylinder engines
  # Load engine dataset
    df <- dataset.engine.2015
  
  # Subset for SC/TC 8 cylinder engines
    stc.8 <- subset(df, grepl("^.C..", df[, 8]) & df[, 3] == 8)
  
  # Parameters
    x <- subset(stc.8, select = 4)
    y <- subset(stc.8, select = 5:7)
    d <- subset(stc.8, select = 2)
  
  # Generate an SOA map
    target.arrival.hdf(x, y, d, 2014, "vrs")
}