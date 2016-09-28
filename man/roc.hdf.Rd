\name{roc.hdf}
\alias{roc.hdf}
\title{Rate of change (RoC) calculation using HDF}
\description{Employs \code{dm.hdf} over time to calculate RoCs.}
\usage{roc.hdf(xdata, ydata, date, t,
        rts="crs", wd=NULL, sg="ssm", ftype="d", cv="convex")}
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
D.-J. Lim,  Internal combustion engine race: naturally aspirated vs turbo/super-charged, \emph{working paper} (2015).
}
\author{Dong-Joon Lim, PhD}
\seealso{
  \code{\link{dm.hdf}} Distance measure using HDF \cr
  \code{\link{roc.hdf}} RoC calculation using HDF \cr
  \code{\link{map.soa.hdf}} SOA mapping using HDF \cr
  \code{\link{target.arrival.hdf}} Arrival target setting using HDF
}
\examples{
# Load engine dataset
  df <- dataset.engine.2015

# Subset for 8 cylinder TC-P engines
  et <- subset(df, df[, 3] == 8 & df[, 8] == "TC-P")
  
# Parameters
  x <- subset(et, select = 4)
  y <- subset(et, select = 5 : 7)
  d <- subset(et, select = 2)
  w <- matrix(c(1, 0, 0), ncol = 3)

# Calc local Roc
  roc.hdf(x, y, d, 2015, "vrs", w, "min")
}