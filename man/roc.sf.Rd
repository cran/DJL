\name{roc.sf}
\alias{roc.sf}
\title{Rate of change (RoC) calculation using SF}
\description{Employs \code{dm.sf} over time to calculate RoCs. This function is valid only when multiplicative form of directional vector is used.}
\usage{roc.sf(xdata, ydata, date, t, 
       rts="crs", g=NULL, wd=NULL, sg="ssm", ftype="d", cv="convex")}
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
  \code{\link{dm.sf}} Distance measure using SF \cr
  \code{\link{roc.sf}} RoC calculation using SF \cr
  \code{\link{map.soa.sf}} SOA mapping using SF \cr
  \code{\link{target.arrival.sf}} Arrival target setting using SF
}
\examples{
# Reproduce Mercedes-Benz CLA45 AMG's local RoC in Table 5 in Lim, D-J. (2015)
  # Load engine dataset
    df <- dataset.engine.2015
  
  # Subset for 4 cylinder engines
    fce <- subset(df, df[, 3] == 4)
    
  # Parameters
    x <- subset(fce, select = 4)
    y <- subset(fce, select = 5 : 7)
    d <- subset(fce, select = 2)
    g <- as.matrix(data.frame(0, y))
    w <- matrix(c(1, 0, 0), ncol = 3)

  # Calc local Roc
    roc.sf(x, y, d, 2014, "crs", g, w, "min")$roc_local[348, ]
}