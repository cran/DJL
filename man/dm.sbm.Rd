\name{dm.sbm}
\alias{dm.sbm}
\title{Distance measure using SBM}
\description{Implements \emph{Tone}'s slack-based model (non-radial & (non-)oriented measure).}
\usage{dm.sbm(xdata, ydata, rts="crs", 
       orientation="n", se=FALSE, sg="ssm", date=NULL, cv="convex")}
\arguments{
  \item{xdata}{Input(s) vector (\emph{n} by \emph{m})}
  \item{ydata}{Output(s) vector (\emph{n} by \emph{s})}
  \item{rts}{Returns to scale assumption \cr
    \code{"crs"  }{Constant RTS (default)} \cr  
    \code{"vrs"  }{Variable RTS} \cr
    \code{"irs"  }{Increasing RTS} \cr
    \code{"drs"  }{Decreasing RTS}
  }
    \item{orientation}{Orientation of the measurement \cr
    \code{"n"  }{Non-orientation (default)} \cr
    \code{"i"  }{Input-orientation} \cr
    \code{"o"  }{Output-orientation}
  }
  \item{se}{Implements super-efficiency model alike \emph{Anderson & Peterson}'s model if \code{TRUE}}
  \item{sg}{Employs second-stage optimization \cr
    \code{"ssm"  }{Slack-sum maximization (default)} \cr
    \code{"max"  }{Date-sum maximization (only if \code{date} is defined)} \cr
    \code{"min"  }{Date-sum minimization (only if \code{date} is defined)}
  }
  \item{date}{Production date (\emph{n} by \emph{1})}
  \item{cv}{Convexity assumption \cr
    \code{"convex" }{Convexity holds (default)} \cr
    \code{"fdh"    }{Free disposal hull (this will override \code{rts})}
  }
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
  X <- data.frame(x1 = c(4, 6, 8, 8, 2),
                  x2 = c(3, 3, 1, 1, 4))
  Y <- data.frame(y1 = c(2, 2, 6, 6, 1),
                  y2 = c(3, 3, 2, 1, 4))
  # go
  dm.sbm(X, Y)

# Reproduce Table 1 in Tone.(2002)
  # Published input slacks are alternate optima (confirmed by Tone)
  # ready
  X <- data.frame(x1 = c(4, 7, 8, 4, 2, 10, 12),
                  x2 = c(3, 3, 1, 2, 4, 1,  1))
  Y <- data.frame(y1 = c(1, 1, 1, 1, 1, 1,  1))
  # go
  dm.sbm(X, Y, se = TRUE)

# Reproduce Table 4 in Tone.(2002)
  # ready
  X <- data.frame(x1 = c(80,  65,  83,  40,   52,  94),
                  x2 = c(600, 200, 400, 1000, 600, 700),
                  x3 = c(54,  97,  72,  75,   20,  36),
                  x4 = c(8,   1,   4,   7,    3,   5))
  Y <- data.frame(y1 = c(90,  58,  60,  80,   72,  96),
                  y2 = c(5,   1,   7,   10,   8,   6))
  # go
  dm.sbm(X, Y, "crs", "i", se = TRUE)
}