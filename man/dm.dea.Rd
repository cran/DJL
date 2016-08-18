\name{dm.dea}
\alias{dm.dea}
\title{Distance measure using DEA}
\description{Implements \emph{Charnes & Cooper}'s data envelopment analysis (radial & oriented measure).}
\usage{dm.dea(xdata, ydata, rts="crs", orientation,
       se=FALSE, sg="ssm", date=NULL, ncv=NULL, env=NULL, cv="convex")}
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
    \code{"i"  }{Input-orientation} \cr
    \code{"o"  }{Output-orientation}
  }
  \item{se}{Implements \emph{Anderson & Peterson}'s super-efficiency model if \code{TRUE}}
  \item{sg}{Employs second-stage optimization \cr
    \code{"ssm"  }{Slack-sum maximization (default)} \cr
    \code{"max"  }{Date-sum maximization (only if \code{date} is defined)} \cr
    \code{"min"  }{Date-sum minimization (only if \code{date} is defined)}
  }
  \item{date}{Production date (\emph{n} by \emph{1})}
  \item{ncv}{Non-controllable variable index(binary) for internal NDF (\emph{1} by \emph{(m+s)})}
  \item{env}{Environment index for external NDF (\emph{n} by \emph{1})}
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
  \item{\code{$vx}}{Input (dual) weight}
  \item{\code{$uy}}{Output (dual) weight}
  \item{\code{$w}}{Free (dual) variable}
}
\references{
Charnes, Abraham, William W. Cooper, and Edwardo Rhodes. "Measuring the efficiency of decision making units." \emph{European journal of operational research} 2.6 (1978): 429~444.

Charnes, Abraham, William W. Cooper, and Edwardo Rhodes. "Evaluating program and managerial efficiency: an application of data envelopment analysis to program follow through." \emph{Management science} 27.6 (1981): 668~697.

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
  X <- data.frame(x1 = c(8, 6,   3,   10, 6,   8,   8, 4),
                  x2 = c(8, 4.6, 1.9, 9,  3.6, 3.6, 9, 1.9))
  Y <- data.frame(y1 = c(8, 5,   2,   9,  4.5, 4.5, 7, 2))
  C <- data.frame(x1 = 0, x2 = 1, y1 = 0)

  # go
  data.frame(ALL_CRS = dm.dea(X, Y, "crs", "i")$eff,
             ALL_VRS = dm.dea(X, Y, "vrs", "i")$eff,
             NDF_CRS = dm.dea(X, Y, "crs", "i", ncv = C)$eff,
             NDF_VRS = dm.dea(X, Y, "vrs", "i", ncv = C)$eff,
             row.names = LETTERS[1 : 8])
}