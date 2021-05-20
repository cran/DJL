\name{dm.network.dea}
\alias{dm.network.dea}
\title{Distance measure using DEA on a two-stage network structure}
\description{Implements \emph{Charnes & Cooper}'s data envelopment analysis (radial & oriented measure) on a two-stage network structure.}
\usage{dm.network.dea(xdata.s1, ydata.s1=NULL, zdata, xdata.s2=NULL, ydata.s2,
               rts="crs", orientation="i", type="nc", leader="1st", ss=10^-4, o=NULL)}
\arguments{
  \item{xdata.s1}{Input(s) vector in Stage 1 (\emph{n} by \emph{m.s1})}
  \item{ydata.s1}{Output(s) vector in Stage 1 (\emph{n} by \emph{s.s1})}
  \item{zdata   }{Intermediate product(s) vector between Stage 1 and Stage 2 (\emph{n} by \emph{p})}
  \item{xdata.s2}{Input(s) vector in Stage 2 (\emph{n} by \emph{m.s2})}
  \item{ydata.s2}{Output(s) vector in Stage 2 (\emph{n} by \emph{s.s2})}
  \item{rts}{Returns to scale assumption \cr
    \code{"crs"  }{Constant RTS (default)} \cr  
    \code{"vrs"  }{Variable RTS} \cr
    \code{"irs"  }{Increasing RTS} \cr
    \code{"drs"  }{Decreasing RTS}
  }
  \item{orientation}{Orientation of the measurement \cr
    \code{"i"  }{Input-orientation (default)} \cr
    \code{"o"  }{Output-orientation}
  }
  \item{type}{Solution method \cr
    \code{"nc"  }{Decentralized model (Stackelberg game approach) (default)} \cr
    \code{"co"  }{Centralized model (cooperative game approach)}
  }
  \item{leader}{Preemptive priority for Decentralized model \cr
    \code{"1st" }{1st stage as the leader (default)} \cr
    \code{"2nd" }{2nd stage as the leader}
  }
  \item{ss}{Step size for heuristic search \code{10^-4} (default)}
  \item{o}{DMU index to operate. \code{NULL} (default) will operate for all}
}
\value{
  \item{\code{$eff.s1}}{Efficiency score of Stage 1}
  \item{\code{$eff.s2}}{Efficiency score of Stage 2}
  \item{\code{$v.s1}}{Weight attached to input in Stage 1}
  \item{\code{$u.s1}}{Weight attached to output in Stage 1}
  \item{\code{$p}}{Weight attached to intermediate product}
  \item{\code{$w.s1}}{Free variable for scaling in Stage 1}
  \item{\code{$v.s2}}{Weight attached to input in Stage 2}
  \item{\code{$u.s2}}{Weight attached to output in Stage 2}
  \item{\code{$w.s2}}{Free variable for scaling in Stage 2}
}
\references{
Kao, Chiang, and Shiuh-Nan Hwang. "Efficiency decomposition in two-stage data envelopment analysis: An application to non-life insurance companies in Taiwan." \emph{European journal of operational research} 185.1 (2008): 418-429.

Cook, Wade D., Liang Liang, and Joe Zhu. "Measuring performance of two-stage network structures by DEA: a review and future perspective." \emph{Omega} 38.6 (2010): 423-430.

Li, Yongjun, Yao Chen, Liang Liang, and Jianhui Xie. "DEA models for extended two-stage network structures." \emph{Omega} 40.5 (2012): 611-618.

Lee, Hsuan-Shih. "Efficiency decomposition of the network DEA in variable returns to scale: An additive dissection in losses." \emph{Omega} 100 (2021): 102212.
}
\author{Dong-Joon Lim, Ph.D.}
\seealso{
  \code{\link{dm.dea}} Distance measure using DEA \cr
  \code{\link{dm.dynamic.bc}} Dynamic DEA for intertemporal budgeting
}
\examples{
# Reproduce Table 2 in W.D. Cook et al.(2010)
  # ready
  X <- data.frame(x1 = c(1178744,1381822,1177494,601320,6699063,2627707,1942833,3789001,
                         1567746,1303249,1962448,2592790,2609941,1396002,2184944,1211716,
                         1453797,757515,159422,145442,84171,15993,54693,163297,1544215),
                  x2 = c(673512,1352755,592790,594259,3531614,668363,1443100,1873530,
                         950432,1298470,672414,650952,1368802,988888,651063,415071,
                         1085019,547997,182338,53518,26224,10502,28408,235094,828963))
  Z <- data.frame(z1 = c(7451757,10020274,4776548,3174851,37392862,9747908,10685457,17267266,
                         11473162,8210389,7222378,9434406,13921464,7396396,10422297,5606013,
                         7695461,3631484,1141950,316829,225888,52063,245910,476419,7832893),
                  z2 = c(856735,1812894,560244,371863,1753794,952326,643412,1134600,
                         546337,504528,643178,1118489,811343,465509,749893,402881,
                         342489,995620,483291,131920,40542,14574,49864,644816,667964))
  Y <- data.frame(y1 = c(984143,1228502,293613,248709,7851229,1713598,2239593,3899530,
                         1043778,1697941,1486014,1574191,3609236,1401200,3355197,854054,
                         3144484,692731,519121,355624,51950,82141,0.1,142370,1602873),
                  y2 = c(681687,834754,658428,177331,3925272,415058,439039,622868,
                         264098,554806,18259,909295,223047,332283,555482,197947,
                         371984,163927,46857,26537,6491,4181,18980,16976,477733))
  
  # go
  res.co    <- dm.network.dea(xdata.s1 = X, zdata = Z, ydata.s2 = Y, type = "co")
  res.nc.LF <- dm.network.dea(xdata.s1 = X, zdata = Z, ydata.s2 = Y, type = "nc", leader = "1st")
  res.nc.FL <- dm.network.dea(xdata.s1 = X, zdata = Z, ydata.s2 = Y, type = "nc", leader = "2nd")

  # print
  data.frame(CO.s1    = res.co$eff.s1,
             CO.s2    = res.co$eff.s2,
             NC.LF.s1 = res.nc.LF$eff.s1,
             NC.LF.s2 = res.nc.LF$eff.s2,
             NC.FL.s1 = res.nc.FL$eff.s1,
             NC.FL.s2 = res.nc.FL$eff.s2)
}