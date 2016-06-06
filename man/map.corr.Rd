\name{map.corr}
\alias{map.corr}
\title{Correlation mapping for reliability test}
\description{Implements a series of correlation analysis by dropping extreme data points one by one using \emph{Mahalanobis} distance measure. Correlation reliability can be investigated with identified anchoring point(s). Correlation map as well as summary table is provided.}
\usage{map.corr(data, from="median", threshold=0.3, r.name=FALSE)}
\arguments{
  \item{data}{Dataframe}
  \item{from}{Datum point from which the distance is measured \cr
    \code{"mean"    }{Mean of each column} \cr  
    \code{"median"  }{Median of each column (default)}
  }
  \item{threshold}{Threshold of correlation change to be noted on the map}
  \item{r.name}{Dropped points are shown in row name when \code{TRUE}}
}
\value{
  \item{\code{$reliability}}{Summary table}
}
\seealso{
  \code{\link{dm.mahalanobis}} Distance measure using Mahalanobis distance
}
\author{Dong-Joon Lim, PhD}
\examples{
# Generate a sample dataframe
df<-data.frame(replicate(2,sample(0:100,50)))

# go
map.corr(df)
}