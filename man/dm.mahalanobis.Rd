\name{dm.mahalanobis}
\alias{dm.mahalanobis}
\title{Distance measure using Mahalanobis distance for outlier detection}
\description{Implements \emph{Mahalanobis} distance measure for outlier detection. In addition to the basic distance measure, boxplots are provided with potential outlier(s) to give an insight into the early stage of data cleansing task.}
\usage{dm.mahalanobis(data, from="mean", p=10, plot=FALSE, v.index=NULL, layout=NULL)}
\arguments{
  \item{data}{Dataframe}
  \item{from}{Datum point from which the distance is measured \cr
    \code{"mean"    }{Mean of each column (default)} \cr  
    \code{"median"  }{Median of each column}
  }
  \item{p}{Percentage to which outlier point(s) is noted (default of 10)}
  \item{plot}{Switch for boxplot(s)}
  \item{v.index}{Numeric vector indicating column(s) to be printed in the boxplot.
        Default value of NULL will present all.}
  \item{layout}{Numeric vector indicating dimension of boxplots.
        Default value of NULL will find an optimal layout.}
}
\value{
  \item{\code{$dist}}{\code{Mahalanobis} distance from \code{from}}
  \item{\code{$order}}{Distance order (decreasing) in row number}
  \item{\code{$suspects}}{Potential outlier(s) in row number}
}
\references{Hair, Joseph F., et al. Multivariate data analysis. Vol. 7. \emph{Upper Saddle River}, NJ: Pearson Prentice Hall, 2006.}
\author{Dong-Joon Lim, PhD}
\examples{
# Generate a sample dataframe
df<-data.frame(replicate(6,sample(0:100,50)))

# go
dm.mahalanobis(df,plot=TRUE)
}