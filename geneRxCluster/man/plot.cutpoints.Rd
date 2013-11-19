\name{plot.cutpoints}
\alias{plot.cutpoints}
\title{plot.cutpoints}
\usage{
  plot.cutpoints(crit, pi.0 = NULL, kvals = NULL, ...)
}
\arguments{
  \item{crit}{- a cutpoint object see
  \code{\link{gRxCluster}}}

  \item{pi.0}{- optional null value to plot}

  \item{kvals}{- which cutpoints to includein the plot}

  \item{...}{passed to barplot}
}
\value{
  list with components of \dQuote{\code{bar.x}} (the value
  of \code{hist()}), \dQuote{\code{kvals}} (window widths
  plotted), and \dQuote{\code{pi.0}} (the input value of
  \code{pi.0})
}
\description{
  Plot a set of cutpoints - Utility
}
\details{
  NOT FOR USERS. Not exported.
}
\author{
  Charles Berry
}

