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
  list with components of \dQuote{bar.x} (the value of
  \code{hist()}), \dQuote{kvals} (window widths plotted),
  and \dQuote{pi.0} (the input vlue of \code{pi.0}
}
\description{
  Plot a set of cutpoints
}
\details{
  not for users. not exported
}
\author{
  Charles Berry
}

