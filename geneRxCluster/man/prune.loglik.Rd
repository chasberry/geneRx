\name{prune.loglik}
\alias{prune.loglik}
\title{prune.loglik}
\usage{
  prune.loglik(x, p.null = 0.5)
}
\arguments{
  \item{x}{a GRanges object}

  \item{p.null}{the probability of category 1 (FALSE)}
}
\description{
  best contiguous region
}
\details{
  prune each end of the region using loglik criterion

  this is to be used as the \code{pruneFun} are of
  \code{gRxCluster}
}
\author{
  Charles Berry
}
\seealso{
  \code{\link{gRxCluster-object}} for details on what this
  function returns.
}

