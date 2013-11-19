\name{noprune}
\alias{noprune}
\title{noprune}
\usage{
  noprune(x, ...)
}
\arguments{
  \item{x}{a GRanges object}

  \item{...}{currently unused}
}
\value{
  same as \code{\link{gRxCluster}} less the \code{metadata}
}
\description{
  join contiguous windows
}
\details{
  return all the candidate sites in a clump without pruning

  this is to be used as the \code{pruneFun} are of
  \code{gRxCluster}
}
\author{
  Charles Berry
}
\seealso{
  \code{\link{gRxCluster-object}} for more details on what
  this function returns.
}

