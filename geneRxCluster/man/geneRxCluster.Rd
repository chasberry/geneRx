\docType{package}
\name{geneRxCluster}
\alias{geneRxCluster}
\alias{geneRxCluster-package}
\title{Differenital Clustering of Integration Sites}
\description{
  geneRxCluster provides the function
  \code{\link{gRxCluster}} and friends.
}
\details{
  Windows defined by \code{k} consecutive integration sites
  are scanned. A two class indicator is tallied to
  determine whether one class dominates. If one does, a
  flag is set and the window is retained. Various values of
  \code{k} are used. Conflicts between overlapping windows
  with the same value of \code{k} can occur --- two windows
  are dominated by the two different classes. In that case,
  the sites of overlap are marked and neither window is
  retained. Conflicts can also arise between windows
  differing in their values of \code{k}. In that case, the
  window having the smaller value of \code{k} is retained
  and the other is discarded.

  Permutation tests and permutation based false discovery
  rates are available.

  Filtering of windows is allowed so that regions which are
  sparsely populated need not be studied.
}

