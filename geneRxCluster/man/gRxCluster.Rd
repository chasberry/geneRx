\name{gRxCluster}
\alias{gRxCluster}
\title{gRxCluster}
\usage{
  gRxCluster(object, starts, group, kvals, nperm = 0L,
    pruneFun = prune.loglik, ..., cutpt.filter.expr,
    cutpt.tail.expr, tmp.env, sample.id, sample.tab)
}
\arguments{
  \item{object}{chromosome names or other grouping of
  starts}

  \item{starts}{ordered chromosome position or ordered
  integer vector}

  \item{group}{logical vector separating two groups}

  \item{kvals}{integer vector of window widths}

  \item{nperm}{number of permutations for FDR calculation}

  \item{pruneFun}{a function like \code{prune.loglik}.}

  \item{...}{other args}

  \item{cutpt.filter.expr}{(optional) R object or call (or
  variable naming a call) with (optional) var x (window
  widths in base pairs) to filter windows. It must evaluate
  to mode "double". If not specified,
  \code{as.double(apply(x,2,median,na.rm=TRUE))} is used.
  If an atomic vector of length one is supplied it is
  expanded to the proper length and coerced to double. If
  this arg is the name of a variable provided in
  \code{tmp.env}, it must be protected with
  \code{quote(\dots)}.}

  \item{cutpt.tail.expr}{R object or call (or variable
  naming a call) with (optional) vars: k,n, and x (as
  above). Returns list like \code{critVal.target}. k is a
  vector of the number of sites in a collection of windows,
  and n is a vector of counts or proportions for the two
  classes of insertion. If not supplied,
  \code{critVal.target(k,n,target=5,posdiff=x)} is used.
  If this arg is the name of a variable provided in
  \code{tmp.env}, it must be protected with
  \code{quote(\dots)}.}

  \item{tmp.env}{(optional) environment to use in
  evaluation of cutpt.* expressions. This is usually needed
  for \code{\link{critVal.power}}, which is first
  calculated and placed in the environment, and the
  supplied object is used in the expression for
  \code{cutpt.filter.expr}.}

  \item{sample.id}{(optional) integer vector indexing cells
  in \code{sample.tab} to be looked up to determine
  \code{group} under permutation. A factor can be used,
  too, but will be coerced to integer.}

  \item{sample.tab}{(optional) integer vector containing 0
  or 1 in each cell. Its length is the same as
  \code{max(sample.id)}. Both or neither \code{sample.id}
  and \code{sample.tab} should be supplied. When supplied
  \code{sample.tab[sample.id]} must equal \code{group}. If
  the arguments are supplied, permutations are of the form
  \code{sample(sample.tab)[sample.id]}. Otherwise they are
  of the form \code{sample(group)}.}
}
\value{
  a GRanges object with a special metadata slot, see
  \code{\link{gRxCluster-object}}
}
\description{
  cluster integration sites - optionally perform the
  permutations needed to estimate the discoveries expected
  under a null hypothesis
}
\author{
  Charles Berry example inst/ex-gRxCluster.R
}

