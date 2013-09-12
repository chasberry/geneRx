\name{gRxCluster}
\alias{gRxCluster}
\title{gRxCluster}
\usage{
  gRxCluster(object, starts, group, kvals, nperm = 0L,
    pruneFun = prune.loglik, ..., cutpt.filter.expr,
    cutpt.tail.expr, tmp.env)
}
\arguments{
  \item{object}{chromosome names or other grouping of
  starts}

  \item{starts}{ordered chromosome position or ordered
  integer vector}

  \item{group}{logical vector separating two groups}

  \item{kvals}{integer vector of window widths}

  \item{nperm}{number of permutations for FDR calculation}

  \item{pruneFun}{a function like prune.loglik}

  \item{...}{other args}

  \item{cutpt.filter.expr}{(optional) R expr with
  (optional) var x (window widths in base pairs) to filter
  windows must eval to mode "double"}

  \item{cutpt.tail.expr}{R expr with (optional) vars: k,n,
  and x (as above). Returns list. k is a vector of the
  number of sites in a collection of windows, and n is a
  vector of counts or proportions for the two classes of
  insertion}

  \item{tmp.env}{(optional) environment to use in
  evaluation of cutpt.* expressions. This is usually needed
  for \code{\link{critVal.power}}, which is first
  calculated and placed in the environment, and the
  supplied object is used in the expression for
  \code{cutpt.filter.expr}.}
}
\value{
  a GRanges object with a special metadata slot, see
  \code{\link{gRxCluster-object}}
}
\description{
  cluster integration sites
}
\examples{

# examples invoking gRxCluster

x.seqnames <- rep(letters[1:3],each=50)
x.starts <- c(seq(1,length=50),seq(1,by=2,length=50),seq(1,by=3,length=50))
x.lens <- rep(c(5,10,15,20,25),each=2)
x.group <- rep(rep(c(TRUE,FALSE),length=length(x.lens)),x.lens)
x.kvals <- as.integer(sort(unique(x.lens)))
x.res <- gRxCluster(x.seqnames,x.starts,x.group,x.kvals)
x.res
rm( x.seqnames, x.starts, x.lens, x.group, x.kvals, x.res)
}
\author{
  Charles Berry
}

