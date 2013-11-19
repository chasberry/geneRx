\name{gRxCluster-object}
\alias{gRxCluster-object}
\title{gRxCluster object }
\description{Overview of the result of \code{gRxCluster(\dots)}
}

\details{

  The object returned is a \code{\link{GRanges}} object.

  If the object is \code{x}, \code{seqnames(x)} and \code{ranges(x)}
  slots demarcate the clusters discovered. There will be one element for
  each cluster (aka \sQuote{clump}) discovered.
  
  Using the default argument \code{pruneFun=prune.loglik} or
  \code{pruneFun=noprune}, \code{mcols(x)} will have these 
  columns:
  
  \describe{
    \item{\code{value1} and
      \code{value2}}{ are the counts of the two classes of
      insertion sites for the clusters of object \code{x}}
    \item{\code{clump.id}}{numbers each cluster.}
  }
  
  If the user supplies a custom \code{pruneFun}, it should return a
  \code{GRanges} with those columns and one element for each unique
  \code{clump.id}. The column \code{target.min} has the smallest nominal
  False Discoveries Expected for each cluster and is added to (or
  replaces) the \code{mcols(x)} produced by the argument supplied as
  \code{pruneFun}.

  
  \code{metadata(x)} will include these components:

  \describe{
    \item{criticalValues}{A list object such as supplied by
      \code{\link{critVal.target}} whose elements each give the cutpoints
      to be used for a window with \code{k} sites.
      \code{attributes(metadata(object)$criticalValues[[i]])} will
      contain elements
      \describe{
        
        \item{fdr}{ with dimension \code{c(k+1,4)} of target false
          discovery expectations and and the one-sided p-values}
        \item{target}{ the target for false discovery which sometimes
          is specified a priori and sometimes results from calculation}
        \item{n}{an upper bound on the number of windows to screen, if
          this number is needed.}  }
      
      In some cases, an attribute is attached to
      \code{metadata(object)$criticalValues}, see
      \code{\link{critVal.power}} for an example.
    }
  \item{kvals}{the number of sites, k, to include in a window}
  \item{perm_cluster_best}{a \code{list} whose canonical element is a vector of
    values like \code{x$target.min} obtained from a permutation of the
    class indicators}
  \item{summary_matrix}{a matrix giving the start, end, depth, and
    counts in each class for every cluster and depth in sequential order }
  \item{call}{the call invoking \code{\link{gRxCluster}}} which may
  include some arguments added by default.            
}
}
\author{Charles Berry \email{ccberry@ucsd.edu}}
\keyword{cluster}
