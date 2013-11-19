\name{gRxSummary}
\alias{gRxSummary}
\title{gRxSummary}
\usage{
  gRxSummary(object, targetFD = NULL)
}
\arguments{
  \item{object}{the result of gRxCluster}

  \item{targetFD}{the critical value target in each tail}
}
\value{
  a list containing the summarized results
}
\description{
  Summarize gRxCluster Results
}
\details{
  Get the FDR and related data for a run of gRxCluster. By
  selecting a value for \code{targetFD} that is smaller
  that what was used in constructing the object, fewer
  clumps will be included in the computation fo the False
  Discovery Rate - akin to what would have been obtained
  from the object if it had been constructed using that
  value.
}
\examples{
x.seqnames <- rep(letters[1:3],each=50)
x.starts <- c(seq(1,length=50),seq(1,by=2,length=50),seq(1,by=3,length=50))
x.lens <- rep(c(5,10,15,20,25),each=2)
x.group <- rep(rep(c(TRUE,FALSE),length=length(x.lens)),x.lens)
x.kvals <- as.integer(sort(unique(x.lens)))
x.res <- gRxCluster(x.seqnames,x.starts,x.group,x.kvals,nperm=100L)
gRxSummary(x.res)
rm( x.seqnames, x.starts, x.lens, x.group, x.kvals, x.res)
}
\author{
  Charles Berry
}

