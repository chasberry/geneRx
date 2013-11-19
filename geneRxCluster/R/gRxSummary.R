##' Summarize gRxCluster Results
##'
##' Get the FDR and related data for a run of gRxCluster. By selecting
##' a value for \code{targetFD} that is smaller that what was used in
##' constructing the object, fewer clumps will be included in the
##' computation fo the False Discovery Rate - akin to what would have
##' been obtained from the object if it had been constructed using
##' that value.
##' 
##' @title gRxSummary
##' @param object the result of gRxCluster
##' @param targetFD the critical value target in each tail
##' @return a list containing the summarized results
##' @export
##' @example inst/ex-gRxSummary.R
##' @include gRxCluster.R
##' @author Charles Berry
gRxSummary<-
  function(object, targetFD=NULL)
{
  if (!(
    inherits(object,"GRanges") &&
    {
      mc <- metadata(object)$call
      isTRUE( mc[[1]] == quote(gRxCluster))
    }))
    stop("This object was not made by gRxCluster.")
  
  
  if (is.null(targetFD)){
      ## get target from $criticalValues
      targetFD <- max(sapply(metadata(object)$criticalValues,
                         function(x) attr(x,"target")))
  }
  
  nperms <- as.list(mc)[["nperm"]]
  
  perms <- metadata(object)$perm_cluster_best
  nd <- sum(object$target.min<=targetFD)
  npd <- sum(unlist(perms)<=targetFD)
  fdr <- sum(npd)/nperms/(1+nd)
  
  res <-
    list(
      Clusters_Discovered=nd,
      FDR=fdr,
      permutations=nperms,
      targetFD=
      if (is.infinite(targetFD)) "UnDetermined" else targetFD,
      call=mc)
  
  res
  
}
