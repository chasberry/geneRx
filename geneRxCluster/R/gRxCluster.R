##' gRxCluster
##' 
##' cluster integration sites - optionally perform the permutations
##' needed to estimate the discoveries expected under a null
##' hypothesis
##' @param object chromosome names or other grouping of starts 
##' @param starts ordered chromosome position or ordered integer
##' vector 
##' @param group logical vector separating two groups
##' @param kvals integer vector of window widths
##' @param nperm number of permutations for FDR calculation
##' @param pruneFun a function like \code{prune.loglik}. 
##' @param ... other args
##' 
##' @param cutpt.filter.expr (optional) R object or call (or variable
##' naming a call) with (optional) var x (window widths in base pairs)
##' to filter windows. It must evaluate to mode "double". If not specified,
##' \code{as.double(apply(x,2,median,na.rm=TRUE))} is used. If an
##' atomic vector of length one is supplied it is expanded to the
##' proper length and coerced to double. If this arg is the name of a
##' variable provided in \code{tmp.env}, it must be protected with
##' \code{quote(\dots)}.
##' 
##' @param cutpt.tail.expr R object or call (or variable naming a
##' call) with (optional) vars: k,n, and x (as above). Returns list
##' like \code{critVal.target}. k is a vector of the number of sites
##' in a collection of windows, and n is a vector of counts or
##' proportions for the two classes of insertion. If not supplied,
##' \code{critVal.target(k,n,target=5,posdiff=x)} is used.  If this
##' arg is the name of a variable provided in \code{tmp.env}, it must
##' be protected with \code{quote(\dots)}.
##'
##' @param tmp.env (optional) environment to use in evaluation of
##' cutpt.* expressions. This is usually needed for
##' \code{\link{critVal.power}}, which is first calculated and placed
##' in the environment, and the supplied object is used in the
##' expression for \code{cutpt.filter.expr}.
##'
##' 
##' @param sample.id (optional) integer vector indexing cells in
##' \code{sample.tab} to be looked up to determine \code{group} under
##' permutation. A factor can be used, too, but will be coerced to
##' integer.
##' 
##' @param sample.tab (optional) integer vector containing 0 or 1 in
##' each cell. Its length is the same as \code{max(sample.id)}. Both
##' or neither \code{sample.id} and \code{sample.tab} should be
##' supplied. When supplied \code{sample.tab[sample.id]} must equal
##' \code{group}. If the arguments are supplied, permutations are of
##' the form \code{sample(sample.tab)[sample.id]}. Otherwise they are
##' of the form \code{sample(group)}.
##' 
##' @return a GRanges object with a special metadata slot, see
##' \code{\link{gRxCluster-object}}
##' 
##' @author Charles Berry
##' example inst/ex-gRxCluster.R
##' @export
##' @import GenomicRanges
##' @import  IRanges
##' @useDynLib geneRxCluster
##' @include dot-gRxCluster.R
gRxCluster <-
    function(object, starts, group, kvals, nperm=0L, pruneFun=prune.loglik, ...,
             cutpt.filter.expr,
             cutpt.tail.expr,  tmp.env, sample.id, sample.tab)
    {
        
        mc <- match.call()

        if (!is.name(mc$cutpt.tail.expr))
            cutpt.tail.expr <-
                if (missing(cutpt.tail.expr)) 
                    quote(critVal.target(k,n,target=5,posdiff=x))
                else
                    mc$cutpt.tail.expr
        
        if (!is.name(mc$cutpt.filter.expr))
            cutpt.filter.expr <-
                if (missing(cutpt.filter.expr)) 
                    quote(as.double(apply(x,2,median,na.rm=TRUE)))
                else
                    mc$cutpt.filter.expr
        
        if (missing(tmp.env)) tmp.env <- new.env()
        
        ## basic sanity checks
        
        if (length(object)!=length(starts) || length(starts) != length(group))
            stop("object, starts, and group must have same lengths")
        
        
        ## object order checked in .gRxCluster(...)
        
        object <- as(object,"Rle")
        chr.lens <- runLength(object)
        chr.starts <- start(object) 
        names(chr.starts) <- as.character(runValue(object))
        chr.ends <- end(object) 
        
        res <- 
            .gRxCluster(chr.starts, chr.ends, starts, group, kvals, nperm,
                        cutpt.filter.expr, cutpt.tail.expr, pruneFun, tmp.env,
                        sample.id, sample.tab)
        
        if (is.null(mc$cutpt.tail.expr)) mc$cutpt.tail.expr <-
            cutpt.tail.expr
        if (is.null(mc$cutpt.filter.expr)) mc$cutpt.filter.expr <-
            cutpt.filter.expr
        
        metadata(res)$call <- mc
        
        res
    }
