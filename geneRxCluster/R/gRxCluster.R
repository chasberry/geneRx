##' gRxCluster
##' 
##' cluster integration sites
##' @param object chromosome names or other grouping of starts 
##' @param starts ordered chromosome position or ordered integer
##' vector 
##' @param group logical vector separating two groups
##' @param kvals integer vector of window widths
##' @param nperm number of permutations for FDR calculation
##' @param pruneFun a function like \code{prune.loglik}. 
##' @param ... other args
##' 
##' @param cutpt.filter.expr (optional) a R call (or variable naming a
##' call) with (optional) var x
##' ##' (window widths in base pairs) to filter windows must eval to mode
##' ##' "double". If not specified,
##' ##' \code{as.double(apply(x,2,median,na.rm=TRUE))} is used. If is
##' an atomic vector of length one is supplied it is expanded to the
##' proper length and coerced to double.
##'
##' 
##' @param cutpt.tail.expr R call (or variable naming a call) with
##' (optional) vars: k,n, and x (as above). Returns list. k is a
##' vector of the number of sites in a collection of windows, and n is
##' a vector of counts or proportions for the two classes of
##' insertion. If not supplied,
##' \code{critVal.target(k,n,target=5,posdiff=x)} is used.
##'
##' @param tmp.env (optional) environment to use in evaluation of
##' cutpt.* expressions. This is usually needed for
##' \code{\link{critVal.power}}, which is first calculated and placed
##' in the environment, and the supplied object is used in the
##' expression for \code{cutpt.filter.expr}.
##' 
##' @return a GRanges object with a special metadata slot, see
##' \code{\link{gRxCluster-object}}
##' 
##' @author Charles Berry
##' @example inst/ex-gRxCluster.R
##' @export
##' @import GenomicRanges
##' @import  IRanges
##' @useDynLib geneRxCluster
##' @include dot-gRxCluster.R
gRxCluster <-
    function(object, starts, group, kvals, nperm=0L, pruneFun=prune.loglik, ...,
             cutpt.filter.expr,
             cutpt.tail.expr,  tmp.env)
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
                        cutpt.filter.expr, cutpt.tail.expr, pruneFun, tmp.env)
        
        if (is.null(mc$cutpt.tail.expr)) mc$cutpt.tail.expr <-
            cutpt.tail.expr
        if (is.null(mc$cutpt.filter.expr)) mc$cutpt.filter.expr <-
            cutpt.filter.expr
        
        metadata(res)$call <- mc
        
        res
    }
