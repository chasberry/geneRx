##' @include grx_clust_Call.R
##' @include Rutils.R
.gRxCluster <-
    function(chr.starts, chr.ends, starts, group, kvals, nperm,
             cutpt.filter.expr = quote(as.double(apply(x,2,median,na.rm=TRUE))),
             cutpt.tail.expr = quote(cutpt.threshhold.binom(k,n,val=0.0005)),
             pruneFun=prune.loglik,
             tmp.env = new.env())
    {
        ## Purpose: Cluster Regions of Relatively High Frequency  
        ## ----------------------------------------------------------------------
        ## Arguments:
        ##   arg               | description                                   
        ##  -------------------+-----------------------------------------------
        ##   chr.starts        | 1-based integer - where chromo starts         
        ##   chr.ends          | 1-based integer - where chromo ends           
        ##   starts            | integer site of attack
        ##   group             | integer 0:1 - same length as starts 
        ##   kvals             | integer window widths                         
        ##   nperm             | integer - number of permutations for FDR      
        ##   cutpt.filter.expr | R expr with (optional) var x to filter windows
        ##                     |   must eval to mode "double"
        ##   cutpt.tail.expr   | R expr with (optional) vars: ccnt,k,n
        ##                     |   list - see help for cutpt.threshhold.val 
        ##   pruneFun          | usually prune.loglik
        ##   tmp.env           | environment in which to eval() cutpt.* expr's 
        ##  -------------------------------------------------------------------
        ##
        ## Author: Charles Berry, Date: 27 Apr 2013, 09:37
    
        ## check chr.starts, chr.ends, starts, kvals, nperm, and (maybe)
        ## cutpt.* expr's


        was.quoted <- function(x){
            is.name(x) ||
            isTRUE(try(is.language(x),TRUE)) &&
                isTRUE(try(is.call(x),TRUE))
        }
        
        nperm <- as.integer(nperm)
        stopifnot( length(chr.starts)==length(chr.ends) )
 
        stopifnot( tail(chr.ends,1)==length(starts) )
        stopifnot( is.integer(kvals) && all(kvals>1) )
        stopifnot( nperm >= 0 )
        if (is.atomic(cutpt.filter.expr))
            {
                cfe <- as.double(cutpt.filter.expr)
                if (length(cfe)!=1 || any(is.na(cfe)))
                    stop("bad cutpt.filter.expr arg")
                cutpt.filter.expr <-
                    bquote(rep(.(atom),length=ncol(x)),
                           list(atom=cfe))
            }
        else
            {
                if (!was.quoted(cutpt.filter.expr))
                    stop("cutpt.filter.expr malformed ??")
            }
        if (!was.quoted(cutpt.tail.expr))
            stop("cutpt.tail.expr not in quote() ??")
        stopifnot( all(as.integer(group) %in% 0:1 ) )
        stopifnot( all( head(kvals,-1) < tail(kvals,-1) ) )
        pruneFun <- match.fun(pruneFun)
        ## coerce to integer
        sn <- names(chr.starts) # hold names
        chr.starts <-  as.integer(chr.starts)
        chr.ends   <-  as.integer(chr.ends  )
        starts     <-  as.integer(starts    )
        kvals <- as.integer(kvals)
        
        ## order of starts must be non-decreasing on each chromo
        starts.dont.decrease <-
            function(i) all( diff( starts[chr.starts[i]:chr.ends[i] ] ) >=0 )
        
        stopifnot( all(sapply(seq_along(chr.starts), starts.dont.decrease )))
        
        ## C code checks is.double( eval( cutpt.filter.expr, tmp.env ) )
        
        res <- grx_clust_Call(
            chr.starts, chr.ends, starts, group, kvals,
            cutpt.filter.expr,
            cutpt.tail.expr,
            tmp.env, nperm)
        
        names(res) <-
            c("kcounts", "ct", "cutptFunRes", "depth", "cluster_id",
              "sitewise_best", "cluster_best", "summary_matrix", "sdiff",
              "cutptSdiff")
        ## prune
        smat <- res[['summary_matrix']][[1]]


        gr <-
            GRanges(seqnames =
                    if (length(sn)) {
                        factor(sn[findInterval( smat[,1] , chr.starts)],sn)
                    } else {
                        as.character(findInterval( smat[,1] , chr.starts))
                    },
                    IRanges(start=starts[smat[,1]],end=starts[smat[,2]]),
                    depth= smat[,3],
                    value1= smat[,4],
                    value2= smat[,5],
                    clump.id= res[["cluster_id"]][smat[,1]])
        
        pruned <-
            if (length(gr))
                pruneFun(gr, prop.table(table(group))[1])
            else
                gr
        
        ## the first cluster_best element is always not in a cluster,
        ## but can be smallish near conflicts
        pruned$target.min <- res[["cluster_best"]][[1]][-1] 
        
        metadata(pruned)[["criticalValues"]] <- res[["cutptFunRes"]]
        metadata(pruned)[[ "kvals" ]] <- kvals
        metadata(pruned)[[ "perm_cluster_best" ]] <-
            lapply( res[["cluster_best"]][-1], '[', -1)
        
        metadata(pruned)[['summary_matrix']] <- res[['summary_matrix']][[1]]
        
        pruned
    }
