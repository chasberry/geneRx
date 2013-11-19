### --- Test setup ---
 
if(FALSE) {
  ## Not really needed, but can be handy when writing tests
  library("RUnit")
  library("geneRxCluster")
}
 
### --- Test functions ---

callCl <- geneRxCluster:::grx_clust_Call

res <- callCl(chromoSts=c(1L,11L), chromoEnds=c(10L,20L), strt=as.integer(1:20),
              grp=as.integer(rep(0:1,each=10)), kvals=c(10L),
              cutptExprs=rep(Inf,1), 
              cutptFunExprs=
              list(structure(as.integer(c(3,8)),
                             fdr=cbind(low=(0:10)/10,up=(10:0)/10))) ,
              tmpEnv=new.env(), nperm=0L) 

names(res) <-
              c("kcounts", "ct", "cutptFunRes", "depth", "cluster_id",
                "sitewise_best", "cluster_best", "summary_matrix", "sdiff",
                "cutptSdiff")

test_simply_grx_clust_Call <- function(){
    checkTrue(all(res$kcounts[c(10,20),1]==c(0,10)))
    checkTrue(all(res$ct==rep(c(-1,1),each=10)))
    checkTrue(all(res$depth==rep(c(-1,1),each=10)))
    checkTrue(all(res$cluster_best[[1]] == c(Inf,0,0)))
    checkTrue(all(res$sdiff[c(10,20),1]==c(9,9)))
    checkTrue(sum(is.na(res$sdiff[,1]))==18)
}

test_crossover_gRxCluster <- function(){

    res <- gRxCluster(rep("a",100),1:100L,rep(c(FALSE,TRUE),each=50),
                      kvals=10L:30L, cutptExprs=rep(Inf,1))
    checkTrue(length(res)==2)
    checkEquals(res$value1,rev(res$value2))
}

test_prune_gRxCluster <- function(){

    res <- gRxCluster(rep("a",100),1:100L,
                      c(rep(TRUE,5),rep(c(FALSE,TRUE),each=45),rep(FALSE,5)),
                      kvals=10L:30L, cutptExprs=rep(Inf,1))
    checkTrue(length(res)==2)
    checkEquals(width(res),rep(45,2))
}

test_bad_args_gRxCluster <- function(){
    checkException(
        gRxCluster(rep("a",100),0:100L,
                   c(rep(TRUE,5),rep(c(FALSE,TRUE),each=45),rep(FALSE,5)),
                   kvals=10L:30L, cutptExprs=rep(Inf,1)),
        "arg lengths not equal")
}
