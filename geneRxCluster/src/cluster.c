#include <R.h>
#include <Rdefines.h>
#include <R_ext/Error.h>


#include <Rmath.h>
#include <R_ext/Utils.h>



/* marks for window counts */

#define NOTAIL 0
#define HITAIL 1
#define LOWTAIL -1
#define HIOUT 2
#define LOWOUT -2
#define CONFLICT 9

void kval_counts( SEXP kcounts, SEXP sdiff, 
                  /* all the usual vars are added */
                  int *grp_ptr, int *kvals_ptr,
                  int strt_n,  int kvals_n)
  {
    int *sdiff_ptr = INTEGER( sdiff );
    int *kcounts_ptr = INTEGER( kcounts );
    
    for(int i = 0; i<kvals_n; i++){
      int k = kvals_ptr[i];
      int count=0;
      for(int j=0; j<k-1; j++){
        count+= grp_ptr[ j ];
        kcounts_ptr[ j + i*strt_n ]=NA_INTEGER;
      };
      for(int j=k-1; j<strt_n; j++){
        count+= grp_ptr[ j ];
        kcounts_ptr[ j + i*strt_n ] = 
          (sdiff_ptr[ j + i*strt_n] == NA_INTEGER) ? NA_INTEGER : count;
        count-= grp_ptr[j-k+1];
      };
    };
  }
 

void cutptApply(
                SEXP cutptFunRes, SEXP kcounts, SEXP ct, 
                int *kvals_ptr, int strt_n, int kvals_n)
{

  
  for (int i = 0; i < kvals_n; i++){
    /* Rprintf("k=%d\n",kvals_ptr[i]); */
    int k_low = REAL(AS_NUMERIC(VECTOR_ELT(cutptFunRes,i)))[0];
    int k_hi = REAL(AS_NUMERIC(VECTOR_ELT(cutptFunRes,i)))[1];
    for (int j = 0; j < strt_n; j++){
      int res; 
      int *kcounts_ptr = INTEGER(kcounts)+i*strt_n;
      int *ct_ptr = INTEGER(ct)+i*strt_n;
      if (kcounts_ptr[j]==NA_INTEGER)
        res = NOTAIL;
      else {
        if (kcounts_ptr[j] < k_low) 
          res = LOWTAIL;
        else 
          if (kcounts_ptr[j] >= k_hi)
            res = HITAIL;
          else
            res = NOTAIL;
      };
      ct_ptr[j] = res;
    };
  };
}

static inline void cutptClean_new_elt(int *pvj, int *ctj, int *pr_low, 
                                      int *pr_hi, int *ct_low, int *ct_hi, 
                                      int *ct_conflict, int kv){
  /* here pvj, ctj reference new elt */
  /* entry prior value */
  if (*pvj!=NOTAIL)
    if (*pvj==LOWTAIL)
      (*pr_low)++;
    else if (*pvj==HITAIL)
      (*pr_hi)++;
    else if (*pvj==CONFLICT){
      (*pr_low)++;
      (*pr_hi)++;
    }

    /* entry window */
  /* ct_confict is incremented sites in entering window */
  if (*ctj!=NOTAIL){
    if (*ctj==HITAIL){ 
      (*ct_hi)++;
      if (*ct_low){
        *ct_conflict=kv;
      }
      else if (*pr_low) {
        *ctj=HIOUT;
      }
    }
    else {
      (*ct_low)++;
      if (*ct_hi){
        *ct_conflict=kv;
      }
      else if (*pr_hi){ 
        *ctj=LOWOUT;
      }
    }
  }
}  

static inline void cutptClean_old_elt(int *pvj, int *ctj, int *pr_low, 
                                      int *pr_hi, int *ct_low, int *ct_hi, 
                                      int *ct_conflict, int kv){
  /* here pvj, ctj reference old elt */

  if (*pvj!=NOTAIL){ 
    if (*pvj==CONFLICT){
      (*pr_hi)--;
      (*pr_low)--;
    }
    else if (*pvj==HITAIL)
      (*pr_hi)--;
    else if (*pvj==LOWTAIL)
      (*pr_low)--;
  }
  

  int cv = *ctj;

  if (*ct_hi && *ct_low){ /* intra-k conflict - possibly alter prior */
    *ctj=NOTAIL;
    if (*pvj==NOTAIL) *pvj=CONFLICT;
    if (cv==HITAIL||cv==HIOUT){
      (*ct_hi)--;
      (*ct_conflict)--;
    }
    else if (cv==LOWTAIL||cv==LOWOUT){
      (*ct_low)--;
      (*ct_conflict)--;
    }
  }
  else if (cv!=NOTAIL){ 
    if (*ct_conflict){ /* intra-k legacy - let prior be */
      (*ct_conflict)--;
      *ctj=NOTAIL;
      if (cv==HITAIL||cv==HIOUT)
        (*ct_hi)--;
      else
        (*ct_low)--;
    }
    else if (cv==HIOUT||cv==LOWOUT){ /* prior already marked */
      *ctj=NOTAIL;
      if (cv==HIOUT)
        (*ct_hi)--;
      else
        (*ct_low)--;
    } 
    else { /* unconflicted window - copy to prior and backfill */ 
      *pvj=cv;
      for (int j=1;j<kv && (pvj[-j]==NOTAIL||pvj[-j]==cv);j++){
        pvj[-j]=cv;
        ctj[-j]=cv;
      }
      if (cv==HITAIL)
        (*ct_hi)--;
      else
        (*ct_low)--;
    }
  } else if (*ct_conflict)
    /* intra-k legacy - let prior and site be */
    (*ct_conflict)--;
}

void cutptClean(
                SEXP cutptFunRes, SEXP kcounts, SEXP ct, 
                int *kvals_ptr, int strt_n, int kvals_n, 
                int *prior_value)
{

  /* At this point ct has the tail flag for each window.  */
  /* Now it will be converted to the flag for each site.  */
  /* Look up tail for each window for one value of k - then classify sites: */
  /* 1. site is in window(s) in lower tail */
  /* 2. site is in window(s) in upper tail */
  /* 3. site is never in a tail or is claimed by at least one window in */
  /*    each tail */
  
  for (int i = 0; i<strt_n; i++) prior_value[i]=0;
  
  /* Accumulate cnt_hi, prev_hi, cnt_low, prev_low on entry; decrement
     on exit. conflict=(cnt_hi+prev_hi)&&(cnt_low+prev_low)
     if conflict, neuter left elt and no backfill else backfill.
     Update prior_value with backfilling. 
     (rethink if no checking is needed and if intra-k corrections are
     needed first)
  */

  for (int i = 0; i<kvals_n; i++){
    /* Rprintf("k=%d",kvals_ptr[i]); */
    
    /* initialize for this k */
    int cnt_hi=0, prev_hi=0, cnt_low=0, prev_low=0, 
      cnt_conflict=0;

    int *ct_ptr = INTEGER(ct)+i*strt_n;
    int kv_offset = kvals_ptr[i]-1;

    for (int j = 0; j<kv_offset; j++)
      cutptClean_new_elt( prior_value+j, ct_ptr+j, &prev_low, 
                          &prev_hi, &cnt_low,  &cnt_hi,  
                          &cnt_conflict, kvals_ptr[i]);
    
    
    for (int j = kv_offset; j<strt_n; j++){
      cutptClean_new_elt( prior_value+j, ct_ptr+j, &prev_low, 
                          &prev_hi, &cnt_low,  &cnt_hi,  
                          &cnt_conflict, kvals_ptr[i]);
      cutptClean_old_elt( prior_value+j-kv_offset, ct_ptr+j-kv_offset, 
                          &prev_low, &prev_hi, &cnt_low,  
                          &cnt_hi,  &cnt_conflict, kvals_ptr[i]);
    };

    /* clear out final kvals_ptr[i] windows */
    for (int j = strt_n-kv_offset; j<strt_n; j++)
      cutptClean_old_elt( prior_value+j, ct_ptr+j, 
                          &prev_low, &prev_hi, &cnt_low,  
                          &cnt_hi,  &cnt_conflict, kvals_ptr[i]);

  };
}

   
 int depthFun(SEXP ct, int *depth, int *cid, int *chromoSts_ptr,
              int *chromoEnds_ptr, int strt_n, int kvals_n, int chromoSts_n)

{
   
   int cid_n=0;
   
   for (int i = 0; i<strt_n; i++) depth[i]=0;
   
   for (int i = 0; i<kvals_n; i++){
     int *ct_ptr = INTEGER(ct)+i*strt_n;
     for (int j = 0; j<strt_n; j++) depth[j] += ct_ptr[j];
   };
   
   /* cluster ID 
    * use sequential numbers for clusters
    * use zero for unassigned sites
    */
   
   for (int m=0; m<chromoSts_n; m++){  
     int lastdepth = 0;
   
     for (int i=chromoSts_ptr[m]-1; i<chromoEnds_ptr[m]; i++){
       if (depth[i] == 0){
         cid[i]=0;
         lastdepth = 0;}
       else 
         {
           if (lastdepth * depth[i] <= 0 )
             cid_n++;
           lastdepth = depth[i];
           cid[i]=cid_n;
         }
     }
   }
   return cid_n;
}

/* get best nominal alpha 
 * for each position, 
 * then for each cluster 
 */

void best_val( 
              SEXP sitewise_best,
              SEXP cluster_best,
              SEXP cutptFunRes,
              SEXP kcounts,
              int *kvals_ptr,
              int *cid,
              int cid_n,
              int strt_n,
              int kvals_n)
{
  
    double *swb = REAL(sitewise_best);
    for (int i = 0; i<strt_n; i++) swb[i]=R_PosInf;
  
    for (int i = 0; i<kvals_n; i++){
  
      double *ctpt_fdr = 
        REAL( getAttrib( VECTOR_ELT( cutptFunRes,i) , install("fdr")));
  
      int *kcounts_ptr = INTEGER(kcounts) + i*strt_n ;
    
      /* recall that ctpt_fdr has 
       * kvals_ptr[i]+1 rows  
       * whose indices are 0:kvals_n[i]
       */
      
      for (int j = 0; j<strt_n; j++){
        int kcp = kcounts_ptr[j]; 
        if ( kcp != NA_INTEGER ){
          double fdr_low = ctpt_fdr[ kcp ];
          double fdr_up =  ctpt_fdr[ kcp + kvals_ptr[i] + 1 ];
          if (fdr_low<swb[j]) swb[j] = fdr_low;
          if (fdr_up<swb[j]) swb[j] = fdr_up;
        };
      };
    };
    
    double *cbest = REAL(cluster_best);
   
    for (int i=0; i <= cid_n; i++) cbest[i]=R_PosInf;
  
    for (int i = 0; i<strt_n; i++)
      if (swb[i] < cbest[ cid[i] ]) cbest[ cid[i] ] = swb[i];
  }

// count the number of cluster/depth combos

int cd_count( SEXP depth_sexp, SEXP cluster_id, int strt_n)
{
  int cdn=0,lastd=0,lastc=0;
  for (int i = 0; i<strt_n; i++){
    int cid = INTEGER(cluster_id)[i];
    if (( cid !=0 ) &&
        ((lastd != INTEGER(depth_sexp)[i] ) ||   
         (lastc != cid))){
      cdn++;
      lastd = INTEGER(depth_sexp)[i];
      lastc = cid;
    };
  };

  return cdn;
}
    
/* get summaries of clusters: 
 * start index, end index , depth, group 1, group 2
 */

void clustsum(SEXP depth_sexp, SEXP cluster_id, SEXP grp,
              SEXP summary_matrix, int cd_combo_n, int strt_n)
{
  int *depth=INTEGER(depth_sexp), *cid = INTEGER(cluster_id), *grp_ptr=INTEGER(grp);
  // columns of results matrix:
  int *start_index = INTEGER(summary_matrix);
  int *end_index =   INTEGER(summary_matrix) +   cd_combo_n;
  int *cl_depth =    INTEGER(summary_matrix) + 2*cd_combo_n;
  int *gr0 =         INTEGER(summary_matrix) + 3*cd_combo_n;
  int *gr1 =         INTEGER(summary_matrix) + 4*cd_combo_n;
  
  int i_last = 0, cd_combo=-1;
  int lastdepth = 0, last_cid=0, ctab[2]={0,0}; 
  for (int i = 0; i<strt_n; i++){
    if (cid[i] != 0) {
      if ( (lastdepth != depth[i]) ||
           (last_cid != cid[i])){ 
        /* finalize last cluster 
           initialize for current cluster
        */
        if (cd_combo > -1){
          end_index[ cd_combo ] = i_last+1; // 1 based index
          gr0[ cd_combo ] = ctab[0];
          gr1[ cd_combo ] = ctab[1];
          cl_depth[ cd_combo ] = depth[i_last];
        };
        cd_combo++;
        start_index[ cd_combo ] = i+1;
        last_cid = cid[i];
        lastdepth = depth[i];
        ctab[0] = ctab[1] = 0;
      };
      i_last = i; 
      ctab[ grp_ptr[i] ]++;
    };
  };
  end_index[ cd_combo ] = i_last+1; // 1 based index
  gr0[ cd_combo ] = ctab[0];
  gr1[ cd_combo ] = ctab[1];
  cl_depth[ cd_combo ] = lastdepth;
  
}

SEXP gRxC_cluster ( SEXP chromoSts, SEXP chromoEnds, SEXP strt, SEXP grp, SEXP kvals, SEXP cutptExprs, SEXP cutptFunExprs, SEXP tmpEnv, SEXP nperm, SEXP sample_id, SEXP sample_tab ) {

   /*
   * Based on this signature: 
   * signature(chromoSts="integer",
   *     chromoEnds="integer",
   *      strt="integer",
   *      grp="integer",
   *      kvals="integer",
   *      cutptExprs="call",
   *      cutptFunExprs="call",
   *      tmpEnv="environment",
   *      nperm="integer",
   *      sample_id="integer",
   *      sample_tab="integer")  
   *   Note that chromoSts, chromoEnds, and strt originate at 1 not 0
   * 
   */
 
 /* set up ptrs */ 
 
 
 int *chromoSts_ptr = INTEGER(chromoSts);
 int *chromoEnds_ptr = INTEGER(chromoEnds);
 int *strt_ptr = INTEGER(strt);
 int *grp_ptr = INTEGER(grp);
 int *kvals_ptr = INTEGER(kvals);
 int *sample_tab_ptr = INTEGER(sample_tab);
 int *sample_id_ptr = INTEGER(sample_id);
 
 /* nperms has length 1. Get the value: 
 */
 
 int perm_n = INTEGER( nperm )[ 0 ];
 
 /* get lengths of objects */
 
 int chromoSts_n = length( chromoSts );
 int strt_n = length( strt );
 int grp_n = length( grp );
 int kvals_n = length( kvals );
 int len_sample = length( sample_tab );
 
 
 /* SEXPs are setup here */
 
 /* SEXPs that can be PROTECTed at onset */
 /* ------------------------------------ */
 
 int nprotect=0;
 
 SEXP grp_orig;
 PROTECT( grp_orig = duplicate( grp ));nprotect++;
 
 SEXP grp_urand;
 PROTECT( grp_urand = allocVector(REALSXP, grp_n));nprotect++;
 
 SEXP pr, n_table;
 PROTECT(pr = allocVector(REALSXP,1));nprotect++;
 PROTECT(n_table = allocVector(REALSXP,2));nprotect++;
 
 SEXP sdiff; 
 PROTECT( sdiff = allocMatrix( INTSXP, strt_n, kvals_n ));nprotect++;
 int *sdiff_ptr = INTEGER( sdiff );
 
 SEXP kcounts;
 PROTECT( kcounts = allocMatrix(INTSXP, strt_n, kvals_n));nprotect++;
 
 SEXP ct;
 PROTECT( ct = allocMatrix(INTSXP, strt_n, kvals_n));nprotect++;
 
 // share the storge here:
 SEXP prior_value_sexp;
 PROTECT( prior_value_sexp = allocVector(INTSXP , strt_n ));nprotect++;
 int *prior_value = INTEGER(prior_value_sexp);
 SEXP depth_sexp = prior_value_sexp;
 int *depth = INTEGER(depth_sexp);
 
 
 SEXP cluster_id, cluster_best_list;
 PROTECT( cluster_id = allocVector(INTSXP,strt_n));nprotect++;
 PROTECT( cluster_best_list = allocVector(VECSXP,perm_n+1));nprotect++;
 int *cid = INTEGER(cluster_id);
 
 SEXP sitewise_best;
 PROTECT( sitewise_best = allocVector(REALSXP, strt_n));nprotect++;
 
 SEXP summary_matrix_list;
 PROTECT( summary_matrix_list = allocVector(VECSXP,perm_n+1));nprotect++;
 
 
  /* SEXPs that must be PROTECTed later */
 /* ----------------------------------- */
 
 
 SEXP final;
 
 // PROTECT(final = allocVector(VECSXP, 8 ));nprotect++;
 
 
 SEXP cutptSdiff;
 
 //  PROTECT(cutptSdiff = eval( cutptExprs, tmpEnv ));nprotect++;
 
 SEXP cutptFunRes;
   
 //  PROTECT(cutptFunRes = eval( cutptFunExprs, tmpEnv ));nprotect++;
 
 SEXP cluster_best;
 
 //  PROTECT( cluster_best = allocVector(REALSXP, 1+cid_n));
 
 SEXP summary_matrix;
 
 //  PROTECT(summary_matrix = allocMatrix(INTSXP,cd_combo_n,5));
 // Rprintf("enter rolling\n");
 /* sdiff is the start difference */
 
  for (int i = 0; i<kvals_n; i++){
     int k = kvals_ptr[i];
     int l = i*strt_n+k-1;
     /* rolling difference of k starts */
     for (int j=k-1; j < strt_n; j++){
       sdiff_ptr[l] = strt_ptr[j]-strt_ptr[j-k+1];
       l++;
     };
     /* omit first k-1 values on each chromo by setting NAs */
     for (int m=0; m<chromoSts_n; m++){
       int max_end = imin2( chromoEnds_ptr[m], chromoSts_ptr[m] + k - 2);
       for (int j=chromoSts_ptr[m]-1; j < max_end; j++) 
         sdiff_ptr[j+strt_n*i] = NA_INTEGER;
     }
   };
 // Rprintf("enter block-broken\n");
 
 for (int i = 0; i<kvals_n; i++){
   int k = kvals_ptr[i];
   /* omit k length groups that only include some position ties */
   for (int m=0; m<chromoSts_n; m++){
     int max_end = chromoEnds_ptr[m];
     for (int j=chromoSts_ptr[m]-1+k; j < max_end; j++){
       if (strt_ptr[j-1]==strt_ptr[j]) {
         sdiff_ptr[j-1+strt_n*i] = NA_INTEGER;
       };
       if (strt_ptr[j-k]==strt_ptr[j-k+1]) 
         sdiff_ptr[j+strt_n*i] = NA_INTEGER;
     }
   }
  };
 
 // Rprintf("enter sdiff\n");
 /* cutptExprs - is the expression to be used in setting up the
    cutpoint for narrow intervals.
 
    tmpEnv - is an environment used for executing R commands
 
  */
 
 if(!isEnvironment(tmpEnv))
   error("tmpEnv should be an environment");
 defineVar( install("x"), sdiff, tmpEnv );
 
 PROTECT(cutptSdiff = eval( cutptExprs, tmpEnv ));nprotect++;
 int cutptSdiff_n = length( cutptSdiff );
 
 if (cutptSdiff_n != kvals_n)
   error( "cutpt.filter.expr returned the wrong length");
 
 if (!isReal(cutptSdiff))
   error("cutpt.filter.expr result must yield double");
 
 int l=0;
 for (int i = 0; i<kvals_n; i++){
   double cutat = REAL(cutptSdiff)[i];
   for (int j=0; j<strt_n; j++){
     if ( sdiff_ptr[l] != NA_INTEGER && (double) sdiff_ptr[l] > cutat ) 
       sdiff_ptr[l] = NA_INTEGER;
     l++;
   };
  };
 
 // Rprintf("enter permute\n");
   
 for (int iperm = perm_n; iperm >= 0; iperm--){
    R_CheckUserInterrupt();
    if (iperm==0) // last time thru, use grp_orig 
      copyVector( grp, grp_orig );
    else 
      { // check if sample ids are tabled
        if (len_sample>0) {
          GetRNGstate();
          for (int i=0; i<len_sample; i++) REAL(grp_urand)[i] = unif_rand();
          PutRNGstate();
          
          rsort_with_index( REAL(grp_urand), sample_tab_ptr, len_sample);
          
          for (int i=0; i<grp_n; i++) grp_ptr[i] = sample_tab_ptr[ sample_id_ptr[ i ]];
 
   }
        else
          {
            /* permute all sites in grp */
            
            GetRNGstate();
            for (int i=0; i<grp_n; i++) REAL(grp_urand)[i] = unif_rand();
            PutRNGstate();
        
            rsort_with_index( REAL(grp_urand), grp_ptr, grp_n);
          }
      }
  
 /* tabulate the group */
 
 int *int_ptr = grp_ptr;
 double ntab[2] = {0.0,0.0};
 for (int i = 0; i < grp_n; i++){
   ntab[ *int_ptr ]++;
   int_ptr++;
  };
 
 REAL(pr)[0] = ntab[1]/(ntab[0]+ntab[1]);
 REAL(n_table)[0] = ntab[0];REAL(n_table)[1] = ntab[1];
 
 // Rprintf("enter counts\n");
 kval_counts(  kcounts,  sdiff, 
               /* all the usual vars are added */
               grp_ptr, kvals_ptr, 
               strt_n, kvals_n);
 
 // Rprintf("enter cutpt\n");
 /* cutptFunExprs - is the expression to be used in setting up the
    cutpoint for kcounts.
    
    tmpEnv - is an environment used for executing R commands
    
 */
 
 if(!isEnvironment(tmpEnv))
   error("tmpEnv should be an environment");
 defineVar( install("k"), kvals, tmpEnv );
 defineVar( install("n"), n_table, tmpEnv );
 
 PROTECT(cutptFunRes = eval( cutptFunExprs, tmpEnv ));nprotect++;
 
 int cutptFunRes_n = length( cutptFunRes );
 
 if (cutptFunRes_n != kvals_n)
   error( "cutptFunExprs returned the wrong length");
 
 /* use the cutpts to classify  kcounts */
 /* resolve conflicts in favor of lesser k values */
 
 cutptApply(  cutptFunRes, kcounts, ct, kvals_ptr,
              strt_n, kvals_n);
 cutptClean(  cutptFunRes, kcounts, ct, kvals_ptr,
              strt_n, kvals_n, prior_value);
 
 // Rprintf("enter depth\n");
 int cid_n = depthFun( ct, depth, cid, chromoSts_ptr, chromoEnds_ptr, strt_n, kvals_n, chromoSts_n);
 // Rprintf("enter miinval\n");
 /* get best nominal alpha for each position, then for each
  * cluster */
 
 
 /* need to protect/unprotect each time */
 
 PROTECT( cluster_best = allocVector(REALSXP, 1+cid_n)); // nprotect++;
     
 best_val( sitewise_best, cluster_best, cutptFunRes, kcounts, 
           kvals_ptr, cid, cid_n, strt_n, kvals_n);
 
 SET_VECTOR_ELT(cluster_best_list, iperm, duplicate(cluster_best));
 
 UNPROTECT(1);
 // Rprintf("enter summary\n");
 int cd_combo_n = cd_count( depth_sexp, cluster_id, strt_n); 
 
 PROTECT(summary_matrix = allocMatrix(INTSXP,cd_combo_n,5)); //nprotect++;
 
 if (cd_combo_n>0)
   clustsum( depth_sexp, cluster_id, grp, summary_matrix, cd_combo_n, strt_n);
 
 SET_VECTOR_ELT( summary_matrix_list, iperm, duplicate(summary_matrix));
 UNPROTECT(1);
 
   
   /* end for (iperm = ... */
   }
 PROTECT(final = allocVector(VECSXP, 10L ));nprotect++;
 SET_VECTOR_ELT(final,0,kcounts);
 SET_VECTOR_ELT(final,1,ct);
 SET_VECTOR_ELT(final,2,cutptFunRes);
 SET_VECTOR_ELT(final,3,depth_sexp);
 SET_VECTOR_ELT(final,4,cluster_id);
 SET_VECTOR_ELT(final,5,sitewise_best);
 SET_VECTOR_ELT(final,6,cluster_best_list);
 SET_VECTOR_ELT(final,7,summary_matrix_list);
 SET_VECTOR_ELT(final,8,sdiff);
 SET_VECTOR_ELT(final,9,cutptSdiff);
 UNPROTECT(nprotect);
   
 return final;
 
  warning("your C program does not return anything!");
  return R_NilValue;
}
