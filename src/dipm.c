
//  ------------------------------------------------------------
//
//    Depth Importance in Precision Medicine (DIPM) Method
//
//      This code creates a tree using the depth 
//      variable importance measure developed by Zhang 
//      for the precision medicine setting.
//
//      Data must have a continuous or right censored
//      survival outcome variable Y and two or more
//      treatment groups.
//
//  ------------------------------------------------------------
//
//    Overview of Functions
//
//      get_bestvar - This function identifies the 
//      "best" variable to split a node.
//
//      get_split - This function identifies the 
//      "best" split of a node.
//
//      maketree - This function creates a tree. Tree
//      information includes: the split variable, the
//      split variable's type, the parent of a node,
//      the depth of a node, and the predicted best 
//      treatment group.
//
//  ------------------------------------------------------------



#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#ifdef _OPENMP
# include <omp.h>        // OpenMP for multiple threads
#endif
#include <R.h>          // for Rprintf
#include <Rmath.h>
#include <Rdefines.h>
#include <Rinternals.h>



const int MAXSIZE=10000;



struct
node

//
//  This struct object represents a single node
//  in a tree. 
//
//  "index" is the integer count of the node
//  starting from 1. "covar" is the integer index 
//  of the split variable starting from 1. "type" 
//  indicates the type of the split variable. 
//  "sign" indicates the sign of the current split.
//  "val" is the current split value. "parent" is 
//  the parent node of the current node. "depth" 
//  is the depth of the current node in an overall 
//  tree starting from 1. "lchild" and "rchild" 
//  are the integer indices of the node's left and 
//  right child nodes respectively. And "G" is the
//  G statistic for the current split (that is
//  calculated when applicable).
//
//  For descriptions about "n", "n0", "n1", "r0",
//  "r1", "p0", "p1" and "predict", see the "cand" 
//  struct below.
//

{
    int index;
    int covar;
    int type;
    int sign;
    double val;
    int parent;
    int depth;
    int lchild;
    int rchild;
    int n;
    int n0;
    int n1;
    double G;
    double r0;
    double r1;
    double p0;
    double p1;
    int predict;
};



struct 
cand

//
//  This struct object holds information about
//  a single candidate node. 
//
//  "n" is the total number of subjects in the 
//  candidate node. "n0" and "n1" are the number 
//  of subjects in treatment groups 0 and 1 
//  respectively. "r0" and "r1" are the response
//  rates to treatments 0 and 1 respectively. 
//  "p0" and "p1" are the proportion of subjects
//  in treatment group 0 or 1 that have Y values
//  greater than "r0" or "r1" respectively.
//  "p" is the predicted best treatment group of 
//  the node. And "val" is a double useful in 
//  calculating a split score.
//

{
    int n;
    int n0;
    int n1;
    double r0;
    double r1;
    double p0;
    double p1;
    int p;
    double val;
};



int
get_rand_int(int n)

//
//  This function accepts as an argument integer
//  "n" and returns a randomly generated integer
//  from 0:(n-1).
//

{
    double u=runif(0,n);
    int rand_int=floor(u);

    return(rand_int);
}



void
sample_without_replace(int n,
                       int nval,
                       int sample[])

//
//  This function returns a random sample of size
//  "nval" that draws integers from 0:(n-1) without
//  replacement. The random sample is saved in 
//  array of integers "sample".
//
//  If "n" equals "nval", then "sample" is a random
//  permutation without replacement of the integers
//  0:(n-1).
//  
//  Note that "nval" should NOT be larger than "n".
//

{
    if ( nval == 1 ) {

        sample[0]=get_rand_int(n);
        return;
    } 

//   initialize sample0 to be the integers 0:(n-1)
    int i;
    int sample0[n];

    for ( i=0; i<n; i++ ) {
        sample0[i]=i;
    }

//   randomly permute sample0
    int j;
    int temp;

    for ( i=0; i<(n-1); i++ ) {

        j=get_rand_int(n-i)+i;  // draw int from i:(n-1)

        temp=sample0[i];
        sample0[i]=sample0[j];
        sample0[j]=temp;
    }

/*
//   check
    Rprintf("\nRandomly permuted sample from 0:%d:\n",n-1);
    for ( i=0; i<n; i++ ) Rprintf("%d  ",sample0[i]);
*/

//   the final sample contains nval values
    for ( i=0; i<nval; i++ ) {
        sample[i]=sample0[i];
    }
}



int
get_number_of_nodes(int n,
                    struct node *tree)

//
//  This function accepts as an argument an array
//  of "n" node structs "tree" and returns the 
//  total number of nodes in "tree" as an integer.
//

{
    int i;
    int max=tree[0].index-1;

    for ( i=0; i<n; i++ ) {

        if ( tree[i].index > max ) {
            max=tree[i].index;
        }
    }

    return(max);
}



int
compare_doubles(const void *a,
                const void *b)

//
//  This function is used to assist the "qsort" 
//  function.
//

{
    const double *da = (const double *) a;
    const double *db = (const double *) b;

    return (*da > *db) - (*da < *db);
}



int
get_nunique(int n,
            double *vector,
            double *unique)

//
//  This function accepts as an argument a vector
//  of "n" doubles called "vector" and returns the 
//  number of unique values as an integer. The 
//  actual unique values are saved in array of
//  doubles "unique".
//

{
    int i,j;
    int nunique=0;

    for ( i=0; i<n; i++ ) {

        for ( j=0; j<(i+1); j++ ) {

            if ( (vector[i] == vector[j]) && (i != j) ) {
                break;
            }

            if ( i == j ) {
                nunique++;
                unique[nunique-1]=vector[i];
            }
        }
    }

/*
//   check
    Rprintf("\n\nUnique values (nunique=%d): ",nunique);
    for ( i=0; i<nunique; i++ ) Rprintf("%f  ",unique[i]);
    Rprintf("\nVector:\n");
    for ( i=0; i<n; i++ ) Rprintf("%f  ",vector[i]);
*/

    return(nunique);
}



double
get_mean(double y[],
         int n)

//
//  This function accepts as an argument an array
//  of "n" doubles "y" and returns its arithmetic 
//  mean.
//

{
    if ( n == 1 ) return(y[0]);

    int i;
    double ybar=0;

    for ( i=0; i<n; i++ ) {
        ybar += y[i];
    }

    ybar /= n;

    return(ybar);
}



void
get_col(double **data,
        int var,
        int n,
        int nc,
        double *col)

//
//  This function accepts as arguments data set
//  of covariates "data" with "n" rows and "nc"
//  columns and the integer index of a variable 
//  "var" (starting from 1 instead of 0).
//
//  This function returns the column of data of
//  the variable denoted by "var" in array of
//  doubles "col".
//

{
    int i;
    for ( i=0; i<n; i++ ) {
        col[i]=data[i][var-1];
    }
}



void
findrows_ord(double *col,
             int n,    // length of "col"; total # rows in data
             int sign, // 0=no sign, 1=LE, 2=GT
             double val,
             int *ifrow)

//
//  This function accepts as arguments an array of 
//  data "col" of length "n" for an ORDINAL variable
//  and one of its values "val" and returns vector of 
//  indicators "ifrow", where "1" means a row of data
//  is either LE (less than or equal to) or GT 
//  (greater than) "val" depending on "sign".
//

{
    if ( sign == 1 ) {

        int i;
        for ( i=0; i<n; i++ ) {

            ifrow[i]=0;
            if ( col[i] <= val ) ifrow[i]=1;
        }

    } else if ( sign == 2 ) {

        int i;
        for ( i=0; i<n; i++ ) {

            ifrow[i]=0;
            if ( col[i] > val ) ifrow[i]=1;
        }
    }
}



void
findrows_bin(double *col,
             int n,  // length of "col"; total # rows in data
             double val,
             int *ifrow)

//
//  This function accepts as arguments an array of
//  data "col" of length "n" for a BINARY variable 
//  and one of its values "val" and returns a vector 
//  of indicators "ifrow", where "1" means a row of 
//  data equals "val".
//

{
    int i;
    for ( i=0; i<n; i++ ) {

        ifrow[i]=0;
        if( col[i] == val ) ifrow[i]=1;
    }
}



void
findrows_nom(double *col,
             int n,  // length of "col"; total # rows in data
             int ncat,
             double val,
             int *ifrow)

//
//  This function accepts as arguments an array of 
//  data "col" of length "n" for a NOMINAL variable
//  and a value "val" representing category membership 
//  within the nominal covariate and returns array of
//  indicators "ifrow", where "1" means a row of
//  data is in the list of categories specified by
//  "val".
//
//  "val" represents category membership by being
//  the decimal version of a number consisting of 0s 
//  and 1s. "val" is the conversion of a binary
//  integer (or bit) to a decimal. The 1s in the bit
//  designate categories that are included at a node.
//

{
    int i;
    for ( i=0; i<n; i++ ) {
        ifrow[i]=0;  // initialize
    }

    int j;
    double val2=val;
    double icat=ncat;

    for ( i=0; i<ncat; i++ ) {

        val2/=2;

//       a remainder of 0 means category (ncat-i) is not in this node
        if ( floor(val2) == val2 ) {
            icat-=1.0;
            continue;
        }

        for ( j=0; j<n; j++ ) {
            if ( col[j] == icat ) ifrow[j]=1;
        }

//       no more categories to consider
        if ( val2 == 0.5 ) break;

        icat-=1.0;
        val2=floor(val2);
    }
}



int
findrows_node(int current_node,
              int n,
              int nc,
              double **data,
              int *ncat,
              struct node *tree,
              int *ifcovar)

//
//  This function accepts as arguments the integer 
//  index of a node and tree information and returns 
//  a vector of indicators "ifcovar", where "1" means 
//  a particular row of data is in the current node.
//

{
    int i;
    for ( i=0; i<n; i++ ) {
        ifcovar[i]=1;  // initialize
    }

//   every subject is in the root node
    if ( current_node == 1 ) return(n);

//   get array of relevant nodes
    int j;
    int current_depth=tree[current_node-1].depth;
    int relevant_nodes[current_depth-1];
    relevant_nodes[0]=current_node;

    if ( current_depth > 2 ) {

        for ( i=1; i<(current_depth-1); i++ ) {        	  

            j=relevant_nodes[i-1];
            relevant_nodes[i]=tree[j-1].parent;
        }
    }

    int k;
    double *col=(double *)malloc(n*sizeof(double));
    int *ifrow=(int *)malloc(n*sizeof(int));

    for ( i=0; i<(current_depth-1); i++ ) {

        j=relevant_nodes[i]-1;

        get_col(data,tree[j].covar,n,nc,col);

        if ( tree[j].type == 1 ) {  // binary

            findrows_bin(col,n,tree[j].val,ifrow);

        } else if ( tree[j].type == 2 ) {  // ordinal

            findrows_ord(col,n,tree[j].sign,tree[j].val,ifrow);

        } else if ( tree[j].type == 3 ) {  // nominal

            findrows_nom(col,n,ncat[tree[j].covar-1],tree[j].val,
                         ifrow);
        }

        for ( k=0; k<n; k++) {
            ifcovar[k]=( ifrow[k] && ifcovar[k] );   
        }
    }

    free(ifrow);
    free(col);

//   nsubj is the total number of subjects in this node
    int nsubj=0;
    for ( i=0; i<n; i++ ) {
        nsubj += ifcovar[i];
    }

    return(nsubj);
}



int
get_subj(int current_node,
         int n,
         int nc,
         double **data,
         int *ncat,
         struct node *tree)

//
//  This function accepts as arguments the integer 
//  index of a node and tree information and returns 
//  the total number of subjects in the node.
//

{
    int *ifcovar=(int *)malloc(n*sizeof(int));
    int nsubj=findrows_node(current_node,n,nc,data,ncat,tree,ifcovar);

    free(ifcovar);

    return(nsubj);
}



int
get_subj_treat(int current_node,
               int treat_group,
               int n,
               int nc,
               double **data,
               int *ncat,
               int *treat,
               struct node *tree)

//
//  This function accepts as arguments the integer 
//  index of a node, the integer index of a treatment 
//  group, and information about the overall tree and
//  returns the total number of subjects in the node
//  who are in the specified treatment group.
//

{
    int *ifcovar=(int *)malloc(n*sizeof(int));
    findrows_node(current_node,n,nc,data,ncat,tree,ifcovar);

    int i;
    int ifmember[n];

    for ( i=0; i<n; i++ ) {

        ifmember[i]=0;  // initialize

//       disregard subjects not in current node
        if ( ifcovar[i] == 0 ) continue;

        if ( treat[i] == treat_group ) {
            ifmember[i]=1;
        }
    }

    free(ifcovar);

    int nsubj_treat=0;
    for ( i=0; i<n; i++ ) {
        nsubj_treat += ifmember[i];
    }

/*
//   check
    Rprintf("\n\n--- current_node=%d ---\n",current_node);
    Rprintf("\nifcovar:\n");
    for ( i=0; i<n; i++ ) Rprintf("%d  ",ifcovar[i]);
    Rprintf("\ntreat:\n");
    for ( i=0; i<n; i++ ) Rprintf("%d  ",treat[i]);
    Rprintf("\nifmember (treatval=%d):\n",treat_group);
    for ( i=0; i<n; i++ ) Rprintf("%d  ",ifmember[i]);
*/

    return(nsubj_treat);
}



int
get_node_predict_y(int current_node,
                   int n,
                   int nc,
                   double *y,
                   double **data,
                   int *ncat,
                   int *treat,
                   double r[],
                   struct node *tree)

//
//  This function accepts as arguments the integer 
//  index of a node as well as information from a
//  tree and returns the predicted best treatment 
//  class of the node as an integer.
//
//  Here, the outcome variable y is used to 
//  determine which treatment is best.
//
//  "-7" means neither treatment is best which
//  occurs when there is a tie or when a treatment 
//  group has no one in it.
//

{
    int *ifcovar=(int *)malloc(n*sizeof(int));
    findrows_node(current_node,n,nc,data,ncat,tree,ifcovar);

    int i;
    double ybar0=0;
    double ybar1=0;
    int counter0=0;
    int counter1=0;

    for ( i=0; i<n; i++ ) {

        if ( ifcovar[i] == 0 ) continue;

        if ( treat[i] == 0 ) {

            ybar0 += y[i];
            counter0++;
        }

        if ( treat[i] == 1 ) {

            ybar1 += y[i]; 
            counter1++;
        }
    }

    if ( (counter0 == 0) || (counter1 == 0) ) {
        free(ifcovar);
        return(-7);
    }

    ybar0 /= counter0;
    ybar1 /= counter1;

    r[0]=ybar0;
    r[1]=ybar1;

    double p0=0;
    double p1=0;

    for ( i=0; i<n; i++ ) {

        if ( ifcovar[i] == 0 ) continue;

        if ( (treat[i] == 0) && (y[i] > ybar0) ) p0 += 1.0;
        if ( (treat[i] == 1) && (y[i] > ybar1) ) p1 += 1.0;
    }

    free(ifcovar);

    r[2]=p0/counter0;
    r[3]=p1/counter1;

/*
//   check
    Rprintf("\n--- node=%d ---\n",current_node);
    Rprintf("\ny values:\n");
    for ( i=0; i<n; i++ ) Rprintf("%f,",y[i]);
    Rprintf("\ntreat:\n");
    for ( i=0; i<n; i++ ) Rprintf("%d ",treat[i]);
    Rprintf("\nifcovar:\n");
    for ( i=0; i<n; i++ ) Rprintf("%d ",ifcovar[i]);

    Rprintf("\nybar0=%f ybar1=%f",ybar0,ybar1);
    Rprintf("\nn0=%d n1=%d\n",counter0,counter1);
*/

    if ( ybar0 > ybar1 ) return(0);
    if ( ybar1 > ybar0 ) return(1);

    return(-7);
}



double
get_diff(int current_node,
         int n,
         int nc,
         double *y,
         double **data,
         int *types,
         int *ncat,
         int *treat,
         struct node *tree)

//
//  This function accepts as arguments the current
//  node (as an integer), data information, and
//  tree information. This function returns the
//  "DIFF" value of the current node in the tree.
//
//  "DIFF" is equal to "the squared difference in
//  the response rates of the two treatments in the"
//  current node (Tsai et al. 2016).
//

{
    int *ifcovar=(int *)malloc(n*sizeof(int));
    findrows_node(current_node,n,nc,data,ncat,tree,ifcovar);

    int i;
    double ybar0=0;
    double ybar1=0;
    int counter0=0;
    int counter1=0;

    for ( i=0; i<n; i++ ) {

        if ( ifcovar[i] == 0 ) continue;

        if ( treat[i] == 0 ) {

            ybar0 += y[i];
            counter0++;
        }

        if ( treat[i] == 1 ) {

            ybar1 += y[i]; 
            counter1++;
        }
    }

    free(ifcovar);

    if ( (counter0 == 0) || (counter1 == 0) ) return(0);

    ybar0 /= counter0;
    ybar1 /= counter1;

//   "DIFF" is the squared difference in response rates by treatment
    double diff=(ybar0-ybar1)*(ybar0-ybar1);

    return(diff);
}



int
invert_matrix3(double m[],
               double inv[])
 
//
//  This function accepts values from a 3 by 3 
//  matrix (where the values are stored in an 
//  array by row) and returns the inverse of the
//  matrix.
//
 
{
    double inv0[9];
 
    inv0[0]=(m[4]*m[8]) - (m[7]*m[5]);
    inv0[1]=(m[2]*m[7]) - (m[1]*m[8]);
    inv0[2]=(m[1]*m[5]) - (m[2]*m[4]);
    inv0[3]=(m[5]*m[6]) - (m[3]*m[8]);
    inv0[4]=(m[0]*m[8]) - (m[2]*m[6]);
    inv0[5]=(m[3]*m[2]) - (m[0]*m[5]);
    inv0[6]=(m[3]*m[7]) - (m[6]*m[4]);
    inv0[7]=(m[6]*m[1]) - (m[0]*m[7]);
    inv0[8]=(m[0]*m[4]) - (m[3]*m[1]);
 
    double det;
    det=(m[0]*inv0[0])+(m[1]*inv0[3])+(m[2]*inv0[6]);
/*
//   check
    Rprintf("\ndet=%f\n",det);
*/
 
    if ( det == 0 ) return(0);
 
    det=1.0/det;
 
    int i;
    for ( i=0; i<9; i++ ) {
        inv[i]=inv0[i]*det;
    }
 
    return(1);
}



double
get_coxph_stat2(int nsubj,
                double **X)
{

//   prepare Cox data for R
    SEXP cox_data;
    PROTECT(cox_data=allocMatrix(REALSXP,nsubj,4));

    int i,j;
    for ( i=0; i<nsubj; i++ ) {
        for ( j=0; j<4; j++ ) {
            REAL(cox_data)[i + nsubj*j]=X[i][j];
        }
    }

//   set up call to R function "coxph_R_to_C"
    SEXP coxph_R;
    PROTECT(coxph_R=lang2(install("coxph_R_to_C"),cox_data));

//   get stat from "coxph" function in R
    SEXP val_R=PROTECT(eval(coxph_R,R_GlobalEnv));
    val_R=coerceVector(val_R,REALSXP);
    double val=REAL(val_R)[0];

/*
//   check
    Rprintf("stat from R=%f\n",val);
*/

    UNPROTECT(3);

    return(val);
}



double
get_lmstat_mc(int nsubj,
              double **X)
{

//   prepare Cox data for R
    SEXP lm_data;
    PROTECT(lm_data=allocMatrix(REALSXP,nsubj,3));

    int i,j;
    for ( i=0; i<nsubj; i++ ) {
        for ( j=0; j<3; j++ ) {
            REAL(lm_data)[i + nsubj*j]=X[i][j];
        }
    }

//   set up call to R function "lm_R_to_C"
    SEXP lm_R;
    PROTECT(lm_R=lang2(install("lm_R_to_C"),lm_data));

//   get stat from "coxph" function in R
    SEXP val_R=PROTECT(eval(lm_R,R_GlobalEnv));
    val_R=coerceVector(val_R,REALSXP);
    double val=REAL(val_R)[0];

/*
//   check
    Rprintf("stat from R=%f\n",val);
*/

    UNPROTECT(3);

    return(val);
}



double
get_coxph_stat2_multi(int nsubj,
                      double **X)
{

//   prepare Cox data for R
    SEXP cox_data;
    PROTECT(cox_data=allocMatrix(REALSXP,nsubj,4));

    int i,j;
    for ( i=0; i<nsubj; i++ ) {
        for ( j=0; j<4; j++ ) {
            REAL(cox_data)[i + nsubj*j]=X[i][j];
        }
    }

//   set up call to R function "coxph_R_to_C_multi"
    SEXP coxph_R;
    PROTECT(coxph_R=lang2(install("coxph_R_to_C_multi"),
            cox_data));

//   get stat from "coxph" function in R
    SEXP val_R=PROTECT(eval(coxph_R,R_GlobalEnv));
    val_R=coerceVector(val_R,REALSXP);
    double val=REAL(val_R)[0];

/*
//   check
    Rprintf("stat from R=%f\n",val);
*/

    UNPROTECT(3);

    return(val);
}



double
get_coxph_stat(int current_node,
               int n,
               int nc,
               double *y,
               double **data,
               int *types,
               int *ncat,
               int *treat,
               int *censor,
               int method,
               struct node *tree)

//
//  This function accepts as arguments the integer
//  index of the current node (starting from 1), 
//  array of survival times "y", the matrix of 
//  covariate data "data" of dimension "n" rows and
//  "nc" columns, array of covariate types, array
//  of the number of categories of each variable
//  "ncat", array of treatment assignments "treat",
//  array of censoring indicators where 1 denotes
//  a non-censored observation, and an array of
//  nodes in a tree.
//
//  For a given tree, this function calculates the
//  z^2 statistic of the interaction term in the
//  treatment by split Cox proportional hazards
//  model as a double.
//

{
//   get data and left split indicator of current node
    int *ifcovar=(int *)malloc(n*sizeof(int));
    findrows_node(current_node,n,nc,data,ncat,tree,ifcovar);

    int *ifleft0=(int *)malloc(n*sizeof(int));
    findrows_node(tree[current_node-1].lchild,n,nc,data,ncat,tree,
                  ifleft0);

    int i;
    int nsubj=tree[current_node-1].n;
    double *yvals=(double *)malloc(nsubj*sizeof(double));
    int *trtvals=(int *)malloc(nsubj*sizeof(int));
    int *cvals=(int *)malloc(nsubj*sizeof(int));
    int *ifleft=(int *)malloc(nsubj*sizeof(int));
    int counter=0;

    for ( i=0; i<n; i++ ) {

        if ( ifcovar[i] == 0 ) continue;

        yvals[counter]=y[i];
        trtvals[counter]=treat[i];
        cvals[counter]=censor[i];
        ifleft[counter]=ifleft0[i];

        counter++;
    }

    free(ifcovar);
    free(ifleft0);

/*
//   check
    Rprintf("vals in C=\n");
    for ( i=0; i<nsubj; i++ ) {
        Rprintf("%d ",ifleft[i]);
    }
    Rprintf("\n");
*/

//   create Cox data matrix for Cox model fitting
    double **cox_data=(double **)malloc(nsubj*sizeof(double *));
    for ( i=0; i<nsubj; i++ ) {
        cox_data[i]=(double *)calloc(4,sizeof(double));
        cox_data[i][0]=yvals[i];
        cox_data[i][1]=cvals[i];
        cox_data[i][2]=trtvals[i];
        cox_data[i][3]=ifleft[i];
    }

    free(yvals);
    free(cvals);
    free(trtvals);
    free(ifleft);

//   get Cox model statistic (z^2)
    double stat=0;
    if ( method >= 20 ) {  // multi-treat methods

        stat=get_coxph_stat2_multi(nsubj,cox_data);

    } else {               // binary treat methods

        stat=get_coxph_stat2(nsubj,cox_data);
    }


    for ( i=0; i<nsubj; i++ ) free(cox_data[i]);
    free(cox_data);

    return(stat);
}



int
get_nunique_int(int n,
                int *vector,
                int *unique)

//
//  This function accepts as an argument a vector
//  of "n" integers called "vector" and returns the 
//  number of unique values as an integer. The 
//  actual unique values are saved in array of
//  integers "unique".
//

{
    int i,j;
    int nunique=0;

    for ( i=0; i<n; i++ ) {

        for ( j=0; j<(i+1); j++ ) {

            if ( (vector[i] == vector[j]) && (i != j) ) {
                break;
            }

            if ( i == j ) {
                nunique++;
                unique[nunique-1]=vector[i];
            }
        }
    }

    return(nunique);
}



int
get_node_predict_multi(int current_node,
                       int n,
                       int nc,
                       double *y,
                       double **data,
                       int *ncat,
                       int *treat,
                       struct node *tree)

//
//  This function accepts as arguments the integer 
//  index of a node as well as information from a
//  tree and returns the predicted best treatment 
//  class of the node as an integer.
//
//  Here, the outcome variable y is used to 
//  determine which treatment is best. The data
//  may contain more than 2 treatments.
//
//  "-7" means no treatment is best.
//

{

//   get unique treatment values
    int unique[n];
    int nunique=get_nunique_int(n,treat,unique);

//   get indicator array of which subjects in current node
    int *ifcovar=(int *)malloc(n*sizeof(int));
    findrows_node(current_node,n,nc,data,ncat,tree,ifcovar);

//   initialize counter
    int i;
    int counter[nunique];

    for ( i=0; i<nunique; i++ ) {
        counter[i]=0;
    }

//   get means for each treatment group using data at current node
    int j;
    double ybar[nunique];

    for ( i=0; i<n; i++ ) {

        if ( ifcovar[i] == 0 ) continue;

        for ( j=0; j<nunique; j++ ) {

            if ( treat[i] != unique[j] ) continue;

            ybar[j] += y[i];
            counter[j]++;
        }
    }

    free(ifcovar);

    for ( i=0; i<nunique; i++ ) {

//       good idea when y always + (like survival times)
        if ( counter[i] == 0 ) ybar[i]=0;  

        ybar[i] /= counter[i];
    }

/*
//   check
    int totaln=0;
    for ( i=0; i<n; i++ ) totaln += ifcovar[i];

    if ( totaln == n ) {

        Rprintf("\n---- sample sizes=\n");
        for ( i=0; i<nunique; i++ ) Rprintf("%d ",counter[i]);

        Rprintf("\n---- treatments=\n");
        for ( i=0; i<nunique; i++ ) Rprintf("%d ",unique[i]);

        Rprintf("\n---- means=\n");
        for ( i=0; i<nunique; i++ ) Rprintf("%f ",ybar[i]);
    }
*/

//   return the treatment group with the largest mean y
    int max_treat=-7;
    double max=0.0; // good idea when y always + (like survival times)

    for ( i=0; i<nunique; i++ ) {

        if ( ybar[i] > max ) {
            max_treat=unique[i];
            max=ybar[i];
        }
    }

    return(max_treat);
}



int
get_min_ntrt(int lll,  // integer index of current node
             int n,    // total number of subjects in data
             int nc,   // total number of covariates in data
             double **data,
             int *ncat,
             int *treat,
             struct node *tree)

//
//  This function accepts as arguments an integer 
//  index for the current candidate node "lll",
//  the integer vector of treatment assignments,
//  and information about a tree "tree" and 
//  returns the number of subjects in the 
//  smallest sized treatment group.
//
//  This function is used to assess the minimum
//  treatment group size at a node in the case 
//  of multiple treatments, i.e., more than 2 
//  treatment groups.
//

{

//   get unique treatment values
    int unique[n];
    int nunique=get_nunique_int(n,treat,unique);

//   find # of subjects in each unique treatment group at current node
    int i;
    int ntrt[nunique];

    for ( i=0; i<nunique; i++ ) {
        ntrt[i]=get_subj_treat(lll,unique[i],n,nc,data,ncat,treat,
                               tree);
    }

//   get smallest sample size value
    int min_trt=ntrt[0];
    for ( i=1; i<nunique; i++ ) {
        if ( ntrt[i] < min_trt ) min_trt=ntrt[i];
    }

    return(min_trt);
}



double
get_diff_mc(int current_node,
            int n,
            int nc,
            double *y,
            double **data,
            int *types,
            int *ncat,
            int *treat,
            struct node *tree)

//
//  This function accepts as arguments the current
//  node (as an integer), data information, and
//  tree information. This function returns the
//  "DIFF" value of the current node in the tree
//  for data with a continuous Y and multiple 
//  treatments.
//

{
//   get unique treatment values
    int unique[n];
    int nunique=get_nunique_int(n,treat,unique);

//   get indicator array of which subjects in current node
    int *ifcovar=(int *)malloc(n*sizeof(int));
    findrows_node(current_node,n,nc,data,ncat,tree,ifcovar);

//   initialize counter
    int i;
    int counter[nunique];

    for ( i=0; i<nunique; i++ ) {
        counter[i]=0;
    }

//   get means for each treatment group using data at current node
    int j;
    double ybar[nunique];

    for ( i=0; i<n; i++ ) {

        if ( ifcovar[i] == 0 ) continue;

        for ( j=0; j<nunique; j++ ) {

            if ( treat[i] != unique[j] ) continue;

            ybar[j] += y[i];
            counter[j]++;
        }
    }

    free(ifcovar);

    for ( i=0; i<nunique; i++ ) {

//       good idea when y always + (like survival times)
        if ( counter[i] == 0 ) ybar[i]=0;  

        ybar[i] /= counter[i];
    }

//   find the largest treatment mean
    double ybar_max=0.0;
    for ( i=0; i<nunique; i++ ) {

        if ( ybar[i] > ybar_max ) {
            ybar_max=ybar[i];
        }
    }

//   find the smallest treatment mean
    double ybar_min=ybar_max;
    for ( i=0; i<nunique; i++ ) {

        if ( ybar[i] < ybar_min ) {
            ybar_min=ybar[i];
        }
    }

//   "DIFF" is the squared difference in response rates by treatment
    double diff=(ybar_max-ybar_min)*(ybar_max-ybar_min);

    return(diff);
}



void
get_split0_d(int i,    // 0=left 1=right candidate node
             int lll,  // integer index of current node
             int n,    // total number of subjects in data
             int nc,   // total number of covariates in data
             double *y,
             double **data,
             int *types,
             int *ncat,
             int *treat,
             int method,
             struct node *tree,
             struct cand *cands)

//
//  This function accepts as arguments integer
//  index "i" for whether the current node is the
//  left or right child node, an integer index for
//  the current candidate node "lll", and 
//  information about a tree "tree" and returns
//  information about the current candidate node 
//  in array of cand structs "cands".
//
//  Here, "val" in "cands" is an intermediate
//  calculation in the overall "DIFF" score.
//

{
    double r[4];
    double diff;

    cands[i].n=get_subj(lll,n,nc,data,ncat,tree);

    if ( method >= 20 ) { // multi-treat methods

        diff=get_diff_mc(lll,n,nc,y,data,types,ncat,treat,tree);

//       set n0 and n1 equal to the smallest sized treatment group
        cands[i].n0=get_min_ntrt(lll,n,nc,data,ncat,treat,tree);
        cands[i].n1=cands[i].n0;

//       get best predicted treatment assignment
        cands[i].p=get_node_predict_multi(lll,n,nc,y,data,ncat,treat,
                                          tree);

//       these values not applicable for multiple treatment case
        cands[i].r0=0;
        cands[i].r1=0;
        cands[i].p0=0;
        cands[i].p1=0;

    } else {

        diff=get_diff(lll,n,nc,y,data,types,ncat,treat,tree);

        cands[i].n0=get_subj_treat(lll,0,n,nc,data,ncat,treat,tree);
        cands[i].n1=get_subj_treat(lll,1,n,nc,data,ncat,treat,tree);
        cands[i].p=get_node_predict_y(lll,n,nc,y,data,ncat,treat,r,
                                      tree);
        cands[i].r0=r[0];
        cands[i].r1=r[1];
        cands[i].p0=r[2];
        cands[i].p1=r[3];
    }

    cands[i].val=cands[i].n*diff;
}



double
get_comb_vals(int current_node,
              int current_var,
              int n,
              int nc,
              double *y,
              double **data,
              int *types,
              int *ncat,
              int *treat,
              int nmin2,
              double temp[],
              double cats[],
              int n_cat,
              int k,
              int ncat_original,
              unsigned long long int n_splits,
              unsigned long long int *counter,
              int cat_index,
              int k_index,
              double max,
              int method,
              struct node temp_tree[],
              struct cand *cands,
              struct node *best_split)

//
//  This function accepts as arguments array 
//  "temp" to temporarily save combinations in
//  array of unique category values "cats", the
//  total number of unique categories "n_cat", the
//  size of a single combination "k", the total 
//  number of categories of the variable originally
//  "ncat_original", the current index in the array
//  of categories "cat_index", and the current
//  index in the combination "k_index".
//
//  This function saves the values of the 
//  combinations to be used in candidate left and 
//  right child nodes in arrays "lvals" and "rvals".
//

{
//   save vals once a combination of size k is complete
    if ( k_index == k ) {

        counter++;
////        Rprintf("counter=%llu\n",*counter);

//       all candidate node values have been assessed
        if ( *counter > n_splits ) {
////            Rprintf("counter=%llu\n",*counter);
            return(max);
        }

//       save categories as doubles for child nodes
        int j;
        int ifleft[n_cat];

        for ( j=0; j<n_cat; j++ ) {
            ifleft[j]=0;  // initialize
        }

        int p;
        double nom_lval=0;

        for ( j=0; j<k; j++ ) {
            for ( p=0; p<n_cat; p++ ) {

                if ( temp[j] == cats[p] ) {
                    nom_lval += pow(2,ncat_original-cats[p]);
                    ifleft[p]=1;
                }
            }
        }

        double nom_rval=0;
        for ( j=0; j<n_cat; j++ ) {
            if ( ifleft[j] == 0 ) {
                nom_rval += pow(2,ncat_original-cats[j]);
            }
        }  
/*
//       check
        Rprintf("lval=%f\n",nom_lval);
*/

//       get values for left child candidate node
        temp_tree[current_node].val=nom_lval;
        get_split0_d(0,current_node+1,n,nc,y,data,types,ncat,
                     treat,method,temp_tree,cands);

//       get values for right child candidate node
        temp_tree[current_node].val=nom_rval;
        get_split0_d(1,current_node+1,n,nc,y,data,types,ncat,
                     treat,method,temp_tree,cands);

//       candidate splits with too few subjects are not useful
        if ( (cands[0].n0 <= nmin2) || 
             (cands[0].n1 <= nmin2) || 
             (cands[1].n0 <= nmin2) ||
             (cands[1].n1 <= nmin2) ) {

            return(max);
        }

        double score=(cands[0].val+cands[1].val)/
                     (cands[0].n+cands[1].n);

        if ( score > max ) {

            max=score;

            best_split[0].covar=current_var;
            best_split[0].type=3;
            best_split[0].sign=0;
            best_split[0].val=nom_lval;
            best_split[0].n=cands[0].n;
            best_split[0].n0=cands[0].n0;
            best_split[0].n1=cands[0].n1;
            best_split[0].r0=cands[0].r0;
            best_split[0].r1=cands[0].r1;
            best_split[0].p0=cands[0].p0;
            best_split[0].p1=cands[0].p1;
            best_split[0].predict=cands[0].p;

            best_split[1].covar=current_var;
            best_split[1].type=3;
            best_split[1].sign=0;
            best_split[1].val=nom_rval;
            best_split[1].n=cands[1].n;
            best_split[1].n0=cands[1].n0;
            best_split[1].n1=cands[1].n1;
            best_split[1].r0=cands[1].r0;
            best_split[1].r1=cands[1].r1;
            best_split[1].p0=cands[1].p0;
            best_split[1].p1=cands[1].p1;
            best_split[1].predict=cands[1].p;
        }

        return(max);
    }

    if ( cat_index >= n_cat ) return(max);

    temp[k_index]=cats[cat_index];

    max=get_comb_vals(current_node,current_var,n,nc,y,data,types,ncat,
                      treat,nmin2,temp,cats,n_cat,k,ncat_original,
                      n_splits,counter,cat_index+1,k_index+1,max,
                      method,temp_tree,cands,best_split);

    max=get_comb_vals(current_node,current_var,n,nc,y,data,types,ncat,
                      treat,nmin2,temp,cats,n_cat,k,ncat_original,
                      n_splits,counter,cat_index+1,k_index,max,
                      method,temp_tree,cands,best_split);

    return(max);
}



int
get_split2_d_all(int current_node,
                 int n,
                 int nc,
                 double *y,
                 double **data,
                 int *types,
                 int *ncat,
                 int *treat,
                 struct node *tree,
                 int nmin2,
                 int method,
                 struct node *best_split)

//  
//  This function accepts as arguments the integer
//  index of the current node "current_node", tree
//  information "tree", and the minimum node size
//  "nmin2" and returns the best split as node 
//  struct array of size 2 "best_split".
//
//  If there is no best split, then the function 
//  returns integer "-7".
//
//  All possible splits are considered, and split
//  scores use the DIFF value.
// 

{
//   initialize arrays in temporary tree
    int i;
    struct node temp_tree[current_node+1];

    for ( i=0; i<current_node; i++ ) {

        temp_tree[i].index=tree[i].index;
        temp_tree[i].covar=tree[i].covar;
        temp_tree[i].type=tree[i].type;
        temp_tree[i].sign=tree[i].sign;
        temp_tree[i].val=tree[i].val;
        temp_tree[i].parent=tree[i].parent;
        temp_tree[i].depth=tree[i].depth;
        temp_tree[i].lchild=tree[i].lchild;
        temp_tree[i].rchild=tree[i].rchild;
        temp_tree[i].n=tree[i].n;
        temp_tree[i].n0=tree[i].n0;
        temp_tree[i].n1=tree[i].n1;
    }

    temp_tree[current_node].index=current_node+1;
    temp_tree[current_node].parent=current_node;
    temp_tree[current_node].depth=tree[current_node-1].depth+1;

//   find row indices of current node
    int *ifcovar=(int *)malloc(n*sizeof(int));
    findrows_node(current_node,n,nc,data,ncat,tree,ifcovar);

//   get "DIFF" value of current node
    double diff;
    if ( method == 22 ) { // multitreat continuous method

        diff=get_diff_mc(current_node,n,nc,y,data,types,ncat,treat,
                         tree);

    } else {

        diff=get_diff(current_node,n,nc,y,data,types,ncat,treat,tree);
    }

//   assess candidate splits by variable type
    int j;
    int col_counter=0;
    int nsubj_node=tree[current_node-1].n;
    int nunique;
    double *col=(double *)malloc(n*sizeof(double));
    double *col_node=(double *)malloc(nsubj_node*sizeof(double));
    double *unique=(double *)malloc(nsubj_node*sizeof(double));
    int lll=current_node+1;
    struct cand *cands=malloc(2*sizeof(struct cand));
    double max=diff;
    double score;

    for ( i=0; i<nc; i++ ) {

//       get data column of current node for ith variable
        get_col(data,i+1,n,nc,col);

        col_counter=0;

        for ( j=0; j<n; j++ ) {

            if ( ifcovar[j] == 1 ) {
                col_node[col_counter]=col[j];
                col_counter++;
            }
        }

//       get and sort unique values of the data column
        nunique=get_nunique(nsubj_node,col_node,unique);

        if ( nunique == 1 ) continue; // if all vals same, no split

        double uniquevals[nunique];
        for ( j=0; j<nunique; j++ ) {
            uniquevals[j]=unique[j];
        }

        qsort(uniquevals,nunique,sizeof(double),compare_doubles);

/*
//       check
        Rprintf("\nuniquevals=\n");
	for ( j=0; j<nunique; j++ ) Rprintf("%f ",uniquevals[j]);
*/

//       set the current covariate in the temporary tree
        temp_tree[current_node].covar=i+1;


        if ( types[i] == 1 ) {  // binary

            temp_tree[current_node].type=1;
            temp_tree[current_node].sign=0;

//           get values for left child candidate node
            temp_tree[current_node].val=0;
            get_split0_d(0,lll,n,nc,y,data,types,ncat,treat,method,
                         temp_tree,cands);

//           get values for right child candidate node
            temp_tree[current_node].val=1;
            get_split0_d(1,lll,n,nc,y,data,types,ncat,treat,method,
                         temp_tree,cands);

//           candidate splits with too few subjects are not useful
            if ( (cands[0].n0 <= nmin2) || (cands[0].n1 <= nmin2) || 
                 (cands[1].n0 <= nmin2) || (cands[1].n1 <= nmin2) ) {

                continue;
            }

            score=(cands[0].val+cands[1].val)/(cands[0].n+cands[1].n);

            if ( score > max ) {

                max=score;

                best_split[0].covar=i+1;
                best_split[0].type=1;
                best_split[0].sign=0;
                best_split[0].val=0;
                best_split[0].n=cands[0].n;
                best_split[0].n0=cands[0].n0;
                best_split[0].n1=cands[0].n1;
                best_split[0].r0=cands[0].r0;
                best_split[0].r1=cands[0].r1;
                best_split[0].p0=cands[0].p0;
                best_split[0].p1=cands[0].p1;
                best_split[0].predict=cands[0].p;

                best_split[1].covar=i+1;
                best_split[1].type=1;
                best_split[1].sign=0;
                best_split[1].val=1;
                best_split[1].n=cands[1].n;
                best_split[1].n0=cands[1].n0;
                best_split[1].n1=cands[1].n1;
                best_split[1].r0=cands[1].r0;
                best_split[1].r1=cands[1].r1;
                best_split[1].p0=cands[1].p0;
                best_split[1].p1=cands[1].p1;
                best_split[1].predict=cands[1].p;
            }

        } else if ( types[i] == 2 ) {  // ordinal

            temp_tree[current_node].type=2;

            for ( j=0; j<(nunique-1); j++ ) {

                temp_tree[current_node].val=uniquevals[j];

//               get values for left child candidate node
                temp_tree[current_node].sign=1; // LE
                get_split0_d(0,lll,n,nc,y,data,types,ncat,treat,
                             method,temp_tree,cands);

//               get values for right child candidate node
                temp_tree[current_node].sign=2; // GT
                get_split0_d(1,lll,n,nc,y,data,types,ncat,treat,
                             method,temp_tree,cands);

//               candidate splits with too few subjects are not useful
                if ( (cands[0].n0 <= nmin2) || 
                     (cands[0].n1 <= nmin2) || 
                     (cands[1].n0 <= nmin2) || 
                     (cands[1].n1 <= nmin2) ) {

                    continue;
                }

                score=(cands[0].val+cands[1].val)/
                      (cands[0].n+cands[1].n);

                if ( score > max ) {

                    max=score;

                    best_split[0].covar=i+1;
                    best_split[0].type=2;
                    best_split[0].sign=1; // LE
                    best_split[0].val=uniquevals[j];
                    best_split[0].n=cands[0].n;
                    best_split[0].n0=cands[0].n0;
                    best_split[0].n1=cands[0].n1;
                    best_split[0].r0=cands[0].r0;
                    best_split[0].r1=cands[0].r1;
                    best_split[0].p0=cands[0].p0;
                    best_split[0].p1=cands[0].p1;
                    best_split[0].predict=cands[0].p;

                    best_split[1].covar=i+1;
                    best_split[1].type=2;
                    best_split[1].sign=2; // GT
                    best_split[1].val=uniquevals[j];
                    best_split[1].n=cands[1].n;
                    best_split[1].n0=cands[1].n0;
                    best_split[1].n1=cands[1].n1;
                    best_split[1].r0=cands[1].r0;
                    best_split[1].r1=cands[1].r1;
                    best_split[1].p0=cands[1].p0;
                    best_split[1].p1=cands[1].p1;
                    best_split[1].predict=cands[1].p;
                }
            }

        } else if ( types[i] == 3 ) {  // nominal

            temp_tree[current_node].type=3;
            temp_tree[current_node].sign=0;

//           calculate total number of possible splits
            unsigned long long int n_splits=pow(2,nunique-1)-1;
////	    Rprintf("\nnunique=%d n_splits=%llu\n",nunique,n_splits);
/*
            unsigned long long int test=pow(2,30-1)-1;
	    Rprintf("\nval=%llu\n",test);
*/

            double nom_lval=0;
            double nom_rval=0;

//           2 categories reduce to binary split
            if ( nunique == 2 ) {  

                nom_lval=pow(2,ncat[i]-uniquevals[0]);
                nom_rval=pow(2,ncat[i]-uniquevals[1]);

//               get values for left child candidate node
                temp_tree[current_node].val=nom_lval;
                get_split0_d(0,lll,n,nc,y,data,types,ncat,treat,
                             method,temp_tree,cands);

//               get values for right child candidate node
                temp_tree[current_node].val=nom_rval;
                get_split0_d(1,lll,n,nc,y,data,types,ncat,treat,
                             method,temp_tree,cands);

//               candidate splits with too few subjects are not useful
                if ( (cands[0].n0 <= nmin2) || 
                     (cands[0].n1 <= nmin2) || 
                     (cands[1].n0 <= nmin2) ||
                     (cands[1].n1 <= nmin2) ) {

                    continue;
                }

                score=(cands[0].val+cands[1].val)/
                      (cands[0].n+cands[1].n);

                if ( score > max ) {

                    max=score;

                    best_split[0].covar=i+1;
                    best_split[0].type=3;
                    best_split[0].sign=0;
                    best_split[0].val=nom_lval;
                    best_split[0].n=cands[0].n;
                    best_split[0].n0=cands[0].n0;
                    best_split[0].n1=cands[0].n1;
                    best_split[0].r0=cands[0].r0;
                    best_split[0].r1=cands[0].r1;
                    best_split[0].p0=cands[0].p0;
                    best_split[0].p1=cands[0].p1;
                    best_split[0].predict=cands[0].p;

                    best_split[1].covar=i+1;
                    best_split[1].type=3;
                    best_split[1].sign=0;
                    best_split[1].val=nom_rval;
                    best_split[1].n=cands[1].n;
                    best_split[1].n0=cands[1].n0;
                    best_split[1].n1=cands[1].n1;
                    best_split[1].r0=cands[1].r0;
                    best_split[1].r1=cands[1].r1;
                    best_split[1].p0=cands[1].p0;
                    best_split[1].p1=cands[1].p1;
                    best_split[1].predict=cands[1].p;
                }

//           for nominal variables with more than 2 categories
            } else {

                int z,k;
                int current_var=i+1;
                int k_index,cat_index;
                int max_k=floor(nunique/2);
                unsigned long long int counter=0;

                for ( z=0; z<max_k; z++ ) {

                    k=z+1;        // size of 1 combination
                    k_index=0;    // current index in a single comb.
                    cat_index=0;  // current index in list of cats.

                    double temp[k]; // array to temporarily store 
                                    // a combination

                    max=get_comb_vals(current_node,current_var,
                                      n,nc,y,data,types,ncat,treat,
                                      nmin2,temp,uniquevals,nunique,
                                      k,ncat[i],n_splits,&counter,
                                      cat_index,k_index,max,method,
                                      temp_tree,cands,best_split);
                }
            }
        }
    }

    free(cands);
    free(ifcovar);
    free(col);
    free(col_node);
    free(unique);

    if ( max == diff ) return(-7);
    return(0);
}



int
get_split2_mtry(int current_node,
                int n,
                int nc,
                double *y,
                double **data,
                int *types,
                int *ncat,
                int *treat,
                struct node *tree,
                int nmin2,
                int mtry,
                int method,
                struct node *best_split)

//  
//  This function accepts as arguments the integer
//  index of the current node "current_node", tree
//  information "tree", the minimum node size 
//  "nmin2", and the number of variables to randomly
//  select "mtry" and returns the best split as node 
//  struct array of size 2 "best_split".
//
//  If there is no best split, then the function 
//  returns integer "-7".
//
//  All possible splits of the mtry selected 
//  variables are considered, and split scores use 
//  the DIFF value.
// 

{
//   select mtry variables
    int mtry_vars[mtry];
    int mtry_types[mtry];

    int i;
    for ( i=0; i<mtry; i++ ) {
        mtry_vars[i]=-7;  // initialize mtry candidate vars
        mtry_types[i]=-7; // initialize mtry candidate types
    }

    int mtry00=mtry;
    if ( mtry > nc ) mtry00=nc;

    int mtry_index[mtry00];
    sample_without_replace(nc,mtry00,mtry_index);

    for ( i=0; i<mtry00; i++ ) {
        mtry_vars[i]=mtry_index[i]+1;
        mtry_types[i]=types[mtry_index[i]];
    }

/*
//   check
    Rprintf("\nRandomly selected mtry=%d variables at node %d:\n",
            mtry,current_node);
    for ( i=0; i<mtry; i++ ) {
        Rprintf("i=%d\tvar=%d\ttype=%d\n",i,mtry_vars[i],
                mtry_types[i]);
    }
*/

//   initialize arrays in temporary tree
    struct node temp_tree[current_node+1];
    for ( i=0; i<current_node; i++ ) {

        temp_tree[i].index=tree[i].index;
        temp_tree[i].covar=tree[i].covar;
        temp_tree[i].type=tree[i].type;
        temp_tree[i].sign=tree[i].sign;
        temp_tree[i].val=tree[i].val;
        temp_tree[i].parent=tree[i].parent;
        temp_tree[i].depth=tree[i].depth;
        temp_tree[i].lchild=tree[i].lchild;
        temp_tree[i].rchild=tree[i].rchild;
        temp_tree[i].n=tree[i].n;
        temp_tree[i].n0=tree[i].n0;
        temp_tree[i].n1=tree[i].n1;
    }

    temp_tree[current_node].index=current_node+1;
    temp_tree[current_node].parent=current_node;
    temp_tree[current_node].depth=tree[current_node-1].depth+1;

//   find row indices of current node
    int *ifcovar=(int *)malloc(n*sizeof(int));
    findrows_node(current_node,n,nc,data,ncat,tree,ifcovar);

//   get "DIFF" value of current node
    double diff;
    if ( method == 23 ) { // multitrt continuous method

        diff=get_diff_mc(current_node,n,nc,y,data,types,ncat,treat,
                         tree);

    } else {

        diff=get_diff(current_node,n,nc,y,data,types,ncat,treat,tree);
    }

//   assess candidate splits by variable type
    int j;
    int col_counter=0;
    int nsubj_node=tree[current_node-1].n;
    int nunique;
    double *col=(double *)malloc(n*sizeof(double));
    double *col_node=(double *)malloc(nsubj_node*sizeof(double));
    double *unique=(double *)malloc(nsubj_node*sizeof(double));
    int lll=current_node+1;
    struct cand *cands=malloc(2*sizeof(struct cand));
    double max=diff;
    double score;

    for ( i=0; i<mtry; i++ ) {

//       type may equal -7 if mtry > nc
        if ( mtry_types[i] == -7 ) break;

//       get data column of current node for ith mtry variable
        get_col(data,mtry_vars[i],n,nc,col);

        col_counter=0;

        for ( j=0; j<n; j++ ) {

            if ( ifcovar[j] == 1 ) {
                col_node[col_counter]=col[j];
                col_counter++;
            }
        }

//       get and sort unique values of the data column
        nunique=get_nunique(nsubj_node,col_node,unique);

        if ( nunique == 1 ) continue; // if all vals same, no split

        double uniquevals[nunique];
        for ( j=0; j<nunique; j++ ) {
            uniquevals[j]=unique[j];
        }

        qsort(uniquevals,nunique,sizeof(double),compare_doubles);

//       set the current covariate in the temporary tree
        temp_tree[current_node].covar=mtry_vars[i];


        if ( mtry_types[i] == 1 ) {  // binary

            temp_tree[current_node].type=1;
            temp_tree[current_node].sign=0;

//           get values for left child candidate node
            temp_tree[current_node].val=0;
            get_split0_d(0,lll,n,nc,y,data,types,ncat,treat,method,
                         temp_tree,cands);

//           get values for right child candidate node
            temp_tree[current_node].val=1;
            get_split0_d(1,lll,n,nc,y,data,types,ncat,treat,method,
                         temp_tree,cands);

//           candidate splits with too few subjects are not useful
            if ( (cands[0].n0 <= nmin2) || (cands[0].n1 <= nmin2) || 
                 (cands[1].n0 <= nmin2) || (cands[1].n1 <= nmin2) ) {

                continue;
            }

            score=(cands[0].val+cands[1].val)/(cands[0].n+cands[1].n);

            if ( score > max ) {

                max=score;

                best_split[0].covar=mtry_vars[i];
                best_split[0].type=1;
                best_split[0].sign=0;
                best_split[0].val=0;
                best_split[0].n=cands[0].n;
                best_split[0].n0=cands[0].n0;
                best_split[0].n1=cands[0].n1;
                best_split[0].r0=cands[0].r0;
                best_split[0].r1=cands[0].r1;
                best_split[0].p0=cands[0].p0;
                best_split[0].p1=cands[0].p1;
                best_split[0].predict=cands[0].p;

                best_split[1].covar=mtry_vars[i];
                best_split[1].type=1;
                best_split[1].sign=0;
                best_split[1].val=1;
                best_split[1].n=cands[1].n;
                best_split[1].n0=cands[1].n0;
                best_split[1].n1=cands[1].n1;
                best_split[1].r0=cands[1].r0;
                best_split[1].r1=cands[1].r1;
                best_split[1].p0=cands[1].p0;
                best_split[1].p1=cands[1].p1;
                best_split[1].predict=cands[1].p;
            }

        } else if ( mtry_types[i] == 2 ) {  // ordinal

            temp_tree[current_node].type=2;

            for ( j=0; j<(nunique-1); j++ ) {

                temp_tree[current_node].val=uniquevals[j];

//               get values for left child candidate node
                temp_tree[current_node].sign=1; // LE
                get_split0_d(0,lll,n,nc,y,data,types,ncat,treat,
                             method,temp_tree,cands);

//               get values for right child candidate node
                temp_tree[current_node].sign=2; // GT
                get_split0_d(1,lll,n,nc,y,data,types,ncat,treat,
                             method,temp_tree,cands);

//               candidate splits with too few subjects are not useful
                if ( (cands[0].n0 <= nmin2) || 
                     (cands[0].n1 <= nmin2) || 
                     (cands[1].n0 <= nmin2) || 
                     (cands[1].n1 <= nmin2) ) {

                    continue;
                }

                score=(cands[0].val+cands[1].val)/
                      (cands[0].n+cands[1].n);

                if ( score > max ) {

                    max=score;

                    best_split[0].covar=mtry_vars[i];
                    best_split[0].type=2;
                    best_split[0].sign=1; // LE
                    best_split[0].val=uniquevals[j];
                    best_split[0].n=cands[0].n;
                    best_split[0].n0=cands[0].n0;
                    best_split[0].n1=cands[0].n1;
                    best_split[0].r0=cands[0].r0;
                    best_split[0].r1=cands[0].r1;
                    best_split[0].p0=cands[0].p0;
                    best_split[0].p1=cands[0].p1;
                    best_split[0].predict=cands[0].p;

                    best_split[1].covar=mtry_vars[i];
                    best_split[1].type=2;
                    best_split[1].sign=2; // GT
                    best_split[1].val=uniquevals[j];
                    best_split[1].n=cands[1].n;
                    best_split[1].n0=cands[1].n0;
                    best_split[1].n1=cands[1].n1;
                    best_split[1].r0=cands[1].r0;
                    best_split[1].r1=cands[1].r1;
                    best_split[1].p0=cands[1].p0;
                    best_split[1].p1=cands[1].p1;
                    best_split[1].predict=cands[1].p;
                }
            }

        } else if ( mtry_types[i] == 3 ) {  // nominal

            temp_tree[current_node].type=3;
            temp_tree[current_node].sign=0;

//           calculate total number of possible splits
            unsigned long long int n_splits=pow(2,nunique-1)-1;
////	    Rprintf("\nnunique=%d n_splits=%llu\n",nunique,n_splits);

//           2 categories reduce to binary split
            if ( nunique == 2 ) {  

                double nom_lval=pow(2,ncat[mtry_vars[i]-1]-
                                    uniquevals[0]);

                double nom_rval=pow(2,ncat[mtry_vars[i]-1]-
                                    uniquevals[1]);

//               get values for left child candidate node
                temp_tree[current_node].val=nom_lval;
                get_split0_d(0,lll,n,nc,y,data,types,ncat,treat,
                             method,temp_tree,cands);

//               get values for right child candidate node
                temp_tree[current_node].val=nom_rval;
                get_split0_d(1,lll,n,nc,y,data,types,ncat,treat,
                             method,temp_tree,cands);

//               candidate splits with too few subjects are not useful
                if ( (cands[0].n0 <= nmin2) || 
                     (cands[0].n1 <= nmin2) || 
                     (cands[1].n0 <= nmin2) ||
                     (cands[1].n1 <= nmin2) ) {

                    continue;
                }

                score=(cands[0].val+cands[1].val)/
                      (cands[0].n+cands[1].n);

                if ( score > max ) {

                    max=score;

                    best_split[0].covar=mtry_vars[i];
                    best_split[0].type=3;
                    best_split[0].sign=0;
                    best_split[0].val=nom_lval;
                    best_split[0].n=cands[0].n;
                    best_split[0].n0=cands[0].n0;
                    best_split[0].n1=cands[0].n1;
                    best_split[0].r0=cands[0].r0;
                    best_split[0].r1=cands[0].r1;
                    best_split[0].p0=cands[0].p0;
                    best_split[0].p1=cands[0].p1;
                    best_split[0].predict=cands[0].p;

                    best_split[1].covar=mtry_vars[i];
                    best_split[1].type=3;
                    best_split[1].sign=0;
                    best_split[1].val=nom_rval;
                    best_split[1].n=cands[1].n;
                    best_split[1].n0=cands[1].n0;
                    best_split[1].n1=cands[1].n1;
                    best_split[1].r0=cands[1].r0;
                    best_split[1].r1=cands[1].r1;
                    best_split[1].p0=cands[1].p0;
                    best_split[1].p1=cands[1].p1;
                    best_split[1].predict=cands[1].p;
                }

//           for nominal variables with more than 2 categories
            } else {

                int z,k;
                int current_var=mtry_vars[i];
                int k_index,cat_index;
                int max_k=floor(nunique/2);
                unsigned long long int counter=0;

                for ( z=0; z<max_k; z++ ) {

                    k=z+1;        // size of 1 combination
                    k_index=0;    // current index in a single comb.
                    cat_index=0;  // current index in list of cats.

                    double temp[k]; // array to temporarily store 
                                    // a combination

                    max=get_comb_vals(current_node,current_var,
                                      n,nc,y,data,types,ncat,treat,
                                      nmin2,temp,uniquevals,nunique,
                                      k,ncat[mtry_vars[i]-1],n_splits,
                                      &counter,cat_index,k_index,max,
                                      method,temp_tree,cands,
                                      best_split);
                }
            }
        }
    }

    free(cands);
    free(ifcovar);
    free(col);
    free(col_node);
    free(unique);

    if ( max == diff ) return(-7);
    return(0);
}



void
maketree_diff(int n,
              int nc,
              double *y,
              double **data,
              int *types,
              int *ncat,
              int *treat,
              struct node *tree,
              int nmin2,
              int maxdepth,
              int method)

//
//  This function returns a tree using the modified 
//  classification tree method developed by Zhang. 
//  This method uses the "DIFF" score to split each
//  node.
//

{
//   initialize root node
    tree[0].index=1;
    tree[0].depth=1;
    tree[0].n=n;

    if ( method >= 20 ) {  // multi-treat methods

        tree[0].n0=get_min_ntrt(1,n,nc,data,ncat,treat,tree);
        tree[0].n1=tree[0].n0;

        tree[0].predict=get_node_predict_multi(1,n,nc,y,data,ncat,
                                               treat,tree);
        tree[0].r0=0;
        tree[0].r1=0;
        tree[0].p0=0;
        tree[0].p1=0;

    } else {               // binary treat methods

        tree[0].n0=get_subj_treat(1,0,n,nc,data,ncat,treat,tree);
        tree[0].n1=get_subj_treat(1,1,n,nc,data,ncat,treat,tree);

        double r[4];
        tree[0].predict=get_node_predict_y(1,n,nc,y,data,ncat,treat,
                                           r,tree);
        tree[0].r0=r[0];
        tree[0].r1=r[1];
        tree[0].p0=r[2];
        tree[0].p1=r[3];
    }

//   build tree
    int i;
    int current_node,isplit,max,d;
    struct node *best_split=malloc(2*sizeof(struct node));

    for ( i=0; i<MAXSIZE; i++ ) {

        current_node=i+1;

        d=tree[current_node-1].depth+1;
        if ( (maxdepth > 0) && (d > maxdepth) ) break;

        isplit=get_split2_d_all(current_node,n,nc,y,data,types,
                                ncat,treat,tree,nmin2,method,
                                best_split);

        if ( isplit != -7 ) {  // "-7" = there should be no split

            max=get_number_of_nodes(MAXSIZE,tree);
            tree[i].lchild=max+1;
            tree[i].rchild=max+2;

//           add left child node to tree
            tree[max].index=tree[i].lchild;
            tree[max].covar=best_split[0].covar;
            tree[max].type=best_split[0].type;
            tree[max].sign=best_split[0].sign;
            tree[max].val=best_split[0].val;
            tree[max].parent=current_node;
            tree[max].depth=d;
            tree[max].n=best_split[0].n;
            tree[max].n0=best_split[0].n0;
            tree[max].n1=best_split[0].n1;
            tree[max].r0=best_split[0].r0;
            tree[max].r1=best_split[0].r1;
            tree[max].p0=best_split[0].p0;
            tree[max].p1=best_split[0].p1;
            tree[max].predict=best_split[0].predict;

//           add right child node to tree
            tree[max+1].index=tree[i].rchild;
            tree[max+1].covar=best_split[1].covar;
            tree[max+1].type=best_split[1].type;
            tree[max+1].sign=best_split[1].sign;
            tree[max+1].val=best_split[1].val;
            tree[max+1].parent=current_node;
            tree[max+1].depth=d;
            tree[max+1].n=best_split[1].n;
            tree[max+1].n0=best_split[1].n0;
            tree[max+1].n1=best_split[1].n1;
            tree[max+1].r0=best_split[1].r0;
            tree[max+1].r1=best_split[1].r1;
            tree[max+1].p0=best_split[1].p0;
            tree[max+1].p1=best_split[1].p1;
            tree[max+1].predict=best_split[1].predict;
/*
	    Rprintf("\nlsplitval=%Lf rsplitval=%Lf\n",
                    tree[max].val,tree[max+1].val);
*/
        }

        max=get_number_of_nodes(MAXSIZE,tree);
        if ( (isplit == -7) && (current_node == max) ) break;
    }

    free(best_split);

/*
//   print tree
    Rprintf("\nTree from maketree_diff (nmin2=%d): \n",nmin2);
    for ( i=0; i<max; i++ ) {

        Rprintf("index=%d covar=%d type=%d sign=%d val=%f lchild=%d rchild=%d n=%d ntrt0=%d ntrt1=%d besttrt=%d\n",
               tree[i].index,tree[i].covar,tree[i].type,
               tree[i].sign,tree[i].val,tree[i].lchild,tree[i].rchild,
               tree[i].n,tree[i].n0,tree[i].n1,tree[i].predict);
    }
*/
}



void
get_split0_cox(int i,    // 0=left 1=right candidate node
               int lll,  // integer index of current node
               int n,    // total number of subjects in data
               int nc,   // total number of covariates in data
               double *y,
               double **data,
               int *types,
               int *ncat,
               int *treat,
               int *censor,
               int method,
               struct node *tree,
               struct cand *cands)

//
//  This function accepts as arguments integer
//  index "i" for whether the current node is the
//  left or right child node, an integer index for
//  the current candidate node "lll", and 
//  information about a tree "tree" and returns
//  information about the current candidate node 
//  in array of cand structs "cands".
//
//  Here, "val" is set to 0.
//

{

    cands[i].n=get_subj(lll,n,nc,data,ncat,tree);

    if ( method >= 20 ) {  // multi-treat methods

//       set n0 and n1 equal to the smallest sized treatment group
        cands[i].n0=get_min_ntrt(lll,n,nc,data,ncat,treat,tree);
        cands[i].n1=cands[i].n0;

//       get best predicted treatment assignment
        cands[i].p=get_node_predict_multi(lll,n,nc,y,data,ncat,treat,
                                          tree);

//       these values not applicable for multiple treatment case
        cands[i].r0=0;
        cands[i].r1=0;
        cands[i].p0=0;
        cands[i].p1=0;

        cands[i].val=0;

    } else {               // binary treat methods

        cands[i].n0=get_subj_treat(lll,0,n,nc,data,ncat,treat,tree);
        cands[i].n1=get_subj_treat(lll,1,n,nc,data,ncat,treat,tree);

        double r[4];
        cands[i].p=get_node_predict_y(lll,n,nc,y,data,ncat,treat,r,
                                      tree);
        cands[i].r0=r[0];
        cands[i].r1=r[1];
        cands[i].p0=r[2];
        cands[i].p1=r[3];

        cands[i].val=0;

/*
        cands[i].val=get_coxph_stat(tree[lll-1].parent,n,nc,y,data,
                                    types,ncat,treat,censor,method,
                                    tree);
*/
    }
}



double
get_comb_vals_surv(int current_node,
                   int current_var,
                   int n,
                   int nc,
                   double *y,
                   double **data,
                   int *types,
                   int *ncat,
                   int *treat,
                   int *censor,
                   int nmin2,
                   int method,
                   double temp[],
                   double cats[],
                   int n_cat,
                   int k,
                   int ncat_original,
                   unsigned long long int n_splits,
                   unsigned long long int *counter,
                   int cat_index,
                   int k_index,
                   double max,
                   struct node temp_tree[],
                   struct cand *cands,
                   struct node *best_split)

//
//  This function accepts as arguments array 
//  "temp" to temporarily save combinations in,
//  array of unique category values "cats", the
//  total number of unique categories "n_cat", the
//  size of a single combination "k", the total 
//  number of categories of the variable originally
//  "ncat_original", the current index in the array
//  of categories "cat_index", and the current
//  index in the combination "k_index".
//
//  This function saves the values of the 
//  combinations to be used in candidate left and 
//  right child nodes in arrays "lvals" and "rvals".
//

{
//   save vals once a combination of size k is complete
    if ( k_index == k ) {

        (*counter)++;
////        Rprintf("counter=%llu\n",*counter);

//       all candidate node values have been assessed
        if ( *counter > n_splits ) {
            return(max);
        }

//       save categories as doubles for child nodes
        int j;
        int ifleft[n_cat];

        for ( j=0; j<n_cat; j++ ) {
            ifleft[j]=0;  // initialize
        }

        int p;
        double nom_lval=0;

        for ( j=0; j<k; j++ ) {
            for ( p=0; p<n_cat; p++ ) {

                if ( temp[j] == cats[p] ) {
                    nom_lval += pow(2,ncat_original-cats[p]);
                    ifleft[p]=1;
                }
            }
        }

        double nom_rval=0;
        for ( j=0; j<n_cat; j++ ) {
            if ( ifleft[j] == 0 ) {
                nom_rval += pow(2,ncat_original-cats[j]);
            }
        }  

//       get values for left child candidate node
        temp_tree[current_node].val=nom_lval;
        get_split0_cox(0,current_node+1,n,nc,y,data,types,ncat,
                       treat,censor,method,temp_tree,cands);

//       get values for right child candidate node
        temp_tree[current_node+1].val=nom_rval;
        get_split0_cox(1,current_node+2,n,nc,y,data,types,ncat,
                       treat,censor,method,temp_tree,cands);

//       candidate splits with too few subjects are not useful
        if ( (cands[0].n0 <= nmin2) || 
             (cands[0].n1 <= nmin2) || 
             (cands[1].n0 <= nmin2) ||
             (cands[1].n1 <= nmin2) ) {

            return(max);
        }

        double score=get_coxph_stat(current_node,n,nc,y,data,types,
                                    ncat,treat,censor,method,
                                    temp_tree);

        if ( score > max ) {

            max=score;
            best_split[0].G=score;

            best_split[0].covar=current_var;
            best_split[0].type=3;
            best_split[0].sign=0;
            best_split[0].val=nom_lval;
            best_split[0].n=cands[0].n;
            best_split[0].n0=cands[0].n0;
            best_split[0].n1=cands[0].n1;
            best_split[0].r0=cands[0].r0;
            best_split[0].r1=cands[0].r1;
            best_split[0].p0=cands[0].p0;
            best_split[0].p1=cands[0].p1;
            best_split[0].predict=cands[0].p;

            best_split[1].covar=current_var;
            best_split[1].type=3;
            best_split[1].sign=0;
            best_split[1].val=nom_rval;
            best_split[1].n=cands[1].n;
            best_split[1].n0=cands[1].n0;
            best_split[1].n1=cands[1].n1;
            best_split[1].r0=cands[1].r0;
            best_split[1].r1=cands[1].r1;
            best_split[1].p0=cands[1].p0;
            best_split[1].p1=cands[1].p1;
            best_split[1].predict=cands[1].p;
        }

        return(max);
    }

    if ( cat_index >= n_cat ) return(max);

    temp[k_index]=cats[cat_index];

    max=get_comb_vals_surv(current_node,current_var,n,nc,y,data,types,
                           ncat,treat,censor,nmin2,method,temp,cats,
                           n_cat,k,ncat_original,n_splits,counter,
                           cat_index+1,k_index+1,max,temp_tree,cands,
                           best_split);

    max=get_comb_vals_surv(current_node,current_var,n,nc,y,data,types,
                           ncat,treat,censor,nmin2,method,temp,cats,
                           n_cat,k,ncat_original,n_splits,counter,
                           cat_index+1,k_index,max,temp_tree,cands,
                           best_split);

    return(max);
}



int
get_split2_cox_all(int current_node,
                   int n,
                   int nc,
                   double *y,
                   double **data,
                   int *types,
                   int *ncat,
                   int *treat,
                   int *censor,
                   struct node *tree,
                   int nmin2,
                   int method,
                   struct node *best_split)

//  
//  This function accepts as arguments the integer
//  index of the current node "current_node", tree
//  information "tree", and the minimum node size
//  "nmin2" and returns the best split as node 
//  struct array of size 2 "best_split".
//
//  If there is no best split, then the function 
//  returns integer "-7".
//
//  All possible splits are considered, and split
//  scores use the z^2 statistic of the interaction
//  term in a Cox proportional hazards model.
// 

{
//   initialize arrays in temporary tree
    int i;
    struct node temp_tree[current_node+2];

    for ( i=0; i<current_node; i++ ) {

        temp_tree[i].index=tree[i].index;
        temp_tree[i].covar=tree[i].covar;
        temp_tree[i].type=tree[i].type;
        temp_tree[i].sign=tree[i].sign;
        temp_tree[i].val=tree[i].val;
        temp_tree[i].parent=tree[i].parent;
        temp_tree[i].depth=tree[i].depth;
        temp_tree[i].lchild=tree[i].lchild;
        temp_tree[i].rchild=tree[i].rchild;
        temp_tree[i].n=tree[i].n;
        temp_tree[i].n0=tree[i].n0;
        temp_tree[i].n1=tree[i].n1;
    }

    temp_tree[current_node-1].lchild=current_node+1;
    temp_tree[current_node].index=current_node+1;
    temp_tree[current_node].parent=current_node;
    temp_tree[current_node].depth=tree[current_node-1].depth+1;

    temp_tree[current_node-1].rchild=current_node+2;
    temp_tree[current_node+1].index=current_node+2;
    temp_tree[current_node+1].parent=current_node;
    temp_tree[current_node+1].depth=tree[current_node-1].depth+1;

//   find row indices of current node
    int *ifcovar=(int *)malloc(n*sizeof(int));
    findrows_node(current_node,n,nc,data,ncat,tree,ifcovar);

//   assess candidate splits by variable type
    int j;
    int col_counter=0;
    int nsubj_node=tree[current_node-1].n;
    int nunique;
    double *col=(double *)malloc(n*sizeof(double));
    double *col_node=(double *)malloc(nsubj_node*sizeof(double));
    double *unique=(double *)malloc(nsubj_node*sizeof(double));
    int lll=current_node+1;
    struct cand *cands=malloc(2*sizeof(struct cand));
    double max=0;
    double score;

    for ( i=0; i<nc; i++ ) {

//       get data column of current node for ith variable
        get_col(data,i+1,n,nc,col);

        col_counter=0;

        for ( j=0; j<n; j++ ) {

            if ( ifcovar[j] == 1 ) {
                col_node[col_counter]=col[j];
                col_counter++;
            }
        }

//       get and sort unique values of the data column
        nunique=get_nunique(nsubj_node,col_node,unique);

        if ( nunique == 1 ) continue; // if all vals same, no split

        double uniquevals[nunique];
        for ( j=0; j<nunique; j++ ) {
            uniquevals[j]=unique[j];
        }

        qsort(uniquevals,nunique,sizeof(double),compare_doubles);

//       set the current covariate in the temporary tree
        temp_tree[current_node].covar=i+1;
        temp_tree[current_node+1].covar=i+1;


        if ( types[i] == 1 ) {  // binary

//           get values for left child candidate node
            temp_tree[current_node].type=1;
            temp_tree[current_node].sign=0;
            temp_tree[current_node].val=0;

            get_split0_cox(0,lll,n,nc,y,data,types,ncat,treat,
                           censor,method,temp_tree,cands);

//           get values for right child candidate node
            temp_tree[current_node+1].type=1;
            temp_tree[current_node+1].sign=0;
            temp_tree[current_node+1].val=1;

            get_split0_cox(1,lll+1,n,nc,y,data,types,ncat,treat,
                           censor,method,temp_tree,cands);

//           candidate splits with too few subjects are not useful
            if ( (cands[0].n0 <= nmin2) || (cands[0].n1 <= nmin2) || 
                 (cands[1].n0 <= nmin2) || (cands[1].n1 <= nmin2) ) {

                continue;
            }


            score=get_coxph_stat(current_node,n,nc,y,data,types,
                                 ncat,treat,censor,method,temp_tree);

            if ( score > max ) {

                max=score;
                best_split[0].G=score;

                best_split[0].covar=i+1;
                best_split[0].type=1;
                best_split[0].sign=0;
                best_split[0].val=0;
                best_split[0].n=cands[0].n;
                best_split[0].n0=cands[0].n0;
                best_split[0].n1=cands[0].n1;
                best_split[0].r0=cands[0].r0;
                best_split[0].r1=cands[0].r1;
                best_split[0].p0=cands[0].p0;
                best_split[0].p1=cands[0].p1;
                best_split[0].predict=cands[0].p;

                best_split[1].covar=i+1;
                best_split[1].type=1;
                best_split[1].sign=0;
                best_split[1].val=1;
                best_split[1].n=cands[1].n;
                best_split[1].n0=cands[1].n0;
                best_split[1].n1=cands[1].n1;
                best_split[1].r0=cands[1].r0;
                best_split[1].r1=cands[1].r1;
                best_split[1].p0=cands[1].p0;
                best_split[1].p1=cands[1].p1;
                best_split[1].predict=cands[1].p;
            }

        } else if ( types[i] == 2 ) {  // ordinal

            temp_tree[current_node].type=2;
            temp_tree[current_node+1].type=2;

            for ( j=0; j<(nunique-1); j++ ) {

//               get values for left child candidate node
                temp_tree[current_node].val=uniquevals[j];
                temp_tree[current_node].sign=1; // LE

                get_split0_cox(0,lll,n,nc,y,data,types,ncat,treat,
                               censor,method,temp_tree,cands);

//               get values for right child candidate node
                temp_tree[current_node+1].val=uniquevals[j];
                temp_tree[current_node+1].sign=2; // GT

                get_split0_cox(1,lll+1,n,nc,y,data,types,ncat,treat,
                               censor,method,temp_tree,cands);

//               candidate splits with too few subjects are not useful
                if ( (cands[0].n0 <= nmin2) || 
                     (cands[0].n1 <= nmin2) || 
                     (cands[1].n0 <= nmin2) || 
                     (cands[1].n1 <= nmin2) ) {

                    continue;
                }

                score=get_coxph_stat(current_node,n,nc,y,data,types,
                                     ncat,treat,censor,method,
                                     temp_tree);

                if ( score > max ) {

                    max=score;
                    best_split[0].G=score;

                    best_split[0].covar=i+1;
                    best_split[0].type=2;
                    best_split[0].sign=1; // LE
                    best_split[0].val=uniquevals[j];
                    best_split[0].n=cands[0].n;
                    best_split[0].n0=cands[0].n0;
                    best_split[0].n1=cands[0].n1;
                    best_split[0].r0=cands[0].r0;
                    best_split[0].r1=cands[0].r1;
                    best_split[0].p0=cands[0].p0;
                    best_split[0].p1=cands[0].p1;
                    best_split[0].predict=cands[0].p;

                    best_split[1].covar=i+1;
                    best_split[1].type=2;
                    best_split[1].sign=2; // GT
                    best_split[1].val=uniquevals[j];
                    best_split[1].n=cands[1].n;
                    best_split[1].n0=cands[1].n0;
                    best_split[1].n1=cands[1].n1;
                    best_split[1].r0=cands[1].r0;
                    best_split[1].r1=cands[1].r1;
                    best_split[1].p0=cands[1].p0;
                    best_split[1].p1=cands[1].p1;
                    best_split[1].predict=cands[1].p;
                }
            }

        } else if ( types[i] == 3 ) {  // nominal

            temp_tree[current_node].type=3;
            temp_tree[current_node].sign=0;

            temp_tree[current_node+1].type=3;
            temp_tree[current_node+1].sign=0;

//           calculate total number of possible splits
            unsigned long long int n_splits=pow(2,nunique-1)-1;
////	    Rprintf("\nnunique=%d n_splits=%llu\n",nunique,n_splits);

            double nom_lval=0;
            double nom_rval=0;

//           2 categories reduce to binary split
            if ( nunique == 2 ) {  

                nom_lval=pow(2,ncat[i]-uniquevals[0]);
                nom_rval=pow(2,ncat[i]-uniquevals[1]);

//               get values for left child candidate node
                temp_tree[current_node].val=nom_lval;
                get_split0_cox(0,lll,n,nc,y,data,types,ncat,treat,
                               censor,method,temp_tree,cands);

//               get values for right child candidate node
                temp_tree[current_node+1].val=nom_rval;
                get_split0_cox(1,lll+1,n,nc,y,data,types,ncat,treat,
                               censor,method,temp_tree,cands);

//               candidate splits with too few subjects are not useful
                if ( (cands[0].n0 <= nmin2) || 
                     (cands[0].n1 <= nmin2) || 
                     (cands[1].n0 <= nmin2) ||
                     (cands[1].n1 <= nmin2) ) {

                    continue;
                }

                score=get_coxph_stat(current_node,n,nc,y,data,types,
                                     ncat,treat,censor,method,
                                     temp_tree);

                if ( score > max ) {

                    max=score;
                    best_split[0].G=score;

                    best_split[0].covar=i+1;
                    best_split[0].type=3;
                    best_split[0].sign=0;
                    best_split[0].val=nom_lval;
                    best_split[0].n=cands[0].n;
                    best_split[0].n0=cands[0].n0;
                    best_split[0].n1=cands[0].n1;
                    best_split[0].r0=cands[0].r0;
                    best_split[0].r1=cands[0].r1;
                    best_split[0].p0=cands[0].p0;
                    best_split[0].p1=cands[0].p1;
                    best_split[0].predict=cands[0].p;

                    best_split[1].covar=i+1;
                    best_split[1].type=3;
                    best_split[1].sign=0;
                    best_split[1].val=nom_rval;
                    best_split[1].n=cands[1].n;
                    best_split[1].n0=cands[1].n0;
                    best_split[1].n1=cands[1].n1;
                    best_split[1].r0=cands[1].r0;
                    best_split[1].r1=cands[1].r1;
                    best_split[1].p0=cands[1].p0;
                    best_split[1].p1=cands[1].p1;
                    best_split[1].predict=cands[1].p;
                }

//           for nominal variables with more than 2 categories
            } else {

                int z,k;
                int current_var=i+1;
                int k_index,cat_index;
                int max_k=floor(nunique/2);
                unsigned long long int counter=0;

                for ( z=0; z<max_k; z++ ) {

                    k=z+1;        // size of 1 combination
                    k_index=0;    // current index in a single comb.
                    cat_index=0;  // current index in list of cats.

                    double temp[k]; // array to temporarily store 
                                    // a combination

                    max=get_comb_vals_surv(current_node,current_var,
                                           n,nc,y,data,types,ncat,
                                           treat,censor,nmin2,method,
                                           temp,uniquevals,nunique,k,
                                           ncat[i],n_splits,&counter,
                                           cat_index,k_index,max,
                                           temp_tree,cands,
                                           best_split);
                }
            }
        }
    }

    free(cands);
    free(ifcovar);
    free(col);
    free(col_node);
    free(unique);

    if ( max <= 0 ) return(-7);
    return(0);
}



int
get_split2_mtry_surv(int current_node,
                     int n,
                     int nc,
                     double *y,
                     double **data,
                     int *types,
                     int *ncat,
                     int *treat,
                     int *censor,
                     struct node *tree,
                     int nmin2,
                     int mtry,
                     int method,
                     struct node *best_split)

//  
//  This function accepts as arguments the integer
//  index of the current node "current_node", tree
//  information "tree", the minimum node size 
//  "nmin2", and the number of variables to randomly
//  select "mtry" and returns the best split as node 
//  struct array of size 2 "best_split".
//
//  If there is no best split, then the function 
//  returns integer "-7".
//
//  All possible splits of the mtry selected 
//  variables are considered, and split scores use 
//  Cox model statistics.
// 

{
//   select mtry variables
    int mtry_vars[mtry];
    int mtry_types[mtry];

    int i;
    for ( i=0; i<mtry; i++ ) {
        mtry_vars[i]=-7;  // initialize mtry candidate vars
        mtry_types[i]=-7; // initialize mtry candidate types
    }

    int mtry00=mtry;
    if ( mtry > nc ) mtry00=nc;

    int mtry_index[mtry00];
    sample_without_replace(nc,mtry00,mtry_index);

    for ( i=0; i<mtry00; i++ ) {
        mtry_vars[i]=mtry_index[i]+1;
        mtry_types[i]=types[mtry_index[i]];
    }

/*
//   check
    Rprintf("\nRandomly selected mtry=%d variables at node %d:\n",
            mtry,current_node);
    for ( i=0; i<mtry; i++ ) {
        Rprintf("i=%d\tvar=%d\ttype=%d\n",i,mtry_vars[i],
                mtry_types[i]);
    }
*/

//   initialize arrays in temporary tree
    struct node temp_tree[current_node+2];
    for ( i=0; i<current_node; i++ ) {

        temp_tree[i].index=tree[i].index;
        temp_tree[i].covar=tree[i].covar;
        temp_tree[i].type=tree[i].type;
        temp_tree[i].sign=tree[i].sign;
        temp_tree[i].val=tree[i].val;
        temp_tree[i].parent=tree[i].parent;
        temp_tree[i].depth=tree[i].depth;
        temp_tree[i].lchild=tree[i].lchild;
        temp_tree[i].rchild=tree[i].rchild;
        temp_tree[i].n=tree[i].n;
        temp_tree[i].n0=tree[i].n0;
        temp_tree[i].n1=tree[i].n1;
    }

    temp_tree[current_node-1].lchild=current_node+1;
    temp_tree[current_node].index=current_node+1;
    temp_tree[current_node].parent=current_node;
    temp_tree[current_node].depth=tree[current_node-1].depth+1;
    temp_tree[current_node].lchild=0;
    temp_tree[current_node].rchild=0;
    temp_tree[current_node].n=0;

    temp_tree[current_node-1].rchild=current_node+2;
    temp_tree[current_node+1].index=current_node+2;
    temp_tree[current_node+1].parent=current_node;
    temp_tree[current_node+1].depth=tree[current_node-1].depth+1;
    temp_tree[current_node+1].lchild=0;
    temp_tree[current_node+1].rchild=0;
    temp_tree[current_node+1].n=0;

/*
//   check
    Rprintf("temp_tree=\n");
    for ( i=0; i<(current_node+2); i++ ) {

      Rprintf("index=%d  lchild=%d  rchild=%d  parent=%d  depth=%d\n",
	      temp_tree[i].index,temp_tree[i].lchild,temp_tree[i].rchild,
              temp_tree[i].parent,temp_tree[i].depth);
    }
*/

//   find row indices of current node
    int *ifcovar=(int *)malloc(n*sizeof(int));
    findrows_node(current_node,n,nc,data,ncat,tree,ifcovar);

//   assess candidate splits by variable type
    int j;
    int col_counter=0;
    int nsubj_node=tree[current_node-1].n;
    int nunique;
    double *col=(double *)malloc(n*sizeof(double));
    double *col_node=(double *)malloc(nsubj_node*sizeof(double));
    double *unique=(double *)malloc(nsubj_node*sizeof(double));
    int lll=current_node+1;
    struct cand *cands=malloc(2*sizeof(struct cand));
    double max=0;
    double score;

    for ( i=0; i<mtry; i++ ) {

//       type may equal -7 if mtry > nc
        if ( mtry_types[i] == -7 ) break;

//       get data column of current node for ith mtry variable
        get_col(data,mtry_vars[i],n,nc,col);

        col_counter=0;

        for ( j=0; j<n; j++ ) {

            if ( ifcovar[j] == 1 ) {
                col_node[col_counter]=col[j];
                col_counter++;
            }
        }

//       get and sort unique values of the data column
        nunique=get_nunique(nsubj_node,col_node,unique);

        if ( nunique == 1 ) continue; // if all vals same, no split

        double uniquevals[nunique];
        for ( j=0; j<nunique; j++ ) {
            uniquevals[j]=unique[j];
        }

        qsort(uniquevals,nunique,sizeof(double),compare_doubles);

//       set the current covariate in the temporary tree
        temp_tree[current_node].covar=mtry_vars[i];
        temp_tree[current_node+1].covar=mtry_vars[i];


        if ( mtry_types[i] == 1 ) {  // binary

//           get values for left child candidate node
            temp_tree[current_node].type=1;
            temp_tree[current_node].sign=0;
            temp_tree[current_node].val=0;

            get_split0_cox(0,lll,n,nc,y,data,types,ncat,treat,
                           censor,method,temp_tree,cands);

//           get values for right child candidate node
            temp_tree[current_node+1].type=1;
            temp_tree[current_node+1].sign=0;
            temp_tree[current_node+1].val=1;

            get_split0_cox(1,lll+1,n,nc,y,data,types,ncat,treat,
                           censor,method,temp_tree,cands);

//           candidate splits with too few subjects are not useful
            if ( (cands[0].n0 <= nmin2) || (cands[0].n1 <= nmin2) || 
                 (cands[1].n0 <= nmin2) || (cands[1].n1 <= nmin2) ) {

                continue;
            }

            score=get_coxph_stat(current_node,n,nc,y,data,types,
                                 ncat,treat,censor,method,temp_tree);

            if ( score > max ) {

                max=score;
                best_split[0].G=score;

                best_split[0].covar=mtry_vars[i];
                best_split[0].type=1;
                best_split[0].sign=0;
                best_split[0].val=0;
                best_split[0].n=cands[0].n;
                best_split[0].n0=cands[0].n0;
                best_split[0].n1=cands[0].n1;
                best_split[0].r0=cands[0].r0;
                best_split[0].r1=cands[0].r1;
                best_split[0].p0=cands[0].p0;
                best_split[0].p1=cands[0].p1;
                best_split[0].predict=cands[0].p;

                best_split[1].covar=mtry_vars[i];
                best_split[1].type=1;
                best_split[1].sign=0;
                best_split[1].val=1;
                best_split[1].n=cands[1].n;
                best_split[1].n0=cands[1].n0;
                best_split[1].n1=cands[1].n1;
                best_split[1].r0=cands[1].r0;
                best_split[1].r1=cands[1].r1;
                best_split[1].p0=cands[1].p0;
                best_split[1].p1=cands[1].p1;
                best_split[1].predict=cands[1].p;
            }

        } else if ( mtry_types[i] == 2 ) {  // ordinal

            temp_tree[current_node].type=2;
            temp_tree[current_node+1].type=2;

            for ( j=0; j<(nunique-1); j++ ) {

//               get values for left child candidate node
                temp_tree[current_node].val=uniquevals[j];
                temp_tree[current_node].sign=1; // LE

                get_split0_cox(0,lll,n,nc,y,data,types,ncat,treat,
                               censor,method,temp_tree,cands);

//               get values for right child candidate node
                temp_tree[current_node+1].val=uniquevals[j];
                temp_tree[current_node+1].sign=2; // GT

                get_split0_cox(1,lll+1,n,nc,y,data,types,ncat,treat,
                               censor,method,temp_tree,cands);

//               candidate splits with too few subjects are not useful
                if ( (cands[0].n0 <= nmin2) || 
                     (cands[0].n1 <= nmin2) || 
                     (cands[1].n0 <= nmin2) || 
                     (cands[1].n1 <= nmin2) ) {

                    continue;
                }

                score=get_coxph_stat(current_node,n,nc,y,data,types,
                                     ncat,treat,censor,method,
                                     temp_tree);

                if ( score > max ) {

                    max=score;
                    best_split[0].G=score;

                    best_split[0].covar=mtry_vars[i];
                    best_split[0].type=2;
                    best_split[0].sign=1; // LE
                    best_split[0].val=uniquevals[j];
                    best_split[0].n=cands[0].n;
                    best_split[0].n0=cands[0].n0;
                    best_split[0].n1=cands[0].n1;
                    best_split[0].r0=cands[0].r0;
                    best_split[0].r1=cands[0].r1;
                    best_split[0].p0=cands[0].p0;
                    best_split[0].p1=cands[0].p1;
                    best_split[0].predict=cands[0].p;

                    best_split[1].covar=mtry_vars[i];
                    best_split[1].type=2;
                    best_split[1].sign=2; // GT
                    best_split[1].val=uniquevals[j];
                    best_split[1].n=cands[1].n;
                    best_split[1].n0=cands[1].n0;
                    best_split[1].n1=cands[1].n1;
                    best_split[1].r0=cands[1].r0;
                    best_split[1].r1=cands[1].r1;
                    best_split[1].p0=cands[1].p0;
                    best_split[1].p1=cands[1].p1;
                    best_split[1].predict=cands[1].p;
                }
            }

        } else if ( mtry_types[i] == 3 ) {  // nominal

            temp_tree[current_node].type=3;
            temp_tree[current_node].sign=0;

            temp_tree[current_node+1].type=3;
            temp_tree[current_node+1].sign=0;

//           calculate total number of possible splits
            unsigned long long int n_splits=pow(2,nunique-1)-1;
////	    Rprintf("\nnunique=%d n_splits=%llu\n",nunique,n_splits);

//           2 categories reduce to binary split
            if ( nunique == 2 ) {  

                double nom_lval=pow(2,ncat[mtry_vars[i]-1]-
                                    uniquevals[0]);

                double nom_rval=pow(2,ncat[mtry_vars[i]-1]-
                                    uniquevals[1]);

//               get values for left child candidate node
                temp_tree[current_node].val=nom_lval;
                get_split0_cox(0,lll,n,nc,y,data,types,ncat,treat,
                               censor,method,temp_tree,cands);

//               get values for right child candidate node
                temp_tree[current_node+1].val=nom_rval;
                get_split0_cox(1,lll+1,n,nc,y,data,types,ncat,treat,
                               censor,method,temp_tree,cands);

//               candidate splits with too few subjects are not useful
                if ( (cands[0].n0 <= nmin2) || 
                     (cands[0].n1 <= nmin2) || 
                     (cands[1].n0 <= nmin2) ||
                     (cands[1].n1 <= nmin2) ) {

                    continue;
                }

                score=get_coxph_stat(current_node,n,nc,y,data,types,
                                     ncat,treat,censor,method,
                                     temp_tree);

                if ( score > max ) {

                    max=score;
                    best_split[0].G=score;

                    best_split[0].covar=mtry_vars[i];
                    best_split[0].type=3;
                    best_split[0].sign=0;
                    best_split[0].val=nom_lval;
                    best_split[0].n=cands[0].n;
                    best_split[0].n0=cands[0].n0;
                    best_split[0].n1=cands[0].n1;
                    best_split[0].r0=cands[0].r0;
                    best_split[0].r1=cands[0].r1;
                    best_split[0].p0=cands[0].p0;
                    best_split[0].p1=cands[0].p1;
                    best_split[0].predict=cands[0].p;

                    best_split[1].covar=mtry_vars[i];
                    best_split[1].type=3;
                    best_split[1].sign=0;
                    best_split[1].val=nom_rval;
                    best_split[1].n=cands[1].n;
                    best_split[1].n0=cands[1].n0;
                    best_split[1].n1=cands[1].n1;
                    best_split[1].r0=cands[1].r0;
                    best_split[1].r1=cands[1].r1;
                    best_split[1].p0=cands[1].p0;
                    best_split[1].p1=cands[1].p1;
                    best_split[1].predict=cands[1].p;
                }

//           for nominal variables with more than 2 categories
            } else {

                int z,k;
                int current_var=mtry_vars[i];
                int k_index,cat_index;
                int max_k=floor(nunique/2);
                unsigned long long int counter=0;

                for ( z=0; z<max_k; z++ ) {

                    k=z+1;        // size of 1 combination
                    k_index=0;    // current index in a single comb.
                    cat_index=0;  // current index in list of cats.

                    double temp[k]; // array to temporarily store 
                                    // a combination

                    max=get_comb_vals_surv(current_node,current_var,
                                           n,nc,y,data,types,ncat,
                                           treat,censor,nmin2,method,
                                           temp,uniquevals,nunique,k,
                                           ncat[mtry_vars[i]-1],
                                           n_splits,&counter,
                                           cat_index,k_index,max,
                                           temp_tree,cands,
                                           best_split);
                }
            }
        }
    }

    free(cands);
    free(ifcovar);
    free(col);
    free(col_node);
    free(unique);
/*
    Rprintf("lchild_nsubj=%d rchild_nsubj=%d\n",best_split[0].n,best_split[1].n);
*/
    if ( max <= 0 ) return(-7);
    return(0);
}



void
maketree_surv(int n,
              int nc,
              double *y,
              double **data,
              int *types,
              int *ncat,
              int *treat,
              int *censor,
              struct node *tree,
              int nmin,
              int maxdepth,
              int method)

//
//  This function returns a tree using the Cox PH
//  model to split each node.
//

{
//   initialize root node
    double r[4];
    tree[0].index=1;
    tree[0].depth=1;
    tree[0].n=n;
    tree[0].n0=get_subj_treat(1,0,n,nc,data,ncat,treat,tree);
    tree[0].n1=get_subj_treat(1,1,n,nc,data,ncat,treat,tree);
    tree[0].predict=get_node_predict_y(1,n,nc,y,data,ncat,treat,r,
                                       tree);
    tree[0].r0=r[0];
    tree[0].r1=r[1];
    tree[0].p0=r[2];
    tree[0].p1=r[3];

//   build tree
    int i;
    int current_node,isplit,max,d;
    struct node *best_split=malloc(2*sizeof(struct node));

    for ( i=0; i<MAXSIZE; i++ ) {

        current_node=i+1;

        d=tree[current_node-1].depth+1;
        if ( (maxdepth > 0) && (d > maxdepth) ) break;

        isplit=get_split2_cox_all(current_node,n,nc,y,data,types,
                                  ncat,treat,censor,tree,nmin,method,
                                  best_split);

        if ( isplit != -7 ) {  // "-7" = there should be no split

            max=get_number_of_nodes(MAXSIZE,tree);
            tree[i].lchild=max+1;
            tree[i].rchild=max+2;

//           add left child node to tree
            tree[max].index=tree[i].lchild;
            tree[max].covar=best_split[0].covar;
            tree[max].type=best_split[0].type;
            tree[max].sign=best_split[0].sign;
            tree[max].val=best_split[0].val;
            tree[max].parent=current_node;
            tree[max].depth=d;
            tree[max].n=best_split[0].n;
            tree[max].n0=best_split[0].n0;
            tree[max].n1=best_split[0].n1;
            tree[max].r0=best_split[0].r0;
            tree[max].r1=best_split[0].r1;
            tree[max].p0=best_split[0].p0;
            tree[max].p1=best_split[0].p1;
            tree[max].predict=best_split[0].predict;

//           add right child node to tree
            tree[max+1].index=tree[i].rchild;
            tree[max+1].covar=best_split[1].covar;
            tree[max+1].type=best_split[1].type;
            tree[max+1].sign=best_split[1].sign;
            tree[max+1].val=best_split[1].val;
            tree[max+1].parent=current_node;
            tree[max+1].depth=d;
            tree[max+1].n=best_split[1].n;
            tree[max+1].n0=best_split[1].n0;
            tree[max+1].n1=best_split[1].n1;
            tree[max+1].r0=best_split[1].r0;
            tree[max+1].r1=best_split[1].r1;
            tree[max+1].p0=best_split[1].p0;
            tree[max+1].p1=best_split[1].p1;
            tree[max+1].predict=best_split[1].predict;
        }

        max=get_number_of_nodes(MAXSIZE,tree);
        if ( (isplit == -7) && (current_node == max) ) break;
    }

    free(best_split);

/*
//   check
    get_coxph_stat(1,n,nc,y,data,types,ncat,treat,censor,method,tree);
*/
}



void
get_OOB(n,boot_n,array,OOB)

int n;       // total sample size
int boot_n;  // size of array
int array[];
int OOB[];

//
//  This function accepts as an argument an array of
//  boot_n integer indices whose possible values are in 
//  0:(n-1). "array" contains the indices of subjects
//  who are in the bootstrap sample.
//
//  This function returns array of indicators "OOB".
//  "1" means an index is NOT in the bootstrap sample
//  and "0" means an index is in the bootstrap sample.
//
//  In other words, 1=OOB, 0=in-bag.
//

{
    int i;
    for ( i=0; i<n; i++ ) {
        OOB[i]=1;  // initialize
    }

    int j;
    for ( i=0; i<boot_n; i++ ) {
        for ( j=0; j<n; j++ ) {

            if ( j == array[i] ) {

                OOB[j]=0; // 0 means subject in bootstrap sample
                break;
            }
        }
    }

/*
//   check
    Rprintf("\nIn-bag indices (n=%d, boot_n=%d):\n",n,boot_n);
    for ( i=0; i<boot_n; i++ ) Rprintf("%d  ",array[i]);
    Rprintf("\nOOB indicator:\n");
    for ( i=0; i<n; i++ ) Rprintf("%d  ",OOB[i]);
*/
}



void
draw_boot(int n,
          int *index)
 
//
//  This function accepts as an argument the total
//  sample size "n" and returns in array of integers 
//  "index" the indices from 0:(n-1) that are in a 
//  bootstrap sample of size "n" sampled with 
//  replacement.
//

{
//   get random row numbers sampled from 0:(n-1)
    int i;
    for ( i=0; i<n; i++ ) {
        index[i]=get_rand_int(n);
    }

/*
//   get array of OOB indicators (1=row index is in OOB sample)
    get_OOB(n,n,index,OOB);
*/
}



int
invert_matrix4(m,inv)

double m[];
double inv[];

//
//  This function accepts values from a 4 by 4 
//  matrix (where the values are stored in an 
//  array by row) and returns the inverse of the
//  matrix.
//
//  Note that this function is the "gluInvertMatrix"
//  function from the "OpenGL Utility Library".
//

{
    double inv0[16];

    inv0[0]= m[5]  * m[10] * m[15] - 
             m[5]  * m[11] * m[14] - 
             m[9]  * m[6]  * m[15] + 
             m[9]  * m[7]  * m[14] +
             m[13] * m[6]  * m[11] - 
             m[13] * m[7]  * m[10];

    inv0[4]= -m[4]  * m[10] * m[15] + 
              m[4]  * m[11] * m[14] + 
              m[8]  * m[6]  * m[15] - 
              m[8]  * m[7]  * m[14] - 
              m[12] * m[6]  * m[11] + 
              m[12] * m[7]  * m[10];

    inv0[8]= m[4]  * m[9] * m[15] - 
             m[4]  * m[11] * m[13] - 
             m[8]  * m[5] * m[15] + 
             m[8]  * m[7] * m[13] + 
             m[12] * m[5] * m[11] - 
             m[12] * m[7] * m[9];

    inv0[12]= -m[4]  * m[9] * m[14] + 
               m[4]  * m[10] * m[13] +
               m[8]  * m[5] * m[14] - 
               m[8]  * m[6] * m[13] - 
               m[12] * m[5] * m[10] + 
               m[12] * m[6] * m[9];

    inv0[1]= -m[1]  * m[10] * m[15] + 
              m[1]  * m[11] * m[14] + 
              m[9]  * m[2] * m[15] - 
              m[9]  * m[3] * m[14] - 
              m[13] * m[2] * m[11] + 
              m[13] * m[3] * m[10];

    inv0[5]= m[0]  * m[10] * m[15] - 
             m[0]  * m[11] * m[14] - 
             m[8]  * m[2] * m[15] + 
             m[8]  * m[3] * m[14] + 
             m[12] * m[2] * m[11] - 
             m[12] * m[3] * m[10];

    inv0[9]= -m[0]  * m[9] * m[15] + 
              m[0]  * m[11] * m[13] + 
              m[8]  * m[1] * m[15] - 
              m[8]  * m[3] * m[13] - 
              m[12] * m[1] * m[11] + 
              m[12] * m[3] * m[9];

    inv0[13]= m[0]  * m[9] * m[14] - 
              m[0]  * m[10] * m[13] - 
              m[8]  * m[1] * m[14] + 
              m[8]  * m[2] * m[13] + 
              m[12] * m[1] * m[10] - 
              m[12] * m[2] * m[9];

    inv0[2]= m[1]  * m[6] * m[15] - 
             m[1]  * m[7] * m[14] - 
             m[5]  * m[2] * m[15] + 
             m[5]  * m[3] * m[14] + 
             m[13] * m[2] * m[7] - 
             m[13] * m[3] * m[6];

    inv0[6]= -m[0]  * m[6] * m[15] + 
              m[0]  * m[7] * m[14] + 
              m[4]  * m[2] * m[15] - 
              m[4]  * m[3] * m[14] - 
              m[12] * m[2] * m[7] + 
              m[12] * m[3] * m[6];

    inv0[10]= m[0]  * m[5] * m[15] - 
              m[0]  * m[7] * m[13] - 
              m[4]  * m[1] * m[15] + 
              m[4]  * m[3] * m[13] + 
              m[12] * m[1] * m[7] - 
              m[12] * m[3] * m[5];

    inv0[14]= -m[0]  * m[5] * m[14] + 
               m[0]  * m[6] * m[13] + 
               m[4]  * m[1] * m[14] - 
               m[4]  * m[2] * m[13] - 
               m[12] * m[1] * m[6] + 
               m[12] * m[2] * m[5];

    inv0[3]= -m[1] * m[6] * m[11] + 
              m[1] * m[7] * m[10] + 
              m[5] * m[2] * m[11] - 
              m[5] * m[3] * m[10] - 
              m[9] * m[2] * m[7] + 
              m[9] * m[3] * m[6];

    inv0[7]= m[0] * m[6] * m[11] - 
             m[0] * m[7] * m[10] - 
             m[4] * m[2] * m[11] + 
             m[4] * m[3] * m[10] + 
             m[8] * m[2] * m[7] - 
             m[8] * m[3] * m[6];

    inv0[11]= -m[0] * m[5] * m[11] + 
               m[0] * m[7] * m[9] + 
               m[4] * m[1] * m[11] - 
               m[4] * m[3] * m[9] - 
               m[8] * m[1] * m[7] + 
               m[8] * m[3] * m[5];

    inv0[15]= m[0] * m[5] * m[10] - 
              m[0] * m[6] * m[9] - 
              m[4] * m[1] * m[10] + 
              m[4] * m[2] * m[9] + 
              m[8] * m[1] * m[6] - 
              m[8] * m[2] * m[5];

    double det;
    det=(m[0]*inv0[0])+(m[1]*inv0[4])+(m[2]*inv0[8])+(m[3]*inv0[12]);
////    Rprintf("\ndet=%f\n",det);

    if ( det == 0 ) return(0);

    det=1.0/det;

    int i;
    for ( i=0; i<16; i++ ) {
        inv[i]=inv0[i]*det;
    }

    return(1);
}



double
get_G(int current_node,
      int n,
      int nc,
      double *y,
      double **data,
      int *types,
      int *ncat,
      int *treat,
      struct node *tree)

//
//  This function calculates the "G" term in the depth
//  variable importance measure developed by Zhang. 
//
//  Since the outcome variable Y is continuous in this
//  setting, the "G" term is the test statistic
//  testing the significance of split by treatment
//  linear regression model with an interaction term.
//

{
//   get y vals, treatments, and left split indicator of current node
//   in embedded tree
    int *ifcovar=(int *)malloc(n*sizeof(int));
    findrows_node(current_node,n,nc,data,ncat,tree,ifcovar);

    int *ifleft0=(int *)malloc(n*sizeof(int));
    findrows_node(tree[current_node-1].lchild,n,nc,data,ncat,tree,
                  ifleft0);

    int i;
    int nsubj=tree[current_node-1].n;
    double yvals[nsubj];
    int trtvals[nsubj];
    int ifleft[nsubj];
    int counter=0;

    for ( i=0; i<n; i++ ) {

        if ( ifcovar[i] == 1 ) {

            yvals[counter]=y[i];
            trtvals[counter]=treat[i];
            ifleft[counter]=ifleft0[i];

            counter++;
        }
    }

    free(ifcovar);
    free(ifleft0);

/*
//   check
    Rprintf("\ncurrent node=%d\n",current_node);
    Rprintf("current lchild=%d\n",tree[current_node-1].lchild);
    Rprintf("subjects in current node=%d\n",nsubj);

    Rprintf("\nNode values:\n");
    for ( i=0; i<nsubj; i++ ) {
        Rprintf("i=%d\tyval=%f\ttrt=%d\tleft=%d\n",
               i,yvals[i],trtvals[i],ifleft[i]);
    }

    Rprintf("\nX matrix:\n");
    for ( i=0; i<nsubj; i++ ) {
        Rprintf("1,%d,%d,%d,",trtvals[i],ifleft[i],
                trtvals[i]*ifleft[i]);
    }

    Rprintf("\ny vector: ");
    for ( i=0; i<nsubj; i++ ) Rprintf("%f,",yvals[i]);
*/

//   create X matrix for regression model
    double X[nsubj][4];
    for ( i=0; i<nsubj; i++ ) {
        X[i][0]=1;
        X[i][1]=trtvals[i];
        X[i][2]=ifleft[i];
        X[i][3]=trtvals[i]*ifleft[i];
    }

//   get transpose of X
    int j;
    double XT[4][nsubj];
    for ( i=0; i<nsubj; i++ ) {
        for ( j=0; j<4; j++ ) {
            XT[j][i]=X[i][j];
        }
    }

//   multiply transpose(X) by X
    int k;
    double XTX[4][4];

    for ( i=0; i<4; i++ ) {
        for ( j=0; j<4; j++ ) {
            XTX[i][j]=0; // initialize
        }
    }

    for ( i=0; i<4; i++ ) {
        for ( j=0; j<4; j++ ) {
            for ( k=0; k<nsubj; k++ ) {

                XTX[i][j] += XT[i][k]*X[k][j];
            }
        }
    }

//   get inverse of XTX
    double XTX0[16];
    double inv0[16];
    counter=0;

    for ( i=0; i<4; i++ ) {
        for ( j=0; j<4; j++ ) {
            XTX0[counter]=XTX[i][j];  // convert matrix to array
            counter++;
        }
    }

    int invertible=invert_matrix4(XTX0,inv0);

    if ( invertible == 0 ) return(0); // XTX not invertible

    double inv[4][4];
    counter=0;

    for ( i=0; i<4; i++ ) {
        for ( j=0; j<4; j++ ) {
            inv[i][j]=inv0[counter];  // convert array to matrix
            counter++;
        }
    }

//   multiply tranpose(X) by Y
    double XTY[4]={0,0,0,0};
    for ( i=0; i<4; i++ ) {
        for ( j=0; j<nsubj; j++ ) {

            XTY[i] += XT[i][j]*yvals[j];
        }
    }

//   calculate beta coefficients of regression model
    double beta[4]={0,0,0,0};
    for ( i=0; i<4; i++ ) {
        for ( j=0; j<4; j++ ) {

            beta[i] += inv[i][j]*XTY[j];
        }
    }

/*
//   check
    Rprintf("Betas:\n");
    for ( i=0; i<4; i++ ) Rprintf("%f  ",beta[i]);
*/

//   when beta1 and beta3 have same sign, both child nodes have
//   the same best treatment
    int signb1=( beta[1] > 0 );
    int signb3=( beta[3] > 0 );
    if ( signb1 == signb3 ) return(0);

//   calculate predicted y values from the model
    double yhats[nsubj];
    for ( i=0; i<nsubj; i++ ) {
        yhats[i]=beta[0] + (beta[1]*X[i][1]) + (beta[2]*X[i][2]) + 
                 (beta[3]*X[i][3]);
    }

//   calculate mean of yvals
    double ybar=get_mean(yvals,nsubj);

//   calculate F statistic
    double SSM=0; // model sum of squares
    double SSE=0; // error sum of squares

    for ( i=0; i<nsubj; i++ ) {
        SSM += ((yhats[i]-ybar)*(yhats[i]-ybar));
        SSE += ((yvals[i]-yhats[i])*(yvals[i]-yhats[i]));
    }

    double F;
    F=(SSM/3) / (SSE/(nsubj-4));

    return(F);
}



double
get_s2(y,ybar,n)

double y[];
double ybar;
int n;

//
//  This function accepts as arguments an array of
//  n doubles "y" as well as its mean "ybar" and
//  returns the sample variance of "y".
//

{
/*
    if ( n == 1 ) return(0);
*/
    int i;
    double s2=0;

    for ( i=0; i<n; i++ ) {
        s2 += ( (y[i]-ybar)*(y[i]-ybar) );
    }

    s2 /= (n-1);

/*
//   check
    Rprintf("\n\ny vector (s2=%f, ybar=%f):\n",s2,ybar);
    for ( i=0; i<n; i++ ) Rprintf("%f,",y[i]);
*/

    return(s2);
}



double
get_G2(int current_node,
       int n,
       int nc,
       double *y,
       double **data,
       int *types,
       int *ncat,
       int *treat,
       struct node *tree)

//
//  This function calculates the "G" term in the depth
//  variable importance measure developed by Zhang. 
//
//  Since the outcome variable Y is continuous in this
//  setting, the "G" term is the test statistic
//  testing the significance of the interaction term
//  in the split by treatment linear regression model.
//

{
//   get sizes of child nodes by treatment group
    int lchild=tree[current_node-1].lchild;
    int rchild=tree[current_node-1].rchild;
    int nodeval[4]={lchild,lchild,rchild,rchild};
    int trtval[4]={1,0,1,0};
    int size[4]={tree[lchild-1].n1,tree[lchild-1].n0,
                 tree[rchild-1].n1,tree[rchild-1].n0};

    int i;
    double w0=0; 

    for ( i=0; i<4; i++ ) {

        if ( (size[i] == 0) || (size[i] == 1) ) {
            return(0);
        }

        w0 += (size[i]-1);  // used to calculate sigma2 below
    }

//   calculate various values by child node and treatment group
    int j;
    int *ifcovar=(int *)malloc(n*sizeof(int));
    int counter;
    double ybar[4];
    double sigma2=0;
    double Gden=0;
    double s2;
    double w;

    for ( i=0; i<4; i++ ) {

//       get y values of current child node and treatment group
        findrows_node(nodeval[i],n,nc,data,ncat,tree,ifcovar);

        double yvals[size[i]];
        counter=0;

        for ( j=0; j<n; j++ ) {

            if ( (ifcovar[j] == 1) && (treat[j] == trtval[i]) ) {
                yvals[counter]=y[j];
                counter++;
            }
        }

//       calculate various values
        ybar[i]=get_mean(yvals,size[i]);
        s2=get_s2(yvals,ybar[i],size[i]);

        w=(size[i] - 1)/w0;
        sigma2 += w*s2;

        Gden += (1/(size[i]+0.0)); // coerce size to double
    }

    free(ifcovar);

//   calculate G
    double Gnum,G;
    Gnum=(ybar[0] - ybar[1]) - (ybar[2] - ybar[3]);
    G=(Gnum*Gnum) / (sigma2*Gden);

    return(G);
}



double
get_G_mc(int current_node,
         int n,
         int nc,
         double *y,
         double **data,
         int *types,
         int *ncat,
         int *treat,
         struct node *tree)

//
//  This function calculates the "G" term in the depth
//  variable importance measure developed by Zhang for
//  data with a continuous outcome variable Y and 
//  multiple treatments.
//

{
//   get y vals, treatments, and left split indicator of current node
//   in embedded tree
    int *ifcovar=(int *)malloc(n*sizeof(int));
    findrows_node(current_node,n,nc,data,ncat,tree,ifcovar);

    int *ifleft0=(int *)malloc(n*sizeof(int));
    findrows_node(tree[current_node-1].lchild,n,nc,data,ncat,tree,
                  ifleft0);

    int i;
    int nsubj=tree[current_node-1].n;
    double *yvals=(double *)malloc(nsubj*sizeof(double));
    int *trtvals=(int *)malloc(nsubj*sizeof(int));
    int *ifleft=(int *)malloc(nsubj*sizeof(int));
    int counter=0;

    for ( i=0; i<n; i++ ) {

        if ( ifcovar[i] == 1 ) {

            yvals[counter]=y[i];
            trtvals[counter]=treat[i];
            ifleft[counter]=ifleft0[i];

            counter++;
        }
    }

    free(ifcovar);
    free(ifleft0);

//   create data matrix for multiple regression model fitting
    double **lm_data=(double **)malloc(nsubj*sizeof(double *));
    for ( i=0; i<nsubj; i++ ) {
        lm_data[i]=(double *)calloc(3,sizeof(double));
        lm_data[i][0]=yvals[i];
        lm_data[i][1]=trtvals[i];
        lm_data[i][2]=ifleft[i];
    }

    free(yvals);
    free(trtvals);
    free(ifleft);

//   get regression model statistic (z^2)
    double stat=0;
    stat=get_lmstat_mc(nsubj,lm_data);


    for ( i=0; i<nsubj; i++ ) free(lm_data[i]);
    free(lm_data);

    return(stat);
}



int
get_child_nodes(int node,
                int lll,
                struct node *tree,
                int child_nodes[])

//
//  This function accepts as arguments the current
//  node of interest as integer "node" (nodes are
//  counted starting at 1), the total number of
//  nodes minus 1 (for the root node) "lll", and 
//  left and right child nodes for each node in a 
//  tree. 
// 
//  This function returns the total number of 
//  descendant nodes of the current node of 
//  interest and saves the indices of the 
//  descendant nodes in the array of integers 
//  "child_nodes".
//

{
    int i;
    for ( i=0; i<lll; i++ ) child_nodes[i]=0;  // initialize

    int counter=0;
    int current_node=node;

    for ( i=0; i<lll; i++ ) {

        if ( current_node == 0 ) break;

        if ( tree[current_node-1].lchild == 0 ) { 
            current_node=child_nodes[i];
            continue;
        }

        child_nodes[counter]=tree[current_node-1].lchild;
        child_nodes[counter+1]=tree[current_node-1].rchild;

        counter += 2;
        current_node=child_nodes[i];
    }

/*
//   check
    Rprintf("\n%d descendants of node %d (lll=%d): ",counter,node,
            lll);
    for ( i=0; i<lll; i++ ) Rprintf("%d  ",child_nodes[i]);
*/

    return(counter);
}



void
get_varimp_replace(double *varimp,
                   int n,
                   int nc,
                   double *y,
                   double **data,
                   int *types,
                   int *ncat,
                   int *treat,
                   int *censor,
                   struct node *tree,
                   int method,
                   int *index)

//
//  This function accepts as arguments data and
//  information from a classification tree and 
//  returns the variable importance measures of 
//  each variable in array of doubles "varimp"
//  using the depth variable importance measure.
//

{
//   lll is the total number of nodes in the tree
    int i;
    int lll=get_number_of_nodes(MAXSIZE,tree);

    if ( lll == 1 ) return;

//   use bootstrap data
    double *boot_y=(double *)malloc(n*sizeof(double));
    int *boot_treat=(int *)malloc(n*sizeof(int));
    int *boot_censor=(int *)malloc(n*sizeof(int));

    double **boot=(double **)malloc(n*sizeof(double *));
    for ( i=0; i<n; i++ ) {
        boot[i]=(double *)malloc(nc*sizeof(double));
    }

    int j;
    for ( i=0; i<n; i++ ) {

        boot_y[i]=y[index[i]];
        boot_treat[i]=treat[index[i]];
        boot_censor[i]=censor[index[i]];

        for ( j=0; j<nc; j++ ) {
            boot[i][j]=data[index[i]][j];
        }
    }

//
//   ... calculate depth term and G term
//
    int current_node;
    int current_lchild;
    int splitvar[lll];
    double depthcoeff[lll];
    double G[lll];

    for ( i=0; i<lll; i++ ) {

        current_node=tree[i].index;
        current_lchild=tree[i].lchild;

//       for terminal nodes, there's no split var, and G=0
        if ( current_lchild == 0 ) {
            G[i]=0;
            splitvar[i]=0;
            depthcoeff[i]=0;
            continue;
        }

//       save index of the variable used to split current node
        splitvar[i]=tree[current_lchild-1].covar;

//       save depth coefficient of current node
        depthcoeff[i]=pow(2,-tree[current_node-1].depth);

//       get G statistic
        G[i]=tree[current_node-1].G;
    }

    free(boot_y);
    free(boot_treat);
    free(boot_censor);
    for ( i=0; i<n; i++ ) free(boot[i]);
    free(boot);

/*
//   check
    Rprintf("\nIn get_varimp, G values:\n");
    for ( i=0; i<lll; i++ ) {
        Rprintf("i=%d\tsplitvar=%d\tG=%f\n",i,splitvar[i],G[i]);
    }
*/

//   replace G values with highest value of any of its descendants
//   if the highest value exceeds G at current node
    int nchild;
    int child_nodes[lll-1];

    for ( i=0; i<lll; i++ ) {

        current_node=tree[i].index;
        current_lchild=tree[i].lchild;

//       terminal nodes have no descendants
        if ( current_lchild == 0 ) continue;

//       find all descendants of current node
        nchild=get_child_nodes(current_node,lll-1,tree,child_nodes);

        for ( j=0; j<nchild; j++ ) {
            if( G[child_nodes[j]-1] > G[i] ) G[i]=G[child_nodes[j]-1];
        }
    }

/*
//   check
    Rprintf("\nAfter replacement, in get_varimp, G values:\n");
    for ( i=0; i<lll; i++ ) {
        Rprintf("i=%d\tsplitvar=%d\tG=%f\n",i,splitvar[i],G[i]);
    }
*/

//
//   ... combine the 2 terms 
//
    for ( i=0; i<lll; i++ ) {
        for ( j=0; j<nc; j++ ) {

            if ( splitvar[i] == (j+1) ) {
                varimp[j] += (G[i]*depthcoeff[i]);
            }
        }
    }
}



void
maketree2(int n,
          int nc,
          double *y,
          double **data,
          int *types,
          int *ncat,
          int *treat,
          struct node *tree,
          int nmin2,
          int mtry,
          int maxdepth2,
          int *index,
          int method)

//
//  This function returns a tree using the modified 
//  classification tree method developed by Zhang. 
//  This method uses the "DIFF" score to split each
//  node.
//

{
//   use bootstrap data
    double *boot_y=(double *)malloc(n*sizeof(double));
    int *boot_treat=(int *)malloc(n*sizeof(int));

    int i;
    double **boot=(double **)malloc(n*sizeof(double *));
    for ( i=0; i<n; i++ ) {
        boot[i]=(double *)malloc(nc*sizeof(double));
    }

    int j;
    for ( i=0; i<n; i++ ) {

        boot_y[i]=y[index[i]];
        boot_treat[i]=treat[index[i]];

        for ( j=0; j<nc; j++ ) {
            boot[i][j]=data[index[i]][j];
        }
    }

//   initialize root node
    tree[0].index=1;
    tree[0].depth=1;
    tree[0].n=n;

    if ( method >= 20 ) {  // multi-treat methods

        tree[0].n0=get_min_ntrt(1,n,nc,boot,ncat,boot_treat,tree);
        tree[0].n1=tree[0].n0;

        tree[0].predict=0;
        tree[0].r0=0;
        tree[0].r1=0;
        tree[0].p0=0;
        tree[0].p1=0;

    } else {               // binary treat methods

        tree[0].n0=get_subj_treat(1,0,n,nc,boot,ncat,boot_treat,tree);
        tree[0].n1=get_subj_treat(1,1,n,nc,boot,ncat,boot_treat,tree);

        double r[4];
        tree[0].predict=get_node_predict_y(1,n,nc,boot_y,boot,ncat,
                                           boot_treat,r,tree);
        tree[0].r0=r[0];
        tree[0].r1=r[1];
        tree[0].p0=r[2];
        tree[0].p1=r[3];
    }

//   build tree
    int current_node,isplit,max,d;
    struct node *best_split=malloc(2*sizeof(struct node));

    for ( i=0; i<MAXSIZE; i++ ) {

        current_node=i+1;

        d=tree[current_node-1].depth+1;
        if ( (maxdepth2 > 0) && (d > maxdepth2) ) break;

        if ( (method == 7) || 
             (method == 8) ||
             (method == 23) ) { // mtry methods

            isplit=get_split2_mtry(current_node,n,nc,boot_y,boot,
                                   types,ncat,boot_treat,tree,nmin2,
                                   mtry,method,best_split);

        } else { // no mtry methods: 5,6,22

            isplit=get_split2_d_all(current_node,n,nc,boot_y,boot,
                                    types,ncat,boot_treat,tree,nmin2,
                                    method,best_split);
        }

        if ( isplit != -7 ) {  // "-7" = there should be no split

            max=get_number_of_nodes(MAXSIZE,tree);
            tree[i].lchild=max+1;
            tree[i].rchild=max+2;

//           add left child node to tree
            tree[max].index=tree[i].lchild;
            tree[max].covar=best_split[0].covar;
            tree[max].type=best_split[0].type;
            tree[max].sign=best_split[0].sign;
            tree[max].val=best_split[0].val;
            tree[max].parent=current_node;
            tree[max].depth=tree[current_node-1].depth+1;
            tree[max].n=best_split[0].n;
            tree[max].n0=best_split[0].n0;
            tree[max].n1=best_split[0].n1;
            tree[max].r0=best_split[0].r0;
            tree[max].r1=best_split[0].r1;
            tree[max].p0=best_split[0].p0;
            tree[max].p1=best_split[0].p1;
            tree[max].predict=best_split[0].predict;

//           add right child node to tree
            tree[max+1].index=tree[i].rchild;
            tree[max+1].covar=best_split[1].covar;
            tree[max+1].type=best_split[1].type;
            tree[max+1].sign=best_split[1].sign;
            tree[max+1].val=best_split[1].val;
            tree[max+1].parent=current_node;
            tree[max+1].depth=tree[current_node-1].depth+1;
            tree[max+1].n=best_split[1].n;
            tree[max+1].n0=best_split[1].n0;
            tree[max+1].n1=best_split[1].n1;
            tree[max+1].r0=best_split[1].r0;
            tree[max+1].r1=best_split[1].r1;
            tree[max+1].p0=best_split[1].p0;
            tree[max+1].p1=best_split[1].p1;
            tree[max+1].predict=best_split[1].predict;

//           get G statistic for current split
            if ( (method == 6) || (method == 7) ) {

                tree[current_node-1].G=get_G2(current_node,n,nc,
                                              boot_y,boot,types,
                                              ncat,boot_treat,tree);

            } else if ( (method == 5) || (method == 8) ) {

                tree[current_node-1].G=get_G(current_node,n,nc,
                                             boot_y,boot,types,
                                             ncat,boot_treat,tree);

            } else {  // method=22 or 23

                tree[current_node-1].G=get_G_mc(current_node,n,nc,
                                                boot_y,boot,types,
                                                ncat,boot_treat,tree);
            }
        }

        max=get_number_of_nodes(MAXSIZE,tree);
        if ( (isplit == -7) && (current_node == max) ) break;
    }

    free(boot_y);
    free(boot_treat);
    for ( i=0; i<n; i++ ) free(boot[i]);
    free(boot);
    free(best_split);

/*
//   print tree
    Rprintf("\nTree (nmin2=%d): \n",nmin2);
    for ( i=0; i<max; i++ ) {

        Rprintf("index=%d covar=%d type=%d sign=%d val=%f lchild=%d rchild=%d n=%d ntrt0=%d ntrt1=%d besttrt=%d\n",
               tree[i].index,tree[i].covar,tree[i].type,
               tree[i].sign,tree[i].val,tree[i].lchild,tree[i].rchild,
               tree[i].n,tree[i].n0,tree[i].n1,tree[i].predict);
    }
*/
}



void
maketree2_surv(int n,
               int nc,
               double *y,
               double **data,
               int *types,
               int *ncat,
               int *treat,
               int *censor,
               struct node *tree,
               int nmin2,
               int mtry,
               int maxdepth2,
               int method,
               int *index)

//
//  This function returns a tree for data with a
//  survival outcome variable "y" and two or more
//  treatment groups. Splits are made using a Cox
//  PH model.
//

{
//   use bootstrap data
    double *boot_y=(double *)malloc(n*sizeof(double));
    int *boot_treat=(int *)malloc(n*sizeof(int));
    int *boot_censor=(int *)malloc(n*sizeof(int));

    int i;
    double **boot=(double **)malloc(n*sizeof(double *));
    for ( i=0; i<n; i++ ) {
        boot[i]=(double *)malloc(nc*sizeof(double));
    }

    int j;
    for ( i=0; i<n; i++ ) {

        boot_y[i]=y[index[i]];
        boot_treat[i]=treat[index[i]];
        boot_censor[i]=censor[index[i]];

        for ( j=0; j<nc; j++ ) {
            boot[i][j]=data[index[i]][j];
        }
    }

//   initialize root node
    tree[0].index=1;
    tree[0].depth=1;
    tree[0].n=n;

    if ( method >= 20 ) {  // multi-treat methods

        tree[0].n0=get_min_ntrt(1,n,nc,boot,ncat,boot_treat,tree);
        tree[0].n1=tree[0].n0;

        tree[0].predict=0;
        tree[0].r0=0;
        tree[0].r1=0;
        tree[0].p0=0;
        tree[0].p1=0;

    } else {               // binary treat methods

        tree[0].n0=get_subj_treat(1,0,n,nc,boot,ncat,boot_treat,tree);
        tree[0].n1=get_subj_treat(1,1,n,nc,boot,ncat,boot_treat,tree);

        double r[4];
        tree[0].predict=get_node_predict_y(1,n,nc,boot_y,boot,ncat,
                                           boot_treat,r,tree);
        tree[0].r0=r[0];
        tree[0].r1=r[1];
        tree[0].p0=r[2];
        tree[0].p1=r[3];
    }

//   build tree
    int current_node,isplit,max,d;
    struct node *best_split=malloc(2*sizeof(struct node));

    for ( i=0; i<MAXSIZE; i++ ) {

        current_node=i+1;

        d=tree[current_node-1].depth+1;
        if ( (maxdepth2 > 0) && (d > maxdepth2) ) break;

        if ( (method == 12) || (method == 20) ) {

            isplit=get_split2_cox_all(current_node,n,nc,boot_y,boot,
                                      types,ncat,boot_treat,
                                      boot_censor,tree,nmin2,method,
                                      best_split);

        } else {  // if ( (method == 13) || (method == 21) )

            isplit=get_split2_mtry_surv(current_node,n,nc,boot_y,boot,
                                        types,ncat,boot_treat,
                                        boot_censor,tree,nmin2,mtry,
                                        method,best_split);
        }

        if ( isplit != -7 ) {  // "-7" = there should be no split

            max=get_number_of_nodes(MAXSIZE,tree);
            tree[i].lchild=max+1;
            tree[i].rchild=max+2;

//           add left child node to tree
            tree[max].index=tree[i].lchild;
            tree[max].covar=best_split[0].covar;
            tree[max].type=best_split[0].type;
            tree[max].sign=best_split[0].sign;
            tree[max].val=best_split[0].val;
            tree[max].parent=current_node;
            tree[max].depth=tree[current_node-1].depth+1;
            tree[max].n=best_split[0].n;
            tree[max].n0=best_split[0].n0;
            tree[max].n1=best_split[0].n1;
            tree[max].r0=best_split[0].r0;
            tree[max].r1=best_split[0].r1;
            tree[max].p0=best_split[0].p0;
            tree[max].p1=best_split[0].p1;
            tree[max].predict=best_split[0].predict;

//           add right child node to tree
            tree[max+1].index=tree[i].rchild;
            tree[max+1].covar=best_split[1].covar;
            tree[max+1].type=best_split[1].type;
            tree[max+1].sign=best_split[1].sign;
            tree[max+1].val=best_split[1].val;
            tree[max+1].parent=current_node;
            tree[max+1].depth=tree[current_node-1].depth+1;
            tree[max+1].n=best_split[1].n;
            tree[max+1].n0=best_split[1].n0;
            tree[max+1].n1=best_split[1].n1;
            tree[max+1].r0=best_split[1].r0;
            tree[max+1].r1=best_split[1].r1;
            tree[max+1].p0=best_split[1].p0;
            tree[max+1].p1=best_split[1].p1;
            tree[max+1].predict=best_split[1].predict;

//           get G statistic for current split
            tree[current_node-1].G=best_split[0].G;
        }

        max=get_number_of_nodes(MAXSIZE,tree);
        if ( (isplit == -7) && (current_node == max) ) break;
    }

    free(boot_y);
    free(boot_treat);
    free(boot_censor);
    for ( i=0; i<n; i++ ) free(boot[i]);
    free(boot);
    free(best_split);

/*
//   print tree
    Rprintf("\nTree from maketree2_surv:\n");
    for ( i=0; i<max; i++ ) {

        Rprintf("index=%d  covar=%d  lchild=%d  rchild=%d  nsubj=%d  depth=%d  G=%f\n",
		tree[i].index,tree[i].covar,tree[i].lchild,tree[i].rchild,
                tree[i].n,tree[i].depth,tree[i].G);
    }
*/
}



int
get_bestvar(int current_node,
            int ntree,
            int n,
            int nc,
            double *y,
            double **data,
            int *types,
            int *ncat,
            int *treat,
            int *censor,
            struct node *tree,
            int nmin2,
            int mtry,
            int maxdepth2,
            int method,
            int ncores)
 
//
//  This function accepts as arguments a data set
//  containing a continuous outcome variable and
//  strictly two treatment groups, data types, the
//  current node, the number of embedded trees to 
//  be constructed, and the minimum node size of 
//  these trees.
//
//  This function returns the integer index of the 
//  variable with the highest variable importance 
//  score.
//
//  Note that variable indices are numbered from 
//  1 to nc.
//

{
//
//   ... get data at current node
//
    int *ifcovar=(int *)malloc(n*sizeof(int));
    int nsubj=findrows_node(current_node,n,nc,data,ncat,tree,ifcovar);

    double *node_y=(double *)malloc(nsubj*sizeof(double));
    int *node_treat=(int *)malloc(nsubj*sizeof(int));
    int *node_censor=(int *)malloc(nsubj*sizeof(int));

    int i,j;
    double **node_data=(double **)malloc(nsubj*sizeof(double *));
    for ( i=0; i<nsubj; i++ ) {
        node_data[i]=(double *)malloc(nc*sizeof(double));
    }

    int counter=0;

    for ( i=0; i<n; i++ ) {

        if ( ifcovar[i] == 1 ) {

            node_y[counter]=y[i];
            node_treat[counter]=treat[i];
            node_censor[counter]=censor[i];

            for ( j=0; j<nc; j++ ) {
                node_data[counter][j]=data[i][j];
            }

            counter++;
            if ( counter == nsubj ) break;
        }
    }

    free(ifcovar);

//
//   ... find the "best" variable based on ntree trees
//
    double **varimp0=(double **)malloc(ntree*sizeof(double *));
    for ( i=0; i<ntree; i++ ) {
        varimp0[i]=(double *)calloc(nc,sizeof(double));
    }

    #ifdef _OPENMP
     omp_set_num_threads(ncores);
    #endif

    #pragma omp parallel private(i,j)
    {

////    Rprintf("Number of threads: %d\n",omp_get_num_threads());

        #pragma omp for
        for ( i=0; i<ntree; i++ ) {

            int *index=(int *)malloc(nsubj*sizeof(int));
            double *varimp00=(double *)calloc(nc,sizeof(double));
            struct node *tree2=calloc(MAXSIZE,sizeof(struct node));

            draw_boot(nsubj,index);

            if ( method > 10 ) {  // survival outcomes

                maketree2_surv(nsubj,nc,node_y,node_data,types,ncat,
                               node_treat,node_censor,tree2,nmin2,
                               mtry,maxdepth2,method,index);

            } else {              // continuous Y

                maketree2(nsubj,nc,node_y,node_data,types,ncat,
                          node_treat,tree2,nmin2,mtry,maxdepth2,
                          index,method);
            }

            get_varimp_replace(varimp00,nsubj,nc,node_y,
                               node_data,types,ncat,node_treat,
                               node_censor,tree2,method,index);

            for ( j=0; j<nc; j++ ) {
                varimp0[i][j] += varimp00[j];
            }

            free(tree2);
            free(index);
            free(varimp00);
        }
    }

    free(node_y);
    free(node_treat);
    free(node_censor);
    for ( i=0; i<nsubj; i++ ) free(node_data[i]);
    free(node_data);

//   calculate final variable importance scores
    double *varimp=(double *)calloc(nc,sizeof(double)); // initialize

    for ( i=0; i<nc; i++ ) {
        for ( j=0; j<ntree; j++ ) {
            varimp[i] += varimp0[j][i];
        }
    }

    for ( i=0; i<ntree; i++ ) free(varimp0[i]);
    free(varimp0);

    for ( i=0; i<nc; i++ ) {
        varimp[i] /= ntree;
    }

/*
//   check 
    if ( (current_node == 1) ) {
        Rprintf("\nVar. importance scores at node %d",current_node); 
        Rprintf(" (get_bestvar, method=%d):\n",method);
        for ( i=0; i<nc; i++ ) {
            Rprintf("covar=%d  varimp=%f \n",i+1,varimp[i]);
        }
    }
*/

//   get integer index of best variable (vars numbered 1:nc)
    int bestvar=-7;
    double max=0;

    for ( i=0; i<nc; i++ ) {

        if ( varimp[i] <= 0 ) continue; 

        if ( varimp[i] > max ) {
            max=varimp[i];
            bestvar=i+1;
        }
    }

    free(varimp);

    return(bestvar);
}



int
get_bestvar_surv(int current_node,
                 int ntree,
                 int n,
                 int nc,
                 double *y,
                 double **data,
                 int *types,
                 int *ncat,
                 int *treat,
                 int *censor,
                 struct node *tree,
                 int nmin2,
                 int mtry,
                 int maxdepth2,
                 int method,
                 int ncores)
 
//
//  This function accepts as arguments a data set
//  containing a right-censored survival outcome
//  variable, data types, the current node, the 
//  number of embedded trees to be constructed, 
//  and the minimum node size of these trees.
//
//  This function returns the integer index of the 
//  variable with the highest variable importance 
//  score.
//
//  Note that variable indices are numbered from 
//  1 to nc.
//

{
//
//   ... get data at current node
//
    int *ifcovar=(int *)malloc(n*sizeof(int));
    int nsubj=findrows_node(current_node,n,nc,data,ncat,tree,ifcovar);

    double *node_y=(double *)malloc(nsubj*sizeof(double));
    int *node_treat=(int *)malloc(nsubj*sizeof(int));
    int *node_censor=(int *)malloc(nsubj*sizeof(int));

    int i,j;
    double **node_data=(double **)malloc(nsubj*sizeof(double *));
    for ( i=0; i<nsubj; i++ ) {
        node_data[i]=(double *)malloc(nc*sizeof(double));
    }

    int counter=0;

    for ( i=0; i<n; i++ ) {

        if ( ifcovar[i] == 1 ) {

            node_y[counter]=y[i];
            node_treat[counter]=treat[i];
            node_censor[counter]=censor[i];

            for ( j=0; j<nc; j++ ) {
                node_data[counter][j]=data[i][j];
            }

            counter++;
            if ( counter == nsubj ) break;
        }
    }

    free(ifcovar);

//
//   ... find the "best" variable based on ntree trees
//
    double **varimp0=(double **)malloc(ntree*sizeof(double *));
    for ( i=0; i<ntree; i++ ) {
        varimp0[i]=(double *)calloc(nc,sizeof(double));
    }

    for ( i=0; i<ntree; i++ ) {

////        Rprintf("\n---- tree = %d -----\n",i+1);

        double *varimp00=(double *)calloc(nc,sizeof(double));
        struct node *tree2=calloc(MAXSIZE,sizeof(struct node));

/*
//       draw bootstrap sample of size nsubj with replacement
        int *index=(int *)malloc(nsubj*sizeof(int));
        draw_boot(nsubj,index);

        maketree2_surv(nsubj,nc,node_y,node_data,types,ncat,
                       node_treat,node_censor,tree2,nmin2,
                       mtry,maxdepth2,method,index);

        get_varimp_replace(varimp00,nsubj,nc,node_y,
                           node_data,types,ncat,node_treat,
                           node_censor,tree2,method,index);
*/

//       bootstrap 80% without replacement only if ntree > 1
        int nval=floor(nsubj*0.8);
        if ( ntree == 1 ) nval=nsubj;
        int *index=(int *)malloc(nval*sizeof(int));

        if ( ntree == 1 ) {

            for ( j=0; j<nsubj; j++ ) {
                index[j]=j;
            }

        } else {

            sample_without_replace(nsubj,nval,index);
        }

/*
//       check
        Rprintf("\nnval=%d\n",nval);
        Rprintf("\nindex=\n");
        for ( j=0; j<10; j++ ) Rprintf("%d  ",index[j]);
*/

        maketree2_surv(nval,nc,node_y,node_data,types,ncat,
                       node_treat,node_censor,tree2,nmin2,
                       mtry,maxdepth2,method,index);

        get_varimp_replace(varimp00,nval,nc,node_y,
                           node_data,types,ncat,node_treat,
                           node_censor,tree2,method,index);

        for ( j=0; j<nc; j++ ) {
            varimp0[i][j] += varimp00[j];
        }

        free(tree2);
        free(index);
        free(varimp00);
    }

    free(node_y);
    free(node_treat);
    free(node_censor);
    for ( i=0; i<nsubj; i++ ) free(node_data[i]);
    free(node_data);

//   calculate final variable importance scores
    double *varimp=(double *)calloc(nc,sizeof(double)); // initialize

    for ( i=0; i<nc; i++ ) {
        for ( j=0; j<ntree; j++ ) {
            varimp[i] += varimp0[j][i];
        }
    }

    for ( i=0; i<ntree; i++ ) free(varimp0[i]);
    free(varimp0);

    for ( i=0; i<nc; i++ ) {
        varimp[i] /= ntree;
    }

/*
//   check 
    if ( (current_node == 1) ) {
        Rprintf("\nVar. importance scores at node %d",current_node); 
        Rprintf(" (get_bestvar, method=%d):\n",method);
        for ( i=0; i<nc; i++ ) {
            Rprintf("covar=%d  varimp=%f \n",i+1,varimp[i]);
        }
    }
*/

//   get integer index of best variable (vars numbered 1:nc)
    int bestvar=-7;
    double max=0;

    for ( i=0; i<nc; i++ ) {

        if ( varimp[i] <= 0 ) continue; 

        if ( varimp[i] > max ) {
            max=varimp[i];
            bestvar=i+1;
        }
    }

    free(varimp);

    return(bestvar);
}



int
get_bestvar_mc(int current_node,
               int ntree,
               int n,
               int nc,
               double *y,
               double **data,
               int *types,
               int *ncat,
               int *treat,
               int *censor,
               struct node *tree,
               int nmin2,
               int mtry,
               int maxdepth2,
               int method)

//
//  This function accepts as arguments a data set
//  containing a continuous outcome variable and
//  more than two treatment groups, data types, 
//  the current node, the number of embedded trees 
//  to be constructed, and the minimum node size of 
//  these trees.
//
//  This function returns the integer index of the 
//  variable with the highest variable importance 
//  score.
//
//  Note that variable indices are numbered from 
//  1 to nc.
//

{
//
//   ... get data at current node
//
    int *ifcovar=(int *)malloc(n*sizeof(int));
    int nsubj=findrows_node(current_node,n,nc,data,ncat,tree,ifcovar);

    double *node_y=(double *)malloc(nsubj*sizeof(double));
    int *node_treat=(int *)malloc(nsubj*sizeof(int));
    int *node_censor=(int *)malloc(nsubj*sizeof(int));

    int i,j;
    double **node_data=(double **)malloc(nsubj*sizeof(double *));
    for ( i=0; i<nsubj; i++ ) {
        node_data[i]=(double *)malloc(nc*sizeof(double));
    }

    int counter=0;

    for ( i=0; i<n; i++ ) {

        if ( ifcovar[i] == 1 ) {

            node_y[counter]=y[i];
            node_treat[counter]=treat[i];
            node_censor[counter]=censor[i];

            for ( j=0; j<nc; j++ ) {
                node_data[counter][j]=data[i][j];
            }

            counter++;
            if ( counter == nsubj ) break;
        }
    }

    free(ifcovar);

//
//   ... find the "best" variable based on ntree trees
//
    double **varimp0=(double **)malloc(ntree*sizeof(double *));
    for ( i=0; i<ntree; i++ ) {
        varimp0[i]=(double *)calloc(nc,sizeof(double));
    }

    for ( i=0; i<ntree; i++ ) {

        int *index=(int *)malloc(nsubj*sizeof(int));
        double *varimp00=(double *)calloc(nc,sizeof(double));
        struct node *tree2=calloc(MAXSIZE,sizeof(struct node));

        draw_boot(nsubj,index);

        maketree2(nsubj,nc,node_y,node_data,types,ncat,
                  node_treat,tree2,nmin2,mtry,maxdepth2,
                  index,method);

        get_varimp_replace(varimp00,nsubj,nc,node_y,
                           node_data,types,ncat,node_treat,
                           node_censor,tree2,method,index);

        for ( j=0; j<nc; j++ ) {
            varimp0[i][j] += varimp00[j];
        }

        free(tree2);
        free(index);
        free(varimp00);
    }

    free(node_y);
    free(node_treat);
    free(node_censor);
    for ( i=0; i<nsubj; i++ ) free(node_data[i]);
    free(node_data);

//   calculate final variable importance scores
    double *varimp=(double *)calloc(nc,sizeof(double)); // initialize

    for ( i=0; i<nc; i++ ) {
        for ( j=0; j<ntree; j++ ) {
            varimp[i] += varimp0[j][i];
        }
    }

    for ( i=0; i<ntree; i++ ) free(varimp0[i]);
    free(varimp0);

    for ( i=0; i<nc; i++ ) {
        varimp[i] /= ntree;
    }

/*
//   check 
    if ( (current_node == 1) ) {
        Rprintf("\nVar. importance scores at node %d",current_node); 
        Rprintf(" (get_bestvar, method=%d):\n",method);
        for ( i=0; i<nc; i++ ) {
            Rprintf("covar=%d  varimp=%f \n",i+1,varimp[i]);
        }
    }
*/

//   get integer index of best variable (vars numbered 1:nc)
    int bestvar=-7;
    double max=0;

    for ( i=0; i<nc; i++ ) {

        if ( varimp[i] <= 0 ) continue; 

        if ( varimp[i] > max ) {
            max=varimp[i];
            bestvar=i+1;
        }
    }

    free(varimp);

    return(bestvar);
}



int
get_split(int current_node,
          int ntree,
          int n,
          int nc,
          double *y,
          double **data,
          int *types,
          int *ncat,
          int *treat,
          int *censor,
          struct node *tree,
          int nmin,
          int nmin2,
          int mtry,
          int maxdepth2,
          int method,
          int ncores,
          struct node *best_split)

//
//  This function accepts as arguments the integer
//  index of the current node "current_node", the
//  number of embedded trees to fit "ntree", tree
//  information "tree", the minimum node size
//  "nmin", and the number of variables to randomly
//  select "mtry" and returns the best split as 
//  node struct array of size 2 "best_split".
//
//  If there is no best split, then the function 
//  returns integer "-7".
//
//  Once the best split variable is found, all
//  possible splits are considered, and split
//  scores use the DIFF value.
//

{
//   get best variable
    int bestvar=-7;

    if ( (method == 22) || (method == 23) ) { // multitrt continuous 
                                              // methods

        bestvar=get_bestvar_mc(current_node,ntree,n,nc,y,data,types,
                               ncat,treat,censor,tree,nmin2,
                               mtry,maxdepth2,method);

    } else {

        bestvar=get_bestvar(current_node,ntree,n,nc,y,data,types,ncat,
                            treat,censor,tree,nmin2,mtry,
                            maxdepth2,method,ncores);
    }

    if ( bestvar == -7 ) return(-7);  // no bestvar, no best split

    best_split[0].covar=bestvar;
    best_split[1].covar=bestvar;

//   col_node is the data column of best variable at current node
    int *ifcovar=(int *)malloc(n*sizeof(int));
    findrows_node(current_node,n,nc,data,ncat,tree,ifcovar);

    double *col=(double *)malloc(n*sizeof(double));
    get_col(data,bestvar,n,nc,col);

    int nsubj_node=tree[current_node-1].n;
    double *col_node=(double *)malloc(nsubj_node*sizeof(double));
    int col_counter=0;

    int i;
    for ( i=0; i<n; i++ ) {

        if ( ifcovar[i] == 1 ) {
            col_node[col_counter]=col[i];
            col_counter++;
        }
    }

    free(ifcovar);
    free(col);

//   get and sort unique values of the data column
    double *unique=(double *)malloc(nsubj_node*sizeof(double));
    int nunique=get_nunique(nsubj_node,col_node,unique);

    free(col_node);

//   if all values are the same, then no split
    if ( nunique == 1 ) {
        free(unique);
        return(-7);
    }

    double uniquevals[nunique];
    for ( i=0; i<nunique; i++ ) {
        uniquevals[i]=unique[i];
    }

    free(unique);

    qsort(uniquevals,nunique,sizeof(double),compare_doubles);

//   initialize arrays in temporary tree
    struct node temp_tree[current_node+1];
    for ( i=0; i<current_node; i++ ) {

        temp_tree[i].index=tree[i].index;
        temp_tree[i].covar=tree[i].covar;
        temp_tree[i].type=tree[i].type;
        temp_tree[i].sign=tree[i].sign;
        temp_tree[i].val=tree[i].val;
        temp_tree[i].parent=tree[i].parent;
        temp_tree[i].depth=tree[i].depth;
        temp_tree[i].lchild=tree[i].lchild;
        temp_tree[i].rchild=tree[i].rchild;
        temp_tree[i].n=tree[i].n;
        temp_tree[i].n0=tree[i].n0;
        temp_tree[i].n1=tree[i].n1;
    }

    temp_tree[current_node].index=current_node+1;
    temp_tree[current_node].covar=bestvar;
    temp_tree[current_node].parent=current_node;
    temp_tree[current_node].depth=tree[current_node-1].depth+1;

//   get "DIFF" value of current node
    double diff;
    if ( (method == 22) || (method == 23) ) { // multitrt continuous 
                                              // methods

        diff=get_diff_mc(current_node,n,nc,y,data,types,ncat,treat,
                         tree);

    } else {

        diff=get_diff(current_node,n,nc,y,data,types,ncat,treat,
                      tree);
    }

//   assess candidate splits by variable type
    int type=types[bestvar-1];  // type of best variable
    int lll=current_node+1;
    struct cand *cands=malloc(2*sizeof(struct cand));
    double max=diff;
    double score;

    if ( type == 1 ) {  // binary

        temp_tree[current_node].type=1;
        temp_tree[current_node].sign=0;

//       get values for left child candidate node
        temp_tree[current_node].val=0;
        get_split0_d(0,lll,n,nc,y,data,types,ncat,treat,method,
                     temp_tree,cands);

//       get values for right child candidate node
        temp_tree[current_node].val=1;
        get_split0_d(1,lll,n,nc,y,data,types,ncat,treat,method,
                     temp_tree,cands);

//       candidate splits with too few subjects are not useful
        if ( (cands[0].n0 <= nmin2) || (cands[0].n1 <= nmin2) || 
             (cands[1].n0 <= nmin2) || (cands[1].n1 <= nmin2) ) {

            free(cands);
            return(-7);
        }

        score=(cands[0].val+cands[1].val) / (cands[0].n+cands[1].n);

        if ( score > max ) {

            max=score;

            best_split[0].type=1;
            best_split[0].sign=0;
            best_split[0].val=0;
            best_split[0].n=cands[0].n;
            best_split[0].n0=cands[0].n0;
            best_split[0].n1=cands[0].n1;
            best_split[0].r0=cands[0].r0;
            best_split[0].r1=cands[0].r1;
            best_split[0].p0=cands[0].p0;
            best_split[0].p1=cands[0].p1;
            best_split[0].predict=cands[0].p;

            best_split[1].type=1;
            best_split[1].sign=0;
            best_split[1].val=1;
            best_split[1].n=cands[1].n;
            best_split[1].n0=cands[1].n0;
            best_split[1].n1=cands[1].n1;
            best_split[1].r0=cands[1].r0;
            best_split[1].r1=cands[1].r1;
            best_split[1].p0=cands[1].p0;
            best_split[1].p1=cands[1].p1;
            best_split[1].predict=cands[1].p;
        }

    } else if ( type == 2 ) {  // ordinal

        temp_tree[current_node].type=2;

        for ( i=0; i<(nunique-1); i++ ) {

            temp_tree[current_node].val=uniquevals[i];

//           get values for left child candidate node
            temp_tree[current_node].sign=1; // LE
            get_split0_d(0,lll,n,nc,y,data,types,ncat,treat,method,
                         temp_tree,cands);

//           get values for right child candidate node
            temp_tree[current_node].sign=2; // GT
            get_split0_d(1,lll,n,nc,y,data,types,ncat,treat,method,
                         temp_tree,cands);

//           candidate splits with too few subjects are not useful
            if ( (cands[0].n0 <= nmin2) || (cands[0].n1 <= nmin2) || 
                 (cands[1].n0 <= nmin2) || (cands[1].n1 <= nmin2) ) {

                continue;
            }

            score=(cands[0].val+cands[1].val)/(cands[0].n+cands[1].n);

            if ( score > max ) {

                max=score;

                best_split[0].type=2;
                best_split[0].sign=1; // LE
                best_split[0].val=uniquevals[i];
                best_split[0].n=cands[0].n;
                best_split[0].n0=cands[0].n0;
                best_split[0].n1=cands[0].n1;
                best_split[0].r0=cands[0].r0;
                best_split[0].r1=cands[0].r1;
                best_split[0].p0=cands[0].p0;
                best_split[0].p1=cands[0].p1;
                best_split[0].predict=cands[0].p;

                best_split[1].type=2;
                best_split[1].sign=2; // GT
                best_split[1].val=uniquevals[i];
                best_split[1].n=cands[1].n;
                best_split[1].n0=cands[1].n0;
                best_split[1].n1=cands[1].n1;
                best_split[1].r0=cands[1].r0;
                best_split[1].r1=cands[1].r1;
                best_split[1].p0=cands[1].p0;
                best_split[1].p1=cands[1].p1;
                best_split[1].predict=cands[1].p;
            }
        }

    } else if ( type == 3 ) {  // nominal

        temp_tree[current_node].type=3;
        temp_tree[current_node].sign=0;

//       calculate total number of possible splits
        unsigned long long int n_splits=pow(2,nunique-1)-1;

//       2 categories reduce to binary split
        if ( nunique == 2 ) {  

            double nom_lval=pow(2,ncat[bestvar-1]-uniquevals[0]);
            double nom_rval=pow(2,ncat[bestvar-1]-uniquevals[1]);

//           get values for left child candidate node
            temp_tree[current_node].val=nom_lval;
            get_split0_d(0,lll,n,nc,y,data,types,ncat,treat,method,
                         temp_tree,cands);

//           get values for right child candidate node
            temp_tree[current_node].val=nom_rval;
            get_split0_d(1,lll,n,nc,y,data,types,ncat,treat,method,
                         temp_tree,cands);

//           candidate splits with too few subjects are not useful
            if ( (cands[0].n0 <= nmin2) || 
                 (cands[0].n1 <= nmin2) || 
                 (cands[1].n0 <= nmin2) ||
                 (cands[1].n1 <= nmin2) ) {

                free(cands);
                return(-7);
            }

            score=(cands[0].val+cands[1].val)/(cands[0].n+cands[1].n);

            if ( score > max ) {

                max=score;

                best_split[0].type=3;
                best_split[0].sign=0;
                best_split[0].val=nom_lval;
                best_split[0].n=cands[0].n;
                best_split[0].n0=cands[0].n0;
                best_split[0].n1=cands[0].n1;
                best_split[0].r0=cands[0].r0;
                best_split[0].r1=cands[0].r1;
                best_split[0].p0=cands[0].p0;
                best_split[0].p1=cands[0].p1;
                best_split[0].predict=cands[0].p;

                best_split[1].type=3;
                best_split[1].sign=0;
                best_split[1].val=nom_rval;
                best_split[1].n=cands[1].n;
                best_split[1].n0=cands[1].n0;
                best_split[1].n1=cands[1].n1;
                best_split[1].r0=cands[1].r0;
                best_split[1].r1=cands[1].r1;
                best_split[1].p0=cands[1].p0;
                best_split[1].p1=cands[1].p1;
                best_split[1].predict=cands[1].p;
            }

//       for nominal variables with more than 2 categories
        } else {

            int z,k;
            int current_var=bestvar;
            int k_index,cat_index;
            int max_k=floor(nunique/2);
            unsigned long long int counter=0;

            for ( z=0; z<max_k; z++ ) {

                k=z+1;        // size of 1 combination
                k_index=0;    // current index in a single comb.
                cat_index=0;  // current index in list of cats.

                double temp[k]; // array to temporarily store 
                                // a combination

                max=get_comb_vals(current_node,current_var,
                                  n,nc,y,data,types,ncat,treat,
                                  nmin2,temp,uniquevals,nunique,
                                  k,ncat[bestvar-1],n_splits,&counter,
                                  cat_index,k_index,max,method,
                                  temp_tree,cands,best_split);
            }
        }
    }

    free(cands);

    if ( max == diff ) return(-7);
    return(0);
}



int
get_split_surv(int current_node,
               int ntree,
               int n,
               int nc,
               double *y,
               double **data,
               int *types,
               int *ncat,
               int *treat,
               int *censor,
               struct node *tree,
               int nmin,
               int nmin2,
               int mtry,
               int maxdepth2,
               int method,
               int ncores,
               struct node *best_split)

//
//  This function accepts as arguments the integer
//  index of the current node "current_node", the
//  number of embedded trees to fit "ntree", tree
//  information "tree", the minimum node size
//  "nmin", and the number of variables to randomly
//  select "mtry" and returns the best split as 
//  node struct array of size 2 "best_split".
//
//  If there is no best split, then the function 
//  returns integer "-7".
//
//  Once the best split variable is found, all
//  possible splits are considered, and split
//  scores are based on a Cox PH model.
//

{
//   get best variable
    int bestvar=get_bestvar_surv(current_node,ntree,n,nc,y,data,types,
                                 ncat,treat,censor,tree,nmin2,
                                 mtry,maxdepth2,method,ncores);

////    Rprintf("bestvar=%d\n",bestvar);

    if ( bestvar == -7 ) return(-7);  // no bestvar, no best split

    best_split[0].covar=bestvar;
    best_split[1].covar=bestvar;

//   col_node is the data column of best variable at current node
    int *ifcovar=(int *)malloc(n*sizeof(int));
    findrows_node(current_node,n,nc,data,ncat,tree,ifcovar);

    double *col=(double *)malloc(n*sizeof(double));
    get_col(data,bestvar,n,nc,col);

    int nsubj_node=tree[current_node-1].n;
    double *col_node=(double *)malloc(nsubj_node*sizeof(double));
    int col_counter=0;

    int i;
    for ( i=0; i<n; i++ ) {

        if ( ifcovar[i] == 1 ) {
            col_node[col_counter]=col[i];
            col_counter++;
        }
    }

    free(ifcovar);
    free(col);

//   get and sort unique values of the data column
    double *unique=(double *)malloc(nsubj_node*sizeof(double));
    int nunique=get_nunique(nsubj_node,col_node,unique);

    free(col_node);

//   if all values are the same, then no split
    if ( nunique == 1 ) {
        free(unique);
        return(-7);
    }

    double uniquevals[nunique];
    for ( i=0; i<nunique; i++ ) {
        uniquevals[i]=unique[i];
    }

    free(unique);

    qsort(uniquevals,nunique,sizeof(double),compare_doubles);

//   initialize arrays in temporary tree
    struct node temp_tree[current_node+2];
    for ( i=0; i<current_node; i++ ) {

        temp_tree[i].index=tree[i].index;
        temp_tree[i].covar=tree[i].covar;
        temp_tree[i].type=tree[i].type;
        temp_tree[i].sign=tree[i].sign;
        temp_tree[i].val=tree[i].val;
        temp_tree[i].parent=tree[i].parent;
        temp_tree[i].depth=tree[i].depth;
        temp_tree[i].lchild=tree[i].lchild;
        temp_tree[i].rchild=tree[i].rchild;
        temp_tree[i].n=tree[i].n;
        temp_tree[i].n0=tree[i].n0;
        temp_tree[i].n1=tree[i].n1;
    }

    temp_tree[current_node-1].lchild=current_node+1;
    temp_tree[current_node].covar=bestvar;
    temp_tree[current_node].index=current_node+1;
    temp_tree[current_node].parent=current_node;
    temp_tree[current_node].depth=tree[current_node-1].depth+1;

    temp_tree[current_node-1].rchild=current_node+2;
    temp_tree[current_node+1].covar=bestvar;
    temp_tree[current_node+1].index=current_node+2;
    temp_tree[current_node+1].parent=current_node;
    temp_tree[current_node+1].depth=tree[current_node-1].depth+1;

//   assess candidate splits by variable type
    int type=types[bestvar-1];  // type of best variable
    int lll=current_node+1;
    struct cand *cands=malloc(2*sizeof(struct cand));
    double max=0;
    double score;

    if ( type == 1 ) {  // binary

//       get values for left child candidate node
        temp_tree[current_node].type=1;
        temp_tree[current_node].sign=0;
        temp_tree[current_node].val=0;

        get_split0_cox(0,lll,n,nc,y,data,types,ncat,treat,
                       censor,method,temp_tree,cands);

//       get values for right child candidate node
        temp_tree[current_node+1].type=1;
        temp_tree[current_node+1].sign=0;
        temp_tree[current_node+1].val=1;

        get_split0_cox(1,lll+1,n,nc,y,data,types,ncat,treat,
                       censor,method,temp_tree,cands);

//       candidate splits with too few subjects are not useful
        if ( (cands[0].n0 <= nmin) || (cands[0].n1 <= nmin) || 
             (cands[1].n0 <= nmin) || (cands[1].n1 <= nmin) ) {

            free(cands);
            return(-7);
        }

        score=get_coxph_stat(current_node,n,nc,y,data,types,
                             ncat,treat,censor,method,temp_tree);

        if ( score > max ) {

            max=score;

            best_split[0].type=1;
            best_split[0].sign=0;
            best_split[0].val=0;
            best_split[0].n=cands[0].n;
            best_split[0].n0=cands[0].n0;
            best_split[0].n1=cands[0].n1;
            best_split[0].r0=cands[0].r0;
            best_split[0].r1=cands[0].r1;
            best_split[0].p0=cands[0].p0;
            best_split[0].p1=cands[0].p1;
            best_split[0].predict=cands[0].p;

            best_split[1].type=1;
            best_split[1].sign=0;
            best_split[1].val=1;
            best_split[1].n=cands[1].n;
            best_split[1].n0=cands[1].n0;
            best_split[1].n1=cands[1].n1;
            best_split[1].r0=cands[1].r0;
            best_split[1].r1=cands[1].r1;
            best_split[1].p0=cands[1].p0;
            best_split[1].p1=cands[1].p1;
            best_split[1].predict=cands[1].p;
        }

    } else if ( type == 2 ) {  // ordinal

        temp_tree[current_node].type=2;
        temp_tree[current_node+1].type=2;

        for ( i=0; i<(nunique-1); i++ ) {

//           get values for left child candidate node
            temp_tree[current_node].val=uniquevals[i];
            temp_tree[current_node].sign=1; // LE

            get_split0_cox(0,lll,n,nc,y,data,types,ncat,treat,
                           censor,method,temp_tree,cands);

//           get values for right child candidate node
            temp_tree[current_node+1].val=uniquevals[i];
            temp_tree[current_node+1].sign=2; // GT

            get_split0_cox(1,lll+1,n,nc,y,data,types,ncat,treat,
                           censor,method,temp_tree,cands);

//           candidate splits with too few subjects are not useful
            if ( (cands[0].n0 <= nmin) || (cands[0].n1 <= nmin) || 
                 (cands[1].n0 <= nmin) || (cands[1].n1 <= nmin) ) {

                continue;
            }

            score=get_coxph_stat(current_node,n,nc,y,data,types,
                                 ncat,treat,censor,method,temp_tree);

            if ( score > max ) {

                max=score;

                best_split[0].type=2;
                best_split[0].sign=1; // LE
                best_split[0].val=uniquevals[i];
                best_split[0].n=cands[0].n;
                best_split[0].n0=cands[0].n0;
                best_split[0].n1=cands[0].n1;
                best_split[0].r0=cands[0].r0;
                best_split[0].r1=cands[0].r1;
                best_split[0].p0=cands[0].p0;
                best_split[0].p1=cands[0].p1;
                best_split[0].predict=cands[0].p;

                best_split[1].type=2;
                best_split[1].sign=2; // GT
                best_split[1].val=uniquevals[i];
                best_split[1].n=cands[1].n;
                best_split[1].n0=cands[1].n0;
                best_split[1].n1=cands[1].n1;
                best_split[1].r0=cands[1].r0;
                best_split[1].r1=cands[1].r1;
                best_split[1].p0=cands[1].p0;
                best_split[1].p1=cands[1].p1;
                best_split[1].predict=cands[1].p;
            }
        }

    } else if ( type == 3 ) {  // nominal

        temp_tree[current_node].type=3;
        temp_tree[current_node].sign=0;

        temp_tree[current_node+1].type=3;
        temp_tree[current_node+1].sign=0;

//       calculate total number of possible splits
        unsigned long long int n_splits=pow(2,nunique-1)-1;

//       2 categories reduce to binary split
        if ( nunique == 2 ) {  

            double nom_lval=pow(2,ncat[bestvar-1]-uniquevals[0]);
            double nom_rval=pow(2,ncat[bestvar-1]-uniquevals[1]);

//           get values for left child candidate node
            temp_tree[current_node].val=nom_lval;
            get_split0_cox(0,lll,n,nc,y,data,types,ncat,treat,
                           censor,method,temp_tree,cands);

//           get values for right child candidate node
            temp_tree[current_node+1].val=nom_rval;
            get_split0_cox(1,lll+1,n,nc,y,data,types,ncat,treat,
                           censor,method,temp_tree,cands);

//           candidate splits with too few subjects are not useful
            if ( (cands[0].n0 <= nmin) || 
                 (cands[0].n1 <= nmin) || 
                 (cands[1].n0 <= nmin) ||
                 (cands[1].n1 <= nmin) ) {

                free(cands);
                return(-7);
            }

            score=get_coxph_stat(current_node,n,nc,y,data,types,
                                 ncat,treat,censor,method,temp_tree);

            if ( score > max ) {

                max=score;

                best_split[0].type=3;
                best_split[0].sign=0;
                best_split[0].val=nom_lval;
                best_split[0].n=cands[0].n;
                best_split[0].n0=cands[0].n0;
                best_split[0].n1=cands[0].n1;
                best_split[0].r0=cands[0].r0;
                best_split[0].r1=cands[0].r1;
                best_split[0].p0=cands[0].p0;
                best_split[0].p1=cands[0].p1;
                best_split[0].predict=cands[0].p;

                best_split[1].type=3;
                best_split[1].sign=0;
                best_split[1].val=nom_rval;
                best_split[1].n=cands[1].n;
                best_split[1].n0=cands[1].n0;
                best_split[1].n1=cands[1].n1;
                best_split[1].r0=cands[1].r0;
                best_split[1].r1=cands[1].r1;
                best_split[1].p0=cands[1].p0;
                best_split[1].p1=cands[1].p1;
                best_split[1].predict=cands[1].p;
            }

//       for nominal variables with more than 2 categories
        } else {

            int z,k;
            int current_var=bestvar;
            int k_index,cat_index;
            int max_k=floor(nunique/2);
            unsigned long long int counter=0;

            for ( z=0; z<max_k; z++ ) {

                k=z+1;        // size of 1 combination
                k_index=0;    // current index in a single comb.
                cat_index=0;  // current index in list of cats.

                double temp[k]; // array to temporarily store 
                                // a combination

                max=get_comb_vals_surv(current_node,current_var,
                                       n,nc,y,data,types,ncat,
                                       treat,censor,nmin2,method,
                                       temp,uniquevals,nunique,k,
                                       ncat[bestvar-1],n_splits,
                                       &counter,cat_index,k_index,max,
                                       temp_tree,cands,best_split);
            }
        }
    }

    free(cands);

    if ( max <= 0 ) return(-7);
    return(0);
}



SEXP
maketree(SEXP R_ntree,
         SEXP R_n,
         SEXP R_nc,
         SEXP R_y,
         SEXP R_data,
         SEXP R_types,
         SEXP R_ncat,
         SEXP R_treat,
         SEXP R_censor,
         SEXP R_nmin,
         SEXP R_nmin2,
         SEXP R_mtry,
         SEXP R_maxdepth,
         SEXP R_maxdepth2,
         SEXP R_method)

//
//  This function accepts as arguments from R: 
//
//      "R_ntree"  - the total number of embedded trees
//                   constructed in the forest at each 
//                   node (when applicable),
//      "R_n"      - the total number of subjects, 
//      "R_nc"     - the total number of candidate split
//                   variables, 
//      "R_y"      - the vector of outcome variable values,
//      "R_data"   - the matrix of split variable values,
//      "R_types"  - the vector of variable types for each 
//                   variable,
//      "R_ncat"   - the vector of the number of categories 
//                   for each candidate split variable, 
//      "R_treat"  - the vector of integer treatment 
//                   assignments,
//      "R_censor" - the vector of censoring indicators
//                   (when applicable),
//      "R_nmin"   - the minimum node size of the tree, 
//      "R_nmin2"  - the minimum node size of embedded
//                   trees (when applicable),
//      "R_mtry"   - the number of variables randomly 
//                   selected at each node of the embedded
//                   trees (when applicable),
//      "R_maxdepth"  - an integer specifying the maximum
//                      depth of the overall tree,
//      "R_maxdepth2" - an integer specifying the maximum
//                      depth of embedded trees (when
//                      applicable),
//      "R_method" - an integer denoting the tree method.
//
//  Overall, this function returns a tree from C to R.
//
//  A key for the different "R_method" values:
//
//                 PM or   Continuous or      Number of    mtry
//      R_method   DIPM    Survival Outcome   Treatments   Yes/No
//      --------   -----   ----------------   ----------   ------
//      -1         PM      Continuous         2            
//       6         DIPM    Continuous         2            No
//       7         DIPM    Continuous         2            Yes
//      11         PM      Survival           2            
//      12         DIPM    Survival           2            No
//      13         DIPM    Survival           2            Yes
//      20         DIPM    Survival           2+           No
//      21         DIPM    Survival           2+           Yes
//      22         DIPM    Continuous         2+           No
//      23         DIPM    Continuous         2+           Yes
//      24         PM      Continuous         2+           
//      25         PM      Survival           2+           
//
//  For candidate split variables, types are: 
//
//      1=binary, 
//      2=ordinal,
//      3=nominal.
//
//  Binary variables must take values of 0 and 1.
//  Nominal variables must take values from 1:ncat.
//

{
//   initialize objects from R into C
    R_ntree=coerceVector(R_ntree,INTSXP);
    R_n=coerceVector(R_n,INTSXP);
    R_nc=coerceVector(R_nc,INTSXP);
    R_nmin=coerceVector(R_nmin,INTSXP);
    R_nmin2=coerceVector(R_nmin2,INTSXP);
    R_mtry=coerceVector(R_mtry,INTSXP);
    R_maxdepth=coerceVector(R_maxdepth,INTSXP);
    R_method=coerceVector(R_method,INTSXP);

    int ntree=INTEGER(R_ntree)[0];
    int n=INTEGER(R_n)[0];
    int nc=INTEGER(R_nc)[0];
    int nmin=INTEGER(R_nmin)[0];
    int nmin2=INTEGER(R_nmin2)[0];
    int mtry=INTEGER(R_mtry)[0];
    int maxdepth=INTEGER(R_maxdepth)[0];
    int maxdepth2=INTEGER(R_maxdepth2)[0];
    int method=INTEGER(R_method)[0];

    double *y=REAL(R_y);
    int *treat=INTEGER(R_treat);
    int *types=INTEGER(R_types);
    int *ncat=INTEGER(R_ncat);
    int *censor=INTEGER(R_censor);

    double **data=(double **)malloc(n*sizeof(double *));
    int i;
    for (i=0; i<n; i++) {
        data[i]=&REAL(R_data)[i*nc];
    }

//   initialize tree object
    struct node *tree=calloc(MAXSIZE,sizeof(struct node));

//   survival tree, split=Cox PH z^2 stat
    if ( (method == 11) || (method == 25) ) {  // 11: 2 treatments
                                               // 25: 2+ treatments

        maketree_surv(n,nc,y,data,types,ncat,treat,censor,tree,nmin,
                      maxdepth,method);

        free(data);

//       return values to R
        SEXP node_index;
        SEXP node_covar;
        SEXP node_types;
        SEXP node_signs;
        SEXP node_vals;
        SEXP node_depth;
        SEXP node_parent;
        SEXP node_lchild;
        SEXP node_rchild;
        SEXP node_n;
        SEXP node_n0;
        SEXP node_n1;
        SEXP node_r0;
        SEXP node_r1;
        SEXP node_p0;
        SEXP node_p1;
        SEXP node_predict;

        int max=get_number_of_nodes(MAXSIZE,tree);
        PROTECT(node_index=allocMatrix(INTSXP,1,max));
        PROTECT(node_covar=allocMatrix(INTSXP,1,max));
        PROTECT(node_types=allocMatrix(INTSXP,1,max));
        PROTECT(node_signs=allocMatrix(INTSXP,1,max));
        PROTECT(node_vals=allocMatrix(REALSXP,1,max));
        PROTECT(node_depth=allocMatrix(INTSXP,1,max));
        PROTECT(node_parent=allocMatrix(INTSXP,1,max));
        PROTECT(node_lchild=allocMatrix(INTSXP,1,max));
        PROTECT(node_rchild=allocMatrix(INTSXP,1,max));
        PROTECT(node_n=allocMatrix(INTSXP,1,max));
        PROTECT(node_n0=allocMatrix(INTSXP,1,max));
        PROTECT(node_n1=allocMatrix(INTSXP,1,max));
        PROTECT(node_r0=allocMatrix(REALSXP,1,max));
        PROTECT(node_r1=allocMatrix(REALSXP,1,max));
        PROTECT(node_p0=allocMatrix(REALSXP,1,max));
        PROTECT(node_p1=allocMatrix(REALSXP,1,max));
        PROTECT(node_predict=allocMatrix(INTSXP,1,max));

        for ( i=0; i<max; i++ ) {
            INTEGER(node_index)[i]=tree[i].index;
            INTEGER(node_covar)[i]=tree[i].covar;
            INTEGER(node_types)[i]=tree[i].type;
            INTEGER(node_signs)[i]=tree[i].sign;
            REAL(node_vals)[i]=tree[i].val;
            INTEGER(node_depth)[i]=tree[i].depth;
            INTEGER(node_parent)[i]=tree[i].parent;
            INTEGER(node_lchild)[i]=tree[i].lchild;
            INTEGER(node_rchild)[i]=tree[i].rchild;
            INTEGER(node_n)[i]=tree[i].n;
            INTEGER(node_n0)[i]=tree[i].n0;
            INTEGER(node_n1)[i]=tree[i].n1;
            REAL(node_r0)[i]=tree[i].r0;
            REAL(node_r1)[i]=tree[i].r1;
            REAL(node_p0)[i]=tree[i].p0;
            REAL(node_p1)[i]=tree[i].p1;
            INTEGER(node_predict)[i]=tree[i].predict;
        }

        free(tree);

        SEXP Rtree;
        PROTECT(Rtree=allocVector(VECSXP,17));

        SET_VECTOR_ELT(Rtree,0,node_index);
        SET_VECTOR_ELT(Rtree,1,node_covar);
        SET_VECTOR_ELT(Rtree,2,node_types);
        SET_VECTOR_ELT(Rtree,3,node_signs);
        SET_VECTOR_ELT(Rtree,4,node_vals);
        SET_VECTOR_ELT(Rtree,5,node_depth);
        SET_VECTOR_ELT(Rtree,6,node_parent);
        SET_VECTOR_ELT(Rtree,7,node_lchild);
        SET_VECTOR_ELT(Rtree,8,node_rchild);
        SET_VECTOR_ELT(Rtree,9,node_n);
        SET_VECTOR_ELT(Rtree,10,node_n0);
        SET_VECTOR_ELT(Rtree,11,node_n1);
        SET_VECTOR_ELT(Rtree,12,node_r0);
        SET_VECTOR_ELT(Rtree,13,node_r1);
        SET_VECTOR_ELT(Rtree,14,node_p0);
        SET_VECTOR_ELT(Rtree,15,node_p1);
        SET_VECTOR_ELT(Rtree,16,node_predict);

        UNPROTECT(18);  // unprotect last 18 protected R objects

        return(Rtree);
    }

//   diff tree for continuous outcomes
    if ( (method == -1) || (method == 24) ) {  // -1: 2 treatments
                                               // 24: 2+ treatments

        maketree_diff(n,nc,y,data,types,ncat,treat,tree,nmin,
                      maxdepth,method);

        free(data);

//       return values to R
        SEXP node_index;
        SEXP node_covar;
        SEXP node_types;
        SEXP node_signs;
        SEXP node_vals;
        SEXP node_depth;
        SEXP node_parent;
        SEXP node_lchild;
        SEXP node_rchild;
        SEXP node_n;
        SEXP node_n0;
        SEXP node_n1;
        SEXP node_r0;
        SEXP node_r1;
        SEXP node_p0;
        SEXP node_p1;
        SEXP node_predict;

        int max=get_number_of_nodes(MAXSIZE,tree);
        PROTECT(node_index=allocMatrix(INTSXP,1,max));
        PROTECT(node_covar=allocMatrix(INTSXP,1,max));
        PROTECT(node_types=allocMatrix(INTSXP,1,max));
        PROTECT(node_signs=allocMatrix(INTSXP,1,max));
        PROTECT(node_vals=allocMatrix(REALSXP,1,max));
        PROTECT(node_depth=allocMatrix(INTSXP,1,max));
        PROTECT(node_parent=allocMatrix(INTSXP,1,max));
        PROTECT(node_lchild=allocMatrix(INTSXP,1,max));
        PROTECT(node_rchild=allocMatrix(INTSXP,1,max));
        PROTECT(node_n=allocMatrix(INTSXP,1,max));
        PROTECT(node_n0=allocMatrix(INTSXP,1,max));
        PROTECT(node_n1=allocMatrix(INTSXP,1,max));
        PROTECT(node_r0=allocMatrix(REALSXP,1,max));
        PROTECT(node_r1=allocMatrix(REALSXP,1,max));
        PROTECT(node_p0=allocMatrix(REALSXP,1,max));
        PROTECT(node_p1=allocMatrix(REALSXP,1,max));
        PROTECT(node_predict=allocMatrix(INTSXP,1,max));

        for ( i=0; i<max; i++ ) {
            INTEGER(node_index)[i]=tree[i].index;
            INTEGER(node_covar)[i]=tree[i].covar;
            INTEGER(node_types)[i]=tree[i].type;
            INTEGER(node_signs)[i]=tree[i].sign;
            REAL(node_vals)[i]=tree[i].val;
            INTEGER(node_depth)[i]=tree[i].depth;
            INTEGER(node_parent)[i]=tree[i].parent;
            INTEGER(node_lchild)[i]=tree[i].lchild;
            INTEGER(node_rchild)[i]=tree[i].rchild;
            INTEGER(node_n)[i]=tree[i].n;
            INTEGER(node_n0)[i]=tree[i].n0;
            INTEGER(node_n1)[i]=tree[i].n1;
            REAL(node_r0)[i]=tree[i].r0;
            REAL(node_r1)[i]=tree[i].r1;
            REAL(node_p0)[i]=tree[i].p0;
            REAL(node_p1)[i]=tree[i].p1;
            INTEGER(node_predict)[i]=tree[i].predict;
        }

        free(tree);

        SEXP Rtree;
        PROTECT(Rtree=allocVector(VECSXP,17));

        SET_VECTOR_ELT(Rtree,0,node_index);
        SET_VECTOR_ELT(Rtree,1,node_covar);
        SET_VECTOR_ELT(Rtree,2,node_types);
        SET_VECTOR_ELT(Rtree,3,node_signs);
        SET_VECTOR_ELT(Rtree,4,node_vals);
        SET_VECTOR_ELT(Rtree,5,node_depth);
        SET_VECTOR_ELT(Rtree,6,node_parent);
        SET_VECTOR_ELT(Rtree,7,node_lchild);
        SET_VECTOR_ELT(Rtree,8,node_rchild);
        SET_VECTOR_ELT(Rtree,9,node_n);
        SET_VECTOR_ELT(Rtree,10,node_n0);
        SET_VECTOR_ELT(Rtree,11,node_n1);
        SET_VECTOR_ELT(Rtree,12,node_r0);
        SET_VECTOR_ELT(Rtree,13,node_r1);
        SET_VECTOR_ELT(Rtree,14,node_p0);
        SET_VECTOR_ELT(Rtree,15,node_p1);
        SET_VECTOR_ELT(Rtree,16,node_predict);

        UNPROTECT(18);  // unprotect last 18 protected R objects

        return(Rtree);
    }

//   set number of cores for multiple-thread processing
    #ifdef _OPENMP
      int ncores=omp_get_num_procs();
    #else
      int ncores=1;
    #endif
    if ( ntree < ncores ) ncores=ntree;

//   initialize root node
    tree[0].depth=1;
    tree[0].index=1;
    tree[0].n=n;

    if ( method >= 20 ) {  // multi-treat methods

        tree[0].n0=get_min_ntrt(1,n,nc,data,ncat,treat,tree);
        tree[0].n1=tree[0].n0;

        tree[0].predict=get_node_predict_multi(1,n,nc,y,data,ncat,
                                               treat,tree);
        tree[0].r0=0;
        tree[0].r1=0;
        tree[0].p0=0;
        tree[0].p1=0;

    } else {               // binary treat methods

        tree[0].n0=get_subj_treat(1,0,n,nc,data,ncat,treat,tree);
        tree[0].n1=get_subj_treat(1,1,n,nc,data,ncat,treat,tree);

        double r[4];
        tree[0].predict=get_node_predict_y(1,n,nc,y,data,ncat,treat,
                                           r,tree);
        tree[0].r0=r[0];
        tree[0].r1=r[1];
        tree[0].p0=r[2];
        tree[0].p1=r[3];
    }

//   build tree
    int current_node,isplit,max,d;
    struct node *best_split=malloc(2*sizeof(struct node));

//   read .Random.seed for R random number generation
    GetRNGstate();

    for ( i=0; i<MAXSIZE; i++ ) {

        current_node=i+1;

        d=tree[current_node-1].depth+1;
        if ( (maxdepth > 0) && (d > maxdepth) ) break;

        if ( (method == 12) || (method == 13) ||
             (method == 20) || (method == 21) ) { // survival depth 
                                                  // methods

            isplit=get_split_surv(current_node,ntree,n,nc,y,data,
                                  types,ncat,treat,censor,
                                  tree,nmin,nmin2,mtry,maxdepth2,
                                  method,ncores,best_split);

        } else {

            isplit=get_split(current_node,ntree,n,nc,y,data,types,
                             ncat,treat,censor,tree,nmin,nmin2,
                             mtry,maxdepth2,method,ncores,best_split);
        }

        if ( isplit != -7 ) {  // "-7" = there should be no split

            max=get_number_of_nodes(MAXSIZE,tree);
            tree[i].lchild=max+1;
            tree[i].rchild=max+2;

//           add left child node to tree
            tree[max].index=tree[i].lchild;
            tree[max].covar=best_split[0].covar;
            tree[max].type=best_split[0].type;
            tree[max].sign=best_split[0].sign;
            tree[max].val=best_split[0].val;
            tree[max].parent=current_node;
            tree[max].depth=d;
            tree[max].n=best_split[0].n;
            tree[max].n0=best_split[0].n0;
            tree[max].n1=best_split[0].n1;
            tree[max].r0=best_split[0].r0;
            tree[max].r1=best_split[0].r1;
            tree[max].p0=best_split[0].p0;
            tree[max].p1=best_split[0].p1;
            tree[max].predict=best_split[0].predict;

//           add right child node to tree
            tree[max+1].index=tree[i].rchild;
            tree[max+1].covar=best_split[1].covar;
            tree[max+1].type=best_split[1].type;
            tree[max+1].sign=best_split[1].sign;
            tree[max+1].val=best_split[1].val;
            tree[max+1].parent=current_node;
            tree[max+1].depth=d;
            tree[max+1].n=best_split[1].n;
            tree[max+1].n0=best_split[1].n0;
            tree[max+1].n1=best_split[1].n1;
            tree[max+1].r0=best_split[1].r0;
            tree[max+1].r1=best_split[1].r1;
            tree[max+1].p0=best_split[1].p0;
            tree[max+1].p1=best_split[1].p1;
            tree[max+1].predict=best_split[1].predict;
        }

        max=get_number_of_nodes(MAXSIZE,tree);
        if ( (isplit == -7) && (current_node == max) ) break;
    }

//   write .Random.seed for R random number generation
    PutRNGstate();

    free(data);
    free(best_split);

/*
//   print tree
    Rprintf("\nTree from maketree:\n");
    for ( i=0; i<max; i++ ) {

        Rprintf("index=%d  covar=%d  type=%d  sign=%d  val=%f  parent=%d  lchild=%d  rchild=%d  depth=%d  nsubj=%d  ntrt0=%d  ntrt1=%d  besttrt=%d\n",
               tree[i].index,tree[i].covar,tree[i].type,
               tree[i].sign,tree[i].val,tree[i].parent,
               tree[i].lchild,tree[i].rchild,tree[i].depth,
               tree[i].n,tree[i].n0,tree[i].n1,tree[i].predict);
    }
*/

//   return values to R
    SEXP node_index;
    SEXP node_covar;
    SEXP node_types;
    SEXP node_signs;
    SEXP node_vals;
    SEXP node_depth;
    SEXP node_parent;
    SEXP node_lchild;
    SEXP node_rchild;
    SEXP node_n;
    SEXP node_n0;
    SEXP node_n1;
    SEXP node_r0;
    SEXP node_r1;
    SEXP node_p0;
    SEXP node_p1;
    SEXP node_predict;

    max=get_number_of_nodes(MAXSIZE,tree);
    PROTECT(node_index=allocMatrix(INTSXP,1,max));
    PROTECT(node_covar=allocMatrix(INTSXP,1,max));
    PROTECT(node_types=allocMatrix(INTSXP,1,max));
    PROTECT(node_signs=allocMatrix(INTSXP,1,max));
    PROTECT(node_vals=allocMatrix(REALSXP,1,max));
    PROTECT(node_depth=allocMatrix(INTSXP,1,max));
    PROTECT(node_parent=allocMatrix(INTSXP,1,max));
    PROTECT(node_lchild=allocMatrix(INTSXP,1,max));
    PROTECT(node_rchild=allocMatrix(INTSXP,1,max));
    PROTECT(node_n=allocMatrix(INTSXP,1,max));
    PROTECT(node_n0=allocMatrix(INTSXP,1,max));
    PROTECT(node_n1=allocMatrix(INTSXP,1,max));
    PROTECT(node_r0=allocMatrix(REALSXP,1,max));
    PROTECT(node_r1=allocMatrix(REALSXP,1,max));
    PROTECT(node_p0=allocMatrix(REALSXP,1,max));
    PROTECT(node_p1=allocMatrix(REALSXP,1,max));
    PROTECT(node_predict=allocMatrix(INTSXP,1,max));

    for ( i=0; i<max; i++ ) {
        INTEGER(node_index)[i]=tree[i].index;
        INTEGER(node_covar)[i]=tree[i].covar;
        INTEGER(node_types)[i]=tree[i].type;
        INTEGER(node_signs)[i]=tree[i].sign;
        REAL(node_vals)[i]=tree[i].val;
        INTEGER(node_depth)[i]=tree[i].depth;
        INTEGER(node_parent)[i]=tree[i].parent;
        INTEGER(node_lchild)[i]=tree[i].lchild;
        INTEGER(node_rchild)[i]=tree[i].rchild;
        INTEGER(node_n)[i]=tree[i].n;
        INTEGER(node_n0)[i]=tree[i].n0;
        INTEGER(node_n1)[i]=tree[i].n1;
        REAL(node_r0)[i]=tree[i].r0;
        REAL(node_r1)[i]=tree[i].r1;
        REAL(node_p0)[i]=tree[i].p0;
        REAL(node_p1)[i]=tree[i].p1;
        INTEGER(node_predict)[i]=tree[i].predict;
    }

    free(tree);

    SEXP Rtree;
    PROTECT(Rtree=allocVector(VECSXP,17));

    SET_VECTOR_ELT(Rtree,0,node_index);
    SET_VECTOR_ELT(Rtree,1,node_covar);
    SET_VECTOR_ELT(Rtree,2,node_types);
    SET_VECTOR_ELT(Rtree,3,node_signs);
    SET_VECTOR_ELT(Rtree,4,node_vals);
    SET_VECTOR_ELT(Rtree,5,node_depth);
    SET_VECTOR_ELT(Rtree,6,node_parent);
    SET_VECTOR_ELT(Rtree,7,node_lchild);
    SET_VECTOR_ELT(Rtree,8,node_rchild);
    SET_VECTOR_ELT(Rtree,9,node_n);
    SET_VECTOR_ELT(Rtree,10,node_n0);
    SET_VECTOR_ELT(Rtree,11,node_n1);
    SET_VECTOR_ELT(Rtree,12,node_r0);
    SET_VECTOR_ELT(Rtree,13,node_r1);
    SET_VECTOR_ELT(Rtree,14,node_p0);
    SET_VECTOR_ELT(Rtree,15,node_p1);
    SET_VECTOR_ELT(Rtree,16,node_predict);

    UNPROTECT(18);  // unprotect last 18 protected R objects

    return(Rtree);
}
