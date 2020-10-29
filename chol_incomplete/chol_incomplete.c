/* -------------------------------------------------------------------
 * CHOL_INCOMPLETE
 *
 * Incomplete Cholesky decomposition for dense symmetric positive
 * definite matrix A (n-by-n).
 * Result is n-by-k lower triangular L, s.t. L*L' can be used as
 * low-rank approximation of P'*A*P, P a permutation matrix. The
 * algorithm is the usual Cholesky decomposition, but in each
 * iteration the next column is chosen to maximize the pivot (diagonal
 * element of L). L is returned in LFACT, P is returned in the index
 * vector PIND, meaning that P maps component I to PIND(I).
 *
 * A is defined via a function computing DIAG(A) and a function
 * computing a column of A. A is selected by the int. code ASEL
 * and an arbitrary number of further input arguments (the number
 * and types depend on the ASEL value). In the moment, the following
 * variants for A are implemented:
 * - 0: Covariance matrix for RBF (Gaussian) covariance function
 * - 1: Covariance matrix for squared-exponential covariance function
 *
 * Input:
 * - N:       Number of rows/columns of A
 * - ASEL:    Selects matrix type A (s.a.)
 * - MAXD:    Maximal number of columns of L
 * - PVTHRES: Algorithm terminates if all remaining pivots are
 *            < PVTHRES [optional, def.: 0]
 * - ...:     Arguments to the functions defining A
 *
 * Return:
 * - LFACT:   S.a. (n-by-d)
 * - PIND:    S.a. (n)
 * -------------------------------------------------------------------
 * Matlab MEX Function
 * Author: Matthias Seeger
 * ------------------------------------------------------------------- */

#include <math.h>
#include "mex.h"
#include "mex_helper.h" /* Helper functions */
#include "blas_headers.h"

/*
 * Variants for matrix A
 *
 * Every matrix variant A is assigned a code in 0,...,'numAVars'-1.
 * It is defined by a function computing diag(A) and a function
 * computing a selected column of A. The corr. function pointers are
 * stored in 'funDiagA', 'funColA'.
 * The last two arguments to these functions are 'nrhs', 'prhs' from
 * the MEX function. Arguments from 'addArgs' are additional arguments
 * to these functions.
 *
 * Arguments to DIAGA functions:
 * - DG:  Diagonal returned here
 * - N:   Number of rows/cols of A
 * - ...: Additional arguments
 *
 * Arguments to COLA functions:
 * - J:   Number of column to be returned
 * - COL: Column returned here
 * - N:   Number of rows/cols of A
 * - ...: Additional arguments
 */

typedef void (*typeDiagA)(double*,int,int,const mxArray**);
typedef void (*typeColA)(int,double*,int,int,const mxArray**);

int addArgs=4;
char errMsg[200];

/*
 * Variant 0: Kernel matrix of RBF (Gaussian) covariance function
 *
 * Covariance function:
 *   K(x,y) = C \exp( -(w/2) \sum_{i=0}^{d-1} (x_i-y_i)^2 ) + v_b,
 * - w:   Inverse squared length scale. Positive
 * - C:   Variance parameter. Positive
 * - v_b: Offset parameter. Positive if given (optional)
 *
 * Additional parameters:
 * - DATA:    Data matrix X (n-by-d)
 * - DVEC:    Vector diag(X*X') (n)
 * - HYPPARS: Hyperparameter vector [w; C; v_b] or [w; C]
 */

void radialbfcf_checkargs(int n,int nrhs,const mxArray** prhs)
{
  if (nrhs<addArgs+3) mexErrMsgTxt("Not enough input arguments");
  if (!mxIsDouble(prhs[addArgs]) || mxGetM(prhs[addArgs])!=n)
    mexErrMsgTxt("Wrong additional arguments");
  if (getVecLen(prhs[addArgs+1],"DVEC")!=n)
    mexErrMsgTxt("Wrong additional arguments");
  if (getVecLen(prhs[addArgs+2],"HYPPARS")<2)
    mexErrMsgTxt("Wrong additional arguments");
}

void radialbfcf_diaga(double* dg,int n,int nrhs,const mxArray** prhs)
{
  int i;
  const double* hypp;
  double val;
  const mxArray* hypO;

  radialbfcf_checkargs(n,nrhs,prhs);
  hypO=prhs[addArgs+2];
  hypp=mxGetPr(hypO);
  val=hypp[1];
  if (mxGetM(hypO)*mxGetN(hypO)>2) val+=hypp[2];
  for (i=0; i<n; i++) dg[i]=val;
}

/*
 * a = X x_j - 0.5 d - 0.5 d_j 1, k = C exp(w a) + v_b 1
 */
void radialbfcf_cola(int j,double* col,int n,int nrhs,const mxArray** prhs)
{
  int i,d,farg1;
  double temp,fdarg1,fdarg2,logcval,wval,vbval;
  const double* hypp,*data,*dvec;
  const mxArray* hypO;

  radialbfcf_checkargs(n,nrhs,prhs);
  d=mxGetN(prhs[addArgs]);
  hypO=prhs[addArgs+2];
  hypp=mxGetPr(hypO);
  logcval=log(hypp[1]); wval=hypp[0];
  vbval=(mxGetM(hypO)*mxGetN(hypO)>2)?hypp[2]:0.0;
  data=mxGetPr(prhs[addArgs]);
  dvec=mxGetPr(prhs[addArgs+1]);
  fdarg1=1.0; fdarg2=0.0; farg1=1;
  BLASFUNC(dgemv) ("N",&n,&d,&fdarg1,data,&n,data+j,&n,&fdarg2,col,&farg1);
  fdarg1=-0.5;
  BLASFUNC(daxpy) (&n,&fdarg1,dvec,&farg1,col,&farg1);
  temp=-0.5*dvec[j];
  for (i=0; i<n; i++)
    col[i]=exp(logcval+wval*(col[i]+temp));
  if (vbval!=0.0)
    for (i=0; i<n; i++) col[i]+=vbval;
}

/*
 * Variant 1: Kernel matrix of squared-exponential covariance function
 *
 * Covariance function:
 *   K(x,y) = C \exp( -(1/2) \sum_{i=0}^{d-1} w_i (x_i-y_i)^2 ) + v_b,
 * - w_i: Inverse squared length scales. Positive
 * - C:   Variance parameter. Positive
 * - v_b: Offset parameter. Positive if given (optional)
 *
 * Additional parameters:
 * - DATA:    Data matrix X (n-by-d)
 * - DVEC:    Vector diag(X*W*X'), W=diag(w_i) (n)
 * - HYPPARS: Hyperparameter vector [w_1; ...; w_d; C; v_b] or
 *            [w_1; ...; w_d; C]
 */

void squaredexpcf_checkargs(int n,int nrhs,const mxArray** prhs)
{
  int d;

  if (nrhs<addArgs+3) mexErrMsgTxt("Not enough input arguments");
  if (!mxIsDouble(prhs[addArgs]) || mxGetM(prhs[addArgs])!=n)
    mexErrMsgTxt("Wrong additional arguments");
  d=mxGetN(prhs[addArgs]);
  if (getVecLen(prhs[addArgs+1],"DVEC")!=n)
    mexErrMsgTxt("Wrong additional arguments");
  if (getVecLen(prhs[addArgs+2],"HYPPARS")<d+1)
    mexErrMsgTxt("Wrong additional arguments");
}

void squaredexpcf_diaga(double* dg,int n,int nrhs,const mxArray** prhs)
{
  int i,d;
  const double* hypp;
  double val;
  const mxArray* hypO;

  squaredexpcf_checkargs(n,nrhs,prhs);
  hypO=prhs[addArgs+2];
  hypp=mxGetPr(hypO);
  d=mxGetN(prhs[addArgs]);
  val=hypp[d];
  if (mxGetM(hypO)*mxGetN(hypO)>d+1) val+=hypp[d+1];
  for (i=0; i<n; i++) dg[i]=val;
}

/*
 * a = X W x_j - 0.5 d - 0.5 d_j 1, k = C exp(a) + v_b 1
 */
void squaredexpcf_cola(int j,double* col,int n,int nrhs,const mxArray** prhs)
{
  int i,d,farg1;
  double temp,fdarg1,fdarg2,logcval,vbval;
  const double* hypp,*data,*dvec;
  double* tvec;
  const mxArray* hypO;

  squaredexpcf_checkargs(n,nrhs,prhs);
  d=mxGetN(prhs[addArgs]);
  hypO=prhs[addArgs+2];
  hypp=mxGetPr(hypO);
  logcval=log(hypp[d]);
  vbval=(mxGetM(hypO)*mxGetN(hypO)>d+1)?hypp[d+1]:0.0;
  data=mxGetPr(prhs[addArgs]);
  dvec=mxGetPr(prhs[addArgs+1]);
  tvec=(double*) mxMalloc(d*sizeof(double));
  for (i=0; i<d; i++) tvec[i]=hypp[i]*data[n*i+j]; /* W x_j */
  fdarg1=1.0; fdarg2=0.0; farg1=1;
  BLASFUNC(dgemv) ("N",&n,&d,&fdarg1,data,&n,tvec,&farg1,&fdarg2,col,&farg1);
  fdarg1=-0.5;
  BLASFUNC(daxpy) (&n,&fdarg1,dvec,&farg1,col,&farg1);
  temp=logcval-0.5*dvec[j];
  for (i=0; i<n; i++)
    col[i]=exp(col[i]+temp);
  if (vbval!=0.0)
    for (i=0; i<n; i++) col[i]+=vbval;
}

int numAVars=2;
typeDiagA funDiagA[] = {&radialbfcf_diaga,&squaredexpcf_diaga};
typeColA  funColA[]  = {&radialbfcf_cola,&squaredexpcf_cola};

/* Main function CHOL_INCOMPLETE */

void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
  int asel,i,j,k,l,n,maxd,farg1,farg2,farg3,deb1,deb2;
  double temp,pvthres=0.0,pivot,fdarg1;
  double* lfact,*pvec,*wkvec,*lnew;
  int* pind;
  /* double* debug; */

  /* Read arguments */
  if (nrhs<3)
    mexErrMsgTxt("Not enough input arguments");
  if (nlhs<2)
    mexErrMsgTxt("Not enough output arguments");
  if ((n=getScalInt(prhs[0],"N"))<=0)
    mexErrMsgTxt("Wrong argument N");
  asel=getScalInt(prhs[1],"ASEL");
  if (asel<0 || asel>=numAVars)
    mexErrMsgTxt("Wrong argument ASEL");
  maxd=getScalInt(prhs[2],"MAXD");
  if (maxd<1 || maxd>n) mexErrMsgTxt("MAXD invalid");
  if (nrhs>3) {
    pvthres=getScalar(prhs[3],"PVTHRES");
    if (pvthres<0.0) mexErrMsgTxt("Wrong argument PVTHRES");
  }

  /* Memory allocations */
  plhs[0]=mxCreateDoubleMatrix(n,maxd,mxREAL);
  lfact=mxGetPr(plhs[0]);
  pind=(int*) mxMalloc(n*sizeof(int));
  pvec=(double*) mxMalloc(n*sizeof(double));
  wkvec=(double*) mxMalloc(n*sizeof(double));

  /* Initialization */
  for (i=0; i<n; i++) pind[i]=i;
  fillVec(lfact,n*maxd,0.0);
  (*funDiagA[asel])(pvec,n,nrhs,prhs); /* Compute diag(A) */

  /* Main loop */
  for (k=0; k<maxd; k++) {
    pivot=pvec[k]; l=k;
    for (i=k+1; i<n; i++)
      if ((temp=pvec[i])>pivot) {
	pivot=temp; l=i;
      }
    if (pivot<=0.0)
      mexErrMsgTxt("Matrix A is not positive definite");
    if (sqrt(pivot)<pvthres) break; /* terminate loop */
    j=pind[l];
    /* Swap entries l,k (and rows in L) */
    pind[l]=pind[k]; pind[k]=j;
    pvec[l]=pvec[k]; pvec[k]=pivot;
    if (k>0)
      BLASFUNC(dswap) (&k,lfact+k,&n,lfact+l,&n);
    lnew=lfact+(k*n); /* New column k for L */
    lnew[k]=pivot=sqrt(pivot);
    (*funColA[asel])(j,wkvec,n,nrhs,prhs); /* Col. j of A into 'wkvec' */
    for (i=k+1; i<n; i++) lnew[i]=wkvec[pind[i]]/pivot;
    if (k>0) {
      /* Matrix-vector multiplication */
      temp=-1.0/pivot;
      farg1=n-k-1; farg2=k; farg3=1; fdarg1=1.0;
      BLASFUNC(dgemv) ("N",&farg1,&farg2,&temp,lfact+(k+1),&n,lfact+k,&n,
		       &fdarg1,lnew+(k+1),&farg3);
    }
    /* Update pivot vector */
    for (i=k+1; i<n; i++) {
      temp=lnew[i];
      pvec[i]-=temp*temp;
    }
  }

  /* Deallocate. Return arguments */
  mxFree((void*) pvec); mxFree((void*) wkvec);
  mxSetN(plhs[0],k);
  plhs[1]=mxCreateDoubleMatrix(n,1,mxREAL);
  pvec=mxGetPr(plhs[1]);
  for (i=0; i<n; i++) pvec[i]=(double) (pind[i]+1);
  mxFree((void*) pind);
}
