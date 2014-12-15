#include "mex.h"
#include "matrix.h"

/* This mex file is equivalent to this Matlab code:
 *
 *    limits             = [pos; pos];
 *    limits(1:np,1)     = limits(1:np,    1) + rad;
 *    limits(np+1:end,1) = limits(np+1:end,1) - rad;
 *    limits(1:np,2)     = limits(1:np,    2) + rad;
 *    limits(np+1:end,2) = limits(np+1:end,2) - rad;
 *
 *  We do it in a MEX file to be faster
 *
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  double   *pos, *posw, *rad, *radw, *limits;
  int      i, np;
  
  if (nrhs<2) {
    mexErrMsgIdAndTxt("pruneAndSweepC:args", "%d is an insufficient number of arguments!!!!", nrhs);
  }
  if (mxGetClassID(prhs[0])!=mxDOUBLE_CLASS) {
    mexErrMsgIdAndTxt("pruneAndSweepC:args", "The first argument must be double!!!!");
  }
  if (mxGetClassID(prhs[1])!=mxDOUBLE_CLASS) {
    mexErrMsgIdAndTxt("pruneAndSweepC:args", "The second argument must be double!!!!");
  }
  if (mxGetN(prhs[0])!=2) {
    mexErrMsgIdAndTxt("pruneAndSweepC:args", "The first argument must have two columns!!!!");
  }
  np             = mxGetM(prhs[0]);
  if (mxGetM(prhs[1])!=np) {
    mexErrMsgIdAndTxt("pruneAndSweepC:args", "The first and second arguments must have the same number of rows, but they have %d and %d!!!!", np, mxGetM(prhs[1]));
  }
  
  pos            = (double *) mxGetPr(prhs[0]);
  rad            = (double *) mxGetPr(prhs[1]);
  plhs[0]        = mxCreateNumericMatrix(np*2, 2, mxDOUBLE_CLASS, mxREAL);
  limits         = (double *) mxGetPr(plhs[0]);

  posw           = pos;
  radw           = rad;
  for (i=0; i<np; i++) {
    *(limits++) = *(posw++) + *(radw++);
  }
  posw           = pos;
  radw           = rad;
  for (i=0; i<np; i++) {
    *(limits++) = *(posw++) - *(radw++);
  }
  posw           = pos+np;
  radw           = rad;
  for (i=0; i<np; i++) {
    *(limits++) = *(posw++) + *(radw++);
  }
  posw           = pos+np;
  radw           = rad;
  for (i=0; i<np; i++) {
    *(limits++) = *(posw++) - *(radw++);
  }
  
}
