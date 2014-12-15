#include "mex.h"
#include "matrix.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  union {double d; unsigned long long l;} id;
  double *values, half, *output;
  int i;
  size_t nin;

  plhs[0] = mxCreateDoubleMatrix((mwSize)mxGetM(prhs[0]), (mwSize)mxGetN(prhs[0]), mxREAL);
  
  output = (double *) mxGetPr(plhs[0]);
  values = (double *) mxGetPr(prhs[0]);
  nin = mxGetNumberOfElements(prhs[0]);
  
  for (i=0;i<nin;++i) {
/*     id.d          = (*values);*/
    id.d          = (*(values++));
    half          = id.d*0.5;
    id.l          = ( 0xBFCDD6A18F6A6F55 - id.l ) >> 1;
    id.d          = id.d*(1.5-half*id.d*id.d); /*newton iterations*/
    id.d          = id.d*(1.5-half*id.d*id.d);
    id.d          = id.d*(1.5-half*id.d*id.d);
/*     (*(values++)) = id.d;*/
    (*(output++)) = id.d;
  }
}
