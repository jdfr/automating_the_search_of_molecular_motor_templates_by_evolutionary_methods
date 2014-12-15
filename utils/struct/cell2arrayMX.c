#include "mex.h"
#include "matrix.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  union {double d; unsigned long long l;} id;
  double *values;
  int i;
  size_t nin;
  mxArray *cell;
  /*THIS FUNCTION ONLY PRODUCES THE EXPECTED BEHAVIOUR WHEN ALL ELEMENTS
    IN THE CELL ARRAY ARE SOUBLE SCALARS!!!! OTHERWISE, IT WILL PRODUCE
    UNEXPECTED RESULTS!!!!!*/

  cell = prhs[0];
  
  plhs[0] = mxCreateDoubleMatrix((mwSize)mxGetM(cell), (mwSize)mxGetN(cell), mxREAL);
  
  values = (double *) mxGetPr(plhs[0]);
  nin = mxGetNumberOfElements(plhs[0]);
  
  for (i=0;i<nin;++i) {
    (*(values++)) = mxGetScalar(mxGetCell(cell, i));
  }
}
