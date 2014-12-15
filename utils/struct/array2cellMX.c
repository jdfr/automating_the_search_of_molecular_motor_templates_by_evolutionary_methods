#include "mex.h"
#include "matrix.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  union {double d; unsigned long long l;} id;
  float *valuesS;
  char *valuesL, *valuesC, mini[2];
  double *valuesD;
  int i;
  size_t nin;
  mxArray *cell, *scalar;
  mwSize dim[2];

  cell = plhs[0] = mxCreateCellMatrix((mwSize)mxGetM(prhs[0]), (mwSize)mxGetN(prhs[0]));
  
  nin = mxGetNumberOfElements(prhs[0]);
  
        valuesD = mxGetPr(prhs[0]);
        for (i=0;i<nin;++i) {
          mxSetCell(cell, i, mxCreateDoubleScalar((*(valuesD++))) );
        }
  
/*  if (nin>0) {
    switch (mxGetClassID(prhs[0])) {
      case mxSINGLE_CLASS:
        valuesS = (float *) mxGetData(prhs[0]);
        for (i=0;i<nin;++i) {
          mxSetCell(cell, i, mxCreateDoubleScalar((double)(*(valuesS++))) );
        }
        break;
      case mxDOUBLE_CLASS:
        valuesD = mxGetPr(prhs[0]);
        for (i=0;i<nin;++i) {
          mxSetCell(cell, i, mxCreateDoubleScalar((*(valuesD++))) );
        }
        break;
      case mxLOGICAL_CLASS:
        valuesL = (char *) mxGetData(prhs[0]);
        for (i=0;i<nin;++i) {
          mxSetCell(cell, i, mxCreateLogicalScalar((*(valuesL++))) );
        }
        break;
      case mxCHAR_CLASS:
        mini[1] = 0;
        valuesC = mxArrayToString(prhs[0]);
        for (i=0;i<nin;++i) {
          *mini = (*(valuesC++));
          mxSetCell(cell, i, mxCreateString(mini) );
        }
        mxFree(valuesC-nin);
        break;
      default:
        mexErrMsgIdAndTxt("array2cellMX:args", "Unsupported class!!!!");
    }
  }*/
  
}
