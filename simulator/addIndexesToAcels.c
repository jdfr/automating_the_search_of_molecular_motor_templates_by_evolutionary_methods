#include "mex.h"
#include "matrix.h"

/* This mex file is equivalent to this Matlab code:
 *
 *  for k=size(units, 1):-1:1
 *    acels(indexes(k,1),:) = acels(indexes(k,1),:)+toadd(k,:);
 *    acels(indexes(k,2),:) = acels(indexes(k,2),:)-toadd(k,:);
 *  end
 *
 *  We do it in a MEX file to be faster
 *
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  double       *acelsX, *acelsY, *acelsZ, *toaddX, *toaddY, *toaddZ;
  unsigned int *indexesA, *indexesB, i, nacels, ni, nd, nd2;
  
  /*if (nrhs<3) {
    mexErrMsgIdAndTxt("addIndexesToAcels:args", "%d is an insufficient number of arguments!!!!", nrhs);
  }
  if (mxGetClassID(prhs[0])!=mxDOUBLE_CLASS) {
    mexErrMsgIdAndTxt("addIndexesToAcels:args", "The first argument must be double!!!!");
  }
  if (mxGetClassID(prhs[1])!=mxDOUBLE_CLASS) {
    mexErrMsgIdAndTxt("addIndexesToAcels:args", "The second argument must be double!!!!");
  }
  if (mxGetClassID(prhs[2])!=mxUINT32_CLASS) {
    mexErrMsgIdAndTxt("addIndexesToAcels:args", "The third argument must be uint32!!!!");
  }
  if (mxGetN(prhs[0])!=2) {
    mexErrMsgIdAndTxt("addIndexesToAcels:args", "The first argument must have two columns!!!!");
  }
  if (mxGetN(prhs[1])!=2) {
    mexErrMsgIdAndTxt("addIndexesToAcels:args", "The second argument must have two columns!!!!");
  }
  if (mxGetN(prhs[2])!=2) {
    mexErrMsgIdAndTxt("addIndexesToAcels:args", "The third argument must have two columns!!!!");
  }*/
  nacels         = (unsigned int)mxGetM(prhs[0]);
  ni             = (unsigned int)mxGetM(prhs[2]);
  nd2            = mxGetN(prhs[0])<3;
  /*if (mxGetM(prhs[1])!=ni) {
    mexErrMsgIdAndTxt("addIndexesToAcels:args", "The third and second arguments must have the same number of rows, but they have %d and %d!!!!", ni, mxGetM(prhs[1]));
  }*/
  
  acelsX         = (double *)       mxGetPr(prhs[0]);
  toaddX         = (double *)       mxGetPr(prhs[1]);
  indexesA       = (unsigned int *) mxGetPr(prhs[2]);
  acelsY         =   acelsX+nacels;
  toaddY         =   toaddX+ni;
  indexesB       = indexesA+ni;
  
  /* This is a trick to use indexes starting at 1 without having to
   * substract 1 every time */
  acelsX--;
  acelsY--;

/*  for k=size(units, 1):-1:1
 *    acels(prox(k,1),:) = acels(prox(k,1),:)+units(k,:);
 *    acels(prox(k,2),:) = acels(prox(k,2),:)-units(k,:);
 *  end
 * */
  if (nd2) {
    for (i=0; i<ni; i++) {
      /* UTTERLY IMPORTANT: we do not check array boundaries when indexing
       *                    acels, so
       *      MAKE DAMN SURE THE BOUNDARIES AREN'T VIOLATED!!!!!!
       * */
      *(acelsX+(*(indexesA)))   += *(toaddX);
      *(acelsY+(*(indexesA++))) += *(toaddY);
      *(acelsX+(*(indexesB)))   -= *(toaddX++);
      *(acelsY+(*(indexesB++))) -= *(toaddY++);
    }
  } else {
    acelsZ  = acelsY+nacels;
    toaddZ  = toaddY+ni;
    for (i=0; i<ni; i++) {
      /* UTTERLY IMPORTANT: we do not check array boundaries when indexing
       *                    acels, so
       *      MAKE DAMN SURE THE BOUNDARIES AREN'T VIOLATED!!!!!!
       * */
      *(acelsX+(*(indexesA)))   += *(toaddX);
      *(acelsY+(*(indexesA)))   += *(toaddY);
      *(acelsZ+(*(indexesA++))) += *(toaddZ);
      *(acelsX+(*(indexesB)))   -= *(toaddX++);
      *(acelsY+(*(indexesB)))   -= *(toaddY++);
      *(acelsZ+(*(indexesB++))) -= *(toaddZ++);
    }
  }
  
}
