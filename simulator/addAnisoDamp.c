#include "mex.h"
#include "matrix.h"

/* This mex file is equivalent to this Matlab code:
 *
 * acels(springEnds(:,1),1) = acels(springEnds(:,1),1) + dampSprings.*au(springEnds(:,1));
 * acels(springEnds(:,1),2) = acels(springEnds(:,1),2) + dampSprings.*au(springEnds(:,1));
 * acels(springEnds(:,2),1) = acels(springEnds(:,2),1) + dampSprings.*au(springEnds(:,2));
 * acels(springEnds(:,2),2) = acels(springEnds(:,2),2) + dampSprings.*au(springEnds(:,2));
 *
 *  We do it in a MEX file to be faster
 *
 */
/*  addAnisoDamp(acels, au, dampSprings, uint32(springEnds)); */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  double       *acelsX, *acelsY, *au, auA, auB, dsX, dsY, *dampSprings;
  unsigned int *springEnds, spA, spB, i, nacels, nsprings;
  
  if (nrhs<3) {
    mexErrMsgIdAndTxt("addAnisoDamp:args", "%d is an insufficient number of arguments!!!!", nrhs);
  }
  if (mxGetClassID(prhs[0])!=mxDOUBLE_CLASS) {
    mexErrMsgIdAndTxt("addAnisoDamp:args", "The first argument must be double!!!!");
  }
  if (mxGetClassID(prhs[1])!=mxDOUBLE_CLASS) {
    mexErrMsgIdAndTxt("addAnisoDamp:args", "The second argument must be double!!!!");
  }
  if (mxGetClassID(prhs[2])!=mxDOUBLE_CLASS) {
    mexErrMsgIdAndTxt("addAnisoDamp:args", "The third argument must be double!!!!");
  }
  if (mxGetClassID(prhs[3])!=mxUINT32_CLASS) {
    mexErrMsgIdAndTxt("addAnisoDamp:args", "The fourth argument must be uint32!!!!");
  }
  if (mxGetN(prhs[0])!=2) {
    mexErrMsgIdAndTxt("addAnisoDamp:args", "The first argument must have two columns!!!!");
  }
  if (mxGetN(prhs[1])!=1) {
    mexErrMsgIdAndTxt("addAnisoDamp:args", "The second argument must have one column!!!!");
  }
  if (mxGetN(prhs[2])!=2) {
    mexErrMsgIdAndTxt("addAnisoDamp:args", "The third argument must have two columns!!!!");
  }
  if (mxGetN(prhs[3])!=2) {
    mexErrMsgIdAndTxt("addAnisoDamp:args", "The fourth argument must have two columns!!!!");
  }
  nacels         = (unsigned int)mxGetM(prhs[0]);
  nsprings       = (unsigned int)mxGetM(prhs[3]);
  if (mxGetM(prhs[1])!=nacels) {
    mexErrMsgIdAndTxt("pruneAndSweepC:args", "The first and second arguments must have the same number of rows, but they have %d and %d!!!!", nacels, mxGetM(prhs[1]));
  }
  if (mxGetM(prhs[2])!=nsprings) {
    mexErrMsgIdAndTxt("pruneAndSweepC:args", "The fourth and third arguments must have the same number of rows, but they have %d and %d!!!!", nsprings, mxGetM(prhs[2]));
  }
  
  acelsX         = (double *)       mxGetPr(prhs[0]);
  au             = (double *)       mxGetPr(prhs[1]);
  dampSprings    = (double *)       mxGetPr(prhs[2]);
  springEnds     = (unsigned int *) mxGetPr(prhs[3]);
  acelsY         =        acelsX+nacels;
  
  /* This is a trick to use indexes starting at 1 without having to
   * substract 1 every time */
  acelsX--;
  acelsY--;
  au--;

/* acels(springEnds(:,1),1) = acels(springEnds(:,1),1) + dampSprings.*au(springEnds(:,1));
 * acels(springEnds(:,1),2) = acels(springEnds(:,1),2) + dampSprings.*au(springEnds(:,1));
 * acels(springEnds(:,2),1) = acels(springEnds(:,2),1) + dampSprings.*au(springEnds(:,2));
 * acels(springEnds(:,2),2) = acels(springEnds(:,2),2) + dampSprings.*au(springEnds(:,2));
 * */
  for (i=0; i<nsprings; i++) {
    /* UTTERLY IMPORTANT: we do not check array boundaries when indexing
     *                    acels, so
     *      MAKE DAMN SURE THE BOUNDARIES AREN'T VIOLATED!!!!!!
     * */
    spB = (*(springEnds+nsprings));
    spA = (*(springEnds++));
    auA = *(au+spA);
    auB = *(au+spB);
    dsY = *(dampSprings+nsprings);
    dsX = *(dampSprings++);
    *(acelsX+spA)   += dsX * auA;
    *(acelsY+spA)   += dsY * auA;
    *(acelsX+spB)   += dsX * auB;
    *(acelsY+spB)   += dsY * auB;
  }
  
}
