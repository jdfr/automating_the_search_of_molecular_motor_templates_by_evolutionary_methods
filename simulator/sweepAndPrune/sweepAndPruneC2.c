#include "mex.h"
#include "matrix.h"

/* Do the heavy part of the sweep and prune algorithm: 
 *   -given a sorted list of opens and ends of objects in two dimensions
 *   -generate a NUMBER_OBJECTS X NUMBER_OBJECTS logical matrix indicating
 *    what objects overlap in BOTH dimensions.
 *   -The opens are indexes 1:NUMBER_OBJECTS.
 *   -The ends  are indexes NUMBER_OBJECTS+1:NUMBER_OBJECTS*2.
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  mxArray  *openmx;
  int      *idxs, idx, idx2, *open, *opentop,N, N2,np, ii, i, j, *openpos, *open2;
  char     *candidatos, *endcand, *cand;
  
  if (nrhs<2) {
    mexErrMsgIdAndTxt("pruneAndSweepC:args", "%d is an insufficient number of arguments!!!!", nrhs);
  }
  if (mxGetClassID(prhs[0])!=mxINT32_CLASS) {
    mexErrMsgIdAndTxt("pruneAndSweepC:args", "The first argument must be int32!!!!");
  }
  
  N              = (int)(mxGetN(prhs[0])*mxGetM(prhs[0]));
  N2             = (int)(mxGetN(prhs[1])*mxGetM(prhs[1]));
  
  if (N!=N2) {
    mexErrMsgIdAndTxt("pruneAndSweepC:args", "Both arguments must have the same number of elements, but first has %d and second has %d!!!!", N, N2);
  }
  
  np             = N>>1; /* (N/2) == (N>>1) */
  plhs[0]        = mxCreateLogicalMatrix(np, np);
  candidatos     = (char *) mxGetPr(plhs[0]);
  
  openmx         = mxCreateNumericMatrix(np, 1, mxINT32_CLASS, mxREAL);
  open           = (int *) mxGetPr(openmx);
  
  for (ii=0; ii<2; ii++) {
    idxs           = (int *)  mxGetPr(prhs[ii]);
    opentop        = open-1;
    for (i=0; i<N; i++) {
      idx = *(idxs++);
      if (idx>np) {
        /* if idx>np, then position(idx) = pos(idx-np)-rad(idx-np), so we 
         * must add it to the open set */
        *(++opentop) = idx-np;
      } else {
        /* if idx<=np, then position(idx) = pos(idx)+rad(idx), so we must
         * erase it from the open set and register as candidates all the
         * remaining indexes in the open set */

        /* first, find its position in the open set */
        openpos = -1;
        for (openpos=open; openpos<=opentop; openpos++) {
          if ( (*openpos)==idx ) {
            break;
          }
        }

        /* check weird error condition */
        if (openpos>opentop) {
          mexErrMsgIdAndTxt("pruneAndSweepC:cagada", "%d was not found in the open set!!!!", idx);
        }

        /*shrink the open set */
        for (open2=openpos; open2<opentop; open2++) {
          *open2 = *(open2+1);
        }

        /*decrease open set number*/
        opentop--;

        /* register as candidates the pairs of this idx and all remaining
         * indexes in the open set */
        for (openpos=open; openpos<=opentop; openpos++) {
          idx2 = *openpos;
          if (idx>idx2) {
            (*(candidatos+(idx -1+(idx2-1)*np)))++;
          } else {
            (*(candidatos+(idx2-1+(idx -1)*np)))++;
          }
        }
      }
    }
    if (opentop>=open) {
      mexErrMsgIdAndTxt("pruneAndSweepC:cagada", "the open set is not empty at exit but has %d items!!!!", (opentop-open)+1);
    }
  }
  
  candidatos     = (char *) mxGetPr(plhs[0]);
  endcand        = candidatos+(np*np);
  for (cand=candidatos; cand<endcand; cand++) {
    /* this is a fast way to make AND: only candidates with both (i.e., 2)
     * marks will be true */
    (*cand) >>= 1;
  }
  
  mxDestroyArray(openmx);
  
  
}
