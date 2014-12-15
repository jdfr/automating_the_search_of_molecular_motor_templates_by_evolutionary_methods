#include <string.h> /* needed for memcpy() */
#include "mex.h"
#include "matrix.h"

/* Do the heavy part of the sweep and prune algorithm: 
 *   -given a sorted list of opens and ends of objects in two dimensions
 *   -generate a NUMBER_OBJECTS X NUMBER_OBJECTS logical matrix indicating
 *    what objects overlap in BOTH dimensions.
 *   -The opens are indexes 1:NUMBER_OBJECTS.
 *   -The ends  are indexes NUMBER_OBJECTS+1:NUMBER_OBJECTS*2.
 */

/*        candidatos          = sweepAndPruneWithSetsCAll(sap.sets, sap.mode, idxsx, idxsy[, idxsz]);
 *        candidatos          = sweepAndPruneWithSetsCAll(idxsx, idxsy, sap.sets, sap.mode);
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  mxArray  *openmx, *candidatosmx, *idxsmx;
  int      *idxs, idx, idx2, *open, *opentop, *openpos, *open2, i1, i2, ii1, ii2;
  unsigned int *istartA, *icandidatosA, *istartB, *icandidatosB, nic, numout, ii, i, np, nd, N;/*, N2,N3;*/
  unsigned char     *sets, notusesets, modesets;
  signed   char     *candidatos;
  
  /*if (nrhs<2) {
    mexErrMsgIdAndTxt("pruneAndSweepC:args", "%d is an insufficient number of arguments!!!!", nrhs);
  }
  if (mxGetClassID(prhs[2])!=mxINT32_CLASS) {
    mexErrMsgIdAndTxt("pruneAndSweepC:args", "The third argument must be int32!!!!");
  }*/
  
  N              = (int)(mxGetN(prhs[3])*mxGetM(prhs[3]));
  /*N2             = (int)(mxGetN(prhs[3])*mxGetM(prhs[3]));*/
  np             = N>>1; /* (N/2) == (N>>1) */
  notusesets     = mxIsEmpty(prhs[0]);/*((nrhs<=2)) || (((unsigned int)(mxGetN(prhs[0])*mxGetM(prhs[0])))!=np);*/
  if (!notusesets) {
    sets         = (char *) mxGetPr(prhs[0]);
    modesets     = (*((char *) mxGetPr(prhs[1])));
  }
  nd = nrhs-3;/*number of dimensions*/
  
  /*if (N!=N2) {
    mexErrMsgIdAndTxt("pruneAndSweepC:args", "Both arguments must have the same number of elements, but first has %d and second has %d!!!!", N, N2);
  }*/
  
  candidatosmx   = prhs[2];/*mxCreateLogicalMatrix(np, np);*/
  if (mxIsEmpty(candidatosmx)) {
    candidatosmx = mxCreateLogicalMatrix(np, np);
  } else {
    candidatosmx = mxDuplicateArray(candidatosmx);
  }
  candidatos     = (signed char *) mxGetPr(candidatosmx);
  
  openmx         = mxCreateNumericMatrix(np+3, 1, mxINT32_CLASS, mxREAL);
  open           = (int *) mxGetPr(openmx);

  i1           = np*10/2; /*(nd<3) ? (np*10/2) : ();*/
  i2           = (np*(np-1)/2);
  numout       = (i1<i2) ? i1 : i2;
/*   plhs[0]      = mxCreateNumericMatrix(numout, 1, mxUINT32_CLASS, mxREAL);*/
/*   plhs[1]      = mxCreateNumericMatrix(numout, 1, mxUINT32_CLASS, mxREAL);*/
  icandidatosA = istartA = (unsigned int *) mxMalloc(numout*sizeof(unsigned int));/*mxGetPr(plhs[0]);*/
  icandidatosB = istartB = (unsigned int *) mxMalloc(numout*sizeof(unsigned int));/*mxGetPr(plhs[1]);*/
  nic            = 0;
  
  for (ii=0; ii<nd; ii++) {
    idxsmx         = prhs[ii+3];
    idxs           = (int *)  mxGetPr(idxsmx);
    if ( (mxGetN(idxsmx)*mxGetM(idxsmx)) != N ) {
      mexErrMsgIdAndTxt("pruneAndSweepC:args", "the arguments after the first two must have the same number of elements, but first has %d and %d has %d!!!!", N, ii+3, mxGetN(idxsmx)*mxGetM(idxsmx));
    }
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
        /*openpos = -1;*/
        for (openpos=open; openpos<=opentop; openpos++) {
          if ( (*openpos)==idx ) {
            break;
          }
        }

        /* check weird error condition */
        if (openpos>opentop) {
          mexErrMsgIdAndTxt("pruneAndSweepC:cagada", "%d was not found in the open set!!!!. Pass: %d, index %d, number of points: %d", idx, ii, i, np);
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
            i1 = (ii1=idx)  - 1;
            i2 = (ii2=idx2) - 1;
          } else {
            i1 = (ii1=idx2) - 1;
            i2 = (ii2=idx)  - 1;
          }
          if ( ((++(*(candidatos+(i1+i2*np))))==nd)) {
            /*if notusesets  => do it*/
            /*if modesets==0 => do it if they are from two different sets*/
            /*if modesets==1 => do it if they are from the same      set*/
            /*/*mexPrintf("Mira %d, %d: %d, %d, %d, %d\n", ii1, ii2, notusesets, modesets, (sets[i1]==sets[i2]), modesets+(sets[i1]==sets[i2]));*/
            if ( notusesets || ((modesets + (sets[i1]==sets[i2]))!=1) ) {
              *(icandidatosA++) = ii1;
              *(icandidatosB++) = ii2;
              if ((++nic)==numout) {
                /*mexPrintf("Double matrix size: %d to %d\n", numout, numout*2);*/
                istartA = mxRealloc(istartA, numout*sizeof(unsigned int)*2);
                istartB = mxRealloc(istartB, numout*sizeof(unsigned int)*2);
                icandidatosA = istartA+nic;
                icandidatosB = istartB+nic;
                numout      *= 2;
              }
            }
          }
        }
      }
    }
    if (opentop>=open) {
      mexErrMsgIdAndTxt("pruneAndSweepC:cagada", "the open set is not empty at exit but has %d items!!!!", (opentop-open)+1);
    }
  }
  
  /*make output matrix*/
  /*/*mexPrintf("make output matrix with %d rows\n", nic);*/
  istartA = mxRealloc(istartA, nic*sizeof(unsigned int)*2);
  memcpy(istartA+nic, istartB, nic*sizeof(unsigned int));
  mxFree(istartB);
  plhs[0] = mxCreateNumericMatrix(0, 2, mxUINT32_CLASS, mxREAL);
  mxSetM(plhs[0], nic);
  mxFree(mxGetData(plhs[0]));
  mxSetData(plhs[0], istartA);
  /*
  nic = icandidatosA-istartA;
  if (nic>numout) {
    mexErrMsgIdAndTxt("pruneAndSweepC:cagada", "We have estimated the number of overlappings to be at most %d, but it was larger (%d)!!!the open set is not empty at exit but has %d items!!!!", numout, nic);
  }
  */
  
  mxDestroyArray(openmx);
  mxDestroyArray(candidatosmx);
  
  plhs[1] = mxCreateNumericMatrix(1, 1, mxUINT32_CLASS, mxREAL);
  (*(unsigned int *) mxGetPr(plhs[1])) = nic;
}
