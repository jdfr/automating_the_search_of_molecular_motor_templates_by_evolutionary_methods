#include "mex.h"
#include "matrix.h"

/*function image = drawLinesAsBoolImage(int32 pos, int32 springEnds, int32 imageDimensions)*/

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  int      *pos, *spEnds, *imDims;
  int       nlines, nfils, ncols, i;
  int       x0, y0, x1, y1, x, y, *p0, *p1, dx, dy, err, ystep;
  char     *img, *loc;
  char      steep, doit;
  int       xdesp, ydesp;
/*   int       npoints, *pp;*/
  
  pos          = (int *)          mxGetPr(prhs[0]);
  spEnds       = (int *)          mxGetPr(prhs[1]);
  imDims       = (int *)          mxGetPr(prhs[2]);
  
/*   npoints      = (int)mxGetN(prhs[0]); //number of columns*/
  nlines       = (int)mxGetN(prhs[1]); /*number of columns*/
  ncols        = imDims[0];
  nfils        = imDims[1];
  
  plhs[0]      = mxCreateLogicalMatrix(nfils, ncols);
  img          = (char *)        mxGetPr(plhs[0]);
  
  for (i=0;i<nlines;++i) {
    p0    = pos + ( (*(spEnds++)) << 1 );
    p1    = pos + ( (*(spEnds++)) << 1 );
    x0    = *(p0++);
    y0    = *(p0);
    x1    = *(p1++);
    y1    = *(p1);
    steep = abs(y1 - y0) > abs(x1 - x0);
    if (steep) {
      x=x0; x0=y0; y0=x;
      x=x1; x1=y1; y1=x;
    }
    if (x0 > x1) {
      x=x0; x0=x1; x1=x;
      x=y0; y0=y1; y1=x;
    }
    dx    = x1 - x0;
    dy    = abs(y1 - y0);
    err   = dx >> 1;
    y     = y0;
    xdesp = (x0<x1) ? 1 : -1;
    ydesp = (y0<y1) ? 1 : -1;
    for (x=x0;x<=x1;x+=xdesp) {
      if ((y>=0) && (x>=0)) {
        if (steep) {
          if ((x<nfils) && (y<ncols)) {
            (*(img+(x+y*nfils))) = 1;
          }
        } else {
          if ((y<nfils) && (x<ncols)) {
            (*(img+(y+x*nfils))) = 1;
          }
        }
      }
      err -= dy;
      if (err<0) {
        err += dx;
        y   += ydesp;
      }
    }
  }
  
}
