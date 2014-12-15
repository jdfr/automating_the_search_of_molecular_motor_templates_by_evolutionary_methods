#include "mex.h"
#include "matrix.h"

/*function image = drawLinesAsImageTrueColor(int32 pos, int32 springEnds, int32 imageDimensions, single springColors)*/

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  int      *pos, *spEnds, *imDims;
  int       nlines, npixels, npixels2, nfils, ncols, i;
  int       x0, y0, x1, y1, x, y, *p0, *p1, dx, dy, err, ystep;
  mwSize    imageDims[3];
  float    *img, *spColors, *loc, cr, cg, cb;
  char      steep, doit;
  int       xdesp, ydesp;
/*   int       npoints, *pp;*/
  
  pos          = (int *)          mxGetPr(prhs[0]);
  spEnds       = (int *)          mxGetPr(prhs[1]);
  imDims       = (int *)          mxGetPr(prhs[2]);
  spColors     = (float *)        mxGetPr(prhs[3]);
  
/*   npoints      = (int)mxGetN(prhs[0]); //number of columns*/
  nlines       = (int)mxGetN(prhs[1]); /*number of columns*/
  
  ncols        = imDims[1];
  nfils        = imDims[2];
  npixels      = nfils*ncols;
  npixels2     = npixels << 1;
  imageDims[0] = nfils;
  imageDims[1] = ncols;
  imageDims[2] = 3;
  
  plhs[0]      = mxCreateNumericArray(3, imageDims, mxSINGLE_CLASS, mxREAL);
  img          = (float *)        mxGetPr(plhs[0]);
  
/*   pp=pos;*/
/*   for (i=0;i<npoints;++i) {*/
/*     x = (*(pp++));*/
/*     y = (*(pp++));*/
/*     if ((y>=0) && (x>=0) && (y<nfils) && (x<ncols)) {*/
/*       loc = img+(y+x*nfils);*/
/*       (*loc)            = 1.0f;*/
/*       (*(loc+npixels))  = 1.0f;*/
/*       (*(loc+npixels2)) = 1.0f;*/
/*     }*/
/*   }*/
  
  for (i=0;i<nlines;++i) {
    cr    = *(spColors++);
    cg    = *(spColors++);
    cb    = *(spColors++);
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
            loc               = img+(x+y*nfils);
            (*loc)            = cr;
            (*(loc+npixels))  = cg;
            (*(loc+npixels2)) = cb;
          }
        } else {
          if ((y<nfils) && (x<ncols)) {
            loc               = img+(y+x*nfils);
            (*loc)            = cr;
            (*(loc+npixels))  = cg;
            (*(loc+npixels2)) = cb;
          }
        }
      }
      err -= dy;
      if (err<0) {
        err += dx;
        y   += ydesp;
      }
    }
/*     if (steep) {*/
/*       if (y0<y1) {*/
/*         if (x0<x1) {*/
/*           for (x=x0;x<=x1;x++) {*/
/*             if ((y>=0) && (x>=0) && (y<=nfils) && (x<=ncols)) {*/
/*               loc               = img+(y+x*nfils);*/
/*               (*loc)            = cr;*/
/*               (*(loc+npixels))  = cg;*/
/*               (*(loc+npixels2)) = cb;*/
/*             }*/
/*             err -= dy;*/
/*             if (err<0) {*/
/*               err += dx;*/
/*               y++;*/
/*             }*/
/*           }*/
/*         } else {*/
/*           for (x=x0;x>=x1;x--) {*/
/*             if ((y>=0) && (x>=0) && (y<=nfils) && (x<=ncols)) {*/
/*               loc               = img+(y+x*nfils);*/
/*               (*loc)            = cr;*/
/*               (*(loc+npixels))  = cg;*/
/*               (*(loc+npixels2)) = cb;*/
/*             }*/
/*             err -= dy;*/
/*             if (err<0) {*/
/*               err += dx;*/
/*               y++;*/
/*             }*/
/*           }*/
/*         }*/
/*       } else {*/
/*         if (x0<x1) {*/
/*           for (x=x0;x<=x1;x++) {*/
/*             if ((y>=0) && (x>=0) && (y<=nfils) && (x<=ncols)) {*/
/*               loc               = img+(y+x*nfils);*/
/*               (*loc)            = cr;*/
/*               (*(loc+npixels))  = cg;*/
/*               (*(loc+npixels2)) = cb;*/
/*             }*/
/*             err -= dy;*/
/*             if (err<0) {*/
/*               err += dx;*/
/*               y--;*/
/*             }*/
/*           }*/
/*         } else {*/
/*           for (x=x0;x>=x1;x--) {*/
/*             if ((y>=0) && (x>=0) && (y<=nfils) && (x<=ncols)) {*/
/*               loc               = img+(y+x*nfils);*/
/*               (*loc)            = cr;*/
/*               (*(loc+npixels))  = cg;*/
/*               (*(loc+npixels2)) = cb;*/
/*             }*/
/*             err -= dy;*/
/*             if (err<0) {*/
/*               err += dx;*/
/*               y--;*/
/*             }*/
/*           }*/
/*         }*/
/*       }*/
/*     } else {*/
/*       if (y0<y1) {*/
/*         if (x0<x1) {*/
/*           for (x=x0;x<=x1;x++) {*/
/*             if ((y>=0) && (x>=0) && (x<=nfils) && (y<=ncols)) {*/
/*               loc               = img+(x+y*nfils);*/
/*               (*loc)            = cr;*/
/*               (*(loc+npixels))  = cg;*/
/*               (*(loc+npixels2)) = cb;*/
/*             }*/
/*             err -= dy;*/
/*             if (err<0) {*/
/*               err += dx;*/
/*               y++;*/
/*             }*/
/*           }*/
/*         } else {*/
/*           for (x=x0;x>=x1;x--) {*/
/*             if ((y>=0) && (x>=0) && (x<=nfils) && (y<=ncols)) {*/
/*               loc               = img+(x+y*nfils);*/
/*               (*loc)            = cr;*/
/*               (*(loc+npixels))  = cg;*/
/*               (*(loc+npixels2)) = cb;*/
/*             }*/
/*             err -= dy;*/
/*             if (err<0) {*/
/*               err += dx;*/
/*               y++;*/
/*             }*/
/*           }*/
/*         }*/
/*       } else {*/
/*         if (x0<x1) {*/
/*           for (x=x0;x<=x1;x++) {*/
/*             if ((y>=0) && (x>=0) && (x<=nfils) && (y<=ncols)) {*/
/*               loc               = img+(x+y*nfils);*/
/*               (*loc)            = cr;*/
/*               (*(loc+npixels))  = cg;*/
/*               (*(loc+npixels2)) = cb;*/
/*             }*/
/*             err -= dy;*/
/*             if (err<0) {*/
/*               err += dx;*/
/*               y--;*/
/*             }*/
/*           }*/
/*         } else {*/
/*           for (x=x0;x>=x1;x--) {*/
/*             if ((y>=0) && (x>=0) && (x<=nfils) && (y<=ncols)) {*/
/*               loc               = img+(x+y*nfils);*/
/*               (*loc)            = cr;*/
/*               (*(loc+npixels))  = cg;*/
/*               (*(loc+npixels2)) = cb;*/
/*             }*/
/*             err -= dy;*/
/*             if (err<0) {*/
/*               err += dx;*/
/*               y--;*/
/*             }*/
/*           }*/
/*         }*/
/*       }*/
/*     }*/
  }
  
  
  
  
  
  
/*  function line(x0, x1, y0, y1)*/
/*      boolean steep := abs(y1 - y0) > abs(x1 - x0)*/
/*      if steep then*/
/*          swap(x0, y0)*/
/*          swap(x1, y1)*/
/*      if x0 > x1 then*/
/*          swap(x0, x1)*/
/*          swap(y0, y1)*/
/*      int deltax := x1 - x0*/
/*      int deltay := abs(y1 - y0)*/
/*      int error := deltax / 2*/
/*      int ystep*/
/*      int y := y0*/
/*      if y0 < y1 then ystep := 1 else ystep := -1*/
/*      for x from x0 to x1*/
/*          if steep then plot(x,y) else plot(y,x)*/
/*          error := error - deltay*/
/*          if error < 0 then*/
/*              y := y + ystep*/
/*              error := error + deltax*/
  
}
