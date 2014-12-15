#include "mex.h"
#include "matrix.h"

/*function [pos vel] = nearersRelaxing(int segs, double pos, double vel, double masses, double ntimes, char noCollisions, double toReduce)*/
/*USE THIS WITH EXTREME CARE!!!!!!!!!!!!!!!*/
/*POS AND VEL WILL BE MODIFIED IN-PLACE!!!!!!*/
/*   So, if you have copied the variables without modifying them,*/
/*   this will change the copies as well!!!!!*/
/*ALSO: SEGS IS EXPECTED TO BE UINT16 INDEXES STARTING FROM 0!!!!!!*/

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  unsigned short int *segsRoot, *segs, p1, p2, *randomRoot, *random, r;
  double *pos, *vel, *toReduceRoot, *toReduce, *mass, mergeVel;
  double *p1x, *p2x, *p1y, *p2y, *v1x, *v2x, *v1y, *v2y;
  double d1x, d1y, d2, d3, d4x, d4y, half, normx, normy, relvel, med, crfac;
  double v1xr, v1yr, v2xr, v2yr, rn1x, rn1y, rn2x, rn2y, m1, m2, m12;
  unsigned long long magic; unsigned char cosa;
  union {double d; unsigned long long l;} id2, isig;
/*   union {float d; unsigned long l;} id2, isig;*/
  int i, j, nsegs, nps, times;
  unsigned char *noCollisionsRoot, *noCollisions;
  
/*   if (sizeof(unsigned short int)!=2) {*/
/*     mexErrMsgTxt("type (unsigned short int) was expected to be 2 bytes wide, but it isn't. You must change the type!!!!");*/
/*   }*/
/*function [pos vel] = collisionRelaxing(int segs, double pos, double vel, double masses, double ntimes, char noCollisions, double toReduce)*/
  segsRoot         = (unsigned short int *) mxGetPr(prhs[0]);
  pos              = (double *)             mxGetPr(prhs[1]);
  vel              = (double *)             mxGetPr(prhs[2]);
  mass             = (double *)             mxGetPr(prhs[3]);
  times            = (int)(*((double *)     mxGetPr(prhs[4])));
  noCollisionsRoot = (unsigned char *)      mxGetPr(prhs[5]);
  toReduceRoot     = (double *)             mxGetPr(prhs[6]);
  randomRoot       = (unsigned short int *) mxGetPr(prhs[7]);
  
  nsegs = (int)mxGetM(prhs[0]); /*number of rows*/
  nps   = (int)mxGetM(prhs[1]); /*number of rows*/
  
  for (j=0;j<times;++j) {
  segs         = segsRoot;
  toReduce     = toReduceRoot;
  noCollisions = noCollisionsRoot;
  
  for (i=0;i<nsegs;++i) {
    /*p1 = segs(i,1); p2 = segs(i,2);*/
    r = *(randomRoot++);
    p1 = *(segsRoot+r);
    p2 = *(segsRoot+(r+nsegs));
/*     p1  = (*segs);*/
/*     p2  = (*((segs++)+nsegs));*/
    
    m1  = mass[p1];
    m2  = mass[p2];
    m12 = m1+m2;
    
    /*pointers to pos([p1, p2],:)*/
    p1x = pos+p1;
    p2x = pos+p2;
    p1y = p1x+nps;
    p2y = p2x+nps;
    /*pointers to vel([p1, p2],:)*/
    v1x = vel+p1;
    v2x = vel+p2;
    v1y = v1x+nps;
    v2y = v2x+nps;
    
    /*d1 = ss.pos(p1,:) - ss.pos(p2,:);*/
    d1x = (*p1x)-(*p2x);
    d1y = (*p1y)-(*p2y);
    
    /*d2 = realsqrt(sum(d1.*d1, 2));*/
    /*compute the inverse of d2*/
    id2.d = d1x*d1x + d1y*d1y;
    half  = id2.d*0.5;
    id2.l = (0xBFCDD6A18F6A6F55 - id2.l) >> 1;
    /*id2.l = magic - id2.l;//magic - id2.l);// >> cosa;*/
    id2.d = id2.d*(1.5-half*id2.d*id2.d); /*two newton iterations*/
    id2.d = id2.d*(1.5-half*id2.d*id2.d);
    if (id2.d>1e6) { /*(d2<1e-6) {*/
      mexErrMsgTxt("Two objects are nearly coincident. This is likely due to a too long timestep");
    }
    
    /*d3 = (d2-minD(m))./d2 == ;*/
    d3 = 1.0 - (*(toReduce+r))*id2.d;
/*     d3 = 1.0 - (*(toReduce++))*id2.d;*/
    /*d4 = d1*d3/2*RELAXINGFACTOR; RELAXINGFACTOR=0.95, 0.9, etc*/
    /*              1.00/2 = .5*/
    /*              0.95/2 = .475*/
    /*              0.9/2  = .45*/
    /*              0.85/2 = .425*/
    /*              0.8/2  = .4*/
    d4x = d1x*d3*0.475;
    d4y = d1y*d3*0.475;
    /*ss.pos(p2,:) = ss.pos(p2,:) + d4;*/
    /*ss.pos(p1,:) = ss.pos(p1,:) - d4;*/
    (*p1x) -= d4x;
    (*p1y) -= d4y;
    (*p2x) += d4x;
    (*p2y) += d4y;
    
    if (!(*(noCollisionsRoot+r))) {
/*     if (!(*(noCollisions++))) {*/

      /*n1 = d1(1)./d2;*/
      /*n2 = d1(2)./d2;*/
      normx = d1x*id2.d;
      normy = d1y*id2.d;
      /*velocities in the direction of collision (d1), calculated as*/
      /*dot product with the properly oriented unit vector*/
      v1xr = (*v1x)*normx + (*v1y)*normy;
      v2xr = (*v2x)*normx + (*v2y)*normy;

  /*     //vars = ss.vel([p1 p2],:); %rows: [rotatedDPos; rotatedBallDPos;]*/
  /*     //rotated = [vars(:,1).*n1+vars(:,2).*n2, vars(:,2).*n1-vars(:,1).*n2];*/
  /*     v1xr = (*v1x)*normx + (*v1y)*normy;*/
  /*     v2xr = (*v2x)*normx + (*v2y)*normy;*/
  /*     v1yr = (*v1y)*normx - (*v1x)*normy;*/
  /*     v2yr = (*v2y)*normx - (*v2x)*normy;*/
      /*relvel = rotated(2,1)-rotated(1,1);*/
      relvel = v2xr-v1xr;
  /*     if relvel>0 %approachingDPos*/
  /*       med                = (masses(m,1)*rotated(1,1)+masses(m,2)*rotated(2,1))/summasses(m);*/
  /*       crfac              = (1-exp(-relvel))*relvel/summasses(m);*/
  /*       rotated([1 2],1)   = med+masses(m,[2 1])'.*[crfac -crfac]';*/
  /*       rn1                = rotated.*n1;*/
  /*       rn2                = rotated.*n2;*/
  /*       derotated          = [rn1(:,1)-rn2(:,2), rn2(:, 1)+rn1(:,2)];*/
  /*       ss.vel([p1 p2],:)  = derotated([1 2],:); */
  /*     end*/
      if (relvel>0.0) {
    /*       med     = (m1*v1xr + m2*v2xr) / m12;*/
    /*       //-we are going to use a sigmoid function to calculate CR (the coefficient of restitution)*/
    /*       //-that is to say, CR will be a function [0..1] of the relative velocity between the objects*/
    /*       //-to keep things simple, we are going to use a very simple function: f(x) = x/sqrt(1+x^2)*/
    /*       //-so, we are going to calculate a fast inverse square root again!!*/
    /* // //       crfac   = relvel/m12; //use this if the coefficient of restitution must always be 1*/
    /* //       relvel *= relvel;*/
    /* //       isig.d  = 1 + relvel; //remember relvel is in fact relvel*relvel*/
    /* //       half    = isig.d*0.5;*/
    /* //       isig.l  = (0xBFCDD6A18F6A6F55 - isig.l) >> 1;*/
    /* //       isig.d  = isig.d*(1.5-half*isig.d*isig.d); //two newton iterations*/
    /* //       isig.d  = isig.d*(1.5-half*isig.d*isig.d);*/
    /* //       crfac   = relvel*isig.d/m12; //remember relvel is in fact relvel*relvel*/
    /* //       //if you want to use CR=0 always, just don't use crfac*/
    /* //       v1xr    = med+m2*crfac;*/
    /* //       v2xr    = med-m1*crfac;*/
    /*       v1xr = v2xr = (m1*v1xr + m2*v2xr) / m12;*/
    /*       //now, de-rotate velocities*/
    /*       (*v1x) = v1xr*normx - v1yr*normy;*/
    /*       (*v2x) = v2xr*normx - v2yr*normy;*/
    /*       (*v1y) = v1yr*normx + v1xr*normy;*/
    /*       (*v2y) = v2yr*normx + v2xr*normy;*/
crfac   = relvel*0.0/m12;
          med  = (m1*v1xr + m2*v2xr) / m12;
          /*calculate the component to be added, such that the speed in that direction is equaled to the desired result*/
          v1xr = med + m2*crfac - v1xr;
          v2xr = med - m1*crfac - v2xr;
          /*now, add those components*/
          (*v1x) += v1xr*normx;
          (*v2x) += v2xr*normx;
          (*v1y) += v1xr*normy;
          (*v2y) += v2xr*normy;
      }
    
    }
    
  }
  
  }
}