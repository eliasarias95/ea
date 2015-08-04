#include <sdw_obj.h>

/***********************************PUBLIC***********************************/
void sdw_obj::init(int k, double smin, double smax, 
    axis *ax1, axis *ax2, axis *ax3) {
  double ds = 1.0/k;
  int ismin = ceil(smin/ds);
  int ismax = floor(smax/ds);
  _axs = new axis(ismin*ds,ds,1+ismax-ismin);
  _ax1 = ax1;
  _ax2 = ax2;
  _ax3 = ax3;
  _r1min = -1.0;
  _r2min = -1.0;
  _r3min = -1.0;
  _r1max = 1.0;
  _r2max = 1.0;
  _r3max = 1.0;
  _k1min = 10;
  _k2min = 10;
  _k3min = 10;
  _epow = 1.0f;
  _esmooth = 1;
}

//Constructors
sdw_obj::sdw_obj(int k, double smin, double smax,
    axis *ax1) {
  init(k,smin,smax,ax1,NULL,NULL);
}

sdw_obj::sdw_obj(int k, double smin, double smax,
    axis *ax1, axis *ax2) {
  init(k,smin,smax,ax1,ax2,NULL);
}

sdw_obj::sdw_obj(int k, double smin, double smax, 
    axis *ax1, axis *ax2, axis *ax3) {
  init(k,smin,smax,ax1,ax2,ax3);
}

//Setter methods
void sdw_obj::setStrainLimits(double r1min, double r1max) {
  setStrainLimits(r1min,r1max,-1.0,1.0,-1.0,1.0);  
}

void sdw_obj::setStrainLimits(double r1min, double r1max,
                              double r2min, double r2max) {
  setStrainLimits(r1min,r1max,r2min,r2max,-1.0,1.0);  
}

void sdw_obj::setStrainLimits(double r1min, double r1max,
                              double r2min, double r2max,
                              double r3min, double r3max) {
  _r1min = r1min; _r1max = r1max;
  _r2min = r2min; _r2max = r2max;
  _r3min = r3min; _r3max = r3max;
}

void sdw_obj::setErrorSmoothing(int esmooth) {
  _esmooth = esmooth;
}

void sdw_obj::setSmoothness(double d1min) {
  double d2 = (_ax2!=NULL)?_ax2->d:1.0;
  setSmoothness(d1min,10.0*d2);
}

void sdw_obj::setSmoothness(double d1min, double d2min) {
  double d3 = (_ax3!=NULL)?_ax3->d:1.0;
  setSmoothness(d1min,d2min,10.0*d3);
}

void sdw_obj::setSmoothness(double d1min, double d2min, double d3min) {
  double d1 = (_ax1!=NULL)?_ax1->d:1.0;
  double d2 = (_ax2!=NULL)?_ax2->d:1.0;
  double d3 = (_ax3!=NULL)?_ax3->d:1.0;
  _k1min = fmax(1,(int)ceil(d1min/d1));
  _k2min = fmax(1,(int)ceil(d2min/d2));
  _k3min = fmax(1,(int)ceil(d3min/d3));
}

//Other public methods
void sdw_obj::findShifts(float **e, float *s) {
  int ns = _axs->n;
  int n1 = _ax1->n;
  double ds = _axs->d;
  double d1 = _ax1->d;
  int k1min = fmin(_k1min,n1-1);
  vector<int> i1k = subsample(n1,k1min);
  int n1k = i1k.size();
  findShiftsFromErrors(_r1min,_r1max,i1k,_axs,_ax1,e,s);
}

void sdw_obj::findShifts(axis *axf, float *f, axis *axg, float *g, float *s) {
  float **e = (float**)mem_alloc2(_axs->n,_ax1->n,sizeof(float));
  computeErrors(axf,f,axg,g,e);
  findShifts(e,s);
}

void sdw_obj::findShifts(
    axis *axf, float **f, axis *axg, float **g, float **s) {
  int ns = _axs->n;
  int n1 = _ax1->n;
  int n2 = _ax2->n;
  vector<int> k1s = subsample(n1,_k1min);
  vector<int> k2s = subsample(n2,_k2min);
  int nk1 = k1s.size();
  int nk2 = k2s.size();
  
  ucsl_printf("findShifts: smoothing in 1st dimension ...\n");
  float ***ek = (float***)mem_alloc3(ns,nk1,n2,sizeof(float));
  float **e1 = (float**)mem_alloc2(ns,n1,sizeof(float));
  float **df = (float**)mem_alloc2(ns,n1,sizeof(float));
  float **dr = (float**)mem_alloc2(ns,n1,sizeof(float));
  for (int i2=0; i2<n2; ++i2) {
    computeErrors(axf,f[i2],axg,g[i2],e1);
    subsampleErrors(_r1min,_r1max,k1s,_axs,_ax1,e1,df,dr,ek[i2]);
  }
  normalizeErrors(ns,nk1,n2,ek);

  ucsl_printf("findShifts: smoothing in 2nd dimension ...\n");
  float ***ekk = (float***)mem_alloc3(ns,nk1,nk2,sizeof(float));
  smoothErrors2(_r2min,_r2max,k2s,_axs,_ax2,nk1,n2,ek,ekk);
  mem_free2((void***)&e1);
  mem_free2((void***)&dr);
  mem_free2((void***)&df);
  mem_free3((void****)&ek);
  normalizeErrors(ns,nk1,nk2,ekk);

  ucsl_printf("findShifts: finding shifts ...\n");
  int **m = (int**)mem_alloc2(ns,nk1,sizeof(int));
  float **d   = (float**)mem_alloc2(ns,nk1,sizeof(float));
  float **skk = (float**)mem_alloc2(nk1,nk2,sizeof(float));
  for (int ik2=0; ik2<nk2; ++ik2) {
    findShiftsFromSubsampledErrors(
        _r1min,_r1max,k1s,_axs,_ax1,m,d,ekk[ik2],skk[ik2]);
  }
  ucsl_printf("findShifts: interpolating shifts ...\n");
  interpolateShifts(_ax1,_ax2,k1s,k2s,skk,s);
  ucsl_printf("findShifts: ... done\n");
  mem_free3((void****)&ekk);
  mem_free2((void***)&skk);
  mem_free2((void***)&d);
  mem_free2((void***)&m);
}

void sdw_obj::findShifts(
    axis *axf, float ***f, axis *axg, float ***g, float *** s) {
  int ns = _axs->n;
  int n1 = _ax1->n;
  int n2 = _ax2->n;
  int n3 = _ax3->n;
  vector<int> k1s = subsample(n1,_k1min);
  vector<int> k2s = subsample(n2,_k2min);
  vector<int> k3s = subsample(n3,_k3min);
  int nk1 = k1s.size();
  int nk2 = k2s.size();
  int nk3 = k3s.size();
  
  ucsl_printf("findShifts: smoothing in 1st dimension ...\n");
  float ****ek = (float****)mem_alloc4(ns,nk1,n2,n3,sizeof(float));
  float **e1 = (float**)mem_alloc2(ns,n1,sizeof(float));
  float **df = (float**)mem_alloc2(ns,n1,sizeof(float));
  float **dr = (float**)mem_alloc2(ns,n1,sizeof(float));
  int n23 = n2*n3;
  for (int i23=0; i23<n23; ++i23) {
    int i2 = i23%n2;
    int i3 = i23/n2;
    computeErrors(axf,f[i3][i2],axg,g[i3][i2],e1);
    subsampleErrors(_r1min,_r1max,k1s,_axs,_ax1,e1,df,dr,ek[i3][i2]);
  }
  normalizeErrors(ns,nk1,n2,n3,ek);

  ucsl_printf("findShifts: smoothing in 2nd dimension ...\n");
  float ****ekk = (float****)mem_alloc4(ns,nk1,nk2,n3,sizeof(float));
  smoothErrors2(_r2min,_r2max,k2s,_axs,_ax2,nk1,n2,n3,ek,ekk);
  mem_free2((void***)&e1);
  mem_free2((void***)&dr);
  mem_free2((void***)&df);
  mem_free4((void*****)&ek);
  normalizeErrors(ns,nk1,nk2,n3,ekk);

  ucsl_printf("findShifts: smoothing in 3rd dimension ...\n");
  float ****ekkk = (float****)mem_alloc4(ns,nk1,nk2,nk3,sizeof(float));
  smoothErrors3(_r3min,_r3max,k3s,_axs,_ax3,nk1,nk2,n3,ekk,ekkk);
  mem_free4((void*****)&ekk);
  normalizeErrors(ns,nk1,nk2,nk3,ekkk);

  for (int is=0; is<_esmooth-1; ++is) {
    smoothSubsampledErrors(_r1min,_r1max,k1s,_r2min,_r2max,k2s,_r3min,_r3max,k3s,
        _axs,_ax1,_ax2,_ax3,ekkk);
    normalizeErrors(ns,nk1,nk2,nk3,ekkk);
  }

  ucsl_printf("findShifts: finding shifts ...\n");
  int **m = (int**)mem_alloc2(ns,nk1,sizeof(int));
  float **d   = (float**)mem_alloc2(ns,nk1,sizeof(float));
  float ***skk = (float***)mem_alloc3(nk1,nk2,nk3,sizeof(float));
  int nk23 = nk2*nk3;
  for (int ik23=0; ik23<nk23; ++ik23) {
    int ik2 = ik23%nk2;
    int ik3 = ik23/nk2;
    findShiftsFromSubsampledErrors(
        _r1min,_r1max,k1s,_axs,_ax1,m,d,ekkk[ik3][ik2],skk[ik3][ik2]);
  }

  ucsl_printf("findShifts: interpolating shifts ...\n");
  interpolateShifts(_ax1,_ax2,_ax3,k1s,k2s,k3s,skk,s);
  ucsl_printf("findShifts: ... done\n");

  mem_free4((void*****)&ekkk);
  mem_free3((void****)&skk);
  mem_free2((void***)&d);
  mem_free2((void***)&m);
}

/**
 * Returns alignment errors computed for specified sequences.
 * @param sf sampling of 1st dimension for the seqeunce f.
 * @param f array of values for sequence f.
 * @param sg sampling of 1st dimension for the seqeunce g.
 * @param g array of values for sequence g.
 * @return array of alignment errors.
 */
void sdw_obj::computeErrors(
    axis *axf, float *f, axis *axg, float *g, float **e) {
  int ns = _axs->n;
  int ne = _ax1->n;
  float *fi = (float*)mem_alloc(ne,sizeof(float));
  float *gi = (float*)mem_alloc(ne,sizeof(float));
  interp(axf,f,_ax1,fi,0.0f);
  for (int is=0; is<ns; ++is) {
    interp(axg,g,_ax1,gi,_axs->get_val(is));
    for (int ie=0; ie<ne; ++ie) {
      e[ie][is] = error(fi[ie],gi[ie]);
    }
  }
  mem_free((void**)&fi);
  mem_free((void**)&gi);
}

//is there overhead with using the fmin and fmax methods from math?
void sdw_obj::normalizeErrors(int n1, int n2, float **e) {
  float emin = e[0][0];
  float emax = e[0][0];
  for (int ie=0; ie<n2; ++ie) {
    for (int is=0; is<n1; ++is) {
      float ei = e[ie][is];
      if (ei<emin) emin = ei;
      if (ei>emax) emax = ei;
    }
  }
  shiftAndScale(emin,emax,n1,n2,e);
}

/***********************************PRIVATE***********************************/
float BIG = 9999999999.9f;

void sdw_obj::interp(axis *axf, float *f, axis *axs, float *fi, float shift) {
  int nf = axf->n;
  int ns = axs->n;
  float df = axf->d;
  float ds = axs->d;
  float w;
  float fsamp;
  int isamp;
  for (int is=0; is<ns; ++is) {
    isamp = fsamp = (is*ds)/df - shift;
    w = fsamp - isamp;
    if (fsamp<0)
      fi[is] = f[0];
    else if (fsamp>=0 && (fsamp+1)<ns)
      fi[is] = (1.0f-w)*f[isamp] + w*f[isamp+1];
    else
      fi[is] = f[ns-1];
  }
}

void sdw_obj::fill(float val, float *x, int nx) {
  for (int i=0; i<nx; ++i)
    x[i] = val;
}

/**
 * Returns an approximately uniformly-sampled subset of indices in [0,n).
 * Indices in the subset are chosen to be approximately uniform, with the
 * difference between consecutive indices not less than the specified
 * minimum increment kmin. Because the first and last indices 0 and n-1 are
 * included in the subset, n must be greater than the minimum increment
 * kmin.
 * @param n number of indices in the set {0,1,2,...,n-1}.
 * @param kmin minimum increment between indices in the subset.
 * @return array of indices in the subset.
 */
vector<int> sdw_obj::subsample(int n, int kmin) {
  if (kmin>=n)
    kmin = n-1;
  int m = 1+(n-1)/kmin;
  double d = static_cast<double>(n-1)/static_cast<double>(m-1);
  vector<int> j(m);
  for (int i=0; i<m; ++i)
    j[i] = static_cast<int>(i*d+0.5);
  return j;
}

float sdw_obj::error(float f, float g) {
  //return pow(abs(f-g),_epow);
  return abs(f-g);
}

void sdw_obj::shiftAndScale(float emin, float emax, int n1, int n2, float **e) {
  float eshift = emin;
  float escale = (emax>emin)?1.0f/(emax-emin):1.0f;
  for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      e[i2][i1] = (e[i2][i1]-eshift)*escale;
    }
  }
}

void sdw_obj::shiftAndScale(
    float emin, float emax, int n1, int n2, int n3, float ***e) {
  float eshift = emin;
  float escale = (emax>emin)?1.0f/(emax-emin):1.0f;
  for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        e[i3][i2][i1] = (e[i3][i2][i1]-eshift)*escale;
      }
    }
  }
}

void sdw_obj::shiftAndScale(
    float emin, float emax, int n1, int n2, int n3, int n4, float ****e) {
  float eshift = emin;
  float escale = (emax>emin)?1.0f/(emax-emin):1.0f;
  for (int i4=0; i4<n4; ++i4) {
    for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          e[i4][i3][i2][i1] = (e[i4][i3][i2][i1]-eshift)*escale;
        }
      }
    }
  }
}

void sdw_obj::smoothErrors2(double r2min, double r2max, vector<int> k2s,
    axis *axs, axis *axe, int nk1, int n2, float ***e, float ***es) {
  int ns = axs->n;
  int nk2 = k2s.size();
  float **e2  = (float**)mem_alloc2(ns,n2,sizeof(float));// errors at index ik1
  float **df  = (float**)mem_alloc2(ns,n2,sizeof(float));// errors at index ik1
  float **dr  = (float**)mem_alloc2(ns,n2,sizeof(float));// errors at index ik1
  float **es2 = (float**)mem_alloc2(ns,nk2,sizeof(float));
  for (int ik1=0; ik1<nk1; ++ik1) {
    for (int i2=0; i2<n2; ++i2) {
      memcpy(e2[i2],e[i2][ik1],ns*sizeof(float));
    }
    subsampleErrors(r2min,r2max,k2s,axs,axe,e2,df,dr,es2);
    for (int ik2=0; ik2<nk2; ++ik2)
      memcpy(es[ik2][ik1],es2[ik2],ns*sizeof(float));
  }
  mem_free2((void***)&dr);
  mem_free2((void***)&df);
  mem_free2((void***)&e2);
  mem_free2((void***)&es2);
}

void sdw_obj::smoothErrors2(
    double r2min, double r2max, vector<int> k2s,
    axis *axs, axis *axe, int nk1, int n2, int n3, float ****e, float ****es) {
  for (int i3=0; i3<n3; ++i3)
    smoothErrors2(r2min,r2max,k2s,axs,axe,nk1,n2,e[i3],es[i3]);
}

void sdw_obj::smoothErrors3(
    double r3min, double r3max, vector<int> k3s,
    axis *axs, axis *axe, int nk1, int nk2, int n3, float ****e, float ****es) {
  int ns = axs->n;
  int nk3 = k3s.size();
  float **e3  = (float**)mem_alloc2(ns,n3,sizeof(float));//smooth errors at i1,i2
  float **df  = (float**)mem_alloc2(ns,n3,sizeof(float));
  float **dr  = (float**)mem_alloc2(ns,n3,sizeof(float));
  float **es3 = (float**)mem_alloc2(ns,nk3,sizeof(float));
  for (int ik1=0; ik1<nk1; ++ik1) {
    for (int ik2=0; ik2<nk2; ++ik2) {
      for (int i3=0; i3<n3; ++i3)
        memcpy(e3[i3],e[i3][ik2][ik1],ns*sizeof(float));
      subsampleErrors(r3min,r3max,k3s,axs,axe,e3,df,dr,es3);
      for (int ik3=0; ik3<nk3; ++ik3)
        memcpy(es[ik3][ik2][ik1],es3[ik3],ns*sizeof(float));
    }
  }
  mem_free2((void***)&dr);
  mem_free2((void***)&df);
  mem_free2((void***)&e3);
  mem_free2((void***)&es3);
}

void sdw_obj::smoothSubsampledErrors(
    double rmin, double rmax, vector<int> kes,
    axis *axs, axis *axe, float **e) {
  int ns = axs->n;
  int nke = kes.size();
  float **ef  = (float**)mem_alloc2(ns,nke,sizeof(float));
  float **er  = (float**)mem_alloc2(ns,nke,sizeof(float));
  accumulateSubsampled( 1,rmin,rmax,kes,axs,axe,e,ef,NULL);
  accumulateSubsampled(-1,rmin,rmax,kes,axs,axe,e,er,NULL);
  float scale = 1.0f/nke;
  for (int ike=0; ike<nke; ++ike) {
    int ke = kes[ike];
    for (int is=0; is<ns; ++is) {
      e[ike][is] = scale*(ef[ike][is]+er[ike][is]-e[ike][is]);
    }
  }
  mem_free2((void***)&er);
  mem_free2((void***)&ef);
}

void sdw_obj::smoothSubsampledErrors(
    double r1min, double r1max, vector<int> k1s,
    double r2min, double r2max, vector<int> k2s,
    axis *axs, axis *ax1, axis *ax2, float ***e) {
  int ns = axs->n;
  int nk1 = k1s.size();
  int nk2 = k2s.size();
  float **e2 = (float**)mem_alloc2(ns,nk2,sizeof(float));
  for (int ik1=0; ik1<nk1; ++ik1) {
    for (int ik2=0; ik2<nk2; ++ik2)
      memcpy(e2[ik2],e[ik2][ik1],ns*sizeof(float));
    smoothSubsampledErrors(r2min,r2max,k2s,axs,ax2,e2);
    for (int ik2=0; ik2<nk2; ++ik2)
      memcpy(e[ik2][ik1],e2[ik2],ns*sizeof(float));
  }
  for (int ik2=0; ik2<nk2; ++ik2) {
    smoothSubsampledErrors(r1min,r1max,k1s,axs,ax1,e[ik2]);
  }
}

void sdw_obj::smoothSubsampledErrors(
    double r1min, double r1max, vector<int> k1s,
    double r2min, double r2max, vector<int> k2s,
    double r3min, double r3max, vector<int> k3s,
    axis *axs, axis *ax1, axis *ax2, axis *ax3, float ****e) {
  int ns = axs->n;
  int nk1 = k1s.size();
  int nk2 = k2s.size();
  int nk3 = k3s.size();
  float **e3 = (float**)mem_alloc2(ns,nk3,sizeof(float));
  for (int ik1=0; ik1<nk1; ++ik1) {
    for (int ik2=0; ik2<nk2; ++ik2) {
      for (int ik3=0; ik3<nk3; ++ik3)
        memcpy(e3[ik3],e[ik3][ik2][ik1],ns*sizeof(float));
      smoothSubsampledErrors(r3min,r3max,k3s,axs,ax3,e3);
      for (int ik3=0; ik3<nk3; ++ik3)
        memcpy(e[ik3][ik2][ik1],e3[ik3],ns*sizeof(float));
    }
  }
  for (int ik3=0; ik3<nk3; ++ik3) {
    smoothSubsampledErrors(
        r1min,r1max,k1s,r2min,r2max,k2s,axs,ax1,ax2,e[ik3]);
  }
}

/**
 * Subsamples alignment errors.
 * @param rmin lower bound on strain.
 * @param rmax upper bound on strain.
 * @param kes array[nke] of indices for quasi-uniform subsampling.
 * @param ss uniform sampling of ns shifts.
 * @param se uniform sampling of ne errors.
 * @param e input array[ne][ns] of alignment errors.
 * @return array[nke][ns] of subsampled errors.
 */
void sdw_obj::subsampleErrors(
    double rmin, double rmax, vector<int> kes, axis *axs, axis *axe,
    float **e, float **df, float **dr, float **d) {
  int ns =  axs->n;
  int ne =  axe->n;
  int nke = kes.size();
  accumulate( 1,rmin,rmax,kes,axs,axe,e,df,NULL);
  accumulate(-1,rmin,rmax,kes,axs,axe,e,dr,NULL);
  //d = df;
  float scale = 1.0f/ne;
  for (int ike=0; ike<nke; ++ike) {
    int ke = kes[ike];
    for (int is=0; is<ns; ++is) {
      d[ike][is] = scale*(df[ike][is]+dr[ike][is]-e[ke][is]);
    }
  }
}

/**
 * Returns shifts found by backtracking with precomputed moves.
 * @param kes array[nke] of indices for quasi-uniform subsampling.
 * @param ss uniform sampling of ns shifts.
 * @param se uniform sampling of ne errors.
 * @param d array[ns] of last forward accumulated errors.
 * @param m array[nke][ns] of moves, changes in shifts.
 * @return array[ne] of shifts.
 */
void sdw_obj::backtrackForShifts(
    vector<int> kes, axis *axs, axis *axe, float *d, int **m, float *ske) {
  int nke = kes.size();
  int ns = axs->n;
  int ne = axe->n;
  int ike = nke-1;
  float dmin = BIG;
  int imin = -1;
  for (int is=0; is<ns; ++is) {
    if (d[is]<dmin) {
      dmin = d[is];
      imin = is;
    }
  }
  int is = imin;
  ske[ike] = axs->get_val(is);
  for (--ike; ike>=0; --ike) {
    is += m[ike+1][is];
    ske[ike] = axs->get_val(is);
  }
}

void sdw_obj::normalizeErrors(int n1, int n2, int n3, float ***e) {
  float emin = BIG;
  float emax = -1.0f*BIG;
  for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        float ei = e[i3][i2][i1];
        if (ei<emin) emin = ei;
        if (ei>emax) emax = ei;
      }
    }
  }
  shiftAndScale(emin,emax,n1,n2,n3,e);
}

void sdw_obj::normalizeErrors(int n1, int n2, int n3, int n4, float ****e) {
  float emin = BIG;
  float emax = -1.0f*BIG;
  for (int i4=0; i4<n4; ++i4) {
    for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          float ei = e[i4][i3][i2][i1];
          if (ei<emin) emin = ei;
          if (ei>emax) emax = ei;
        }
      }
    }
  }
  shiftAndScale(emin,emax,n1,n2,n3,n4,e);
}

/**
 * Accumulates alignment errors in forward or reverse direction.
 * Does not subsample the accumulated errors.
 * @param dir direction, 1 for forward, -1 for reverse.
 * @param rmin lower bound on strain.
 * @param rmax upper bound on strain.
 * @param me number of errors e summed per accumulated error d
 * @param ss uniform sampling of ns shifts.
 * @param se uniform sampling of ne errors.
 * @param e input array[ne][ns] of alignment errors.
 * @param d output array[ne][ns] of accumulated errors.
 */
void sdw_obj::accumulate(
    int dir, double rmin, double rmax, int me, 
    axis *axs, axis *axe, float **e, float **d) {
  int ns = axs->n;
  int ne = axe->n;
  double ds = axs->d;
  double de = axe->d;
  int ied = dir>0?1:-1;
  int ieb = dir>0?0:ne-1;
  int iee = dir>0?ne:-1;
  for (int is=0; is<ns; ++is)
    d[ieb][is] = e[ieb][is];
  static float *dprev = NULL;
  if (!dprev) dprev = (float*)mem_alloc(ns,sizeof(float));
  
  for (int ie=ieb+ied; ie!=iee; ie+=ied) {
    int je = fmax(0,fmin(ne-1,ie-dir*me));
    int ke = ie-je;
    int msmin,msmax;
    if (ke>0) {
      msmin = static_cast<int>( ceil(-rmax*ke*de/ds));
      msmax = static_cast<int>(floor(-rmin*ke*de/ds));
    } else {
      msmin = static_cast<int>( ceil(-rmin*ke*de/ds));
      msmax = static_cast<int>(floor(-rmax*ke*de/ds));
    }
    fill(BIG,d[ie],ns);
    for (int ms=msmin; ms<=msmax; ++ms) {
      int islo = fmax(0,-ms);
      int ishi = fmin(ns,ns-ms);
      for (int is=islo; is<ishi; ++is)
        dprev[is] = d[ie-ied][is+ms];
      updateSumsOfErrors(ie,je,ms,e,dprev,d[ie],NULL,ns);
    }
  }
}

/**
 * Accumulates alignment errors in forward or reverse direction.
 * @param dir direction, 1 for forward, -1 for reverse.
 * @param rmin lower bound on strain.
 * @param rmax upper bound on strain.
 * @param kes array[nke] of indices for quasi-uniform subsampling.
 * @param ss uniform sampling of ns shifts.
 * @param se uniform sampling of ne errors.
 * @param e input array[ne][ns] of alignment errors.
 * @param d output array[nke][ns] of accumulated errors.
 * @param m output array[nke][ns] of minimizing moves; or null.
 */
void sdw_obj::accumulate(
    int dir, double rmin, double rmax, vector<int> kes, 
    axis *axs, axis *axe, float **e, float **d, int **m) {
  int ns = axs->n;
  double ds = axs->d;
  double de = axe->d;
  int nke = kes.size();
  int iked = dir>0?1:-1;
  int ikeb = dir>0?0:nke-1;
  int ikee = dir>0?nke:-1;
  for (int is=0; is<ns; ++is)
    d[ikeb][is] = e[kes[ikeb]][is];
  static float *dprev = NULL;
  if (!dprev) dprev = (float*)mem_alloc(ns,sizeof(float));

  for (int ike=ikeb+iked; ike!=ikee; ike+=iked) {
    int je = kes[ike-iked];
    int ie = kes[ike];
    int me = ie-je;
    int msmin,msmax;
    if (me>0) {
      msmin = static_cast<int>( ceil(-rmax*me*de/ds));
      msmax = static_cast<int>(floor(-rmin*me*de/ds));
    } else {
      msmin = static_cast<int>( ceil(-rmin*me*de/ds));
      msmax = static_cast<int>(floor(-rmax*me*de/ds));
    }
    if (msmin>msmax) {
      ucsl_printf("ie=%d",ie," je=%d",je," me=%d",me,
          " msmin=%d",msmin," msmax=%d",msmax);
    }
    assert(msmin<=msmax);
    fill(BIG,d[ike],ns);
    for (int ms=msmin; ms<=msmax; ++ms) {
      int islo = fmax(0,-ms);
      int ishi = fmin(ns,ns-ms);
      for (int is=islo; is<ishi; ++is)
        dprev[is] = d[ike-iked][is+ms];
      if (m!=NULL) {
        updateSumsOfErrors(ie,je,ms,e,dprev,d[ike],m[ike],ns);
        //cout << "m!=NULL\n";
      }
      else {
        updateSumsOfErrors(ie,je,ms,e,dprev,d[ike],NULL,ns);
        //cout << "m==NULL\n";
      }
    }
  }
}

/**
 * Accumulates subsampled errors in forward or reverse direction.
 * @param dir direction, 1 for forward, -1 for reverse.
 * @param rmin lower bound on strain.
 * @param rmax upper bound on strain.
 * @param kes array[nke] of indices for quasi-uniform subsampling.
 * @param ss uniform sampling of ns shifts.
 * @param se uniform sampling of ne errors.
 * @param e input array[nke][ns] of subsampled errors.
 * @param d output array[nke][ns] of accumulated errors.
 * @param m output array[nke][ns] of minimizing moves; or NULL.
 */
void sdw_obj::accumulateSubsampled(
    int dir, double rmin, double rmax, vector<int> kes,
    axis *axs, axis *axe, float **e, float **d, int **m) {
  int ns = axs->n;
  double ds = axs->d;
  double de = axe->d;
  int nke = kes.size();
  int iked = dir>0?1:-1;
  int ikeb = dir>0?0:nke-1;
  int ikee = dir>0?nke:-1;
  for (int is=0; is<ns; ++is)
    d[ikeb][is] = e[ikeb][is];
  static float *dprev = NULL;
  if (!dprev) dprev = (float*)mem_alloc(ns,sizeof(float));  

  for (int ike=ikeb+iked; ike!=ikee; ike+=iked) {
    memcpy(dprev,d[ike-iked],ns*sizeof(float));
    int me = kes[ike]-kes[ike-iked];
    int msmin,msmax;
    if (me>0) {
      msmin = (int) ceil(-rmax*me*de/ds);
      msmax = (int)floor(-rmin*me*de/ds);
    } else {
      msmin = (int) ceil(-rmin*me*de/ds);
      msmax = (int)floor(-rmax*me*de/ds);
    }
    for (int is=0; is<ns; ++is) {
      float dmin = BIG;
      int mmin = -1;
      for (int ms=msmin; ms<=msmax; ++ms) {
        int js = is+ms;
        if (0<=js && js<ns) {
          float dj = dprev[js];
          if (dj<dmin) {
            dmin = dj;
            mmin = ms;
          }
        }
      }
      d[ike][is] = dmin+e[ike][is];
      if (m!=NULL) {
        m[ike][is] = mmin;
      }
    }
  }
}

void sdw_obj::findShiftsFromErrors(double rmin, double rmax, 
    vector<int> kes, axis *axs, axis *axe, float **e, float *s) {
  int nke = kes.size();
  int ns = axs->n;
  float** d = (float**)mem_alloc2(ns,nke,sizeof(float));
  int** m = (int**)mem_alloc2(ns,nke,sizeof(int));
  accumulate(1,rmin,rmax,kes,axs,axe,e,d,m);
  float *ske = (float*)mem_alloc(nke,sizeof(float));
  backtrackForShifts(kes,axs,axe,d[nke-1],m,ske);
  interpolateShifts(axe,kes,ske,s);
}

void sdw_obj::findShiftsFromSubsampledErrors(
    double rmin, double rmax, vector<int> kes, axis *axs, axis *axe, int **m, 
    float **d, float **e, float *s) {
  int nke = kes.size();
  int ns = axs->n;
  accumulateSubsampled(1,rmin,rmax,kes,axs,axe,e,d,m);
  backtrackForShifts(kes,axs,axe,d[nke-1],m,s);
}

void sdw_obj::interpolateShifts(
    axis *ax1, vector<int> k1s, float *sk, float *s) {
  int n1 = ax1->n;
  int nk1 = k1s.size();
  double d1 = ax1->d;
  double dk1 = k1s[1] - k1s[0];
  float w1;
  float x1;
  int isamp1, isamp1p;
  for (int i1=0; i1<n1; ++i1) {
    isamp1 = x1 = i1*d1/dk1;
    if (x1<0 || (x1+1)>=nk1) {isamp1p = isamp1;}
    else {isamp1p = isamp1+1;}
    w1 = x1 - isamp1;
    s[i1] = (1.0f-w1)*sk[isamp1] + w1*sk[isamp1p];
  }
}

void sdw_obj::interpolateShifts(axis *ax1, axis *ax2, vector<int> k1s, 
    vector<int> k2s, float **skk, float **s) {
  int n1 = ax1->n;
  int n2 = ax2->n;
  double d1 = ax1->d;
  double d2 = ax2->d;
  int nk1 = k1s.size();
  int nk2 = k2s.size();
  double dk1 = k1s[1] - k1s[0];
  double dk2 = k2s[1] - k2s[0];

  // Interpolate.
  float w1, w2;
  float x1, x2;
  int isamp1, isamp1p, isamp2, isamp2p;
  for (int i2=0; i2<n2; ++i2) {
    isamp2 = x2 = i2*d2/dk2;
    if (x2<0 || (x2+1)>=nk2) {isamp2p = isamp2;}
    else {isamp2p = isamp2+1;}
    w2 = x2 - isamp2;
    for (int i1=0; i1<n1; ++i1) {
      isamp1 = x1 = i1*d1/dk1;
      if (x1<0 || (x1+1)>=nk1) {isamp1p = isamp1;}
      else {isamp1p = isamp1+1;}
      w1 = x1 - isamp1;
      s[i2][i1] = (1.0f-w2)*(1.0f-w1)*skk[isamp2 ][isamp1 ] + 
                         w2*(1.0f-w1)*skk[isamp2p][isamp1 ] +
                         (1.0f-w2)*w1*skk[isamp2 ][isamp1p] +
                                w2*w1*skk[isamp2p][isamp1p];
    }
  }
}

void sdw_obj::interpolateShifts(axis *ax1, axis *ax2, axis *ax3, 
    vector<int> k1s, vector<int> k2s, vector<int> k3s,
    float ***skk, float ***s) {
  int n2 = ax2->n;
  int n3 = ax3->n;
  int nk1 = k1s.size();
  int nk2 = k2s.size();
  int nk3 = k3s.size();

  float **sk23 = (float**)mem_alloc2(nk2,nk3,sizeof(float));
  float ***s23 = (float***)mem_alloc3(nk1,n2,n3,sizeof(float));
  for (int ik1=0; ik1<nk1; ++ik1) {
    for (int ik3=0; ik3<nk3; ++ik3) {
      for (int ik2=0; ik2<nk2; ++ik2) {
        sk23[ik3][ik2] = skk[ik3][ik2][ik1];
      }
    }
    interpolateShifts(ax2,ax3,k2s,k3s,sk23,s23[ik1]);
  }

  for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
      interpolateShifts(ax1,k1s,s23[i3][i2],s[i3][i2]);
    }
  }
  mem_free2((void***)&sk23);
  mem_free3((void****)&s23);
}

/**
 * Updates sums of errors for one shift between two error sample indices.
 * Computes sums of errors along linear trajectories according to:
 * <pre>
 * d[is] += sum from ke=ie to ke!=je of e[ke][is+(ie-ke)*ms/(ie-je)]
 * </pre>
 * Where the complicated last subscript (typically) is not an integer, this
 * method uses linear interpolation of the alignment errors e. After the
 * sums of errors have been computed for all shift indices is, this method
 * updates the minimum sum of errors and the corresponding change in shift.
 * @param ie error sample index at which to begin sum.
 * @param je error sample index at which to end (not) sum.
 * @param ms change in shift at error sample index je, not in sum.
 * @param e[ne][ns] input array of alignment errors.
 * @param d[ns] input/output array in which to accumulate errors.
 * @param dmin[ns] input/output array of minimum accumulated errors.
 * @param mmin[ns] input/output array of minimizing moves; or null.
 */
void sdw_obj::updateSumsOfErrors(int ie, int je, int ms, float **e, 
    float *d, float *dmin, int *mmin, int ns) {
  int islo = fmax(0,-ms); // update only for is >= islo
  int ishi = fmin(ns,ns-ms); // update only for is < ishi
  for (int is=islo; is<ishi; ++is)
    d[is] += e[ie][is];
  int me = ie-je;
  int de = me>0?-1:1;
  if (ms==0) { // if no shift, no interpolation required
    for (int ke=ie+de; ke!=je; ke+=de) {
      for (int is=islo; is<ishi; ++is) {
        d[is] += e[ke][is];
      }
    }
  } else { // else, use linear interpolation of errors
    float r = static_cast<float>(ms)/static_cast<float>(me); // strain
    for (int ke=ie+de; ke!=je; ke+=de) {
      float sk = r*(ie-ke);
      int ks = (int)sk;
      if (sk<0.0f) --ks;
      int ksa = ks+1;
      int ksb = ks;
      float wsa = sk-ks;
      float wsb = 1.0f-wsa;
      for (int is=islo; is<ishi; ++is)
        d[is] += wsa*e[ke][is+ksa]+wsb*e[ke][is+ksb];
    }
  }
  for (int is=islo; is<ishi; ++is) {
    if (d[is]<dmin[is]) {
      dmin[is] = d[is];
      if (mmin!=NULL)
        mmin[is] = ms;
    }
  }
}
