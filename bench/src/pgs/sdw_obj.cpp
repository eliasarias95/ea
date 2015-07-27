#include <sdw_obj.h>

/***********************************PUBLIC***********************************/
void sdw_obj::init(int k, double smin, double smax, 
    axis *ax1, axis *ax2, axis *ax3) {
  double ds = (ax1->d)/k;
  int ismin = ceil(smin/ds);
  int ismax = floor(smin/ds);
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
  //_si = new sinc_interp();
  //_si.setExtrapolation(sinc_interp::Extrapolation::CONSTANT);
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
  for (int i2=0; i2<n2; ++i2) {
    computeErrors(axf,f[i2],axg,g[i2],e1);
    subsampleErrors(_r1min,_r1max,k1s,_axs,_ax1,e1,ek[i2]);
  }
  normalizeErrors(ns,nk1,n2,ek);

  ucsl_printf("findShifts: smoothing in 2nd dimension ...\n");
  float ***ekk = (float***)mem_alloc3(ns,nk1,nk2,sizeof(float));
  smoothErrors2(_r2min,_r2max,k2s,_axs,_ax2,nk1,n2,ek,ekk);
  normalizeErrors(ns,nk1,nk2,ekk);

  ucsl_printf("findShifts: finding shifts ...\n");
  float **skk = (float**)mem_alloc2(nk1,nk2,sizeof(float));
  for (int ik2=0; ik2<nk2; ++ik2) {
    findShiftsFromSubsampledErrors(
        _r1min,_r1max,k1s,_axs,_ax1,ekk[ik2],skk[ik2]);
  }
  ucsl_printf("findShifts: interpolating shifts ...\n");
  interpolateShifts(_ax1,_ax2,k1s,k2s,skk,s);
  ucsl_printf("findShifts: ... done\n");
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
  //ucsl_printf("in compute errors\n");
  int ns = _axs->n;
  int ne = _ax1->n;
  int nf =  axf->n;
  int ng =  axg->n;
  float *fi = (float*)mem_alloc(ne,sizeof(float));
  float *gi = (float*)mem_alloc(ne,sizeof(float));
  //_si.interpolate(axf,f,_ax1,fi);
  float sum = 0.0f;
  for (int is=0; is<ns; ++is) {
    //_si.interpolate(
    //    ng,axg->d,axg->o,g,
    //    ne,_ax1->d,_ax1->o+_axs->get_val(is),gi);
    for (int ie=0; ie<ne; ++ie) {
      e[ie][is] = error(f[ie],g[ie]);
    }
  }
}

//is there overhead with using the fmin and fmax methods from math?
void sdw_obj::normalizeErrors(int n1, int n2, float **e) {
  ucsl_printf("in normalize errors\n");
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

void sdw_obj::fill(double val, float *x, int nx) {
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
  return pow(abs(f-g),_epow);
}

void sdw_obj::shiftAndScale(float emin, float emax, int n1, int n2, float **e) {
  ucsl_printf("in shiftAndScale errors\n");
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
  ucsl_printf("in shiftAndScale errors\n");
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
  //ucsl_printf("in shiftAndScale errors\n");
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
  //float[][][] es = new float[nk2][nk1][ns]; // smoothed errors
  for (int ik1=0; ik1<nk1; ++ik1) {
    float **e2 = (float**)mem_alloc2(ns,n2,sizeof(float));// errors at index ik1
    for (int i2=0; i2<n2; ++i2)
      memcpy(e2[i2],e[i2][ik1],ns*sizeof(float));
    float **es2 = (float**)mem_alloc2(ns,nk2,sizeof(float));
    subsampleErrors(r2min,r2max,k2s,axs,axe,e2,es2);
    for (int ik2=0; ik2<nk2; ++ik2)
      memcpy(es[ik2][ik1],es2[ik2],ns*sizeof(float));
  }
}

void sdw_obj::smoothErrors2(
    double r2min, double r2max, vector<int> k2s,
    axis *axs, axis *axe, float ****e, float ****es) {

}

void sdw_obj::smoothErrors3(
    double r3min, double r3max, vector<int> k3s,
    axis *axs, axis *axe, float ****e, float ****es) {

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
void sdw_obj::subsampleErrors(double rmin, double rmax, 
    vector<int> kes, axis *axs, axis *axe, float **e, float **d) {
  //ucsl_printf("in subsample errors\n");
  int ns = axs->n;
  int ne = axe->n;
  int nke = kes.size();
  float **df = (float**)mem_alloc2(ns,nke,sizeof(float));
  float **dr = (float**)mem_alloc2(ns,nke,sizeof(float));
  accumulate( 1,rmin,rmax,kes,axs,axe,e,df,NULL);
  accumulate(-1,rmin,rmax,kes,axs,axe,e,dr,NULL);
  d = df; // is this correct??
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
  float dmin = (float)HUGE_VAL;
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
  float emin = (float)HUGE_VAL;
  float emax = (float)HUGE_VAL*-1.0f;
  ucsl_printf("in normalize errors\n");
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
  float emin = (float)HUGE_VAL;
  float emax = (float)HUGE_VAL*-1.0f;
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
  vector<float> dprev(ns);
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
    fill(HUGE_VAL,d[ie],ns);
    for (int ms=msmin; ms<=msmax; ++ms) {
      int islo = fmax(0,-ms);
      int ishi = fmin(ns,ns-ms);
      for (int is=islo; is<ishi; ++is)
        dprev[is] = d[ie-ied][is+ms];
      updateSumsOfErrors(ie,je,ms,e,dprev,d[ie],NULL);
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
  vector<float> dprev(ns);
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
      ucsl_printf("In here.\n");
      //ucsl_printf("ie="+ie+" je="+je+" me="+me+" msmin="+msmin+" msmax="+msmax);
    }
    assert(msmin<=msmax);
    fill(HUGE_VAL,d[ike],ns);
    for (int ms=msmin; ms<=msmax; ++ms) {
      int islo = fmax(0,-ms);
      int ishi = fmin(ns,ns-ms);
      for (int is=islo; is<ishi; ++is)
        dprev[is] = d[ike-iked][is+ms];
      if (m!=NULL)
        updateSumsOfErrors(ie,je,ms,e,dprev,d[ike],m[ike]);
      else
        updateSumsOfErrors(ie,je,ms,e,dprev,d[ike],NULL);
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
 * @param m output array[nke][ns] of minimizing moves; or null.
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
  float *dprev = (float*)mem_alloc(ns,sizeof(float));  
  for (int is=0; is<ns; ++is)
    d[ikeb][is] = e[ikeb][is];
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
      float dmin = (float)HUGE_VAL;
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
        d[ike][is] = dmin+e[ike][is];
        if (m!=NULL)
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

void sdw_obj::findShiftsFromSubsampledErrors(double rmin, double rmax,
    vector<int> kes, axis *axs, axis *axe, float **e, float *s) {
  //ucsl_printf("in findShiftsFromSubsampledErrors\n");
  int nke = kes.size();
  int ns = axs->n;
  float** d = (float**)mem_alloc2(ns,nke,sizeof(float));
  int** m = (int**)mem_alloc2(ns,nke,sizeof(int));
  accumulateSubsampled(1,rmin,rmax,kes,axs,axe,e,d,m);
  backtrackForShifts(kes,axs,axe,d[nke-1],m,s);
}

void sdw_obj::interpolateShifts(
    axis *ax1, vector<int> k1s, float *sk, float *s) {
  int n1 = ax1->n;
  int nk1 = k1s.size();
  float* xk1 = (float*)mem_alloc(nk1,sizeof(float));
  for (int jk1=0; jk1<nk1; ++jk1)
    xk1[jk1] = ax1->get_val(k1s[jk1]);
  //Some interpolation routine
  //CubicInterpolator ci = makeInterpolator1(xk1,uk);
  for (int j1=0; j1<n1; ++j1) {
    float x1 = ax1->get_val(j1);
    //Some interpolation routine
    //s[j1] = ci.interpolate(x1);
  }
}

void sdw_obj::interpolateShifts(axis *ax1, axis *ax2, vector<int> k1s, 
    vector<int> k2s, float **skk, float **s) {
  int n1 = ax1->n;
  int n2 = ax2->n;
  int nk1 = k1s.size();
  int nk2 = k2s.size();

  // Coarse sampling of 1st and 2nd dimensions.
  float* xk1 = (float*)mem_alloc(nk1,sizeof(float));
  for (int jk1=0; jk1<nk1; ++jk1)
    xk1[jk1] = ax1->get_val(k1s[jk1]);
  float* xk2 = (float*)mem_alloc(nk2,sizeof(float));
  for (int jk2=0; jk2<nk2; ++jk2)
    xk2[jk2] = ax2->get_val(k2s[jk2]);

  // Interpolate.
  //Some interpolation routine
  //BicubicInterpolator2 bc = new BicubicInterpolator2(
  //  BicubicInterpolator2.Method.MONOTONIC,
  //  BicubicInterpolator2.Method.SPLINE,
  //  xk1,xk2,skk);
  for (int j2=0; j2<n2; ++j2) {
    float x2 = ax2->get_val(j2);
    for (int j1=0; j1<n1; ++j1) {
      float x1 = ax1->get_val(j1);
      //Some interpolation routine
      //s[j2][j1] = bc.interpolate(x1,x2);
    }
  }
}

void sdw_obj::interpolateShifts(axis *ax1, axis *ax2, axis *ax3, 
    vector<int> k1s, vector<int> k2s, vector<int> k3s,
    float ***skk, float ***s) {
  
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
    vector<float> d, float *dmin, int *mmin) {
  int ns = d.size();
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
