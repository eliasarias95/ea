#include <sdw_obj.h>

/*****************PRIVATE*****************/
//Need to fix this method.
void sdw_obj::trace(std::string s) {
  std::cout << s << "\n";
}

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
static std::vector<int> sdw_obj::subsample(int n, int kmin) {
  if (kmin>=n)
    kmin = n-1;
  int m = 1+(n-1)/kmin;
  double d = static_cast<double>(n-1)/static_cast<double>(m-1);
  std::vector<int> j(m);
  for (int i=0; i<m; ++i)
    j[i] = static_cast<int>(i*d+0.5);
  return j;
}

float sdw_obj::error(float f, float g) {
  return pow(std::abs(f-g),_epow);
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
static void sdw_obj::subsampleErrors(double rmin, double rmax, 
    std::vector<int> *kes, axis axs, axis axe, float **e, float **d) {
  int ns = axs->n;
  int ne = axe->n;
  int nke = kes.size();
  float **df = static_cast<float**>(mem_alloc2(ns,nke,sizeof(float)));
  float **dr = static_cast<float**>(mem_alloc2(ns,nke,sizeof(float)));
  accumulate( 1,rmin,rmax,kes,ss,se,e,df,NULL);
  accumulate(-1,rmin,rmax,kes,ss,se,e,dr,NULL);
  float **d = df; // is this correct??
  float scale = 1.0f/ne;
  for (int ike=0; ike<nke; ++ike) {
    int ke = kes[ike];
    for (int is=0; is<ns; ++is) {
      d[ike][is] = scale*(df[ike][is]+dr[ike][is]-e[ke][is]);
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
static void sdw_obj::accumulate(
    int dir, double rmin, double rmax, std::vector<int> kes, 
    axis axs, axis axe, float **e, float **d, int **m) {
  int ns = axs->n;
  double ds = axs->d;
  double de = axe->d;
  int nke = kes.size();
  int iked = dir>0?1:-1;
  int ikeb = dir>0?0:nke-1;
  int ikee = dir>0?nke:-1;
  for (int is=0; is<ns; ++is)
    d[ikeb][is] = e[kes[ikeb]][is];
  std::vector<float> dprev(ns);
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
      trace("ie="+ie+" je="+je+" me="+me+" msmin="+msmin+" msmax="+msmax);
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
static void sdw_obj::accumulate(
    int dir, double rmin, double rmax, int me, 
    axis axs, axis axe, float **e, float **d) {
  int ns = axs->n;
  int ne = axe->n;
  double ds = axs->d;
  double de = axe->d;
  int ied = dir>0?1:-1;
  int ieb = dir>0?0:ne-1;
  int iee = dir>0?ne:-1;
  for (int is=0; is<ns; ++is)
    d[ieb][is] = e[ieb][is];
  std::vector<float> dprev(ns);
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
static void sdw_obj::updateSumsOfErrors(int ie, int je, int ms, float **e, 
    std::vector<float> d, float *dmin, int *mmin) {
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

/*****************PUBLIC*****************/
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

void sdw_obj::init(int k, double smin, double smax, 
    axis *ax1, axis *ax2, axis *ax3) {
  double ds = (ax1->d)/k;
  int ismin = static_cast<int>(ceil(smin/ds));
  int ismin = static_cast<int>(floor(smin/ds));
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
  //_si = new sinc_interp();
  //_si.setExtrapolation(sinc_interp::Extrapolation::CONSTANT);
}

void sdw_obj::findShifts(axis axf, float **f, axis axg, float **g, float **s) {
  int n1 = _ax1->n;
  int n2 = _ax2->n;
  std::vector<int> k1s = subsample(n1,_k1min);
  std::vector<int> k2s = subsample(n2,_k2min);
  int nk1 = k1s.size();
  int nk2 = k2s.size();

  trace("finsShifts: smoothing in 1st dimension ...");
  float ***ek = new float[n2][][];
  float **e1 = static_cast<float**>(mem_alloc2(_axs->n,n1,sizeof(float)));
  for (int i2=0; i2<n2; ++i2) {
    computeErrors(axf,f[i2],axg,g[i2],e1);
    ek[i2] = subsampleErrors(_r1min,_r1max,k1s,_axs,_ax1,e1);
  }
  /*
  normalizeErrors(ek);

  trace("finsShifts: smoothing in 2nd dimension ...");
  float ***ekk = smoothErrors(_r2min,_r2max,k2s,_axs,_ax2,ek);
  normalizeErrors(ekk);

  trace("findShifts: finding shifts ...");
  float **ukk = new float[nk2][];
  for (int ik2=0; ik2<nk2; ++ik2) {
    ukk[ik2] = findShiftsFromSubsampledErrors(
        _r1min,_r1max,k1s,_axs,_ax1,ekk[ik2]);
  }
  trace("findShifts: interpolating shifts ...");
  float **u = interpolateShiftsBl(_ax1,_ax2,k1s,k2s,ukk);
  trace("findShifts: ... done");
  */
}

/**
 * Returns alignment errors computed for specified sequences.
 * @param sf sampling of 1st dimension for the seqeunce f.
 * @param f array of values for sequence f.
 * @param sg sampling of 1st dimension for the seqeunce g.
 * @param g array of values for sequence g.
 * @return array of alignment errors.
 */
void computeErrors(axis axf, float *f, axis axg, float *g, float **e) {
  int ns = _axs->n;
  int ne = _ax1->n;
  int nf =  axf->n;
  int ng =  axg->n;
  float *fi = static_cast<float*>(mem_alloc(ne,sizeof(float)));
  float *gi = static_cast<float*>(mem_alloc(ne,sizeof(float)));
  //_si.interpolate(axf,f,_ax1,fi);
  float sum = 0.0f;
  for (int is=0; is<ns; ++is) {
    //_si.interpolate(
    //    ng,axg->d,axg->o,g,
    //    ne,_ax1->d,_ax1->o+_axs.get_val(is),gi);
    for (int ie=0; ie<ne; ++ie) {
      e[ie][is] = error(fi[ie],gi[ie]);
    }
  }
  return e;
}
