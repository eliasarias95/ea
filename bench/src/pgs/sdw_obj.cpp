#include <sdw_obj.h>

/*****************PRIVATE*****************/
sdw_obj::trace(std::string s) {
  std::cout << s << "\n";
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
static std::vector<int> subsample(int n, int kmin) {
  if (kmin>=n)
    kmin = n-1;
  int m = 1+(n-1)/kmin;
  double d = static_cast<double>(n-1)/static_cast<double>(m-1);
  std::vector<int> j(m);
  for (int i=0; i<m; ++i)
    j[i] = static_cast<int>(i*d+0.5);
  return j;
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
  for (int i2=0; i2<n2; ++i2) {
    float **e1 = computeErrors(axf,f[i2],axg,g[i2]);
    ek[i2] = subsampleErrors(_r1min,_r1max,k1s,_axs,_ax1,e1);
  }
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
  Sampling ss = _ss;
  Sampling se = _s1;
  int ns = _axs->n;
  int ne = _ax1->n;
  int nf = axf->;
  int ng = axg->;
  float[] fi = new float[ne];
  float[] gi = new float[ne];
  _si.interpolate(sf,f,se,fi);
  float sum = 0;
  for (int is=0; is<ns; ++is) {
    _si.interpolate(
        ng,sg.getDelta(),sg.getFirst(),g,
        ne,se.getDelta(),se.getFirst()+ss.getValue(is),gi);
    for (int ie=0; ie<ne; ++ie) {
      e[ie][is] = error(fi[ie],gi[ie]);
      //if (e[ie][is]>= 2.0f)
      //  ++sum;
    }
  }
  //trace("e length = "+e[0].length*e.length);
  //trace("sum= "+sum);
  //trace("max error= "+max(e));
  return e;
}
