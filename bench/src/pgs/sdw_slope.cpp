/********************************************************************
 * Slope estimation using a modified smooth dynamic warping algorithm.
 * For 2D inputs, 1 file of output slope values is produced.
 * For 3D inputs, 2 file of output slope values are produced (1 for xline
 * slope estimates and 1 for subline slope estimates). Computational speed,
 * memory, and slope estimate smoothness are all highly dependent on the
 * parameters chosen.
 *
 * Author: Elias Arias
 * Version: 09.09.2015
 * Adapted from Java code written by Elias Arias.
 ********************************************************************/

#include <sdw_slope.h>

/***********************************PUBLIC***********************************/

/**
 * Initializes private members of this class.
 * 
 * @param k slope sampling interval becomes 1/k.
 * @param pmax maximum slope to be computed.
 * @param h1 subsampling interval for which slopes are computed in 1st dim.
 * @param h2 subsampling interval for which slopes are computed in 2nd dim.
 * @param h3 subsampling interval for which slopes are computed in 3rd dim.
 * @param r1 strain limit in 1st dimension.
 * @param r2 strain limit in 2nd dimension.
 * @param r3 strain limit in 3rd dimension.
 * @param ax1 unit axis object of slopes for 1st dimension.
 * @param ax2 unit axis object of slopes for 2nd dimension.
 * @param ax3 unit axis object of slopes for 3rd dimension.
 */
void sdw_slope::init(int k, double pmax, double h1, double h2, double h3, 
    double r1, double r2, double r3, axis *ax1, axis *ax2, axis *ax3) {
  _k = k;
  _pmax = pmax;
  _h1 = h1;
  _h2 = h2;
  _h3 = h3;
  _r1 = r1;
  _r2 = r2;
  _r3 = r3;
  _ax1 = ax1;
  _ax2 = ax2;
  _ax3 = ax3;
  _sdw = new sdw_obj(k,-pmax,pmax,ax1,ax2,ax3);
  _sdw->setSmoothness(h1,h2,h3);
  _sdw->setStrainLimits(-r1,r1,-r2,r2,-r3,r3);
  if (_ax3 == NULL) _sdw->getMemoryCost2();
  else _sdw->getMemoryCost3();
}

/*******************constructors*******************/

/**
 * Constructs an sdw_slope object.
 * 
 * @param pmax maximum slope to be computed.
 * @param ax1 unit axis object of slopes for 1st dimension.
 * @param ax2 unit axis object of slopes for 2nd dimension.
 */
sdw_slope::sdw_slope(double pmax, axis *ax1, axis *ax2) {
  init(1,pmax,1.0,1.0,1.0,1.0,1.0,1.0,ax1,ax2,NULL);
}

/**
 * Constructs an sdw_slope object.
 * 
 * @param k slope sampling interval becomes 1/k.
 * @param pmax maximum slope to be computed.
 * @param h1 subsampling interval for which slopes are computed in 1st dim.
 * @param h2 subsampling interval for which slopes are computed in 2nd dim.
 * @param r1 strain limit in 1st dimension.
 * @param r2 strain limit in 2nd dimension.
 * @param ax1 unit axis object of slopes for 1st dimension.
 * @param ax2 unit axis object of slopes for 2nd dimension.
 */
sdw_slope::sdw_slope(int k, double pmax, double h1, double h2, 
    double r1, double r2, axis *ax1, axis *ax2) {
  init(k,pmax,h1,h2,1.0,r1,r2,1.0,ax1,ax2,NULL);
}

/**
 * Constructs an sdw_slope object.
 * 
 * @param k slope sampling interval becomes 1/k.
 * @param pmax maximum slope to be computed.
 * @param h1 subsampling interval for which slopes are computed in 1st dim.
 * @param h2 subsampling interval for which slopes are computed in 2nd dim.
 * @param h3 subsampling interval for which slopes are computed in 3rd dim.
 * @param r1 strain limit in 1st dimension.
 * @param r2 strain limit in 2nd dimension.
 * @param r3 strain limit in 3rd dimension.
 * @param ax1 unit axis object of slopes for 1st dimension.
 * @param ax2 unit axis object of slopes for 2nd dimension.
 * @param ax3 unit axis object of slopes for 3rd dimension.
 */
sdw_slope::sdw_slope(int k, double pmax, double h1, double h2, double h3, 
    double r1, double r2, double r3, axis *ax1, axis *ax2, axis *ax3) {
  init(k,pmax,h1,h2,h3,r1,r2,r3,ax1,ax2,ax3);
}

/*******************other functions*******************/

/**
 * Set the slope sampling interval.
 * E.g. for k=10, slope accuracy will be up to 1/k or 0.1 samples/trace.
 *
 * @param k slope sampling interval becomes 1/k.
 */
void sdw_slope::setK(int k) {
  _k = k;
}

/**
 * Sets the number of nonlinear smoothings of alignment errors.
 * In dynamic warping, alignment errors are smoothed the specified 
 * number of times, along all dimensions (in order 1, 2, ...), 
 * before estimating slopes by accumulating and backtracking along 
 * only the 1st dimension. 
 *
 * The default number of smoothings is zero, which is best for 1D
 * sequences. For 2D and 3D images, two smoothings are recommended.
 *
 * @param esmooth number of nonlinear smoothings.
 */
void sdw_slope::setErrorSmoothing(int esmooth) {
  _sdw->setErrorSmoothing(esmooth);
}

/**
 * Returns slopes computed for specified 2D images.
 *
 * @param axf unit axis object of 1st dimension for the image f.
 * @param f array of values for image f.
 * @param p array of slopes.
 */
void sdw_slope::findSlopes(axis *axf, float **f, float **p) {
  int n1 = _ax1->n;
  int n2 = _ax2->n;
  float **fm   = (float**)mem_alloc2(n1,n2,sizeof(float));

  memcpy(p[0],f[0],n1*sizeof(float));
  memcpy(p[n2-1],f[n2-2],n1*sizeof(float));
  for (int i2=1; i2<n2-1; ++i2) {
    memcpy(p[i2],f[i2-1],n1*sizeof(float));
  }

  _sdw->findShifts(axf,p,axf,f,fm);
  interpolateSlopes(fm,p);
  mem_free2((void***)&fm);
}

/**
 * Returns slopes computed for specified 3D images.
 *
 * @param axf unit axis object of 1st dimension for the image f.
 * @param f  array of values for image f.
 * @param p2 array of xline slopes.
 * @param p3 array of subline slopes.
 */
void sdw_slope::findSlopes(axis *axf, float ***f, float ***p2, float ***p3) {
  size_t n1 = _ax1->n;
  size_t n2 = _ax2->n;
  size_t n3 = _ax3->n;

  float ***fm = (float***)mem_alloc3(n1,n2,n3,sizeof(float));
  for (int i3=0; i3<n3; ++i3) {
    memcpy(p2[i3][0],f[i3][0],n1*sizeof(float));
    memcpy(p2[i3][n2-1],f[i3][n2-2],n1*sizeof(float));
    for (int i2=1; i2<n2-1; ++i2) {
      memcpy(p2[i3][i2],f[i3][i2-1],n1*sizeof(float));
    }
  }

  _sdw->findShifts(axf,p2,axf,f,fm);
  interpolateSlopes2(fm,p2);
  memset(fm[0][0],0,n1*n2*n3*sizeof(float));

  memcpy(p3[0][0],f[0][0],n1*n2*sizeof(float));
  memcpy(p3[n3-1][0],f[n3-2][0],n1*n2*sizeof(float));
  for (int i3=1; i3<n3-1; ++i3) {
    memcpy(p3[i3][0],f[i3-1][0],n1*n2*sizeof(float));
  }

  _sdw->findShifts(axf,p3,axf,f,fm);
  interpolateSlopes3(fm,p3);
  mem_free3((void****)&fm);
}

/***********************************PRIVATE***********************************/

/**
 * Using smooth dynamic warping, slopes are estimated on a grid that is
 * half way between image samples. This method interpolates slopes back
 * onto the desired grid.
 *
 * @param p array of slopes on incorrect grid.
 * @param pi array of slopes on correct grid.
 */
void sdw_slope::interpolateSlopes(float **p, float **pi) {
  int n1 = _ax1->n;
  int n2 = _ax2->n;
  // Interpolate.
#ifdef _USE_OMP
#pragma omp parallel for
#endif
  for (int i2=0; i2<n2; ++i2) {
    int i2p;
    if ((i2+1)>=n2) {i2p = n2-1;}
    else {i2p = i2+1;}
    float w2 = 0.5f;
    for (int i1=0; i1<n1; ++i1) {
      pi[i2][i1] = w2*p[i2][i1] + w2*p[i2p][i1];
    }
  }
}

/**
 * Using smooth dynamic warping, slopes are estimated on a grid that is
 * half way between image samples. This method interpolates slopes back
 * onto the desired grid.
 *
 * @param p array of slopes on incorrect grid.
 * @param pi array of slopes on correct grid.
 */
void sdw_slope::interpolateSlopes2(float ***p, float ***pi) {
  int n1 = _ax1->n;
  int n2 = _ax2->n;
  int n3 = _ax3->n;
  // Interpolate.
#ifdef _USE_OMP
#pragma omp parallel for
#endif
  for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
      int i2p;
      if ((i2+1)>=n2) {i2p = n2-1;}
      else {i2p = i2+1;}
      float w2 = 0.5f;
      for (int i1=0; i1<n1; ++i1) {
        pi[i3][i2][i1] = w2*p[i3][i2][i1] + w2*p[i3][i2p][i1];
      }
    }
  }
}

/**
 * Using smooth dynamic warping, slopes are estimated on a grid that is
 * half way between image samples. This method interpolates slopes back
 * onto the desired grid.
 *
 * @param p array of slopes on incorrect grid.
 * @param pi array of slopes on correct grid.
 */
void sdw_slope::interpolateSlopes3(float ***p, float ***pi) {
  int n1 = _ax1->n;
  int n2 = _ax2->n;
  int n3 = _ax3->n;
  // Interpolate.
#ifdef _USE_OMP
#pragma omp parallel for
#endif
  for (int i3=0; i3<n3; ++i3) {
    int i3p;
    if ((i3+1)>=n3) {i3p = n3-1;}
    else {i3p = i3+1;}
    float w3 = 0.5f;
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        pi[i3][i2][i1] = w3*p[i3][i2][i1] + w3*p[i3p][i2][i1];
      }
    }
  }
}
