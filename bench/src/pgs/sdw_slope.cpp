#include <sdw_slope.h>

/*****************PRIVATE*****************/
void sdw_slope::trace(std::string s) {
  std::cout << s << "\n";
}

void sdw_slope::interpolateSlopes(float **p) {
  int n1 = _ax1->n;
  int n2 = _ax2->n;
  //float **pi   = (float**)mem_alloc2(n1,n2,sizeof(float));
  float *x1   = static_cast<float*>(mem_alloc(n1,sizeof(float)));
  float *x2   = static_cast<float*>(mem_alloc(n2,sizeof(float)));
  for (int i1=0; i1<n1; ++i1)
    x1[i1] = i1;
  for (int i2=0; i2<n2; ++i2)
    x2[i2] = i2-0.5f;

  //BilinearInterpolator2 bc = new BilinearInterpolator2(x1,x2,p);
  //p = bc.interpolate00(_ax1,_ax2);
  mem_free((void**)&x1);
  mem_free((void**)&x2);
}

/*****************PUBLIC*****************/
//Constructors 
sdw_slope::sdw_slope(double pmax, axis *ax1, axis *ax2) {
  init(1,pmax,1.0,1.0,1.0,1.0,1.0,1.0,ax1,ax2,NULL);
}

sdw_slope::sdw_slope(int k, double pmax, double h1, double h2, 
    double r1, double r2, axis *ax1, axis *ax2) {
  init(k,pmax,h1,h2,1.0,r1,r2,1.0,ax1,ax2,NULL);
}

sdw_slope::sdw_slope(int k, double pmax, double h1, double h2, double h3, 
    double r1, double r2, double r3, axis *ax1, axis *ax2, axis *ax3) {
  init(k,pmax,h1,h2,h3,r1,r2,r3,ax1,ax2,ax3);
}

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
  _sdw.setSmoothness(h1,h2,h3);
  _sdw.setStrainLimits(-r1,r1,-r2,r2,-r3,r3);
}

//Functions
void sdw_slope::setK(int k) {
  _k = k;
}

void sdw_slope::setErrorSmoothing(int esmooth) {
  _sdw.setErrorSmoothing(esmooth);
}

void sdw_slope::findSlopes(axis *axf, float **f, float **p) {
  int n1 = _ax1->n;
  int n2 = _ax2->n;
  float **fm   = static_cast<float**>(mem_alloc2(n1,n2,sizeof(float)));
  float **temp = static_cast<float**>(mem_alloc2(n1,n2,sizeof(float)));

  memcpy(fm[0],f[0],n1*sizeof(float));
  memcpy(fm[n2-1],f[n2-2],n1*sizeof(float));
  for (int i2=1; i2<n2-1; ++i2) {
    memcpy(fm[i2],f[i2-1],n1*sizeof(float));
  }

  //Fix to be done by address/reference so no need for temp
  p = _sdw.findShifts(axf,fm,axf,f);
  interpolateSlopes(p);

  mem_free2((void***)&fm);
  mem_free2((void***)&temp);
}

void sdw_slope::findSlopes(axis *axf, float ***f, float ***p2, float ***p3) {
  int n1 = _ax1->n;
  int n2 = _ax2->n;
  int n3 = _ax3->n;
  float ***f2m   = static_cast<float***>(mem_alloc3(n1,n2,n3,sizeof(float)));
  float ***f3m   = static_cast<float***>(mem_alloc3(n1,n2,n3,sizeof(float)));
  float ***temp2 = static_cast<float***>(mem_alloc3(n1,n2,n3,sizeof(float)));
  float ***temp3 = static_cast<float***>(mem_alloc3(n1,n2,n3,sizeof(float)));

  for (int i3=0; i3<n3; ++i3) {
    memcpy(f2m[i3][0],f[i3][0],n1*sizeof(float));
    memcpy(f2m[i3][n2-1],f[i3][n2-2],n1*sizeof(float));
    for (int i2=1; i2<n2-1; ++i2) {
      memcpy(f2m[i3][i2],f[i3][i2-1],n1*sizeof(float));
    }
  }

  memcpy(f3m[0],f[0],n1*n2*sizeof(float));
  memcpy(f3m[n3-1],f[n3-2],n1*n2*sizeof(float));
  for (int i3=1; i3<n3-1; ++i3) {
    memcpy(f3m[i3],f[i3-1],n1*n2*sizeof(float));
  }

  temp2 = _dwk.findShifts(axf,f2m,axf,f);
  temp2 = interpolateSlopes(temp2,2);
  temp3 = _dwk.findShifts(axf,f3m,axf,f);
  temp3 = interpolateSlopes(temp3,3);

  mem_free3((void****)&f2m);
  mem_free3((void****)&f3m);
  mem_free3((void****)&temp2);
  mem_free3((void****)&temp3);
}
