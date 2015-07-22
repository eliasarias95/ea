#ifndef SDW_H
#define SDW_H
#include <ucsl.h>
#include <iostream>
#include <vector>

class sdw_obj {
  private:
    axis *_ax1, *_ax2, *_ax3, *_axs;
    double _r1min, _r2min, _r3min;
    double _r1max, _r2max, _r3max;
    float _epow = 1.00f;
    int _k1min, _k2min, _k3min;
    int _esmooth = 1;
    void trace(std::string s);
    static int *subsample(int n, int kmin);

  public:
    void init(int k, double smin, double smax, 
        axis *ax1, axis *ax2, axis *ax3);

    //Constructors
    sdw_obj(int k, double smin, double smax, axis *ax1);
    sdw_obj(int k, double smin, double smax, axis *ax1, axis *ax2);
    sdw_obj(int k, double smin, double smax, axis *ax1, axis *ax2, axis *ax3);

    void setStrainLimits(double r1min, double r1max);
    void setStrainLimits(double r1min, double r1max, 
                         double r2min, double r2max);
    void setStrainLimits(double r1min, double r1max, 
                         double r2min, double r2max, 
                         double r3min, double r3max);
    void setErrorSmoothing(int esmooth);
    void setSmoothness(double d1min);
    void setSmoothness(double d1min, double d2min);
    void setSmoothness(double d1min, double d2min, double d3min);
    void findShifts(axis axf, float **f, asix axg, float **g, float **s);
    void computeErrors(axis axf, float *f, axis axg, float *g, float *e);
};

#endif
