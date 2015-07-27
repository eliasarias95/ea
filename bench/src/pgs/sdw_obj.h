#ifndef SDW_H
#define SDW_H
#include <ucsl.h>
#include <math.h>
#include <assert.h>
//#include <sinc_interp.h>
#include <cmath>
#include <vector>
#include <iostream>

class sdw_obj {
  public:
    void init(int k, double smin, double smax, 
        axis *ax1, axis *ax2, axis *ax3);

    //Constructors
    sdw_obj(int k, double smin, double smax, axis *ax1);
    sdw_obj(int k, double smin, double smax, axis *ax1, axis *ax2);
    sdw_obj(int k, double smin, double smax, axis *ax1, axis *ax2, axis *ax3);

    //Setter methods
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

    //Other public methods
    void findShifts(float **e, float *s);    
    void findShifts(axis *axf, float *f, axis *axg, float *g, float *s);
    void findShifts(axis *axf, float **f, axis *axg, float **g, float **s);

    void computeErrors(axis *axf, float *f, axis *axg, float *g, float **e);
    static void normalizeErrors(int n1, int n2, float **e);

  private:
    axis *_ax1, *_ax2, *_ax3, *_axs;
    //sinc_interp *_si;
    double _r1min, _r2min, _r3min;
    double _r1max, _r2max, _r3max;
    float _epow;
    int _k1min, _k2min, _k3min;
    int _esmooth;

    static void fill(double val, float *x, int nx);
    vector<int> subsample(int n, int kmin);
    float error(float f, float g);

    static void shiftAndScale(
        float emin, float emax, int n1, int n2, float **e);
    static void shiftAndScale(
        float emin, float emax, int n1, int n2, int n3, float ***e);
    static void shiftAndScale(
        float emin, float emax, int n1, int n2, int n3, int n4, float ****e);

    static void smoothErrors2(double r2min, double r2max, vector<int> k2s,
        axis *axs, axis *axe, int nk1, int n2, float ***e, float ***es);
    static void smoothErrors2(double r2min, double r2max, vector<int> k2s,
        axis *axs, axis *axe, float ****e, float ****es);
    static void smoothErrors3(double r3min, double r3max, vector<int> k3s,
        axis *axs, axis *axe, float ****e, float ****es);

    static void smoothSubsampledErrors(
        double rmin, double rmax, vector<int> kes,
        axis *axs, axis *axe, float **e);
    static void smoothSubsampledErrors(
        double r1min, double r1max, vector<int> k1s,
        double r2min, double r2max, vector<int> k2s,
        axis *axs, axis *ax1, axis *ax2, float ***e);
    static void smoothSubsampledErrors(
        double r1min, double r1max, vector<int> k1s,
        double r2min, double r2max, vector<int> k2s,
        double r3min, double r3max, vector<int> k3s,
        axis *axs, axis *ax1, axis *ax2, axis *ax3, float ****e);

    static void subsampleErrors(double rmin, double rmax, 
        vector<int> kes, axis *axs, axis *axe, float **e, float **d);
    static void backtrackForShifts(vector<int> kes, axis *axs, axis *axe,
        float *d, int **m, float* ske);    

    static void normalizeErrors(int n1, int n2, int n3, float ***e);
    static void normalizeErrors(int n1, int n2, int n3, int n4, float ****e);

    static void accumulate(
        int dir, double rmin, double rmax, int me, 
        axis *axs, axis *axe, float **e, float **d);
    static void accumulate(
        int dir, double rmin, double rmax, vector<int> kes, 
        axis *axs, axis *axe, float **e, float **d, int **m);
    static void accumulateSubsampled(
        int dir, double rmin, double rmax, vector<int> kes,
        axis *axs, axis *axe, float **e, float **d, int **m);

    static void findShiftsFromErrors(double rmin, double rmax,
        vector<int> kes, axis *axs, axis *axe, float **e, float *s);
    static void findShiftsFromSubsampledErrors(double rmin, double rmax,
        vector<int> kes, axis *axs, axis *axe, float **e, float *s);

    static void interpolateShifts(
        axis *ax1, vector<int> k1s, float *sk, float *s);
    static void interpolateShifts(
        axis *ax1, axis *ax2, vector<int> k1s, vector<int> k2s, 
        float **skk, float **s);
    static void interpolateShifts(
        axis *ax1, axis *ax2, axis *ax3, 
        vector<int> k1s, vector<int> k2s, vector<int> k3s, 
        float ***skk, float ***s);

    static void updateSumsOfErrors(int ie, int je, int ms, float **e, 
        vector<float> d, float *dmin, int *mmin);
};
#endif
