/********************************************************************
 * Dynamic warping to find shifts between two sequences or images.
 * For 1D sequences f(x) and g(x), dynamic warping computes shifts s(x) such
 * that f(x) ~ g(x+s(x)). For 2D images f(x1,x2) and g(x1,x2), it finds shifts
 * u(x1,x2) such that f(x1,x2) ~ g(x1+s(x1,x2),x2). Note that the shifts
 * u(x1,x2,...) are computed for only the first dimension of multi-dimensional
 * images. For example, if the 1st dimension of an image is time, then only
 * time shifts are computed.
 *
 * Constraints are placed on strains, the rates at which shifts change in any
 * dimension. For example, strain r1 = du/d1 is the derivative of shift
 * s(x1,x2,...) with respect to the x1 coordinate, and is constrained to lie
 * between lower and upper bounds r1min and r1max. Default bounds are 
 * r1min = -1.0 and r1max = 1.0.
 *
 * In many applications of warping, strains derived from estimated shifts may
 * be an important by-product. However, when computing shifts, only a finite
 * number of strains are permitted, and this quantization may yield
 * unrealistic estimates of strain. To address this problem, this dynamic
 * warping provides control over the sampling of strains, in the form of a
 * smoothness value. The number of strains sampled is proportional to this
 * value.
 *
 * Smoothness values represent approximate intervals for a quasi-uniform
 * subsampling grid on which shifts are computed. For example, for a time
 * sampling interval of 4 ms, one might specify a smoothness of 200 ms, so
 * that shifts are computed on a grid that is 50 times more coarse than the 4
 * ms sampling grid. However, strains on this coarse grid would be sampled 50
 * times more finely. After initially computing shifts on the coarse grid,
 * shifts are interpolated onto the finer grid.
 *
 * Author: Elias Arias and Sergey Frolov
 * Version: 09.09.2015
 * Adapted from Java code written by Dave Hale and Elias Arias.
 ********************************************************************/

#ifndef SDW_H
#define SDW_H
#include <ucsl.h>
#include <math.h>
#include <assert.h>
#include <cmath>
#include <vector>
#include <iostream>
#include <omp.h>
#include <time.h>

class sdw_obj {
  public:
    void init(int k, double smin, double smax, 
        axis *ax1, axis *ax2, axis *ax3);

/*******************constructors*******************/
    sdw_obj(int k, double smin, double smax, axis *ax1);
    sdw_obj(int k, double smin, double smax, axis *ax1, axis *ax2);
    sdw_obj(int k, double smin, double smax, axis *ax1, axis *ax2, axis *ax3);

/*******************setter methods*******************/
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

/*****************other public methods*****************/
    void findShifts(float **e, float *s);    
    void findShifts(axis *axf, float *f, axis *axg, float *g, float *s);
    void findShifts(axis *axf, float **f, axis *axg, float **g, float **s);
    void findShifts(axis *axf, float ***f, axis *axg, float ***g, float ***s);

    void computeErrors(axis *axf, float *f, axis *axg, float *g, float **e);
    static void normalizeErrors(int n1, int n2, float **e);
    void getMemoryCost2();
    void getMemoryCost3();

/***********************************PRIVATE***********************************/
  private:
    axis *_ax1, *_ax2, *_ax3, *_axs;
    double _r1min, _r2min, _r3min;
    double _r1max, _r2max, _r3max;
    float _epow;
    int _k1min, _k2min, _k3min;
    int _esmooth;

/*******************utility methods*******************/
    static void interp(axis *axf, float *f, axis *axs, float *fi, float shift);
    static void fill(float val, float *x, int nx);
    vector<int> subsample(int n, int kmin);
    float error(float f, float g);

/*******************shift and scale*******************/
    static void shiftAndScale(
        float emin, float emax, int n1, int n2, float **e);
    static void shiftAndScale(
        float emin, float emax, int n1, int n2, int n3, float ***e);
    static void shiftAndScale(
        float emin, float emax, int n1, int n2, int n3, int n4, float ****e);

/*******************smooth errors*******************/
    static float ***smoothErrors2(
        double r2min, double r2max, vector<int> k2s,axis *axs, axis *axe, 
        int nk1, int n2, float ***e);
    static float ****smoothErrors2(
        double r2min, double r2max, vector<int> k2s, axis *axs, axis *axe, 
        int nk1, int n2, int n3, float ****e);
    static float ****smoothErrors3(
        double r3min, double r3max, vector<int> k3s, axis *axs, axis *axe, 
        int nk1, int nk2, int n3, float ****e);

/*******************smooth subsampled errors*******************/
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

/*******************normalization*******************/
    static void normalizeErrors(int n1, int n2, int n3, float ***e);
    static void normalizeErrors(int n1, int n2, int n3, int n4, float ****e);

/*******************accumulation*******************/
    static void accumulate(
        int dir, double rmin, double rmax, int me, 
        axis *axs, axis *axe, float **e, float **d);
    static void accumulate(
        int dir, double rmin, double rmax, vector<int> kes, 
        axis *axs, axis *axe, float **e, float **d, int **m);
    static void accumulateSubsampled(
        int dir, double rmin, double rmax, vector<int> kes,
        axis *axs, axis *axe, float **e, float **d, int **m);

/*******************find shifts from errors*******************/
    static void findShiftsFromErrors(double rmin, double rmax,
        vector<int> kes, axis *axs, axis *axe, float **e, float *s);
    static void findShiftsFromSubsampledErrors(
        double rmin, double rmax, vector<int> kes, axis *axs, axis *axe,
        int **m, float **d, float **e, float *s);

/*******************interpolation*******************/
    static void interpolateShifts(
        axis *ax1, vector<int> k1s, float *sk, float *s);
    static void interpolateShifts(
        axis *ax1, axis *ax2, vector<int> k1s, vector<int> k2s, 
        float **skk, float **s);
    static void interpolateShifts(
        axis *ax1, axis *ax2, axis *ax3, 
        vector<int> k1s, vector<int> k2s, vector<int> k3s, 
        float ***skk, float ***s);

/*******************other private methods*******************/
    static void subsampleErrors(double rmin, double rmax, vector<int> kes,
        axis *axs, axis *axe, float **e, float **d);
    static void backtrackForShifts(vector<int> kes, axis *axs, axis *axe,
        float *d, int **m, float* ske);    
    static void updateSumsOfErrors(int ie, int je, int ms, float **e, 
        float *d, float *dmin, int *mmin, int ns);
};
#endif
