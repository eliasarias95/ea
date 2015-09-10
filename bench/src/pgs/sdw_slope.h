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

#ifndef SDW_SLOPE_H
#define SDW_SLOPE_H
#include <ucsl.h>
#include <sdw_obj.h>

class sdw_slope {
  public:
    void init(int k, double pmax, double h1, double h2, double h3, 
        double r1, double r2, double r3, axis *ax1, axis *ax2, axis *ax3);

/*******************constructors*******************/
    sdw_slope(double pmax, axis *ax1, axis *ax2);
    sdw_slope(int k, double pmax, double h1, double h2, 
        double r1, double r2, axis *ax1, axis *ax2);
    sdw_slope(int k, double pmax, double h1, double h2, double h3, 
        double r1, double r2, double r3, axis *ax1, axis *ax2, axis *ax3);

/*******************other functions*******************/
    void setK(int k);
    void setErrorSmoothing(int esmooth);

    void findSlopes(axis *axf, float **f, float **p);
    void findSlopes(axis *axf, float ***f, float ***p2, float ***p3);

/***********************************PRIVATE***********************************/
  private:
    int _k;
    double _pmax,_h1,_h2,_h3,_r1,_r2,_r3;
    axis *_ax1, *_ax2, *_ax3;
    sdw_obj *_sdw;
    void interpolateSlopes(float **p, float **pi);
    void interpolateSlopes2(float ***p, float ***pi);
    void interpolateSlopes3(float ***p, float ***pi);
};
#endif
