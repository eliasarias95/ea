#ifndef SDW_SLOPE_H
#define SDW_SLOPE_H
#include <ucsl.h>
#include <sdw_obj.h>
#include <iostream>

class sdw_slope {
  private:
    int _k;
    double _pmax,_h1,_h2,_h3,_r1_r2,_r3;
    axis *_ax1, *_ax2, *ax3;
    sdw_obj *_sdw;
    void trace(std::string s);
    void setValueAt(int i1, int i2, float val);
    void setValueAt(int i1, int i2, int i3, float val);
    float getValueAt(int i1, int i2);
    float getValueAt(int i1, int i2, int i3);
    void interpolateSlopes(float **p);

  public:
    void init(int k, double pmax, double h1, double h2, double h3, 
        double r1, double r2, double r3, axis *ax1, axis *ax2, axis *ax3);

    //Constructors 
    sdw_slope(double pmax, axis *ax1, axis *ax2);
    sdw_slope(int k, double pmax, double h1, double h2, 
        double r1, double r2, axis *ax1, axis *ax2);
    sdw_slope(int k, double pmax, double h1, double h2, double h3, 
        double r1, double r2, double r3, axis *ax1, axis *ax2, axis *ax3);

    void setK(int k);
    void setErrorSmoothing(int esmooth);

    void findSlopes(axis *axf, float **f, float **p);
    void findSlopes(axis *axf, float ***f, float ***p2, float ***p3);
};
#endif
