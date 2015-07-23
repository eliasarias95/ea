#ifndef SINC_INTERP_H
#define SINC_INTERP_H
#include <cmath>
#include <unordered_map>
#include <math.h>
#include <axis.h>
#include <kaiser_window.h>

/**************************************************************************
  Copyright (c) 2012, Colorado School of Mines and others. All rights reserved.
  This program and accompanying materials are made available under the terms of
  the Common Public License - v1.0, which accompanies this distribution, and is
  available at http://www.eclipse.org/legal/cpl-v10.html
 ***************************************************************************/

/**
 *@author Dave Hale, Colorado School of Mines
 *@author Bill Harlan, Landmark Graphics
 *C++ translation by Elias Arias, Colorado School of Mines 
 *@version 2012.12.21
 *@version 2015.07.23
 */
class sinc_interp {

  public:
    enum class Extrapolation {
      ZERO,
      CONSTANT,
    };

  public:
    static sinc_interp *fromErrorAndLength(double emax, int lmax);
    static sinc_interp *fromErrorAndFrequency(double emax, double fmax);
    static sinc_interp *fromFrequencyAndLength(double fmax, int lmax);

    sinc_interp() : sinc_interp(0.0,0.3,8) {}

    virtual double getMaximumError();
    virtual double getMaximumFrequency();
    virtual int getMaximumLength();
    virtual long long getTableBytes();
    virtual Extrapolation getExtrapolation();
    virtual void setExtrapolation(Extrapolation extrap);

    virtual float interpolate(
        int nxu, double dxu, double fxu, float yu[], double xi);
    virtual void interpolate(
        int nxu, double dxu, double fxu, float yu[],
        int nxi, float xi[], float yi[]);
    virtual void interpolate(
        int nxu, double dxu, double fxu, float yu[], 
        int nxi, double dxi, double fxi, float yi[]);
    virtual float interpolate(axis *sxu, float yu[], double xi);
    virtual void interpolate(axis *sxu, float yu[], axis *sxi, float yi[]);

    virtual void accumulate(
        double xa, float ya, int nxu, double dxu, double fxu, float yu[]);
    virtual void accumulate(
        int nxa, float xa[], float ya[], 
        int nxu, double dxu, double fxu, float yu[]);

    virtual float **getTable();
    virtual int getNumberInTable();
    virtual int getLengthInTable();

  private:
    static const double EWIN_FRAC = 0.9;
    // Maximum table size, when maximum error and frequency are specified.
    static const int NTAB_MAX = 16385;
    // Extrapolation method.
    Extrapolation _extrap = Extrapolation::ZERO;
    // Table of sinc interpolation coefficients.
    Table *_table; // with all fields cached below
    int _lsinc, _nsinc, _ishift; // length of sinc approximations
    double _dsinc, _nsincm1; // sampling interval in table
    //need to delete somewhere below (where appropriate).
    float **_asinc; // array[nsinc][lsinc] of sinc approximations
    sinc_interp(double emax, double fmax, int lmax);

  private:
    class Design {
      public:
        double emax;
        double fmax;
        int lmax;
        Design(double emax, double fmax, int lmax) {
          this->emax = emax;
          this->fmax = fmax;
          this->lmax = lmax;
        }
    };

  private:
    class Table {
      public:
        Design *design; // here, all three design parameters are non-zero
        int lsinc ,nsinc, nsincm1, ishift;
        double dsinc;
        //need to delete somewhere below (where appropriate).
        float **asinc;
    };

  private:
    static Table *makeTable(Design *design);
    static Table *makeTable(int nsinc, int lsinc, kaiser_window *kwin);
    static double sinc(double x);
    //May be a future source of error
    static const std::unordered_map<Design*, Table*> _tables = 
      std::unordered_map<Design*, Table*>();
    static Table *getTable(double emax, double fmax, int lmax);

    float interpolate(
        double xscale, double xshift, int nxum, int nxu, float yu[], double x);

    void shift(
        int nxu, double dxu, double fxu, float yu[],
        int nxi,             double fxi, float yi[]);

    void accumulate(
        double xscale, double xshift, int nxum, 
        double x, float y, int nxu, float yu[]);
};
#endif
