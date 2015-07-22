/****************************************************************************
  Copyright (c) 2012, Colorado School of Mines and others. All rights reserved.
  This program and accompanying materials are made available under the terms of
  the Common Public License - v1.0, which accompanies this distribution, and is
  available at http://www.eclipse.org/legal/cpl-v10.html
 ****************************************************************************/
package slopes;

import java.util.HashMap;

import static edu.mines.jtk.util.ArrayMath.*;
import edu.mines.jtk.dsp.Sampling;

/**
 * @author Dave Hale, Colorado School of Mines
 * @author Bill Harlan, Landmark Graphics
 * @version 2012.12.21
 */
public class SincInterpolator {

  /**
   * The method used to extrapolate samples when interpolating.
   * Sampled functions are defined implicitly by extrapolation outside the 
   * domain for which samples are specified explicitly, with either zero or 
   * constant values. If constant, the extrapolated values are the first and 
   * last uniform sample values. The default is extrapolation with zeros.
   */
  public enum Extrapolation {
    ZERO, 
      CONSTANT,
  }

  /**
   * Returns a sinc interpolator with specified maximum error and length.
   * Computes the maximum frequency fmax. Note that for some parameters
   * emax and lmax, the maximum frequency fmax may be zero. In this case,
   * the returned interpolator is useless.
   * @param emax the maximum error for frequencies less than fmax; e.g., 
   *  0.01 for 1% percent error. 0.0 &lt; emax &lt;= 0.1 is required.
   * @param lmax the maximum interpolator length, in samples. 
   *  Must be an even integer not less than 8.
   * @return the sinc interpolator.
   */
  public static SincInterpolator fromErrorAndLength(
      double emax, int lmax)
  {
    return new SincInterpolator(emax,0.0,lmax);
  }

  /**
   * Returns a sinc interpolator with specified maximum error and frequency.
   * Computes the maximum length lmax.
   * @param emax the maximum error for frequencies less than fmax; e.g., 
   *  0.01 for 1% percent error. Must be greater than 0.0 and less than 1.0.
   * @param fmax the maximum frequency, in cycles per sample. 
   *  Must be greater than 0.0 and less than 0.5.
   * @return the sinc interpolator.
   */
  public static SincInterpolator fromErrorAndFrequency(
      double emax, double fmax)
  {
    return new SincInterpolator(emax,fmax,0);
  }

  /**
   * Returns a sinc interpolator with specified maximum frequency and length.
   * Computes the maximum error emax.
   * <p>
   * The product (1-2*fmax)*lmax must be greater than one. For when this 
   * product is less than one, a useful upper bound on interpolation error 
   * cannot be computed.
   * @param fmax the maximum frequency, in cycles per sample. 
   *  Must be greater than 0.0 and less than 0.5*(1.0-1.0/lmax).
   * @param lmax the maximum interpolator length, in samples. 
   *  Must be an even integer not less than 8 and greater than 
   *  1.0/(1.0-2.0*fmax).
   * @return the sinc interpolator.
   */
  public static SincInterpolator fromFrequencyAndLength(
      double fmax, int lmax)
  {
    return new SincInterpolator(0.0,fmax,lmax);
  }

  /**
   * Constructs a default sinc interpolator. The default design parameters 
   * are fmax = 0.3 cycles/sample (60% of Nyquist) and lmax = 8 samples.
   * For these parameters, the computed maximum error is less than 0.007
   * (0.7%). In testing, observed maximum error is less than 0.005 (0.5%).
   */
  public SincInterpolator() {
    this(0.0,0.3,8);
  }

  /**
   * Gets the maximum error for this interpolator.
   * @return the maximum error.
   */
  public double getMaximumError() {
    return _table.design.emax;
  }

  /**
   * Gets the maximum frequency for this interpolator.
   * @return the maximum frequency.
   */
  public double getMaximumFrequency() {
    return _table.design.fmax;
  }

  /**
   * Gets the maximum length for this interpolator.
   * @return the maximum length.
   */
  public int getMaximumLength() {
    return _table.design.lmax;
  }

  /**
   * Gets the number of bytes consumed by the table of interpolators.
   * The use of interpolators with small emax and large lmax may require 
   * the computation of large tables. This method can be used to determine 
   * how much memory is consumed by the table for an interpolator.
   * @return the number of bytes.
   */
  public long getTableBytes() {
    long nbytes = 4L;
    nbytes *= _table.lsinc;
    nbytes *= _table.nsinc;
    return nbytes;
  }

  /**
   * Gets the extrapolation method for this interpolator.
   * @return the extrapolation method.
   */
  public Extrapolation getExtrapolation() {
    return _extrap;
  }

  /**
   * Sets the extrapolation method for this interpolator.
   * The default extrapolation method is extrapolation with zeros.
   * @param extrap the extrapolation method.
   */
  public void setExtrapolation(Extrapolation extrap) {
    _extrap = extrap;
  }

  /**
   * Interpolates one real value y(x).
   * @param nxu number of input samples.
   * @param dxu input sampling interval.
   * @param fxu first input sampled x value.
   * @param yu input array of sampled values y(x).
   * @param xi value x at which to interpolate.
   * @return interpolated value y(x).
   */
  public float interpolate(
      int nxu, double dxu, double fxu, float[] yu, double xi)
  {
    double xscale = 1.0/dxu;
    double xshift = _lsinc-fxu*xscale;
    int nxum = nxu-_lsinc;
    return interpolate(xscale,xshift,nxum,nxu,yu,xi);
  }

  /**
   * Interpolates multiple real values y(x).
   * @param nxu number of input samples.
   * @param dxu input sampling interval.
   * @param fxu first input sampled x value.
   * @param yu input array of sampled values y(x).
   * @param nxi number of output samples.
   * @param xi input array of x values at which to interpolate.
   * @param yi output array of interpolated values y(x).
   */
  public void interpolate(
      int nxu, double dxu, double fxu, float[] yu, 
      int nxi, float[] xi, float[] yi)
  {
    double xscale = 1.0/dxu;
    double xshift = _lsinc-fxu*xscale;
    int nxum = nxu-_lsinc;
    for (int ixi=0; ixi<nxi; ++ixi)
      yi[ixi] = interpolate(xscale,xshift,nxum,nxu,yu,xi[ixi]);
  }

  /**
   * Interpolates multiple real values y(x).
   * @param nxu number of input samples.
   * @param dxu input sampling interval.
   * @param fxu first input sampled x value.
   * @param yu input array of sampled values y(x).
   * @param nxi number of output samples.
   * @param dxi output sampling interval.
   * @param fxi first output sampled x value.
   * @param yi output array of interpolated values y(x).
   */
  public void interpolate(
      int nxu, double dxu, double fxu, float[] yu, 
      int nxi, double dxi, double fxi, float[] yi)
  {
    if (dxu==dxi) {
      shift(nxu,dxu,fxu,yu,nxi,fxi,yi);
    } else {
      double xscale = 1.0/dxu;
      double xshift = _lsinc-fxu*xscale;
      int nxum = nxu-_lsinc;
      for (int ixi=0; ixi<nxi; ++ixi) {
        double xi = fxi+ixi*dxi;
        yi[ixi] = interpolate(xscale,xshift,nxum,nxu,yu,xi);
      }
    }
  }

  /**
   * Interpolates one real value y(x).
   * @param sxu sampling of input samples.
   * @param yu input array of uniformly sampled values y(x).
   * @param xi value x at which to interpolate.
   * @return interpolated value y(x).
   */
  public float interpolate(Sampling sxu, float[] yu, double xi) {
    return interpolate(sxu.getCount(),sxu.getDelta(),sxu.getFirst(),yu,xi);
  }

  /**
   * Interpolates multiple real values y(x).
   * @param sxu sampling of input samples.
   * @param yu input array of uniformly sampled values y(x).
   * @param sxi sampling of output samples.
   * @param yi output array of interpolated values y(x).
   */
  public void interpolate(
      Sampling sxu, float[] yu, 
      Sampling sxi, float[] yi) 
  {
    if (sxi.isUniform()) {
      interpolate(sxu.getCount(),sxu.getDelta(),sxu.getFirst(),yu,
          sxi.getCount(),sxi.getDelta(),sxi.getFirst(),yi);
    } else {
      int nxu = sxu.getCount();
      int nxi = sxi.getCount();
      double xscale = 1.0/sxu.getDelta();
      double xshift = _lsinc-sxu.getFirst()*xscale;
      int nxum = nxu-_lsinc;
      for (int ixi=0; ixi<nxi; ++ixi) {
        double xi = sxi.getValue(ixi);
        yi[ixi] = interpolate(xscale,xshift,nxum,nxu,yu,xi);
      }
    }
  }

  /**
   * Accumulates a specified real value y(x) into uniformly-sampled yu.
   * Accumulation is the transpose (not the inverse) of interpolation.
   * Whereas interpolation gathers from uniformly sampled values.
   * accumulation scatters into uniformly sampled values.
   * @param xa value x at which to accumulate.
   * @param ya value y(x) to accumulate.
   * @param nxu number of input/output samples.
   * @param dxu input/output sampling interval.
   * @param fxu first input/output sampled x value.
   * @param yu input/output array of sampled values y(x).
   */
  public void accumulate(
      double xa, float ya,
      int nxu, double dxu, double fxu, float[] yu) 
  {
    double xscale = 1.0/dxu;
    double xshift = _lsinc-fxu*xscale;
    int nxum = nxu-_lsinc;
    accumulate(xscale,xshift,nxum,xa,ya,nxu,yu);
  }

  /**
   * Accumulates a specified real value y(x) into uniformly-sampled yu.
   * Accumulation is the transpose (not the inverse) of interpolation.
   * Whereas interpolation gathers from uniformly sampled values.
   * accumulation scatters into uniformly sampled values.
   * @param nxa number of values to accumulate.
   * @param xa input array of values x at which to accumulate.
   * @param ya input array of values y(x) to accumulate.
   * @param nxu number of input/output samples.
   * @param dxu input/output sampling interval.
   * @param fxu first input/output sampled x value.
   * @param yu input/output array of sampled values y(x).
   */
  public void accumulate(
      int nxa, float[] xa, float[] ya,
      int nxu, double dxu, double fxu, float[] yu) 
  {
    double xscale = 1.0/dxu;
    double xshift = _lsinc-fxu*xscale;
    int nxum = nxu-_lsinc;
    for (int ixa=0; ixa<nxa; ++ixa)
      accumulate(xscale,xshift,nxum,xa[ixa],ya[ixa],nxu,yu);
  }

  /**
   * Get a copy of the interpolation table.  Returns a copy of this
   * interpolator's table of sinc interpolation coefficients.
   * @return A copy of the table.
   */
  public float[][] getTable() {
    return copy(_table.asinc);
  }

  /**
   * Gets the number of sampled sinc approximations in the table.
   * @return the number of tabulated sampled sinc approximations.
   */
  public int getNumberInTable() {
    return _table.asinc.length;
  }

  /**
   * Gets the length of sampled sinc approximations in the table.
   * @return the length of tabulated sampled sinc approximations.
   */
  public int getLengthInTable() {
    return _table.asinc[0].length;
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  // Fraction of error due to windowing; remainder is due to table lookup.
  private static final double EWIN_FRAC = 0.9;

  // Maximum table size, when maximum error and frequency are specified.
  private static final int NTAB_MAX = 16385;

  // Extrapolation method.
  private Extrapolation _extrap = Extrapolation.ZERO;

  // Table of sinc interpolation coefficients.
  private Table _table; // with all fields cached below
  private int _lsinc; // length of sinc approximations
  private int _nsinc; // number of sinc approximations
  private double _dsinc; // sampling interval in table
  private float[][] _asinc; // array[nsinc][lsinc] of sinc approximations
  private double _nsincm1; // nsinc-1
  private int _ishift; // -lsinc-lsinc/2+1

  /**
   * Constructs a sinc interpolator with specified parameters.
   * Exactly one of the parameters must be zero, and is computed here.
   */
  private SincInterpolator(double emax,double fmax,int lmax) {
    _table = getTable(emax,fmax,lmax);
    _lsinc = _table.lsinc;
    _nsinc = _table.nsinc;
    _nsincm1 = _table.nsincm1;
    _ishift = _table.ishift;
    _dsinc = _table.dsinc;
    _asinc = _table.asinc;
  }

  /**
   * Design parameters (can be turned into a struct cpp).
   */
  private static class Design {
    double emax;
    double fmax; 
    int lmax;
    Design(double emax, double fmax, int lmax) {
      this.emax = emax;
      this.fmax = fmax;
      this.lmax = lmax;
    }
  }

  /**
   * Table of sinc interpolator coefficients.
   */
  private static class Table {
    Design design; // here, all three design parameters are non-zero
    int lsinc,nsinc,nsincm1,ishift;
    double dsinc;
    float[][] asinc;
  }

  /**
   * Builds a table of interpolators for specified design parameters.
   * Exactly one of the design parameters must be zero.
   */
  private static Table makeTable(Design design) {
    double emax = design.emax;
    double fmax = design.fmax;
    int lmax = design.lmax;

    // The Kaiser window transition width is twice the difference 
    // between the Nyquist frequency 0.5 and the maximum frequency.
    double wwin = 2.0*(0.5-fmax);

    // The Kaiser window accounts for a hard-wired fraction of the maximum 
    // interpolation error. The other error will be due to table lookup.
    double ewin = emax*EWIN_FRAC;
    KaiserWindow kwin = KaiserWindow.fromWidthAndLength(wwin,lmax);
    ewin = 3.0*kwin.getError();
    emax = ewin/EWIN_FRAC;
    double etabMin = 1.1*PI*fmax/(NTAB_MAX-1);
    double emaxMin = etabMin/(1.0-EWIN_FRAC);
    if (emax<emaxMin) {
      emax = emaxMin;
      ewin = emax*EWIN_FRAC;
    }

    double etab = emax-ewin;
    double dsinc = (fmax>0.0)?etab/(PI*fmax):1.0;
    int nsincMin = 1+(int)ceil(1.0/dsinc);
    int nsinc = 2;
    while (nsinc<nsincMin)
      nsinc *= 2;
    ++nsinc;
    int lsinc = lmax;
    Table table = makeTable(nsinc,lsinc,kwin);
    table.design = new Design(emax,fmax,lmax);
    _tables.put(design,table); // key is design with one zero parameter
    return table;
  }

  /**
   * Builds a table of interpolators for a specified Kaiser window.
   */
  private static Table makeTable(int nsinc, int lsinc, KaiserWindow kwin) {
    float[][] asinc = new float[nsinc][lsinc];
    int nsincm1 = nsinc-1;
    int ishift = -lsinc-lsinc/2+1;
    double dsinc = 1.0/(nsinc-1);

    // The first and last interpolators are shifted unit impulses.
    // Handle these two cases exactly, with no rounding errors.
    for (int j=0; j<lsinc; ++j) {
      asinc[      0][j] = 0.0f;
      asinc[nsinc-1][j] = 0.0f;
    }
    asinc[      0][lsinc/2-1] = 1.0f;
    asinc[nsinc-1][lsinc/2  ] = 1.0f;

    // Other interpolators are sampled Kaiser-windowed sinc functions.
    for (int isinc=1; isinc<nsinc-1; ++isinc) {
      double x = -lsinc/2+1-dsinc*isinc;
      for (int i=0; i<lsinc; ++i,x+=1.0) {
        asinc[isinc][i] = (float)(sinc(x)*kwin.evaluate(x));
      }
    }
    Table table = new Table();
    table.lsinc = lsinc;
    table.nsinc = nsinc;
    table.nsincm1 = nsincm1;
    table.ishift = ishift;
    table.dsinc = dsinc;
    table.asinc = asinc;
    return table;
  }
  private static double sinc(double x) {
    return (x!=0.0)?sin(PI*x)/(PI*x):1.0;
  }

  /**
   * Map from design parameters to tables of coefficients.
   * This map saves both time and space required to compute the tables.
   */
  private final static HashMap<Design,Table> _tables = new HashMap<Design,Table>();
  private static Table getTable(double emax, double fmax, int lmax) {
    Design design = new Design(emax,fmax,lmax);
    synchronized(_tables) {
      Table table = _tables.get(design);
      if (table==null)
        table = makeTable(design);
      return table;
    }
  }

  private float interpolate(
      double xscale, double xshift, int nxum, int nxu, 
      float[] yu, double x)
  {
    // Which uniform samples?
    double xn = xshift+x*xscale;
    int ixn = (int)xn;
    int kyu = _ishift+ixn;

    // Which sinc approximation?
    double frac = xn-ixn;
    if (frac<0.0)
      frac += 1.0;
    int ksinc = (int)(frac*_nsincm1+0.5);
    float[] asinc = _asinc[ksinc];

    // If no extrapolation is necessary, use a fast loop.
    // Otherwise, extrapolate uniform samples, as necessary.
    float yr = 0.0f;
    if (kyu>=0 && kyu<=nxum) {
      for (int isinc=0; isinc<_lsinc; ++isinc,++kyu)
        yr += yu[kyu]*asinc[isinc];
    } else if (_extrap==Extrapolation.ZERO) {
      for (int isinc=0; isinc<_lsinc; ++isinc,++kyu) {
        if (0<=kyu && kyu<nxu)
          yr += yu[kyu]*asinc[isinc];
      }
    } else if (_extrap==Extrapolation.CONSTANT) {
      for (int isinc=0; isinc<_lsinc; ++isinc,++kyu) {
        int jyu = (kyu<0)?0:(nxu<=kyu)?nxu-1:kyu;
        yr += yu[jyu]*asinc[isinc];
      }
    }
    return yr;
  }

  private void shift(
      int nxu, double dxu, double fxu, float[] yu,
      int nxi,             double fxi, float[] yi)
  {
    double lxu = fxu+(nxu-1)*dxu;
    double xscale = 1.0/dxu;
    double xshift = _lsinc-fxu*xscale;
    int nxum = nxu-_lsinc;

    // Which output samples are near beginning and end of uniform sequence?
    double dx = dxu;
    double x1 = fxu+dxu*_lsinc/2;
    double x2 = lxu-dxu*_lsinc/2;
    double x1n = (x1-fxi)/dx;
    double x2n = (x2-fxi)/dx;
    int ix1 = max(0,min(nxi,(int)x1n+1));
    int ix2 = max(0,min(nxi,(int)x2n-1));

    // Interpolate output samples near beginning of uniform sequence.
    for (int ixi=0; ixi<ix1; ++ixi) {
      double xi = fxi+ixi*dx;
      yi[ixi] = interpolate(xscale,xshift,nxum,nxu,yu,xi);
    }

    // Interpolate output samples near end of uniform sequence.
    for (int ixi=ix2; ixi<nxi; ++ixi) {
      double xi = fxi+ixi*dx;
      yi[ixi] = interpolate(xscale,xshift,nxum,nxu,yu,xi);
    }

    // Now we ignore the ends, and use a single sinc approximation.

    // Which uniform samples?
    double xn = xshift+(fxi+ix1*dx)*xscale;
    int ixn = (int)xn;
    int kyu = _ishift+ixn;

    // Which sinc approximation?
    double frac = xn-ixn;
    if (frac<0.0)
      frac += 1.0;
    int ksinc = (int)(frac*_nsincm1+0.5);
    float[] asinc = _asinc[ksinc];

    // Interpolate for output indices ix1 <= ix <= ix2.
    for (int ix=ix1; ix<ix2; ++ix,++kyu) {
      float yr = 0.0f;
      for (int isinc=0,jyu=kyu; isinc<_lsinc; ++isinc,++jyu)
        yr += yu[jyu]*asinc[isinc];
      yi[ix] = yr;
    }
  }

  private void accumulate(
      double xscale, double xshift, int nxum,
      double x, float y, int nxu, float[] yu) 
  {
    // Which uniform samples?
    double xn = xshift+x*xscale;
    int ixn = (int)xn;
    int kyu = _ishift+ixn;

    // Which sinc approximation?
    double frac = xn-ixn;
    if (frac<0.0)
      frac += 1.0;
    int ksinc = (int)(frac*_nsincm1+0.5);
    float[] asinc = _asinc[ksinc];

    // If no extrapolation is necessary, use a fast loop.
    // Otherwise, extrapolate uniform samples, as necessary.
    if (kyu>=0 && kyu<=nxum) {
      for (int isinc=0; isinc<_lsinc; ++isinc,++kyu)
        yu[kyu] += y*asinc[isinc];
    } else if (_extrap==Extrapolation.ZERO) {
      for (int isinc=0; isinc<_lsinc; ++isinc,++kyu) {
        if (0<=kyu && kyu<nxu)
          yu[kyu] += y*asinc[isinc];
      }
    } else if (_extrap==Extrapolation.CONSTANT) {
      for (int isinc=0; isinc<_lsinc; ++isinc,++kyu) {
        int jyu = (kyu<0)?0:(nxu<=kyu)?nxu-1:kyu;
        yu[jyu] += y*asinc[isinc];
      }
    }
  }
}
