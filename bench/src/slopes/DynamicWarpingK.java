/****************************************************************************
Copyright (c) 2012, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package slopes;

import edu.mines.jtk.dsp.Sampling;
import edu.mines.jtk.dsp.Tensors2;
import edu.mines.jtk.dsp.EigenTensors2;
import edu.mines.jtk.interp.*;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.ArrayMath.*;

import dnp.*;
import interp.*;
import warp.*;

// FOR DEVELOPMENT ONLY
import edu.mines.jtk.mosaic.*;

/**
 * Dynamic warping to find shifts between two sequences or images.
 * For 1D sequences f(x) and g(x), dynamic warping computes shifts u(x) such
 * that f(x) ~ g(x+u(x)). For 2D images f(x1,x2) and g(x1,x2), it finds shifts
 * u(x1,x2) such that f(x1,x2) ~ g(x1+u(x1,x2),x2). Note that the shifts
 * u(x1,x2,...) are computed for only the first dimension of multi-dimensional
 * images. For example, if the 1st dimension of an image is time, then only
 * time shifts are computed.
 * <p>
 * Constraints are placed on strains, the rates at which shifts change in any
 * dimension. For example, strain r1 = du/d1 is the derivative of shift
 * u(x1,x2,...) with respect to the x1 coordinate, and is constrained to lie
 * between lower and upper bounds r1min and r1max. Default bounds are 
 * r1min = -1.0 and r1max = 1.0.
 * <p>
 * In many applications of warping, strains derived from estimated shifts may
 * be an important by-product. However, when computing shifts, only a finite
 * number of strains are permitted, and this quantization may yield
 * unrealistic estimates of strain. To address this problem, this dynamic
 * warping provides control over the sampling of strains, in the form of a
 * smoothness value. The number of strains sampled is proportional to this
 * value.
 * <p>
 * Smoothness values represent approximate intervals for a quasi-uniform
 * subsampling grid on which shifts are computed. For example, for a time
 * sampling interval of 4 ms, one might specify a smoothness of 200 ms, so
 * that shifts are computed on a grid that is 50 times more coarse than the 4
 * ms sampling grid. However, strains on this coarse grid would be sampled 50
 * times more finely. After initially computing shifts on the coarse grid,
 * shifts are interpolated onto the finer grid.
 *
 * @author Dave Hale, Colorado School of Mines
 * @author Elias Arias, Colorado School of Mines
 * @version 25.3.2015
 */
public class DynamicWarpingK {

  /**
   * Constructs a dynamic warping.
   * If this warping is used for 2D or 3D images, then default unit samplings
   * are assumed for 2nd and 3rd dimensions.
   * @param k discretization value for shift sampling interval.
   * @param smin lower bound on shift.
   * @param smax upper bound on shift.
   * @param s1 sampling of shifts for 1st dimension.
   */
  public DynamicWarpingK(int k, double smin, double smax, Sampling s1) {
    this(k,smin,smax,s1,null,null);
  }

  /**
   * Constructs a dynamic warping.
   * If this warping is used for 3D images, then default unit samplings
   * are assumed for the 3rd dimensions.
   * @param k discretization value for shift sampling interval.
   * @param smin lower bound on shift.
   * @param smax upper bound on shift.
   * @param s1 sampling of shifts for 1st dimension.
   * @param s2 sampling of shifts for 2nd dimension.
   */
  public DynamicWarpingK(int k, double smin, double smax, 
                         Sampling s1, Sampling s2) {
    this(k,smin,smax,s1,s2,null);
  }

  /**
   * Constructs a dynamic warping.
   * @param k discretization value for shift sampling interval.
   * @param smin lower bound on shift.
   * @param smax upper bound on shift.
   * @param s1 sampling of shifts for 1st dimension.
   * @param s2 sampling of shifts for 2nd dimension.
   * @param s3 sampling of shifts for 3rd dimension.
   */
  public DynamicWarpingK(
    int k, double smin, double smax, 
    Sampling s1, Sampling s2, Sampling s3) 
  {
    double ds = s1.getDelta()/k; // shift sampling interval
    int ismin = (int) ceil(smin/ds);
    int ismax = (int)floor(smax/ds);
    _ss = new Sampling(1+ismax-ismin,ds,ismin*ds);
    _s1 = s1;
    _s2 = s2;
    _s3 = s3;
    _r1min = -1.0;
    _r2min = -1.0;
    _r3min = -1.0;
    _r1max =  1.0;
    _r2max =  1.0;
    _r3max =  1.0;
    _k1min = 10;
    _k2min = 10;
    _k3min = 10;
    _si = new SincInterpolator();
    _si.setExtrapolation(SincInterpolator.Extrapolation.CONSTANT);
  }

  /**
   * Gets the sampling of shifts s used to compute alignment errors.
   * @return the sampling.
   */
  public Sampling getSamplingS() {
    return _ss;
  }

  /**
   * Gets the sampling in the 1st dimension.
   * @return the sampling.
   */
  public Sampling getSampling1() {
    return _s1;
  }

  /**
   * Gets the sampling in the 2nd dimension.
   * @return the sampling.
   */
  public Sampling getSampling2() {
    return _s2;
  }

  /**
   * Gets the sampling in the 3rd dimension.
   * @return the sampling.
   */
  public Sampling getSampling3() {
    return _s3;
  }

  /**
   * Sets bounds on strains for this dynamic warping.
   * Default lower and upper bounds are -1.0 and 1.0, respectively.
   * @param r1min lower bound on strain in 1st dimension.
   * @param r1max upper bound on strain in 1st dimension.
   */
  public void setStrainLimits(double r1min, double r1max) {
    setStrainLimits(r1min,r1max,-1.0,1.0,-1.0,1.0);
  }

  /**
   * Sets bounds on strains for this dynamic warping.
   * Default lower and upper bounds are -1.0 and 1.0, respectively.
   * @param r1min lower bound on strain in 1st dimension.
   * @param r1max upper bound on strain in 1st dimension.
   * @param r2min lower bound on strain in 2nd dimension.
   * @param r2max upper bound on strain in 2nd dimension.
   */
  public void setStrainLimits(
    double r1min, double r1max,
    double r2min, double r2max)
  {
    setStrainLimits(r1min,r1max,r2min,r2max,-1.0,1.0);
  }

  /**
   * Sets bounds on strains for this dynamic warping.
   * Default lower and upper bounds are -1.0 and 1.0, respectively.
   * @param r1min lower bound on strain in 1st dimension.
   * @param r1max upper bound on strain in 1st dimension.
   * @param r2min lower bound on strain in 2nd dimension.
   * @param r2max upper bound on strain in 2nd dimension.
   * @param r3min lower bound on strain in 3rd dimension.
   * @param r3max upper bound on strain in 3rd dimension.
   */
  public void setStrainLimits(
    double r1min, double r1max,
    double r2min, double r2max,
    double r3min, double r3max)
  {
    _r1min = r1min; _r1max = r1max;
    _r2min = r2min; _r2max = r2max;
    _r3min = r3min; _r3max = r3max;
  }

  /**
   * Sets the number of nonlinear smoothings of alignment errors.
   * In dynamic warping, alignment errors are smoothed the specified 
   * number of times, along all dimensions (in order 1, 2, ...), 
   * before estimating shifts by accumulating and backtracking along 
   * only the 1st dimension. 
   * <p> 
   * The default number of smoothings is zero, which is best for 1D
   * sequences. For 2D and 3D images, two smoothings are recommended.
   * @param esmooth number of nonlinear smoothings.
   */
  public void setErrorSmoothing(int esmooth) {
    _esmooth = esmooth;
  }

  /**
   * Sets the smoothness of shifts computed by this dynamic warping.
   * <em>
   * Units of smoothness are the same as those used in the corresponding
   * sampling of the 1st dimension.
   * </em>
   * Default smoothness values are 10 times the sampling intervals for
   * corresponding dimensions.
   * @param d1min smoothness in 1st dimension.
   */
  public void setSmoothness(double d1min) {
    double d2 = (_s2!=null)?_s2.getDelta():1.0;
    setSmoothness(d1min,10.0*d2);
  }

  /**
   * Sets the smoothness of shifts computed by this dynamic warping.
   * <em>
   * Units of smoothness are the same as those used in corresponding
   * samplings of 1st and 2nd dimensions.
   * </em>
   * Default smoothness values are 10 times the sampling intervals for
   * corresponding dimensions.
   * @param d1min smoothness in 1st dimension.
   * @param d2min smoothness in 2nd dimension.
   */
  public void setSmoothness(double d1min, double d2min) {
    double d3 = (_s3!=null)?_s3.getDelta():1.0;
    setSmoothness(d1min,d2min,10.0*d3);
  }

  /**
   * Sets the smoothness of shifts computed by this dynamic warping.
   * <em>
   * Units of smoothness are the same as those used in corresponding
   * samplings of 1st, 2nd and 3rd dimensions.
   * </em>
   * Default smoothness values are 10 times the sampling intervals for
   * corresponding dimensions.
   * @param d1min smoothness in 1st dimension.
   * @param d2min smoothness in 2nd dimension.
   * @param d3min smoothness in 3rd dimension.
   */
  public void setSmoothness(double d1min, double d2min, double d3min) {
    double d1 = (_s1!=null)?_s1.getDelta():1.0;
    double d2 = (_s2!=null)?_s2.getDelta():1.0;
    double d3 = (_s3!=null)?_s3.getDelta():1.0;
    _k1min = max(1,(int)ceil(d1min/d1));
    _k2min = max(1,(int)ceil(d2min/d2));
    _k3min = max(1,(int)ceil(d3min/d3));
  }

  /**
   * Returns shifts computed for specified 1D sequences.
   * @param sf sampling of 1st dimension for the seqeunce f.
   * @param f array of values for sequence f.
   * @param sg sampling of 1st dimension for the seqeunce g.
   * @param g array of values for sequence g.
   * @return array of shifts.
   */
  public float[] findShifts(
    Sampling sf, float[] f,
    Sampling sg, float[] g)
  {
    float[][] e = computeErrors(sf,f,sg,g);
    return findShifts(e);
  }

  /**
   * Returns shifts computed for specified 2D images.
   * @param sf sampling of 1st dimension for the image f.
   * @param f array of values for image f.
   * @param sg sampling of 1st dimension for the image g.
   * @param g array of values for image g.
   * @return array of shifts.
   */
  public float[][] findShifts(
    final Sampling sf, final float[][] f,
    final Sampling sg, final float[][] g)
  {
    // Samplings.
    final Sampling ss = _ss;
    final Sampling s1 = _s1;
    final Sampling s2 = sampling2(f);
    final int n1 = s1.getCount();
    final int n2 = s2.getCount();

    // Quasi-uniform subsamplings.
    final int[] k1s = subsample(n1,_k1min);
    final int[] k2s = subsample(n2,_k2min);
    final int nk1 = k1s.length;
    final int nk2 = k2s.length;

    trace("findShifts: smoothing in 1st dimension ...");
    final float[][][] ek = new float[n2][][];
    Parallel.loop(n2,new Parallel.LoopInt() {
    public void compute(int i2) {
      float[][] e1 = computeErrors(sf,f[i2],sg,g[i2]);
      ek[i2] = subsampleErrors(_r1min,_r1max,k1s,ss,s1,e1);
    }});
    normalizeErrors(ek);

    trace("findShifts: smoothing in 2nd dimension ...");
    float[][][] ekk = smoothErrors2(_r2min,_r2max,k2s,ss,s2,ek);
    normalizeErrors(ekk);

    /*
    float[][][] es = new float[nk2][][];
    for (int is=0; is<_esmooth-1; ++is) {

    }
    */

    trace("findShifts: finding shifts ...");
    float[][] ukk = new float[nk2][];
    for (int ik2=0; ik2<nk2; ++ik2) {
      ukk[ik2] = findShiftsFromSubsampledErrors(
        _r1min,_r1max,k1s,ss,s1,ekk[ik2]);
    }
    trace("findShifts: interpolating shifts ...");
    float[][] u = interpolateShifts(s1,s2,k1s,k2s,ukk);
    trace("findShifts: ... done");
    return u;
  }

  /**
   * Returns shifts computed for specified 2D images.
   * @param sf sampling of 1st dimension for the image f.
   * @param f array of values for image f.
   * @param sg sampling of 1st dimension for the image g.
   * @param g array of values for image g.
   * @return array of shifts.
   */
  public float[][][] findShifts(
    final Sampling sf, final float[][][] f,
    final Sampling sg, final float[][][] g)
  {
    // Samplings.
    final Sampling ss = _ss;
    final Sampling s1 = _s1;
    final Sampling s2 = sampling2(f);
    final Sampling s3 = sampling3(f); 
    final int n1 = s1.getCount();
    final int n2 = s2.getCount();
    final int n3 = s3.getCount();

    // Quasi-uniform subsamplings.
    final int[] k1s = subsample(n1,_k1min);
    final int[] k2s = subsample(n2,_k2min);
    final int[] k3s = subsample(n3,_k3min);
    final int nk1 = k1s.length;
    final int nk2 = k2s.length;
    final int nk3 = k3s.length;

    trace("findShifts: smoothing in 1st dimension ...");
    final float[][][][] ek = new float[n3][n2][][]; 
    int n23 = n2*n3;
    Parallel.loop(n23,new Parallel.LoopInt() {
    public void compute(int i23) {
      int i2 = i23%n2;
      int i3 = i23/n2;
      float[][] e1 = computeErrors(sf,f[i3][i2],sg,g[i3][i2]);
      ek[i3][i2] = subsampleErrors(_r1min,_r1max,k1s,ss,s1,e1);
    }});
    normalizeErrors(ek);

    trace("findShifts: smoothing in 2nd dimension ...");
    float[][][][] ekk = smoothErrors2(_r2min,_r2max,k2s,ss,s2,ek);
    normalizeErrors(ekk);

    trace("findShifts: smoothing in 3rd dimension ...");
    final float[][][][] ekkk = smoothErrors3(_r3min,_r3max,k3s,ss,s3,ekk);
    normalizeErrors(ekkk);

    float[][][][] es = ekkk;
    for (int is=0; is<_esmooth-1; ++is) {
      smoothSubsampledErrors(_r1min,_r1max,k1s,
                                   _r2min,_r2max,k2s,
                                   _r3min,_r3max,k3s,ss,s1,s2,s3,es);
      normalizeErrors(es);
    }
    final float[][][][] e = es;

    trace("findShifts: finding shifts ...");
    final float[][][] ukk = new float[nk3][nk2][];
    int nk23 = nk2*nk3;
    Parallel.loop(nk23,new Parallel.LoopInt() {
    public void compute(int ik23) {
      int ik2 = ik23%nk2;
      int ik3 = ik23/nk2;
      ukk[ik3][ik2] = findShiftsFromSubsampledErrors(
        _r1min,_r1max,k1s,ss,s1,e[ik3][ik2]);
    }});

    trace("findShifts: interpolating shifts ...");
    float[][][] u = interpolateShifts(s1,s2,s3,k1s,k2s,k3s,ukk);
    trace("findShifts: ... done");
    return u;
  }

  /**
   * Returns alignment errors computed for specified sequences.
   * @param sf sampling of 1st dimension for the seqeunce f.
   * @param f array of values for sequence f.
   * @param sg sampling of 1st dimension for the seqeunce g.
   * @param g array of values for sequence g.
   * @return array of alignment errors.
   */
  public float[][] computeErrors(
    Sampling sf, float[] f,
    Sampling sg, float[] g)
  {
    //Sampling ss = _ss;
    //Sampling se = _s1;
    int ns = _ss.getCount();
    int ne = _s1.getCount();
    int nf = sf.getCount();
    int ng = sg.getCount();
    float[][] e = new float[ne][ns];
    float[] fi = new float[ne];
    float[] gi = new float[ne];
    _si.interpolate(sf,f,_s1,fi);
    float sum = 0;
    for (int is=0; is<ns; ++is) {
      _si.interpolate(
        ng,sg.getDelta(),sg.getFirst(),g,
        ne,_s1.getDelta(),_s1.getFirst()+_ss.getValue(is),gi);
      for (int ie=0; ie<ne; ++ie) {
        e[ie][is] = error(f[ie],gi[ie]);
      }
    }
    return e;
  }

  /**
   * Normalizes alignment errors, in place.
   * After normalizing, minimum error = 0 and maximum error = 1.
   */
  public static void normalizeErrors(float[][] e) {
    int ns = e[0].length;
    int ne = e.length;
    float emin = e[0][0];
    float emax = e[0][0];
    for (int ie=0; ie<ne; ++ie) {
      for (int is=0; is<ns; ++is) {
        float ei = e[ie][is];
        if (ei<emin) emin = ei;
        if (ei>emax) emax = ei;
      }
    }
    shiftAndScale(emin,emax,e);
  }

  /**
   * Returns shifts estimated for specified alignment errors.
   * @param e array[n1][ns] of alignment errors.
   * @return array[n1] of shifts.
   */
  public float[] findShifts(float[][] e) {
    int ns = _ss.getCount();
    int n1 = _s1.getCount();
    double ds = _ss.getDelta();
    double d1 = _s1.getDelta();
    int k1min = min(_k1min,n1-1);
    int[] i1k = subsample(n1,k1min);
    //dump(i1k);
    int n1k = i1k.length;
    return findShiftsFromErrors(_r1min,_r1max,i1k,_ss,_s1,e);
  }

  /**
   * Returns uniformly sampled warped sequence h(x1) = g(x1+u(x1)).
   * @param sg sampling of the sequence g to be warped.
   * @param g array for the sequence g to be warped.
   * @param u array of shifts.
   * @return array for the warped sequence h.
   */
  public float[] applyShifts(Sampling sg, float[] g, float[] u) {
    Sampling s1 = _s1;
    int ng = sg.getCount();
    int n1 = s1.getCount();
    float[] h = new float[n1];
    for (int i1=0; i1<n1; ++i1) {
      double x1 = s1.getValue(i1)+u[i1];
      h[i1] = _si.interpolate(sg,g,x1);
    }
    return h;
  }

  /**
   * Returns uniformly sampled warped image h(x1,x2) = g(x1+u(x1,x2),x2).
   * @param sg sampling of the sequence g to be warped.
   * @param g array for the sequence g to be warped.
   * @param u array of shifts.
   * @return array for the warped sequence h.
   */
  public float[][] applyShifts(Sampling sg, float[][] g, float[][] u) {
    int n2 = g.length;
    float[][] h = new float[n2][];
    for (int i2=0; i2<n2; ++i2)
      h[i2] = applyShifts(sg,g[i2],u[i2]);
    return h;
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private Sampling _ss,_s1,_s2,_s3;
  private double _r1min,_r2min,_r3min;
  private double _r1max,_r2max,_r3max;
  private int _k1min,_k2min,_k3min;
  private int _esmooth = 1;
  private SincInterpolator _si;
  private float _epow = 1.00f;

  private static CubicInterpolator makeInterpolator1(
    float[] x, float[] y) 
  {
    //return new CubicInterpolator(CubicInterpolator.Method.LINEAR,x,y);
    return new CubicInterpolator(CubicInterpolator.Method.SPLINE,x,y);
  }
  private static CubicInterpolator makeInterpolator2(
    float[] x, float[] y) 
  {
    return new CubicInterpolator(CubicInterpolator.Method.SPLINE,x,y);
  }
  private static CubicInterpolator makeCubicInterpolator(
    float[] x, float[] y, float[] yd) 
  {
    return new CubicInterpolator(x,y,yd);
  }

  private static void trace(String s) {
    System.out.println(s);
  }

  private float error(float f, float g) {
    return pow(abs(f-g),_epow);
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
  private static int[] subsample(int n, int kmin) {
    if (kmin>=n)
      kmin = n-1;
    int m = 1+(n-1)/kmin;
    double d = (double)(n-1)/(double)(m-1);
    int[] j = new int[m];
    for (int i=0; i<m; ++i)
      j[i] = (int)(i*d+0.5);
    return j;
  }
  private static void subsampleTest() {
    int kmin = 3;
    for (int n=kmin+1; n<5*kmin+1; ++n) {
      int[] j = subsample(n,kmin);
      System.out.println("n="+n);
      dump(j);
    }
  }
  private static void main(String[] args) {
    subsampleTest();
  }

  /**
   * Shifts and scales alignment errors to be in range [0,1].
   * @param emin minimum alignment error before normalizing.
   * @param emax maximum alignment error before normalizing.
   * @param e input/output array of alignment errors.
   */
  private static void shiftAndScale(float emin, float emax, float[][] e) {
    int nl = e[0].length;
    int n1 = e.length;
    float eshift = emin;
    float escale = (emax>emin)?1.0f/(emax-emin):1.0f;
    for (int i1=0; i1<n1; ++i1) {
      for (int il=0; il<nl; ++il) {
        e[i1][il] = (e[i1][il]-eshift)*escale;
      }
    }
  }

  /**
   * Subsamples alignment errors.
   * @param rmin lower bound on strain.
   * @param rmax upper bound on strain.
   * @param kes array[nke] of indices for quasi-uniform subsampling.
   * @param ss uniform sampling of ns shifts.
   * @param se uniform sampling of ne errors.
   * @param e input array[ne][ns] of alignment errors.
   * @return array[nke][ns] of subsampled errors.
   */
  private static float[][] subsampleErrors(
    double rmin, double rmax, int[] kes, 
    Sampling ss, Sampling se, float[][] e) 
  {
    int ns = ss.getCount();
    int ne = se.getCount();
    int nke = kes.length;
    float[][] df = new float[nke][ns];
    float[][] dr = new float[nke][ns];
    accumulate( 1,rmin,rmax,kes,ss,se,e,df,null);
    accumulate(-1,rmin,rmax,kes,ss,se,e,dr,null);
    float[][] d = df;
    float scale = 1.0f/ne;
    for (int ike=0; ike<nke; ++ike) {
      int ke = kes[ike];
      for (int is=0; is<ns; ++is) {
        d[ike][is] = scale*(df[ike][is]+dr[ike][is]-e[ke][is]);
      }
    }
    return d;
  }

  /**
   * Returns alignment errors smoothed in the second dimension.
   * Returned errors are sparse in the second dimension, and
   * unchanged in the first dimension.
   * @param r2min minimum strain the second dimension.
   * @param r2max maximum strain in the second dimension.
   * @param k2s second dimension sparse grid indices.
   * @param ss uniform sampling of ns shifts.
   * @param se uniform sampling of ne errors.
   * @param e alignment errors.
   * @return smoothed alignment errors with size
   *  [k2s.length][e[0].length][e[0][0].length].
   */
  private static float[][][] smoothErrors2(
      final double r2min, final double r2max, final int[] k2s, 
      final Sampling ss, final Sampling se, final float[][][] e) {
    final int ns = ss.getCount();
    final int n2 = e.length;
    final int nk1 = e[0].length;
    final int nk2 = k2s.length;
    final float[][][] es = new float[nk2][nk1][ns]; // smoothed errors
    Parallel.loop(nk1,new Parallel.LoopInt() {
    public void compute(int ik1) {
      float[][]  e2 = new float[n2][ns]; // errors at index ik1
      for (int i2=0; i2<n2; ++i2)
        e2[i2] = e[i2][ik1];
      float[][] es2 = subsampleErrors(r2min,r2max,k2s,ss,se,e2);
      for (int ik2=0; ik2<nk2; ++ik2)
        es[ik2][ik1] = es2[ik2];
    }});
    return es;
  }

  /**
   * Returns alignment errors smoothed in the second dimension.
   * Returned errors are sparse in the second dimension, and
   * unchanged in the first dimension and third dimension.
   * @param r2Min minimum strain in the second dimension.
   * @param r2Max maximum strain in the second dimension.
   * @param k2s second dimension sparse grid indices.
   * @param ss uniform sampling of ns shifts.
   * @param se uniform sampling of ne errors.
   * @param e alignment errors.
   * @return smoothed alignment errors with size
   *  [e.length][k2s.length][e[0][0].length][e[0][0][0].length].
   */
  private static float[][][][] smoothErrors2(
      final double r2min, final double r2max, final int[] k2s, 
      final Sampling ss, final Sampling se, final float[][][][] e) {
    final int n3 = e.length;
    final float[][][][] es = new float[n3][][][]; // smoothed errors
    for (int i3=0; i3<n3; ++i3)
      es[i3] = smoothErrors2(r2min,r2max,k2s,ss,se,e[i3]);
    return es;
  }

  /**
   * Returns alignment errors smoothed in the third dimension.
   * Returned errors are sparse in the third dimension, and
   * unchanged in the first and second dimension.
   * @param r3Min minimum strain in the third dimension.
   * @param r3Max maximum strain in the third dimension.
   * @param k3s third dimension sparse grid indices.
   * @param ss uniform sampling of ns shifts.
   * @param se uniform sampling of ne errors.
   * @param e alignment errors.
   * @return smoothed alignment errors with size
   *  [k3s.length][e[0].length][e[0][0].length][e[0][0][0].length].
   */
  private static float[][][][] smoothErrors3(
      final double r3min, final double r3max, final int[] k3s, 
      final Sampling ss, final Sampling se, final float[][][][] e) {
    final int ns = e[0][0][0].length;
    final int nk1 = e[0][0].length;
    final int nk2 = e[0].length;
    final int n3 = e.length;
    final int nk3 = k3s.length;
    final float[][][][] es = new float[nk3][nk2][nk1][ns]; // smoothed errors
    Parallel.loop(nk1,new Parallel.LoopInt() {
    public void compute(int ik1) {
      for (int ik2=0; ik2<nk2; ++ik2) {
        float[][] e3 = new float[n3][ns]; // smooth errors at index i1,i2
        for (int i3=0; i3<n3; ++i3)
          e3[i3] = e[i3][ik2][ik1];
        float[][] es3 = subsampleErrors(r3min,r3max,k3s,ss,se,e3);
        for (int ik3=0; ik3<nk3; ++ik3)
          es[ik3][ik2][ik1] = es3[ik3];
      }
    }});
    return es;
  }

  private static void smoothSubsampledErrors(
    double rmin, double rmax, int[] kes, 
    Sampling ss, Sampling se, float[][] e) 
  {
    int ns = ss.getCount();
    int nke = kes.length;
    float[][] ef = new float[nke][ns];
    float[][] er = new float[nke][ns];
    accumulateSubsampled( 1,rmin,rmax,kes,ss,se,e,ef,null);
    accumulateSubsampled(-1,rmin,rmax,kes,ss,se,e,er,null);
    float scale = 1.0f/nke;
    for (int ike=0; ike<nke; ++ike) {
      int ke = kes[ike];
      for (int is=0; is<ns; ++is) {
        e[ike][is] = scale*(ef[ike][is]+er[ike][is]-e[ike][is]);
      }
    }
  }

  private static void smoothSubsampledErrors(
    double r1min, double r1max, int[] k1s,
    double r2min, double r2max, int[] k2s,
    Sampling ss, Sampling s1, Sampling s2, float[][][] e) 
  {
    int ns = ss.getCount();
    int nk1 = k1s.length;
    int nk2 = k2s.length;
    float[][] e2 = new float[nk2][];
    for (int ik1=0; ik1<nk1; ++ik1) {
      for (int ik2=0; ik2<nk2; ++ik2)
        e2[ik2] = e[ik2][ik1];
      smoothSubsampledErrors(r2min,r2max,k2s,ss,s2,e2);
      for (int ik2=0; ik2<nk2; ++ik2)
        e[ik2][ik1] = e2[ik2];
    }
    for (int ik2=0; ik2<nk2; ++ik2) {
      smoothSubsampledErrors(r1min,r1max,k1s,ss,s1,e[ik2]);
    }
  }

  private static void smoothSubsampledErrors(
    double r1min, double r1max, int[] k1s,
    double r2min, double r2max, int[] k2s,
    double r3min, double r3max, int[] k3s,
    Sampling ss, Sampling s1, Sampling s2, Sampling s3, float[][][][] e) 
  {
    int ns = ss.getCount();
    int nk1 = k1s.length;
    int nk2 = k2s.length;
    int nk3 = k3s.length;
    float[][] e3 = new float[nk3][];
    for (int ik1=0; ik1<nk1; ++ik1) {
      for (int ik2=0; ik2<nk2; ++ik2) {
        for (int ik3=0; ik3<nk3; ++ik3)
          e3[ik3] = e[ik3][ik2][ik1];
        smoothSubsampledErrors(r3min,r3max,k3s,ss,s3,e3);
        for (int ik3=0; ik3<nk3; ++ik3)
          e[ik3][ik2][ik1] = e3[ik3];
      }
    }
    for (int ik3=0; ik3<nk3; ++ik3) {
      smoothSubsampledErrors(r1min,r1max,k1s,r2min,r2max,k2s,ss,s1,s2,e[ik3]);
    }
  }

  /**
   * Finds shifts from alignment errors.
   * @param rmin lower bound on strain.
   * @param rmax upper bound on strain.
   * @param kes array[nke] of indices for quasi-uniform subsampling.
   * @param ss uniform sampling of ns shifts.
   * @param se uniform sampling of ne errors.
   * @param e input array[ne][ns] of alignment errors.
   * @return array[ne] of shifts.
   */
  private static float[] findShiftsFromErrors(
    double rmin, double rmax, int[] kes,
    Sampling ss, Sampling se, float[][] e) 
  {
    int nke = kes.length;
    int ns = ss.getCount();
    float[][] d = new float[nke][ns];
    int[][] m = new int[nke][ns];
    accumulate(1,rmin,rmax,kes,ss,se,e,d,m);
    float[] uke = backtrackForShifts(kes,ss,se,d[nke-1],m);
    return interpolateShifts(se,kes,uke);
  }

  /**
   * Finds shifts from subsampled alignment errors.
   * @param rmin lower bound on strain.
   * @param rmax upper bound on strain.
   * @param kes array[nke] of indices for quasi-uniform subsampling.
   * @param ss uniform sampling of ns shifts.
   * @param se uniform sampling of ne errors.
   * @param e array[nke][ns] of subsampled alignment errors.
   * @return array[ne] of shifts.
   */
  private static float[] findShiftsFromSubsampledErrors(
    double rmin, double rmax, int[] kes,
    Sampling ss, Sampling se, float[][] e) 
  {
    int nke = kes.length;
    int ns = ss.getCount();
    float[][] d = new float[nke][ns];
    int[][] m = new int[nke][ns];
    accumulateSubsampled(1,rmin,rmax,kes,ss,se,e,d,m);
    return backtrackForShifts(kes,ss,se,d[nke-1],m);
  }

  /**
   * Accumulates alignment errors in forward or reverse direction.
   * @param dir direction, 1 for forward, -1 for reverse.
   * @param rmin lower bound on strain.
   * @param rmax upper bound on strain.
   * @param kes array[nke] of indices for quasi-uniform subsampling.
   * @param ss uniform sampling of ns shifts.
   * @param se uniform sampling of ne errors.
   * @param e input array[ne][ns] of alignment errors.
   * @param d output array[nke][ns] of accumulated errors.
   * @param m output array[nke][ns] of minimizing moves; or null.
   */
  private static void accumulate(
    int dir, double rmin, double rmax, int[] kes, 
    Sampling ss, Sampling se, float[][] e, float[][] d, int[][] m) 
  {
    int ns = ss.getCount();
    double ds = ss.getDelta();
    double de = se.getDelta();
    int nke = kes.length;
    int iked = dir>0?1:-1;
    int ikeb = dir>0?0:nke-1;
    int ikee = dir>0?nke:-1;
    for (int is=0; is<ns; ++is)
      d[ikeb][is] = e[kes[ikeb]][is];
    float[] dprev = new float[ns];
    for (int ike=ikeb+iked; ike!=ikee; ike+=iked) {
      int je = kes[ike-iked];
      int ie = kes[ike];
      int me = ie-je;
      int msmin,msmax;
      if (me>0) {
        msmin = (int) ceil(-rmax*me*de/ds);
        msmax = (int)floor(-rmin*me*de/ds);
      } else {
        msmin = (int) ceil(-rmin*me*de/ds);
        msmax = (int)floor(-rmax*me*de/ds);
      }
      if (msmin>msmax) {
        trace("ie="+ie+" je="+je+" me="+me+" msmin="+msmin+" msmax="+msmax);
      }
      assert msmin<=msmax:"msmin<=msmax";
      fill(Float.MAX_VALUE,d[ike]);
      for (int ms=msmin; ms<=msmax; ++ms) {
        int islo = max(0,-ms);
        int ishi = min(ns,ns-ms);
        for (int is=islo; is<ishi; ++is)
          dprev[is] = d[ike-iked][is+ms];
        if (m!=null)
          updateSumsOfErrors(ie,je,ms,e,dprev,d[ike],m[ike]);
        else
          updateSumsOfErrors(ie,je,ms,e,dprev,d[ike],null);
      }
    }
  }

  /**
   * Accumulates alignment errors in forward or reverse direction.
   * Does not subsample the accumulated errors.
   * @param dir direction, 1 for forward, -1 for reverse.
   * @param rmin lower bound on strain.
   * @param rmax upper bound on strain.
   * @param me number of errors e summed per accumulated error d
   * @param ss uniform sampling of ns shifts.
   * @param se uniform sampling of ne errors.
   * @param e input array[ne][ns] of alignment errors.
   * @param d output array[ne][ns] of accumulated errors.
   */
  private static void accumulate(
    int dir, double rmin, double rmax, int me, 
    Sampling ss, Sampling se, float[][] e, float[][] d) 
  {
    int ns = ss.getCount();
    int ne = se.getCount();
    double ds = ss.getDelta();
    double de = se.getDelta();
    int ied = dir>0?1:-1;
    int ieb = dir>0?0:ne-1;
    int iee = dir>0?ne:-1;
    for (int is=0; is<ns; ++is)
      d[ieb][is] = e[ieb][is];
    float[] dprev = new float[ns];
    for (int ie=ieb+ied; ie!=iee; ie+=ied) {
      int je = max(0,min(ne-1,ie-dir*me));
      int ke = ie-je;
      int msmin,msmax;
      if (ke>0) {
        msmin = (int) ceil(-rmax*ke*de/ds);
        msmax = (int)floor(-rmin*ke*de/ds);
      } else {
        msmin = (int) ceil(-rmin*ke*de/ds);
        msmax = (int)floor(-rmax*ke*de/ds);
      }
      fill(Float.MAX_VALUE,d[ie]);
      for (int ms=msmin; ms<=msmax; ++ms) {
        int islo = max(0,-ms);
        int ishi = min(ns,ns-ms);
        for (int is=islo; is<ishi; ++is)
          dprev[is] = d[ie-ied][is+ms];
        updateSumsOfErrors(ie,je,ms,e,dprev,d[ie],null);
      }
    }
  }

  /**
   * Updates sums of errors for one shift between two error sample indices.
   * Computes sums of errors along linear trajectories according to:
   * <pre>
   * d[is] += sum from ke=ie to ke!=je of e[ke][is+(ie-ke)*ms/(ie-je)]
   * </pre>
   * Where the complicated last subscript (typically) is not an integer, this
   * method uses linear interpolation of the alignment errors e. After the
   * sums of errors have been computed for all shift indices is, this method
   * updates the minimum sum of errors and the corresponding change in shift.
   * @param ie error sample index at which to begin sum.
   * @param je error sample index at which to end (not) sum.
   * @param ms change in shift at error sample index je, not in sum.
   * @param e[ne][ns] input array of alignment errors.
   * @param d[ns] input/output array in which to accumulate errors.
   * @param dmin[ns] input/output array of minimum accumulated errors.
   * @param mmin[ns] input/output array of minimizing moves; or null.
   */
  private static void updateSumsOfErrors(
    int ie, int je, int ms, float[][] e, float[] d, float[] dmin, int[] mmin) 
  {
    int ns = d.length;
    int islo = max(0,-ms); // update only for is >= islo
    int ishi = min(ns,ns-ms); // update only for is < ishi
    for (int is=islo; is<ishi; ++is)
      d[is] += e[ie][is];
    int me = ie-je;
    int de = me>0?-1:1;
    if (ms==0) { // if no shift, no interpolation required
      for (int ke=ie+de; ke!=je; ke+=de) {
        for (int is=islo; is<ishi; ++is) {
          d[is] += e[ke][is];
        }
      }
    } else { // else, use linear interpolation of errors
      float r = (float)ms/(float)me; // strain
      for (int ke=ie+de; ke!=je; ke+=de) {
        float sk = r*(ie-ke);
        int ks = (int)sk;
        if (sk<0.0f) --ks;
        int ksa = ks+1;
        int ksb = ks;
        float wsa = sk-ks;
        float wsb = 1.0f-wsa;
        for (int is=islo; is<ishi; ++is)
          d[is] += wsa*e[ke][is+ksa]+wsb*e[ke][is+ksb];
      }
    }
    for (int is=islo; is<ishi; ++is) {
      if (d[is]<dmin[is]) {
        dmin[is] = d[is];
        if (mmin!=null)
          mmin[is] = ms;
      }
    }
  }

  /**
   * Accumulates subsampled errors in forward or reverse direction.
   * @param dir direction, 1 for forward, -1 for reverse.
   * @param rmin lower bound on strain.
   * @param rmax upper bound on strain.
   * @param kes array[nke] of indices for quasi-uniform subsampling.
   * @param ss uniform sampling of ns shifts.
   * @param se uniform sampling of ne errors.
   * @param e input array[nke][ns] of subsampled errors.
   * @param d output array[nke][ns] of accumulated errors.
   * @param m output array[nke][ns] of minimizing moves; or null.
   */
  private static void accumulateSubsampled(
    int dir, double rmin, double rmax, int[] kes, 
    Sampling ss, Sampling se, float[][] e, float[][] d, int[][] m) 
  {
    int ns = ss.getCount();
    double ds = ss.getDelta();
    double de = se.getDelta();
    int nke = kes.length;
    int iked = dir>0?1:-1;
    int ikeb = dir>0?0:nke-1;
    int ikee = dir>0?nke:-1;
    for (int is=0; is<ns; ++is)
      d[ikeb][is] = e[ikeb][is];
    for (int ike=ikeb+iked; ike!=ikee; ike+=iked) {
      float[] dprev = d[ike-iked];
      int me = kes[ike]-kes[ike-iked];
      int msmin,msmax;
      if (me>0) {
        msmin = (int) ceil(-rmax*me*de/ds);
        msmax = (int)floor(-rmin*me*de/ds);
      } else {
        msmin = (int) ceil(-rmin*me*de/ds);
        msmax = (int)floor(-rmax*me*de/ds);
      }
      for (int is=0; is<ns; ++is) {
        float dmin = Float.MAX_VALUE;
        int mmin = -1;
        for (int ms=msmin; ms<=msmax; ++ms) {
          int js = is+ms;
          if (0<=js && js<ns) {
            float dj = dprev[js];
            if (dj<dmin) {
              dmin = dj;
              mmin = ms;
            }
          }
          d[ike][is] = dmin+e[ike][is];
          if (m!=null)
            m[ike][is] = mmin;
        }
      }
    }
  }

  /**
   * Returns shifts found by backtracking with precomputed moves.
   * @param kes array[nke] of indices for quasi-uniform subsampling.
   * @param ss uniform sampling of ns shifts.
   * @param se uniform sampling of ne errors.
   * @param d array[ns] of last forward accumulated errors.
   * @param m array[nke][ns] of moves, changes in shifts.
   * @return array[ne] of shifts.
   */
  private static float[] backtrackForShifts(
    int[] kes, Sampling ss, Sampling se, float[] d, int[][] m) 
  {
    int nke = kes.length;
    int ns = ss.getCount();
    int ne = se.getCount();
    int ike = nke-1;
    float dmin = Float.MAX_VALUE;
    int imin = -1;
    for (int is=0; is<ns; ++is) {
      if (d[is]<dmin) {
        dmin = d[is];
        imin = is;
      }
    }
    int is = imin;
    float[] uke = new float[nke];
    uke[ike] = (float)ss.getValue(is);
    for (--ike; ike>=0; --ike) {
      is += m[ike+1][is];
      uke[ike] = (float)ss.getValue(is);
    }
    return uke;
  }

  private static float[] interpolateShifts(
    Sampling s1, int[] k1s, float[] uk) 
  {
    int n1 = s1.getCount();
    int nk1 = k1s.length;
    float[] xk1 = new float[nk1];
    for (int jk1=0; jk1<nk1; ++jk1)
      xk1[jk1] = (float)s1.getValue(k1s[jk1]);
    CubicInterpolator ci = makeInterpolator1(xk1,uk);
    float[] u = new float[n1];
    for (int j1=0; j1<n1; ++j1) {
      float x1 = (float)s1.getValue(j1);
      u[j1] = ci.interpolate(x1);
    }
    return u;
  }
  
  private static float[][] interpolateShifts(
    Sampling s1, Sampling s2, int[] k1s, int[] k2s, float[][] ukk) 
  {
    int n1 = s1.getCount();
    int n2 = s2.getCount();
    int nk1 = k1s.length;
    int nk2 = k2s.length;

    // Coarse sampling of 1st and 2nd dimensions.
    float[] xk1 = new float[nk1];
    for (int jk1=0; jk1<nk1; ++jk1)
      xk1[jk1] = (float)s1.getValue(k1s[jk1]);
    float[] xk2 = new float[nk2];
    for (int jk2=0; jk2<nk2; ++jk2)
      xk2[jk2] = (float)s2.getValue(k2s[jk2]);

    // Interpolate.
    BicubicInterpolator2 bc = new BicubicInterpolator2(
      BicubicInterpolator2.Method.MONOTONIC,
      BicubicInterpolator2.Method.SPLINE,
      xk1,xk2,ukk);
    float[][] u = new float[n2][n1];
    for (int j2=0; j2<n2; ++j2) {
      float x2 = (float)s2.getValue(j2);
      for (int j1=0; j1<n1; ++j1) {
        float x1 = (float)s1.getValue(j1);
        u[j2][j1] = bc.interpolate(x1,x2);
      }
    }
    return u;
  }

  /**
   * Added by Elias Arias, slightly modified from DynamicWarpingC 
   * by Stefan Compton.
   */
  private static float[][][] interpolateShifts(
      Sampling s1, Sampling s2, Sampling s3,
      int[] k1s, int[] k2s, int[] k3s, float[][][] ukk)
  {
    int n1 = s1.getCount();
    int n2 = s2.getCount();
    int n3 = s3.getCount();
    int nk1 = k1s.length;
    int nk2 = k2s.length;
    int nk3 = k3s.length;
    double d1 = s1.getDelta();
    double d2 = s2.getDelta();
    double d3 = s3.getDelta();
    double f1 = s1.getFirst();
    double f2 = s2.getFirst();
    double f3 = s3.getFirst();

     // Coarse sampling of 1st and 2nd dimensions.
    float[] xk1 = new float[nk1];
    for (int jk1=0; jk1<nk1; ++jk1)
      xk1[jk1] = (float)s1.getValue(k1s[jk1]);
    float[] xk2 = new float[nk2];
    for (int jk2=0; jk2<nk2; ++jk2)
      xk2[jk2] = (float)s2.getValue(k2s[jk2]);
    float[] xk3 = new float[nk3];
    for (int jk3=0; jk3<nk3; ++jk3)
      xk3[jk3] = (float)s3.getValue(k3s[jk3]);

    // Interpolate along the second and third dimension.
    float[][] uk23 = new float[nk3][nk2];
    float[][][] u23 = new float[n3][n2][nk1];
    for (int ik1=0; ik1<nk1; ++ik1) {
      for (int ik3=0; ik3<nk3; ++ik3)
        for (int ik2=0; ik2<nk2; ++ik2)
          uk23[ik3][ik2] = ukk[ik3][ik2][ik1];
      BicubicInterpolator2 bc = new BicubicInterpolator2(
        BicubicInterpolator2.Method.MONOTONIC,
        BicubicInterpolator2.Method.MONOTONIC,
        xk2,xk3,uk23);
      double v3 = f3;
      for (int i3=0; i3<n3; ++i3, v3=f3+i3*d3) {
        float v3f = (float)v3;
        double v2 = f2;
        for (int i2=0; i2<n2; ++i2, v2=f2+i2*d2)
          u23[i3][i2][ik1] = bc.interpolate((float)v2,v3f);
      }
    }

    // Interpolate along the first dimension.
    float[][][] u = new float[n3][n2][n1];
    for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        CubicInterpolator ci = 
          new CubicInterpolator(CubicInterpolator.Method.MONOTONIC,
                xk1,u23[i3][i2]);
        double v1 = f1;
        for (int i1=0; i1<n1; ++i1, v1=f1+i1*d1)
          u[i3][i2][i1] = ci.interpolate((float)v1);
      }
    }
    return u;
  }

  private Sampling sampling2(float[][] f) {
    if (_s2!=null) {
      Check.argument(f.length==_s2.getCount(),"valid sampling2");
      return _s2;
    } else {
      return new Sampling(f.length,1.0,0.0);
    }
  }
  private Sampling sampling2(float[][][] f) {
    if (_s2!=null) {
      Check.argument(f[0].length==_s2.getCount(),"valid sampling2");
      return _s2;
    } else {
      return new Sampling(f[0].length,1.0,0.0);
    }
  }
  private Sampling sampling3(float[][][] f) {
    if (_s3!=null) {
      Check.argument(f.length==_s3.getCount(),"valid sampling3");
      return _s3;
    } else {
      return new Sampling(f.length,1.0,0.0);
    }
  }

  /**
   * Normalizes alignment errors to be in range [0,1].
   * @param e input/output array of alignment errors.
   */
  private static void normalizeErrors(float[][][] e) {
    final int ns = e[0][0].length;
    final int n1 = e[0].length;
    final int n2 = e.length;
    final float[][][] ef = e;
    MinMax mm = Parallel.reduce(n2,new Parallel.ReduceInt<MinMax>() {
    public MinMax compute(int i2) {
      float emin =  Float.MAX_VALUE;
      float emax = -Float.MAX_VALUE;
      for (int i1=0; i1<n1; ++i1) {
        for (int is=0; is<ns; ++is) {
          float ei = ef[i2][i1][is];
          if (ei<emin) emin = ei;
          if (ei>emax) emax = ei;
        }
      }
      return new MinMax(emin,emax);
    }
    public MinMax combine(MinMax mm1, MinMax mm2) {
      return new MinMax(min(mm1.emin,mm2.emin),max(mm1.emax,mm2.emax));
    }});
    shiftAndScale(mm.emin,mm.emax,e);
  }

  /**
   * Normalizes alignment errors to be in range [0,1].
   * @param e input/output array of alignment errors.
   */
  private static void normalizeErrors(float[][][][] e) {
    final int ns = e[0][0][0].length;
    final int n1 = e[0][0].length;
    final int n2 = e[0].length;
    final int n3 = e.length;
    final float[][][][] ef = e;
    MinMax mm = Parallel.reduce(n3,new Parallel.ReduceInt<MinMax>() {
      public MinMax compute(int i3) {
        float emin =  Float.MAX_VALUE;
        float emax = -Float.MAX_VALUE;
        for (int i2=0; i2<n2; ++i2) {
          for (int i1=0; i1<n1; ++i1) {
            for (int is=0; is<ns; ++is) {
              float ei = ef[i3][i2][i1][is];
              if (ei<emin) emin = ei;
              if (ei>emax) emax = ei;
            }
          }
        }
        return new MinMax(emin,emax);
      }
      public MinMax combine(MinMax mm1, MinMax mm2) {
        return new MinMax(min(mm1.emin,mm2.emin),max(mm1.emax,mm2.emax));
      }});
      shiftAndScale(mm.emin,mm.emax,e);
    }

    private static class MinMax {
      float emin,emax;
      MinMax(float emin, float emax) {
        this.emin = emin;
        this.emax = emax;
      }
    }

  /**
   * Shifts and scales alignment errors to be in range [0,1].
   * @param emin minimum alignment error before normalizing.
   * @param emax maximum alignment error before normalizing.
   * @param e input/output array of alignment errors.
   */
  private static void shiftAndScale(float emin, float emax, float[][][] e) {
    trace("shiftAndScale: emin="+emin+" emax="+emax);
    final int ns = e[0][0].length;
    final int n1 = e[0].length;
    final int n2 = e.length;
    final float eshift = emin;
    final float escale = (emax>emin)?1.0f/(emax-emin):1.0f;
    final float[][][] ef = e;
    Parallel.loop(n2,new Parallel.LoopInt() {
    public void compute(int i2) {
      for (int i1=0; i1<n1; ++i1) {
        for (int is=0; is<ns; ++is) {
          ef[i2][i1][is] = (ef[i2][i1][is]-eshift)*escale;
        }
      }
    }});
  }

  /**
   * Shifts and scales alignment errors to be in range [0,1].
   * @param emin minimum alignment error before normalizing.
   * @param emax maximum alignment error before normalizing.
   * @param e input/output array of alignment errors.
   */
  private static void shiftAndScale(float emin, float emax, float[][][][] e) {
    trace("shiftAndScale: emin="+emin+" emax="+emax);
    final int ns = e[0][0][0].length;
    final int n1 = e[0][0].length;
    final int n2 = e[0].length;
    final int n3 = e.length;
    final float eshift = emin;
    final float escale = (emax>emin)?1.0f/(emax-emin):1.0f;
    final float[][][][] ef = e;
    Parallel.loop(n3,new Parallel.LoopInt() {
    public void compute(int i3) {
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          for (int is=0; is<ns; ++is) {
            ef[i3][i2][i1][is] = (ef[i3][i2][i1][is]-eshift)*escale;
          }
        }
      }
    }});
  }
}
