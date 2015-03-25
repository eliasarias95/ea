package slopes;

import util.Plot;
import util.Util;
import edu.mines.jtk.dsp.*;
import edu.mines.jtk.interp.*;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * @author Elias Arias, Colorado School of Mines, CWP
 * @version 19.2.2015
 */
public class DynamicWarpingSlopes {

  public DynamicWarpingSlopes(double pmax, Sampling s1, Sampling s2) {
    this(1,pmax,1.0,1.0,1.0,1.0,s1,s2);
  }

  public DynamicWarpingSlopes(int k, double pmax, double h1, double h2,
      double r1, double r2, Sampling s1, Sampling s2) {
    this(k,pmax,h1,h2,10,r1,r2,1,s1,s2,null);
  }

  public DynamicWarpingSlopes(int k, double pmax, 
      double h1, double h2, double h3,
      double r1, double r2, double r3,
      Sampling s1, Sampling s2, Sampling s3) {
    _k = k;
    _pmax = pmax;
    _h1 = h1;
    _h2 = h2;
    _h3 = h3;
    _r1 = r1;
    _r2 = r2;
    _r3 = r3;
    _s1 = s1;
    _s2 = s2;
    _s3 = s3;
    _dwk = new DynamicWarpingK(k,-pmax,pmax,s1,s2,s3);
    _dwk.setSmoothness(h1,h2,h3);
    _dwk.setStrainLimits(-r1,r1,-r2,r2,-r3,r3);
    _si = new SincInterpolator();
  }

  public void setK(int k) {
    _k = k;
  }

  public float[][] findSlopes(float[][] f) {
    int n2 = f.length;
    int n1 = f[0].length;

    float[][] fp = new float[n2][n1];
    float[][] fm = new float[n2][n1];

    fp[0]    = f[1];
    fm[0]    = f[0];
    fp[n2-1] = f[n2-1];
    fm[n2-1] = f[n2-2];

    for (int i2=1; i2<n2-1; ++i2) {
      fp[i2] = f[i2+1];
      fm[i2] = f[i2-1];
    }
    
    float[][] g = superSample(f);
    float[][] gp = superSample(fp);
    float[][] gm = superSample(fm);
    float[][] pp = _dwk.findShifts(_s1,g,_s1,gp);
    float[][] pm = _dwk.findShifts(_s1,g,_s1,gm);
    float[][] pa = sub(pp,pm);
    pa = mul(pa,0.5f);
    pa = mul(pa,1.0f/(float)_k);
    float[][] qa = subSample(pa);
    return qa;
  }

  public void findSmoothSlopes(Sampling sf, float[][] f, float[][] p) {
    int n2 = f.length;
    int n1 = f[0].length;
    float[][] fm   = new float[n2][n1];
    float[][] temp = new float[n2][n1];

    fm[0]    = f[0];
    fm[n2-1] = f[n2-2];
    for (int i2=1; i2<n2-1; ++i2) {
      fm[i2] = f[i2-1];
    }

    temp = _dwk.findShifts(sf,f,sf,fm);
    temp = mul(temp,-1.0f);
    temp = interpolateSlopes(temp);
    for (int i2=0; i2<n2; ++i2)
      for (int i1=0; i1<n1; ++i1)
        p[i2][i1] = temp[i2][i1];
  }

  public void findSmoothSlopes(Sampling sf, float[][][] f, 
                            float[][][] p2, float[][][] p3) {
    int n3 = f.length;
    int n2 = f[0].length;
    int n1 = f[0][0].length;
    float[][][] f2m  = new float[n3][n2][n1];
    float[][][] f3m  = new float[n3][n2][n1];
    float[][][] temp2  = new float[n3][n2][n1];
    float[][][] temp3  = new float[n3][n2][n1];

    for (int i3=0; i3<n3; ++i3) {
      f2m[i3][0]    = f[i3][0];
      f2m[i3][n2-1] = f[i3][n2-2];
      for (int i2=1; i2<n2-1; ++i2) {
        f2m[i3][i2] = f[i3][i2-1];
      }
    }

    f3m[0]    = f[0];
    f3m[n3-1] = f[n3-2];
    for (int i3=1; i3<n3-1; ++i3) {
      f3m[i3] = f[i3-1];
    }

    temp2 = _dwk.findShifts(sf,f,sf,f2m);
    temp2 = mul(temp2,-1.0f);
    temp2 = interpolateSlopes(temp2);
    temp3 = _dwk.findShifts(sf,f,sf,f3m);
    temp3 = mul(temp3,-1.0f);
    temp3 = interpolateSlopes(temp3);

    for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          p2[i3][i2][i1] = temp2[i3][i2][i1];
          p3[i3][i2][i1] = temp3[i3][i2][i1];
        }
      }
    }
  }

/////////////////////////PRIVATE/////////////////////////

  private static final float pi = FLT_PI;
  private int _k;
  private double _pmax,_h1,_h2,_h3,_r1,_r2,_r3;
  private Sampling _s1,_s2,_s3;
  private DynamicWarpingK _dwk;
  private SincInterpolator _si;

  private static void trace(String s) {
    System.out.println(s);
  }

  private float[][] subSample(float[][] f) {
    int n2 = f.length;
    int n1f = f[0].length;
    int n1g = (n1f-1)/_k+1;
    return copy(n1g,n2,0,0,_k,1,f);
  }

  private float[][] superSample(float[][] f) {
    int n2 = f.length;
    int n1f = f[0].length;
    float d2f = 1.0f;
    float d1f = 1.0f;
    float f2f = 0.0f;
    float f1f = 0.0f;
    Sampling s1f = new Sampling(n1f,d1f,f1f);

    int n1g = (n1f-1)*_k+1;
    float d2g = d2f;
    float d1g = 1.0f/(float)_k;
    float f2g = f2f;
    float f1g = 0.0f;
    float[][] g = new float[n2][n1g];

    int n1h = n1f-1;
    float d1h = d1f;
    float f1h = 0.0f;
    float[] h = new float[n1h];

    for (int i=0; i<_k; ++i, f1h+=d1g) {
      Sampling s1h = new Sampling(n1h,d1h,f1h);
      for (int i2=0; i2<n2; ++i2) {
        _si.interpolate(s1f,f[i2],s1h,h);
        copy(n1h,0,1,h,i,_k,g[i2]);
      }
    }
    for (int i2=0; i2<n2; ++i2)
      for (int i1=0; i1<n1f; ++i1)
        g[i2][i1*_k] = f[i2][i1];
    return g;
  }

  private float[][] interpolateSlopes(float[][] p) {
    int n2 = p.length;
    int n1 = p[0].length;
    float[][] pi = new float[n2][n1];
    float[] x1 = new float[n1];
    float[] x2 = new float[n2];
    for (int i1=0; i1<n1; ++i1)
      x1[i1] = i1;
    for (int i2=0; i2<n2; ++i2)
      x2[i2] = i2-0.5f;

    //BicubicInterpolator2 bc = new BicubicInterpolator2(
    //  BicubicInterpolator2.Method.MONOTONIC,
    //  BicubicInterpolator2.Method.SPLINE,
    //  x1,x2,p);
    BilinearInterpolator2 bc = new BilinearInterpolator2(x1,x2,p);
    pi = bc.interpolate00(_s1,_s2);
    return pi;
  }

  private float[][][] interpolateSlopes(float[][][] p) {
    int n3 = p.length;
    int n2 = p[0].length;
    int n1 = p[0][0].length;
    float[][][] pi = new float[n3][][];
    float[] x1 = new float[n1];
    float[] x2 = new float[n2];
    float[] x3 = new float[n3];
    for (int i1=0; i1<n1; ++i1)
      x1[i1] = i1;
    for (int i2=0; i2<n2; ++i2)
      x2[i2] = i2-0.5f;
    for (int i3=0; i3<n3; ++i3)
      x3[i3] = i3-0.5f;

    //TricubicInterpolator3 tc = new TricubicInterpolator3(
    //  TricubicInterpolator3.Method.MONOTONIC,
    //  TricubicInterpolator3.Method.SPLINE,
    //  TricubicInterpolator3.Method.SPLINE,
    //  x1,x2,x3,p);
    TrilinearInterpolator3 tc = new TrilinearInterpolator3(x1,x2,x3,p);
    pi = tc.interpolate000(_s1,_s2,_s3);
    return pi;
  }
}
