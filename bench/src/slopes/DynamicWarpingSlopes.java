package slopes;

import util.Plot;
import util.Util;
import edu.mines.jtk.dsp.*;
import edu.mines.jtk.interp.*;
import static edu.mines.jtk.util.ArrayMath.*;

public class DynamicWarpingSlopes {

  public DynamicWarpingSlopes(double pmax, Sampling s1, Sampling s2) {
    this(1,pmax,1.0,1.0,-1.0,1.0,-1.0,1.0,s1,s2);
  }

  public DynamicWarpingSlopes(int k, double pmax, double h1, double h2,
      double r1min, double r1max, double r2min, double r2max, 
      Sampling s1, Sampling s2) {
    _k = k;
    _pmax = pmax;
    _h1 = h1;
    _h2 = h2;
    _r1min = r1min;
    _r1max = r1max;
    _r2min = r2min;
    _r2max = r2max;
    _s1 = s1;
    _s2 = s2;
    _dwk = new DynamicWarpingK(k,-pmax,pmax,s1,s2);
    _dwk.setSmoothness(h1,h2);
    _dwk.setStrainLimits(r1min,r1max,r2min,r2max);
    _si = new SincInterpolator();
  }

  public void setK(int k) {
    _k = k;
    _dwk = new DynamicWarpingK(k,-_pmax,_pmax,_s1,_s2);
    _dwk.setSmoothness(_h1,_h2);
    _dwk.setStrainLimits(_r1min,_r1max,_r2min,_r2max);
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
    
    trace("n1*k= "+n1*_k);
    Sampling sg = new Sampling(n1*_k);
    float[][] g = superSample(f);
    float[][] gp = superSample(fp);
    float[][] gm = superSample(fm);
    float[][] pp = _dwk.findShifts(sg,g,sg,gp);
    float[][] pm = _dwk.findShifts(sg,g,sg,gm);
    float[][] pa = sub(pp,pm);
    pa = mul(pa,0.5f);
    pa = mul(pa,1.0f/(float)_k);
    float[][] qa = subSample(pa);
    return qa;
  }

  public float[][] findSmoothSlopes(float[][] f) {
    int n2 = f.length;
    int n1 = f[0].length;
    float[][] pm = new float[n2][n1];
    float[][] fm = new float[n2][n1];

    fm[0]    = f[0];
    fm[n2-1] = f[n2-2];
    for (int i2=1; i2<n2-1; ++i2) {
      fm[i2] = f[i2-1];
    }

    pm = _dwk.findShifts(_s1,f,_s1,fm);
    pm = mul(pm,-1.0f);
    pm = interpolateSlopes(pm);
    return pm;
  }

  public static float[][][] makeConstantSlope(float freq, float p2, 
                                              int n1, int n2) {
    float[][][] fandp = new float[2][n2][n1];
    float[][] f = new float[n2][n1];
    float[][] p = new float[n2][n1];
    for (int i2=0; i2<n2; ++i2){
      for (int i1=0; i1<n1; ++i1){
        f[i2][i1] = cos(2.0f*pi*(p2*freq*i2-freq*i1));
        p[i2][i1] = p2;
      }
    }
    fandp[0] = copy(f);
    fandp[1] = copy(p);
    return fandp;
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
      x1[i1] = i1-0.5f;
    for (int i2=0; i2<n2; ++i2)
      x2[i2] = i2-0.5f;

    BicubicInterpolator2 bc = new BicubicInterpolator2(
      BicubicInterpolator2.Method.MONOTONIC,
      BicubicInterpolator2.Method.SPLINE,
      x1,x2,p);
    pi = bc.interpolate00(_s1,_s2);
    return pi;
  }

  private static void trace(String s) {
    System.out.println(s);
  }

  private static void goTestSuperAndSub() {
    int k = 5;
    int n = 10;
    Sampling s = new Sampling(n);
    DynamicWarpingSlopes dws = new DynamicWarpingSlopes(1,s,s);
    dws.setK(k);
    float[][] x = new float[n][n];
    for (int i2=0; i2<n; ++i2) {
      for (int i1=0; i1<n; ++i1) {
        x[i2][i1] = i1;
        x[i2][0] = 3;
        System.out.println("x["+i2+"]["+i1+"]= "+x[i2][i1]);
      }
    }

    float[][] y = dws.superSample(x);
    float[][] z = dws.subSample(y);
    int n1 = y[0].length;
    int n2 = y.length;
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        System.out.println("y["+i2+"]["+i1+"]= "+y[i2][i1]);
      }
    }

    n1 = z[0].length;
    n2 = z.length;
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        System.out.println("z["+i2+"]["+i1+"]= "+z[i2][i1]);
      }
    }
  }

  private static void goTestDWS() {
    int k = 10;
    int pmax = 4;
    int m1 = 501;
    int m2 = 501;
    Sampling s1 = new Sampling(m1);
    Sampling s2 = new Sampling(m2);
    DynamicWarpingSlopes dws = new DynamicWarpingSlopes(pmax*k,s1,s2);
    dws.setK(k);
    float freq = 0.1f;
    float p2 = -3.7f;
    float[][][] fandp = makeConstantSlope(freq,p2,m1,m2);
    float[][] f = fandp[0];
    float[][] g = dws.superSample(f); //super sampled image
    float[][] p = fandp[1];
    float[][] q = dws.superSample(p); //super sampled slopes 

    float[][] pe = dws.findSlopes(f);
    System.out.println("k= "+k+" slope= "+pe[230][230]);
    Sampling sg = new Sampling(m1*k,1.0f/k,0.0f);
    System.out.println("max slope= "+max(p));
    // clip, interp, title, paint, colorbar, color
    String hl = "Traces"; //horizontal label
    String vl = "Samples"; //vertical label
    String cbl = ""; //colorbar label
    Plot.plot(sg,s2,g,"Super Synthetic",hl,vl,cbl,
        0.75f,0.9f,-4,4,
        false,false,true,false,true,false);
    Plot.plot(sg,s2,f,"Synthetic",hl,vl,cbl,
        0.75f,0.9f,-4,4,
        false,false,true,false,true,false);

    cbl = "slope (samples/trace)"; //colorbar label
    Plot.plot(sg,s2,q,"Known slopes",hl,vl,cbl,
        0.75f,0.9f,-4,4,
        false,false,true,false,true,true);

    // interp, title, paint, colorbar, color
    Plot.plot(s1,s2,pe,"DW",hl,vl,cbl,
        0.75f,0.9f,-4,4,
        false,false,true,false,true,true);
    float error = Util.rmsError(pe,p,1.0f/(float)k,1.0f,true);
  }

  public static void main(String[] args) {
    goTestDWS();
    //goTestSuperAndSub();
  }

  private static final float pi = FLT_PI;
  private int _k;
  private double _pmax,_h1,_h2,_r1min,_r1max,_r2min,_r2max;
  private Sampling _s1,_s2;
  private DynamicWarpingK _dwk;
  private SincInterpolator _si;
}
