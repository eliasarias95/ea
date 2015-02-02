package slopes;

import util.Plot;
import util.Util;
import edu.mines.jtk.dsp.*;
import static edu.mines.jtk.util.ArrayMath.*;

public class DynamicWarpingSlopes {

  public DynamicWarpingSlopes(int pmax, float usmooth1, float usmooth2,
      float strainMax1, float strainMax2, int esmooth) {
    _dw = new DynamicWarping(-pmax,pmax);
    _dw.setShiftSmoothing(usmooth1,usmooth2);
    _dw.setStrainMax(strainMax1,strainMax2);
    _dw.setErrorSmoothing(esmooth);
    _si = new SincInterpolator();
  }

  public float[][] findSlopes(int k, float[][] f) {
    int n2 = f.length;
    int n1 = f[0].length;
    float[][] p = new float[n2][n1];

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
    
    float[][] g = superSample(k,f);
    float[][] gp = superSample(k,fp);
    float[][] gm = superSample(k,fm);
    float[][] pp = _dw.findShifts(g,gp);
    float[][] pm = _dw.findShifts(g,gm);
    float[][] pa = sub(pp,pm);
    pa = mul(pa,0.5f);
    pa = mul(pa,1.0f/(float)k);
    float[][] qa = subSample(k,pa);
    return qa;
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

  private float[][] subSample(int k, float[][] f) {
    int n2 = f.length;
    int n1f = f[0].length;
    int n1g = (n1f-1)/k+1;
    return copy(n1g,n2,0,0,k,1,f);
  }

  private float[][] superSample(int k, float[][] f) {
    int n2 = f.length;
    int n1f = f[0].length;
    float d2f = 1.0f;
    float d1f = 1.0f;
    float f2f = 0.0f;
    float f1f = 0.0f;
    Sampling s1f = new Sampling(n1f,d1f,f1f);

    int n1g = (n1f-1)*k+1;
    float d2g = d2f;
    float d1g = 1.0f/(float)k;
    float f2g = f2f;
    float f1g = 0.0f;
    float[][] g = new float[n2][n1g];

    int n1h = n1f-1;
    float d1h = d1f;
    float f1h = 0.0f;
    float[] h = new float[n1h];

    for (int i=0; i<k; ++i, f1h+=d1g) {
      Sampling s1h = new Sampling(n1h,d1h,f1h);
      for (int i2=0; i2<n2; ++i2) {
        _si.interpolate(s1f,f[i2],s1h,h);
        copy(n1h,0,1,h,i,k,g[i2]);
      }
    }
    for (int i2=0; i2<n2; ++i2)
      for (int i1=0; i1<n1f; ++i1)
        g[i2][i1*k] = f[i2][i1];
    return g;
  }

  private static void goTestSuperAndSub() {
    DynamicWarpingSlopes dws = new DynamicWarpingSlopes(5,18.0f,
                                                        1.0f,1.0f,1.0f,2);
    int n = 10;
    float[][] x = new float[n][n];
    for (int i2=0; i2<n; ++i2) {
      for (int i1=0; i1<n; ++i1) {
        x[i2][i1] = i1;
        x[i2][0] = 3;
        System.out.println("x["+i2+"]["+i1+"]= "+x[i2][i1]);
      }
    }

    int k = 5;
    float[][] y = dws.superSample(k,x);
    float[][] z = dws.subSample(k,y);
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
    int pmax = 6;
    DynamicWarpingSlopes dws = new DynamicWarpingSlopes(pmax*k,k*18.0f,
                                                        2.0f,0.2f,0.2f,2);
    int m1 = 501;
    int m2 = 501;
    float freq = 0.1f;
    float p2 = -5.7f;
    float[][][] fandp = makeConstantSlope(freq,p2,m1,m2);
    float[][] f = fandp[0];
    float[][] g = dws.superSample(k,f); //super sampled image
    float[][] p = fandp[1];
    float[][] q = dws.superSample(k,p); //super sampled slopes 

    float[][] pe = dws.findSlopes(k,f);
    System.out.println("k= "+k+" slope= "+pe[230][230]);
    Sampling s1 = new Sampling(m1);
    Sampling sg = new Sampling(m1*k,1.0f/k,0.0f);
    Sampling s2 = new Sampling(m2);
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
  private DynamicWarping _dw;
  private SincInterpolator _si;
}
