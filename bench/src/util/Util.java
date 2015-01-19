package util;

import edu.mines.jtk.io.*;
import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.Check;
import edu.mines.jtk.util.RandomFloat;
import edu.mines.jtk.mosaic.SimplePlot;

import dnp.LocalSlopeFinder;
import dnp.PlaneWaveDestructor;

import slopes.*;

import static edu.mines.jtk.util.ArrayMath.*;

import java.io.*;
import java.nio.*;
import java.util.Random;

public class Util {

  public static float rmsError(float[][] pe, float[][] pk, float d1, float d2,
      boolean print_error) {
    int n1 = pe[0].length;
    int n2 = pe.length;
    pe = mul(pe,d1/d2);
    pe = sub(pe,pk);
    pe = pow(pe,2);
    float rmserror = sqrt(sum(pe)/(n2*n1));
    if (print_error == true)
      System.out.println("Error = "+rmserror);
    return rmserror;
  }

  public static float[][] DWSlopesAvg(DynamicWarping dw, 
      float[][] f) {
    int n1 = f[0].length;
    int n2 = f.length;

    float[][] pp = new float[n2][n1];
    float[][] pm = new float[n2][n1];
    float[][] pa = new float[n2][n1];
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

    pp = dw.findShifts(f,fp);
    pm = dw.findShifts(f,fm);
    pa = sub(pp,pm);
    pa = mul(pa,0.5f);
    return pa;
  }

  public static float[][] DWSlopes1D(DynamicWarping dw, 
      float[][] f) {
    int n1 = f[0].length;
    int n2 = f.length;
    float[][] p = new float[n2][n1];
    for (int i2=1; i2<n2-1; ++i2)
      p[i2] = dw.findShifts(f[i2-1],f[i2]);
    p[0] = p[1];
    return p;
  }

  /**
   * Tests for the parameters that produce the smallest RMS error.
   * The m variable requires a number 1,2,3 or 4.
   * For optimal structure tensor parameters, method = 1
   * For optimal Madagascar PWD parameters,   method = 2
   * For optimal dynamic warping parameters,  method = 3
   * For optimal Dave's PWD parameters,       method = 4
   */
  public static void testOptimalParameters(Sampling s1_param, Sampling s2_param, 
      Sampling s1, Sampling s2, float[] param1, float[] param2, 
      int m, String method) {
    Check.argument(m==1 || m==2 || m==3 || m== 4,"valid method");
    
    int n1_param = s1_param.getCount();
    int n2_param = s2_param.getCount();
    int n1 = s1.getCount();
    int n2 = s2.getCount();

    float d1 = (float)s1.getDelta();
    float d2 = (float)s2.getDelta();
    
    float noise = 0.5f;
    float[][][] fandpk = FakeData.seismicAndSlopes2d2014B(noise,false);
    float[][] f = fandpk[0];
    float[][] pk = fandpk[1];
    float[][] pe = new float[n2][n1];
    int pmax = 5;

    if (m== 1) {
      float[][] rmserror = new float[n2_param][n1_param]; // RMS error
      LocalSlopeFinder lsf;
      for(int i2=0; i2<n2_param; ++i2) {
        for(int i1=0; i1<n1_param; ++i1) {
          lsf = new LocalSlopeFinder(param1[i1],param2[i2],pmax);
          lsf.findSlopes(f,pe);
          rmserror[i2][i1] = rmsError(pe,pk,d1,d2,false);
        }
      }
      writeBinary(rmserror,
        "/Users/earias/Home/git/ea/bench/src/slopes/data/"+method+"_errors.dat");
    }

    else if (m== 2) {
      float[][] rmserror = new float[n2_param][n1_param]; // RMS error
      Sfdip sd = new Sfdip(-pmax,pmax);
      sd.setOrder(2);
      for(int i2=0; i2<n2_param; ++i2) {
        for(int i1=0; i1<n1_param; ++i1) {
          sd.setRect((int)param1[i1],(int)param2[i2]);
          sd.findSlopes(s1,s2,f,pe);
          rmserror[i2][i1] = rmsError(pe,pk,d1,d2,false);
        }
      }
      writeBinary(rmserror,
        "/Users/earias/Home/git/ea/bench/src/slopes/data/"+method+"_errors.dat");
    }

    else if (m== 3) {
      float[][] rmserror = new float[n2_param][n1_param]; // RMS error
      DynamicWarping dw = new DynamicWarping(-pmax,pmax);
      dw.setStrainMax(0.5,0.5);
      dw.setErrorSmoothing(2);
      for(int i2=0; i2<n2_param; ++i2) {
        for(int i1=0; i1<n1_param; ++i1) {
          dw.setShiftSmoothing(param1[i1],param2[i2]);
          pe = DWSlopesAvg(dw,f);
          rmserror[i2][i1] = rmsError(pe,pk,d1,d2,false);
        }
      }
      writeBinary(rmserror,
        "/Users/earias/Home/git/ea/bench/src/slopes/data/"+method+"_errors.dat");
    }

    else if (m== 4) {
      float[][] rmserror = new float[n2_param][n1_param]; // RMS error
      PlaneWaveDestructor pwd = new PlaneWaveDestructor(-pmax,pmax);
      for(int i2=0; i2<n2_param; ++i2) {
        for(int i1=0; i1<n1_param; ++i1) {
          pwd.setSmoothness(param1[i1],param2[i2]);
          pe = pwd.findSlopes(f);
          rmserror[i2][i1] = rmsError(pe,pk,d1,d2,false);
        }
      }
      writeBinary(rmserror,
        "/Users/earias/Home/git/ea/bench/src/slopes/data/"+method+"_errors.dat");
    }
  }

  public static void plotOptimalParameters(Sampling s1_param, Sampling s2_param,
      Sampling s1, Sampling s2, float[] param1, float[] param2, String method) {
    int n1_param = s1_param.getCount();
    int n2_param = s2_param.getCount();

    int[] error_index = new int[2];    
    float[][] rmserror = readImage(n1_param,n2_param,
      "/Users/earias/Home/git/ea/bench/src/slopes/data/"+method+"_errors.dat");
    float min_error = min(rmserror,error_index);
    System.out.println("Parameter1= "+param1[error_index[0]]+
                      " Parameter2= "+param2[error_index[1]]+
                      " Minimum Error Value= "+min_error);    
    float fw = 0.75f; //fraction width for slide
    float fh = 0.9f; //fraction height for slide
    // sigma, interp, title, paint, colorbar, color
    Plot.plot(s1_param,s2_param,rmserror,"RMS Error "+method,fw,fh,true,true,false,
        true,true,true);
  }

  public static float[] addNoise(double nrms, float[] f) {
    int n1 = f.length;
    Random r = new Random(1);
    float[] g = mul(2.0f,sub(randfloat(r,n1),0.5f));
    RecursiveGaussianFilter rgf = new RecursiveGaussianFilter(2.0);
    rgf.apply1(g,g); // 1st derivative enhances high-frequencies
    g = mul(g,(float)nrms*rms(f)/rms(g));
    return add(f,g);
  }

  public static float rms(float[] f) {
    int n1 = f.length;
    double sum = 0.0;
    for (int i1=0; i1<n1; ++i1) {
      float fi = f[i1];
      sum += fi*fi;
    }
    return (float)sqrt(sum/n1);
  }

  public static void newCurve(int n) {
    RandomFloat rf = new RandomFloat();
    float[] a = new float[n];
    float[] b = new float[n];
    for (int i=0; i<n; ++i) {
      a[i] = rf.normal();
    }
    writeBinary(a,
        "/Users/earias/Home/git/ea/bench/src/slopes/data/dtw_test.dat");
  }

  /**
   * Ungains the data input.
   * @param x array[n2][n1] of data to be ungained
   * @return array[n2][n1] of ungained data
   */
  public static float[][] slog(float[][] x) {
    return mul(sgn(x),sub(exp(abs(x)),1.0f));
  }

  /**
   * Gains the data input.
   * @param x array[n2][n1] of data to be gained
   * @return array[n2][n1] of gained data
   */
  public static float[][] sexp(float[][] x) {
    return mul(sgn(x),log(add(abs(x),1.0f)));
  }

  /**
   * Reads a binary file.
   * @param n1 the length of floats in the 1st-dimension
   * @param n2 the length of floats in the 2nd-dimension
   * @param fileName the name of the file to be read
   * @return array[n2][n1] of floats read from file
   */
  public static float[][] readImage(int n1, int n2, String fileName) {
    try {
      ArrayInputStream ais = new ArrayInputStream(fileName);
      float[][] x = new float[n2][n1];
      ais.readFloats(x);
      ais.close();
      return x;
    } catch (IOException e) {
      throw new RuntimeException(e);
    }
  }

  public static float[] readImage(int n1, String fileName) {
    try {
      ArrayInputStream ais = new ArrayInputStream(fileName);
      float[] x = new float[n1];
      ais.readFloats(x);
      ais.close();
      return x;
    } catch (IOException e) {
      throw new RuntimeException(e);
    }
  }

  /**
   * Writes seismic data to binary file.
   * @param x array[n2][n1] of data to write to the binary file
   * @param fileName name of output binary file
   */
  public static void writeBinary(float[][] x, String fileName) {
    ByteOrder byteOrder = ByteOrder.BIG_ENDIAN;
    try {
      ArrayOutputStream aos = new ArrayOutputStream(fileName,byteOrder);
      aos.writeFloats(x);
      aos.close();
    } catch (IOException e) {
      throw new RuntimeException(e);
    }
  }

  /**
   * Writes seismic data to binary file.
   * @param x array[n2][n1] of data to write to the binary file
   * @param fileName name of output binary file
   */
  public static void writeBinary(float[] x, String fileName) {
    ByteOrder byteOrder = ByteOrder.BIG_ENDIAN;
    try {
      ArrayOutputStream aos = new ArrayOutputStream(fileName,byteOrder);
      aos.writeFloats(x);
      aos.close();
    } catch (IOException e) {
      throw new RuntimeException(e);
    }
  }
}
