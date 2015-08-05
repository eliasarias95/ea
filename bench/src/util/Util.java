package util;

import edu.mines.jtk.io.*;
import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.RandomFloat;
import edu.mines.jtk.util.SimpleFloat3;

import static edu.mines.jtk.util.ArrayMath.*;

import java.io.*;
import java.nio.*;
import java.util.Random;

/**
 *  
 * @author Elias Arias, Colorado School of Mines, CWP
 * @version 19.1.2015
 */

public class Util {

  /**
   * Converts a 1D array of doubles into an array of floats.
   * @param d array[n1] of doubles to be converted.
   * @return array[n1] of floats.
   */
  public static float[] f(double[] x) {
    int n1 = x.length;
    float[] f = new float[n1];
    for (int i1=0; i1<n1; ++i1)
      f[i1] = (float)x[i1];
    return f;
  }

  /**
   * Converts a 1D array of doubles into an array of ints.
   * @param d array[n1] of doubles to be converted.
   * @return array[n1] of ints.
   */
  public static int[] i(double[] x) {
    int n1 = x.length;
    int[] i = new int[n1];
    for (int i1=0; i1<n1; ++i1)
      i[i1] = (int)x[i1];
    return i;
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

  public static float[][] addNoise(double nrms, float[][] f) {
    int n1 = f[0].length;
    int n2 = f.length;
    Random r = new Random(1);
    //nrms *= max(abs(f));
    float[][] g = mul(2.0f,sub(randfloat(r,n1,n2),0.5f));
    RecursiveGaussianFilter rgf = new RecursiveGaussianFilter(2.0);
    rgf.apply10(g,g); // 1st derivative enhances high-frequencies
    float frms = sqrt(sum(mul(f,f))/n1/n2);
    float grms = sqrt(sum(mul(g,g))/n1/n2);
    g = mul(g,(float)nrms*frms/grms);
    return add(f,g);
  }

  public static float[][][] addNoise(double nrms, float[][][] f) {
    int n1 = f[0][0].length;
    int n2 = f[0].length;
    int n3 = f.length;
    Random r = new Random(1); // 31415
    float[][][] g = mul(2.0f,sub(randfloat(r,n1,n2,n3),0.5f));
    RecursiveGaussianFilter rgf = new RecursiveGaussianFilter(1.0);
    rgf.apply100(g,g); // 1st derivative enhances high-frequencies
    g = mul(g,(float)nrms*rms(f)/rms(g));
    return add(f,g);
  }

  /*
   * A parallel transpose method with a block size of 1.
   */
  public static void transpose23(float[][][] x, float[][][] y) {
    int n3 = x.length;
    int n2 = x[0].length;
    int n1 = x[0][0].length;
    for (int i2=0; i2<n2; ++i2)
      for (int i3=0; i3<n3; ++i3)
        for (int i1=0; i1<n1; ++i1)
          y[i2][i3][i1] = x[i3][i2][i1];
  }

  public static float[][] flip2(float[][] x) {
    int n2 = x.length;
    int n1 = x[0].length;
    float[][] xf = new float[n2][n1];
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        xf[i2][i1] = x[n2-1-i2][i1];
      }
    }
    return xf;
  }

  public static float[][] rotateR(float[][] x) {
    int n2 = x.length;
    int n1 = x[0].length;
    float[][] xf = new float[n1][n2];
    for (int i1=0; i1<n1; ++i1) {
      for (int i2=0; i2<n2; ++i2) {
        xf[i1][i2] = x[i2][n1-1-i1];
      }
    }
    return xf;
  }

  public static float[][] rotateL(float[][] x) {
    int n2 = x.length;
    int n1 = x[0].length;
    float[][] xf = new float[n1][n2];
    for (int i1=0; i1<n1; ++i1) {
      for (int i2=0; i2<n2; ++i2) {
        xf[i1][i2] = x[n2-1-i2][i1];
      }
    }
    return xf;
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

  public static float rms(float[][] f) {
    int n1 = f[0].length;
    int n2 = f.length;
    double sum = 0.0;
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        float fi = f[i2][i1];
        sum += fi*fi;
      }
    }
    return (float)sqrt(sum/n1/n2);
  }

  public static float rms(float[][][] f) {
    int n1 = f[0][0].length;
    int n2 = f[0].length;
    int n3 = f.length;
    double sum = 0.0;
    for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          float fi = f[i3][i2][i1];
          sum += fi*fi;
        }
      }
    }
    return (float)sqrt(sum/n1/n2/n3);
  }

  public static float rmsError(
      float[][] pe, float[][] p, boolean print_error) {
    int n1 = pe[0].length;
    int n2 = pe.length;
    pe = sub(pe,p);
    pe = pow(pe,2);
    float rmserror = sqrt(sum(pe)/(float)(n2*n1));
    if (print_error)
      System.out.println("Error = "+rmserror);
    return rmserror;
  }

  public static float rmsError(
      float[][][] pe, float[][][] p, boolean print_error) {
    int n1 = pe[0][0].length;
    int n2 = pe[0].length;
    int n3 = pe.length;
    pe = sub(pe,p);
    pe = pow(pe,2);
    float rmserror = sqrt(sum(pe)/(float)(n3*n2*n1));
    if (print_error)
      System.out.println("Error = "+rmserror);
    return rmserror;
  }

  public static float errorStatistic(
      float[][][] p2, float[][][] p3, 
      float[] e, float[] p, float[] pe, float[] pbar) {
    int n = 20;
    float[] m  = new float[n];
    float[] am = new float[n];

    pe[0]  = p3[55 ][147][108];
    pe[1]  = p2[52 ][215][105];
    pe[2]  = p2[38 ][197][49 ];
    pe[3]  = p2[49 ][199][108];
    pe[4]  = p2[58 ][194][195];
    pe[5]  = p2[83 ][124][43 ];
    pe[6]  = p3[65 ][0  ][69 ];
    pe[7]  = p2[81 ][56 ][65 ]; // +-
    pe[8]  = p3[77 ][205][158]; // +-
    pe[9]  = p2[139][58 ][82 ]; // +-
    pe[10] = p2[109][246][123];
    pe[11] = p3[124][134][152];
    pe[12] = p3[132][82 ][124];
    pe[13] = p3[91 ][227][108];
    pe[14] = p3[111][235][81 ];
    pe[15] = p3[118][257][24 ];
    pe[16] = p2[81 ][335][75 ];
    pe[17] = p3[93 ][317][166];
    pe[18] = p3[61 ][278][182];
    pe[19] = p2[160][70 ][242];
    //dump(pe);
    
    p[0]  = -0.885f;
    p[1]  = -0.427f;
    p[2]  = -0.369f;
    p[3]  = -0.357f;
    p[4]  = -0.321f;
    p[5]  = -0.149f;
    p[6]  = -0.119f;
    p[7]  = -0.045f;// +-
    p[8]  =  0.011f;// +-
    p[9]  =  0.014f;// +-
    p[10] =  0.263f;
    p[11] =  0.324f;
    p[12] =  0.328f;
    p[13] =  0.409f;
    p[14] =  0.447f;
    p[15] =  0.527f;
    p[16] =  0.560f;
    p[17] =  0.600f;
    p[18] =  0.605f;
    p[19] =  0.847f;
    //dump(p);
    for (int i=0; i<n; ++i)
      e[i] = pe[i] - pbar[i];

    m  = copy(e);
    quickSort(m);
    float med  = ( m[9]+ m[10])/2.0f;
    return med;
  }

  public static float[] reSampleSinc(Sampling sf, float[] f, Sampling sg) {
    float[] g = new float[sg.getCount()];
    _si.interpolate(sf,f,sg,g);
    return g;
  }

  public static void newCurve(int n) {
    RandomFloat rf = new RandomFloat();
    float[] a = new float[n];
    float[] b = new float[n];
    for (int i=0; i<n; ++i) {
      a[i] = rf.normal();
    }
    writeBinary(a,
        "/users/elias.arias/Home/git/ea/bench/src/util/data/dtw_test.dat");
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
    return mul(sgn(x),log(add(abs(x),0.5f)));
  }

  public static float[][] slice12(int k3, float[][][] f) {
    int n1 = f[0][0].length;
    int n2 = f[0].length;
    float[][] s = new float[n2][n1];
    new SimpleFloat3(f).get12(n1,n2,0,0,k3,s);
    return s;
  }

  public static float[][] slice13(int k2, float[][][] f) {
    int n1 = f[0][0].length;
    int n3 = f.length;
    float[][] s = new float[n3][n1];
    new SimpleFloat3(f).get13(n1,n3,0,k2,0,s);
    return s;
  }

  public static float[][] slice23(int k1, float[][][] f) {
    int n2 = f[0].length;
    int n3 = f.length;
    float[][] s = new float[n3][n2];
    new SimpleFloat3(f).get23(n2,n3,k1,0,0,s);
    return s;
  }

  /**
   * Reads a binary file.
   * @param n1 the length of floats in the 1st-dimension
   * @param fileName the name of the file to be read
   * @return array[n1] of floats read from file
   */
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

  /**
   * Reads a binary file.
   * @param n1 the length of floats in the 1st-dimension
   * @param n2 the length of floats in the 2nd-dimension
   * @param fileName the name of the file to be read
   * @return array[n2][n1] of floats read from file
   */
  public static float[][] readImageL(int n1, int n2, String fileName) {
    ByteOrder byteOrder = ByteOrder.LITTLE_ENDIAN;
    try {
      ArrayInputStream ais = new ArrayInputStream(fileName,byteOrder);
      float[][] x = new float[n2][n1];
      ais.readFloats(x);
      ais.close();
      return x;
    } catch (IOException e) {
      throw new RuntimeException(e);
    }
  }

  /**
   * Reads a binary file.
   * @param n1 the length of floats in the 1st-dimension
   * @param n2 the length of floats in the 2nd-dimension
   * @param n3 the length of floats in the 3rd-dimension
   * @param fileName the name of the file to be read
   * @return array[n3][n2][n1] of floats read from file
   */
  public static float[][][] readImage(int n1, int n2, int n3, String fileName) {
    try {
      ArrayInputStream ais = new ArrayInputStream(fileName);
      float[][][] x = new float[n3][n2][n1];
      ais.readFloats(x);
      ais.close();
      return x;
    } catch (IOException e) {
      throw new RuntimeException(e);
    }
  }

  /**
   * Reads a binary file.
   * @param n1 the length of floats in the 1st-dimension
   * @param n2 the length of floats in the 2nd-dimension
   * @param n3 the length of floats in the 3rd-dimension
   * @param fileName the name of the file to be read
   * @return array[n3][n2][n1] of floats read from file
   */
  public static float[][][] readImageL(int n1, int n2, int n3, String fileName) {
    ByteOrder byteOrder = ByteOrder.LITTLE_ENDIAN;
    try {
      ArrayInputStream ais = new ArrayInputStream(fileName,byteOrder);
      float[][][] x = new float[n3][n2][n1];
      ais.readFloats(x);
      ais.close();
      return x;
    } catch (IOException e) {
      throw new RuntimeException(e);
    }
  }

  /**
   * Writes seismic data to binary file.
   * @param x array[n1] of data to write to the binary file
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
  public static void writeBinaryL(float[][] x, String fileName) {
    ByteOrder byteOrder = ByteOrder.LITTLE_ENDIAN;
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
   * Java default byteorder is BIG_ENDIAN.
   * @param x array[n3][n2][n1] of data to write to the binary file
   * @param fileName name of output binary file
   */
  public static void writeBinary(float[][][] x, String fileName) {
    ByteOrder byteOrder = ByteOrder.BIG_ENDIAN;
    try {
      ArrayOutputStream aos = new ArrayOutputStream(fileName,byteOrder);
      aos.writeFloats(x);
      aos.close();
    } catch (IOException e) {
      throw new RuntimeException(e);
    }
  }

  private static void trace(String s) {
    System.out.println(s);
  }

  private static SincInterpolator _s = new SincInterpolator();
  private static SincInterpolator _si = _s.fromErrorAndFrequency(0.003,0.4);
}
