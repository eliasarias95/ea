package util;

import edu.mines.jtk.io.*;
import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.RandomFloat;
import edu.mines.jtk.mosaic.SimplePlot;

import dnp.LocalSlopeFinder;
import dnp.PlaneWaveDestructor;

import slopes.*;

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

  public static float rmsError(float[][] pe, float[][] pk, float d1, float d2,
      boolean print_error) {
    int n1 = pe[0].length;
    int n2 = pe.length;
    pe = sub(pe,pk);
    pe = pow(pe,2);
    float rmserror = sqrt(sum(pe)/(float)(n2*n1));
    if (print_error == true)
      System.out.println("Error = "+rmserror);
    return rmserror;
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

  private static SincInterpolator _s = new SincInterpolator();
  private static SincInterpolator _si = _s.fromErrorAndFrequency(0.003,0.4);
}
