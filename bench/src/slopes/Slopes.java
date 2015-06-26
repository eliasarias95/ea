package slopes;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.sgl.*;
import edu.mines.jtk.util.Check;
import edu.mines.jtk.util.Stopwatch;


import static edu.mines.jtk.util.ArrayMath.*;

import util.*;
import dnp.PlaneWaveDestructor;

import java.awt.Color;
import java.util.Random;
import javax.swing.*;

/**
 *  
 * @author Elias Arias, Colorado School of Mines, CWP
 * @version 25.3.2015
 */

public class Slopes{
  
  public Slopes(float noise, float pmax, Sampling s1, Sampling s2) {
    this(noise,pmax,s1,s2,null);
  }

  public Slopes(
      float noise, float pmax, Sampling s1, Sampling s2, Sampling s3){
    if (noise==0.0f)      num = 1;
    else if (noise==0.5f) num = 2;
    else if (noise==1.0f) num = 3;
    else                  num = 4;
    _noise = noise;
    _s1 = s1;
    _s2 = s2;
    _s3 = s3;
    _pmax = pmax;
  }

  /**
   * Creates a 2D synthetic seismic image with complex geologic structures and
   * known slope values.
   * @param noise RMS noise to RMS signal ratio for the generated image.
   * @param f array[501][501] for the seismic image.
   * @param p array[501][501] for the known slope values.
   * @param r array[501][501] for the reflectivity.
   */
  public static void makeSyntheticComplex(
      float noise, float[][] f, float[][] p, float[][] r) {
    int n2 = f.length;
    int n1 = f[0].length;
    float[][][] fandp = FakeData.seismicAndSlopes2d2014A(noise,false);
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        f[i2][i1] = fandp[0][i2][i1];
        p[i2][i1] = fandp[1][i2][i1];
        r[i2][i1] = fandp[2][i2][i1];
      }
    }
  }

  /**
   * Creates a 2D synthetic seismic image with constant, known slope values.
   * @param freq frequency of the 
   * @param p2 constant slope value for the image.
   * @param f array[n2][n1] for the seismic image.
   * @param p array[n2][n1] for the know slope values.
   */
  public static void makeSyntheticConstant(
      float freq, float p2, float[][] f, float[][] p) {
    int n2 = f.length;
    int n1 = f[0].length;
    for (int i2=0; i2<n2; ++i2){
      for (int i1=0; i1<n1; ++i1){
        f[i2][i1] = cos(2.0f*pi*(p2*freq*i2-freq*i1));
        p[i2][i1] = p2;
      }
    }
  }

  /**
   * Creates a 2D synthetic seismic image with constant, known slope values.
   * @param freq frequency of the 
   * @param p2 constant slope value for the image.
   * @param f array[n2][n1] for the seismic image.
   * @param p array[n2][n1] for the know slope values.
   */
  public static void makeSyntheticConstant(float freq, float pc2, float pc3, 
                            float[][][] f, float[][][] p2, float[][][] p3) {
    int n3 = f.length;
    int n2 = f[0].length;
    int n1 = f[0][0].length;
    for (int i3=0; i3<n3; ++i3){
      for (int i2=0; i2<n2; ++i2){
        for (int i1=0; i1<n1; ++i1){
          f[i3][i2][i1] = cos(2.0f*pi*(pc3*freq*i3+pc2*freq*i2-freq*i1));
          p2[i3][i2][i1] = pc2;
          p3[i3][i2][i1] = pc3;
        }
      }
    }
  }

  /**
   * Creates a 3D synthetic seismic image with complex geologic structures and
   * known slope values.
   * @param noise RMS noise to RMS signal ratio for the generated image.
   * @param f array[103][102][101] for the seismic image.
   * @param p array[103][102][101] for the know slope values.
   */
  public static void makeSyntheticComplex(
      float noise, float[][][] f, float[][][] p2, float[][][] p3) {
    float[][][][] fandp = FakeData.seismicAndSlopes3d2014A(noise,false);
    int n3 = f.length;
    int n2 = f[0].length;
    int n1 = f[0][0].length;
    for (int i3=0; i3<n3; ++i3){
      for (int i2=0; i2<n2; ++i2){
        for (int i1=0; i1<n1; ++i1){
          f[i3][i2][i1]  = fandp[0][i3][i2][i1];
          p2[i3][i2][i1] = fandp[1][i3][i2][i1];
          p3[i3][i2][i1] = fandp[2][i3][i2][i1];
        }
      }
    }
  }

  /**
   * Takes a real seismic image from the GOM and places it in the array f.
   * @param f array[920][301] that holds the data for the GOM seismic image.
   */
  public static void makeRealGOM(float[][] f) {
    int n1 = f[0].length;
    int n2 = f.length;
    float[][] temp = (Util.readImage(n1,n2,PATH+"data/gom.dat"));
    for (int i2=0; i2<n2; ++i2){
      for (int i1=0; i1<n1; ++i1){
        f[i2][i1] = temp[i2][i1];
      }
    }
  }

  /**
   * Takes a real seismic image from the GOM and places it in the array f.
   * @param f array[920][301] that holds the data for the GOM seismic image.
   */
  public static void makeRealTp(float[][] f) {
    int n1 = f[0].length;
    int n2 = f.length;
    float[][] temp = (Util.readImage(n1,n2,PATH+"data/tp/tp73.dat"));
    //temp = Util.addNoise(0.0,temp);
    for (int i2=0; i2<n2; ++i2){
      for (int i1=0; i1<n1; ++i1){
        f[i2][i1] = temp[i2][i1];
      }
    }
  }

  /**
   * Takes a real seismic image from the GOM and places it in the array f.
   * @param f array[920][301] that holds the data for the GOM seismic image.
   */
  public static void makeRealTp(float[][][] f) {
    int n1 = f[0][0].length;
    int n2 = f[0].length;
    int n3 = f.length;
    //float[][] g = new float[n3][n1];
    float[][][] temp = Util.readImage(n1,n2,n3,PATH+
          //"data/tp/tpsz_subz_51_4_1400.dat");
          "data/tp/tpsz_subz_401_4_400.dat");
    //temp = Util.addNoise(0.0001,temp);
    for (int i3=0; i3<n3; ++i3){
      for (int i2=0; i2<n2; ++i2){
        for (int i1=0; i1<n1; ++i1){
          f[i3][i2][i1] = temp[i3][i2][i1];
        }
      }
    }

    /*
    for (int i3=0; i3<n3; ++i3){
      for (int i1=0; i1<n1; ++i1){
        g[i3][i1] = f[i3][317][i1];
      }
    }
    Util.writeBinary(g,PATH+"data/tp/tp317.dat");
    */
  }

/****************************ESTIMATING SLOPES****************************/

  /**
   * Structure tensor: plots the estimated slopes, RMS error, and time.
   */
  public void estimateLSF(float[][] f, float[][] p, String title) {
    int n1 = _s1.getCount();
    int n2 = _s2.getCount();
    LocalSlopeFinder lsf = new LocalSlopeFinder(23.0f,1.0f,_pmax);
    float[][] pe = new float[n2][n1];
    ZeroMask zm = new ZeroMask(f);
    _sw.restart();
    //lsf.findSlopesE(f,pe,null);
    lsf.findSlopes(f,pe);
    _sw.stop();
    zm.apply(0.0f,pe);
    trace("Structure tensor time = "+_sw.time());    
    Util.writeBinary(pe,PATH+"data/"+title+"_p.dat");
    trace("Structure tensor:");
    float error;
    if (p!=null) {
      error = Util.rmsError(pe,p,T);
      trace("mean slope= "+sum(p)/(n1*n2));
    }
  }

  /**
   * Structure tensor: plots the estimated slopes.
   */
  public void estimateLSF(
      float[][][] f, float[][][] p2, float[][][] p3, String title) {
    int n1 = _s1.getCount();
    int n2 = _s2.getCount();
    int n3 = _s3.getCount();
    LocalSlopeFinder lsf = new LocalSlopeFinder(8.0f,2.0f,2.0f,_pmax);
    float[][][] p2e = new float[n3][n2][n1];
    float[][][] p3e = new float[n3][n2][n1];
    ZeroMask zm = new ZeroMask(f);
    _sw.restart();
    //lsf.findSlopesE(f,p2e,p3e,null);
    lsf.findSlopes(f,p2e,p3e,null);
    _sw.stop();
    zm.apply(0.0f,p2e);
    zm.apply(0.0f,p3e);
    trace("Structure tensor time = "+_sw.time());    
    Util.writeBinary(p2e,PATH+"data/"+title+"_p2.dat");
    Util.writeBinary(p3e,PATH+"data/"+title+"_p3.dat");
    trace("Structure tensor:");
    float error2,error3;
    if (p2!=null && p3!=null) {
      trace("p2:");
      error2 = Util.rmsError(p2e,p2,T);
      trace("p3:");
      error3 = Util.rmsError(p3e,p3,T);
      trace("p2:max= "+max(p2)+" min= "+min(p2)+" mean= "+sum(p2)/(n1*n2*n3));
      trace("p3:max= "+max(p3)+" min= "+min(p3)+" mean= "+sum(p3)/(n1*n2*n3));
    }
  }

  /**
   * Structure tensor: plots the estimated slopes.
   */
  public void estimateTransLSF(float[][][] f, String title) {
    int n1 = _s1.getCount();
    int n2 = _s2.getCount();
    int n3 = _s3.getCount();
    LocalSlopeFinder lsf = new LocalSlopeFinder(8.0f,2.0f,2.0f,_pmax);
    float[][][] p2 = new float[n2][n3][n1];
    float[][][] p3 = new float[n2][n3][n1];
    ZeroMask zm = new ZeroMask(f);
    _sw.restart();
    //lsf.findSlopesE(f,p2,p3,null);
    lsf.findSlopes(f,p2,p3,null);
    _sw.stop();
    zm.apply(0.0f,p2);
    zm.apply(0.0f,p3);
    trace("Structure tensor time = "+_sw.time());    
    Util.writeBinary(p2,PATH+"data/"+title+"_p2.dat");
    Util.writeBinary(p3,PATH+"data/"+title+"_p3.dat");
  }

  /**
   * PWD Madagascar: plots the estimated slopes, RMS error, and time.
   */
  public void estimatePWDM(float[][] f, float[][] p, String title) {
    int n1 = _s1.getCount();
    int n2 = _s2.getCount();
    Sfdip sd = new Sfdip(-_pmax,_pmax);
    sd.setRect(40,6);
    sd.setOrder(4);
    sd.setNiter(_niter);
    sd.setNj(1);
    sd.setBoth("n");
    float[][] pe = new float[n2][n1]; //pwd w/initial p=0 (Madagascar)
    _sw.restart();
    sd.findSlopes(_s1,_s2,f,pe);
    _sw.stop();
    trace("Madagascar PWD time = "+_sw.time());    
    Util.writeBinary(pe,PATH+"data/"+title+"_p.dat");
    trace("Madagascar PWD:");
    float error;
    if (p!=null) error = Util.rmsError(pe,p,T);
  }

  /**
   * PWD Madagascar: plots the estimated slopes.
   */
  public void estimatePWDM(
      float[][][] f, float[][][] p2, float[][][] p3, String title) {
    int n1 = _s1.getCount();
    int n2 = _s2.getCount();
    int n3 = _s3.getCount();
    Sfdip sd = new Sfdip(-_pmax,_pmax);
    sd.setRect(10,9,9);
    sd.setOrder(4);
    sd.setNiter(_niter);
    sd.setN4(2);
    sd.setBoth("n");
    float[][][] p2e = new float[n3][n2][n1];
    float[][][] p3e = new float[n3][n2][n1];
    _sw.restart();
    sd.findSlopes(_s1,_s2,_s3,f,p2e,p3e);
    //sd.ffindSlopes(_s1,_s2,_s3,f,p2e,p3e); //fast 3D find slopes
    _sw.stop();
    trace("Madagascar PWD time = "+_sw.time());    
    Util.writeBinary(p2e,PATH+"data/"+title+"_p2.dat");
    Util.writeBinary(p3e,PATH+"data/"+title+"_p3.dat");
    trace("Madagascar PWD:");
    float error2,error3;
    if (p2!=null && p3!=null) {
      trace("p2:");
      error2 = Util.rmsError(p2e,p2,T);
      trace("p3:");
      error3 = Util.rmsError(p3e,p3,T);
    }
  }

  /**
   * PWD Madagascar: plots the estimated slopes.
   */
  public void estimateTransPWDM(float[][][] f, String title) {
    int n1 = _s1.getCount();
    int n2 = _s2.getCount();
    int n3 = _s3.getCount();
    Sfdip sd = new Sfdip(-_pmax,_pmax);
    sd.setRect(30,9,9);
    sd.setOrder(2);
    sd.setNiter(_niter);
    sd.setN4(2);
    sd.setBoth("n");
    float[][][] p2 = new float[n2][n3][n1];
    float[][][] p3 = new float[n2][n3][n1];
    _sw.restart();
    sd.findSlopes(_s1,_s3,_s2,f,p2,p3);
    //sd.ffindSlopes(_s1,_s2,_s3,f,p2,p3); //fast 3D find slopes
    _sw.stop();
    trace("Madagascar PWD time = "+_sw.time());    
    Util.writeBinary(p2,PATH+"data/"+title+"_p2.dat");
    Util.writeBinary(p3,PATH+"data/"+title+"_p3.dat");
  }

  /**
   * PWD Dave: plots the estimated slopes, RMS error, and time.
   */
  public void estimatePWDD(float[][] f, float[][] p, String title) {
    boolean one = true;
    int n1 = _s1.getCount();
    int n2 = _s2.getCount();
    PlaneWaveDestructor pwd = new PlaneWaveDestructor(-_pmax,_pmax);
    pwd.setSmoothness(6,1);
    pwd.setLateralBias(0.0);
    pwd.setOuterIterations(_niter); // default is 5
    pwd.setInnerIterations(20); // default is 20
    float[][] pe = new float[n2][n1]; //pwd w/initial p=0 (Dave)
    _sw.restart();
    pe = pwd.findSlopes(f);
    pwd.updateSlopes(f,pe);
    _sw.stop();
    trace("Dave's PWD time = "+_sw.time());    
    Util.writeBinary(pe,PATH+"data/"+title+"_p.dat");
    trace("Dave's PWD:");
    float error;
    if (p!=null) error = Util.rmsError(pe,p,T);
  }

  /**
   * Dynamic warping: plots the estimated slopes, RMS error, and time.
   */
  public void estimateDW(int k, float[][] f, float[][] p, String title) {
    boolean one = true;
    int n1 = _s1.getCount();
    int n2 = _s2.getCount();
    Sampling ss1 = new Sampling((n1-1)*k+1);
    Sampling ss2 = new Sampling(n2);
    DynamicWarpingSlopes dws = new DynamicWarpingSlopes((int)_pmax*k,ss1,ss2);
    dws.setK(k);
    _sw.restart();
    float[][] pe = dws.findSlopes(f);
    _sw.stop();
    trace("Dynamic warping time = "+_sw.time());    
    Util.writeBinary(pe,PATH+"data/"+title+"_p.dat");
    trace("Dynamic warping:");
    float error;
    if (p!=null) error = Util.rmsError(pe,p,T);
  }

  /**
   * Smooth dynamic warping: plots the estimated slopes, RMS error, and time.
   */
  public void estimateSDW(int k, float[][] f, float[][] p, String title) {
    int n1 = _s1.getCount();
    int n2 = _s2.getCount();
    double r1 = 0.1;
    double r2 = 0.3;
    double h1 = 72.0;
    double h2 =  6.0;
    Sampling ss1 = new Sampling(n1);
    Sampling ss2 = new Sampling(n2);
    DynamicWarpingSlopes dws = new DynamicWarpingSlopes(k,_pmax,h1,h2,r1,r2,
                                                        ss1,ss2);
    dws.setErrorSmoothing(2);
    _sw.restart();
    float[][] pe = new float[n2][n1];
    ZeroMask zm = new ZeroMask(f);
    dws.findSmoothSlopes(ss1,f,pe);
    _sw.stop();
    zm.apply(0.0f,pe);
    trace("Smooth dynamic warping time = "+_sw.time());    
    Util.writeBinary(pe,PATH+"data/"+title+"_p.dat");
    trace("Smooth dynamic warping:");
    float error;
    if (p!=null) error = Util.rmsError(pe,p,T);
  }

   /**
   * Smooth dynamic warping: plots the estimated slopes, RMS error, and time.
   */
  public void estimateSDW(
      int k, float[][][] f, float[][][] p2, float[][][] p3, String title) {
    int n1 = _s1.getCount();
    int n2 = _s2.getCount();
    int n3 = _s3.getCount();
    double r1 = 0.1;
    double r2 = 0.7;
    double r3 = 0.7;
    double h1 = 20.0;
    double h2 =  9.0;
    double h3 =  9.0;
    Sampling ss1 = new Sampling(n1);
    Sampling ss2 = new Sampling(n2);
    Sampling ss3 = new Sampling(n3);
    DynamicWarpingSlopes dws = new DynamicWarpingSlopes(k,_pmax,h1,h2,h3,
                                       r1,r2,r3,ss1,ss2,ss3);
    dws.setErrorSmoothing(2);
    float[][][] p2e = new float[n3][n2][n1];
    float[][][] p3e = new float[n3][n2][n1];
    ZeroMask zm = new ZeroMask(f);
    _sw.restart();
    dws.findSmoothSlopes(_s1,f,p2e,p3e);
    _sw.stop();
    zm.apply(0.0f,p2e);
    zm.apply(0.0f,p3e);
    trace("Smooth dynamic warping time = "+_sw.time());    
    Util.writeBinary(p2e,PATH+"data/"+title+"_p2.dat");
    Util.writeBinary(p3e,PATH+"data/"+title+"_p3.dat");
    trace("Smooth dynamic warping:");
    float error2,error3;
    if (p2!=null && p3!=null) {
      trace("p2:");
      error2 = Util.rmsError(p2e,p2,T);
      trace("p3:");
      error3 = Util.rmsError(p3e,p3,T);
    }
  }

  public void estimateTransSDW(int k, float[][][] f, String title) {
    int n1 = _s1.getCount();
    int n2 = _s2.getCount();
    int n3 = _s3.getCount();
    double r1 = 0.1;
    double r2 = 0.2;
    double r3 = 0.2;
    double h1 = 20.0;
    double h2 =  9.0;
    double h3 =  9.0;
    Sampling ss1 = new Sampling(n1);
    Sampling ss2 = new Sampling(n2);
    Sampling ss3 = new Sampling(n3);
    DynamicWarpingSlopes dws = new DynamicWarpingSlopes(k,_pmax,h1,h3,h2,
                                       r1,r3,r2,ss1,ss3,ss2);
    dws.setErrorSmoothing(2);
    _sw.restart();
    float[][][] p2 = new float[n2][n3][n1];
    float[][][] p3 = new float[n2][n3][n1];
    ZeroMask zm = new ZeroMask(f);
    dws.findSmoothSlopes(_s1,f,p2,p3);
    _sw.stop();
    zm.apply(0.0f,p2);
    zm.apply(0.0f,p3);
    trace("Smooth dynamic warping time = "+_sw.time());    
    Util.writeBinary(p2,PATH+"data/"+title+"_p2.dat");
    Util.writeBinary(p3,PATH+"data/"+title+"_p3.dat");
  }

  public void retranspose(String title) {
    int n1 = _s1.getCount();
    int n2 = _s2.getCount();
    int n3 = _s3.getCount();
    float[][][] p2 = Util.readImage(n1,n3,n2,PATH+"data/"+title+"_p2.dat");
    float[][][] p3 = Util.readImage(n1,n3,n2,PATH+"data/"+title+"_p3.dat");
    float[][][] p2t = new float[n3][n2][n1];
    float[][][] p3t = new float[n3][n2][n1];
    Util.transpose23(p2,p3t);
    Util.transpose23(p3,p2t);
    Util.writeBinary(p2t,PATH+"data/"+title+"_p2.dat");
    Util.writeBinary(p3t,PATH+"data/"+title+"_p3.dat");
  }

/**********************CHOOSING OPTIMAL PARAMETERS**********************/

  /**
   * Tests for the parameters that produce the smallest RMS error.
   * for noise = 0.0   14,1
   * for noise = 0.5   23,1
   */
  public void testOptimalSmoothLSF(String fileName, float[][] f, float[][] p, 
                                   Sampling sp1, Sampling sp2) {
    int n1  = _s1.getCount();
    int n2  = _s2.getCount();
    int np1 = sp1.getCount();
    int np2 = sp2.getCount();
    float[] param1 = Util.f(sp1.getValues());
    float[] param2 = Util.f(sp2.getValues());

    float[][] pe = new float[n2][n1];
    float[][] rmserror = new float[np2][np1]; // RMS error
    LocalSlopeFinder lsf;
    for(int i2=0; i2<np2; ++i2) {
      for(int i1=0; i1<np1; ++i1) {
        lsf = new LocalSlopeFinder(param1[i1],param2[i2],_pmax);
        lsf.findSlopes(f,pe);
        rmserror[i2][i1] = Util.rmsError(pe,p,F);
      }
    }
    Util.writeBinary(rmserror,PATH+"data/"+fileName);
  }

  /**
   * Tests for the parameters that produce the smallest RMS error.
   * for noise = 0.0   34,2
   * for noise = 0.5   75,6
   */
  public void testOptimalSmoothPWDM(String fileName, float[][] f, float[][] p, 
                                    Sampling sp1, Sampling sp2) {
    int n1  = _s1.getCount();
    int n2  = _s2.getCount();
    int np1 = sp1.getCount();
    int np2 = sp2.getCount();
    int[] param1 = Util.i(sp1.getValues());
    int[] param2 = Util.i(sp2.getValues());

    float[][] pe = new float[n2][n1];
    float[][] rmserror = new float[np2][np1]; // RMS error
    Sfdip sd = new Sfdip(-_pmax,_pmax);
    sd.setOrder(2);
    for(int i2=0; i2<np2; ++i2) {
      for(int i1=0; i1<np1; ++i1) {
        sd.setRect(param1[i1],param2[i2]);
        sd.findSlopes(_s1,_s2,f,pe);
        rmserror[i2][i1] = Util.rmsError(pe,p,F);
      }
    }
    Util.writeBinary(rmserror,PATH+"data/"+fileName);
  }

  /**
   * Tests for the parameters that produce the smallest RMS error.
   * for noise = 0.0    5,1
   * for noise = 0.5    6,1
   */
  public void testOptimalSmoothPWDD(String fileName, float[][] f, float[][] p, 
                                    Sampling sp1, Sampling sp2) {
    int n1 = _s1.getCount();
    int n2 = _s2.getCount();
    int np1 = sp1.getCount();
    int np2 = sp2.getCount();
    float[] param1 = Util.f(sp1.getValues());
    float[] param2 = Util.f(sp2.getValues());

    float[][] pe = new float[n2][n1];
    float[][] rmserror = new float[np2][np1]; // RMS error
    PlaneWaveDestructor pwd = new PlaneWaveDestructor(-_pmax,_pmax);
    for(int i2=0; i2<np2; ++i2) {
      for(int i1=0; i1<np1; ++i1) {
        pwd.setSmoothness(param1[i1],param2[i2]);
        pe = pwd.findSlopes(f);
        rmserror[i2][i1] = Util.rmsError(pe,p,F);
      }
    }
    Util.writeBinary(rmserror,PATH+"data/"+fileName);
  }

  /**
   * Tests for the parameters that produce the smallest RMS error.
   * for noise = 0.0   72,12
   * for noise = 0.5   30,9
   */
  public void testOptimalSmoothSDW(String fileName, int k, 
              float[][] f, float[][] p, Sampling sp1, Sampling sp2) {
    int n1  = _s1.getCount();
    int n2  = _s2.getCount();
    int np1 = sp1.getCount();
    int np2 = sp2.getCount();
    double[] param1 = sp1.getValues();
    double[] param2 = sp2.getValues();
    double r1 = 0.1;
    double r2 = 0.4;
    Sampling ss1 = new Sampling(n1);
    Sampling ss2 = new Sampling(n2);

    float[][] pe = new float[n2][n1];
    float[][] rmserror = new float[np2][np1]; // RMS error
    DynamicWarpingSlopes dws;
    for(int i2=0; i2<np2; ++i2) {
      for(int i1=0; i1<np1; ++i1) {
        dws = new DynamicWarpingSlopes(k,_pmax,param1[i1],param2[i2],
                                       r1,r2,ss1,ss2);
        dws.findSmoothSlopes(_s1,f,pe);
        rmserror[i2][i1] = Util.rmsError(pe,p,F);
      }
    }
    Util.writeBinary(rmserror,PATH+"data/"+fileName);
  }

  /**
   * Tests for the parameters that produce the smallest RMS error.
   * for noise = 0.0  0.1,0.6
   * for noise = 0.5  0.1,0.4
   */
  public void testOptimalStrainSDW(String fileName, int k, 
              float[][] f, float[][] p, Sampling sp1, Sampling sp2) {
    int n1  = _s1.getCount();
    int n2  = _s2.getCount();
    int np1 = sp1.getCount();
    int np2 = sp2.getCount();
    float[] param1 = Util.f(sp1.getValues());
    float[] param2 = Util.f(sp2.getValues());
    double h1 = 72.0;
    double h2 = 12.0;
    Sampling ss1 = new Sampling(n1);
    Sampling ss2 = new Sampling(n2);

    float[][] pe = new float[n2][n1];
    float[][] rmserror = new float[np2][np1]; // RMS error
    DynamicWarpingSlopes dws;
    for(int i2=0; i2<np2; ++i2) {
      for(int i1=0; i1<np1; ++i1) {
        dws = new DynamicWarpingSlopes(k,_pmax,h1,h2,
                                       param1[i1],param2[i2],ss1,ss2);
        dws.findSmoothSlopes(_s1,f,pe);
        rmserror[i2][i1] = Util.rmsError(pe,p,F);
      }
    }
    Util.writeBinary(rmserror,PATH+"data/"+fileName);
  }

/**********************OTHER EVALUATION METHODS**********************/

  public void testOrderVsTime(Sampling s, String fileName, 
      float[][] f, float[][] p) {
    int norder = s.getCount();
    int n1 = _s1.getCount();
    int n2 = _s2.getCount();
    float[] ovt = new float[norder];
    Sfdip sd = new Sfdip(-_pmax,_pmax);
    sd.setRect(75,6);
    sd.setNiter(_niter);
    float[][] p_pwdm = new float[n2][n1]; //pwd w/initial p=0 (Madagascar)
    for (int i=1; i<norder; ++i) {
        _sw.restart();
        sd.setOrder(i);
        sd.findSlopes(_s1,_s2,f,p_pwdm);
        ovt[i] = (float)_sw.time();
        trace("ovt["+i+"]= "+ovt[i]);
    }
    Util.writeBinary(ovt,PATH+fileName);
  }

  public void testRmsErrorCurveLSF(Sampling sn, String fileName, 
      float sigma1, float sigma2) {
    trace("Structure tensor:");
    int n1 = _s1.getCount();
    int n2 = _s2.getCount();
    int nrms = sn.getCount();
    float[] nsratio = Util.f(sn.getValues());

    float[][][] fandp = new float[2][n2][n1];
    float[][] f = new float[n2][n1];
    float[][] p = new float[n2][n1];

    //Structure tensor
    LocalSlopeFinder lsf = new LocalSlopeFinder(sigma1,sigma2,_pmax);
    float[][] pe = new float[n2][n1];
    float[] rms_error = new float[nrms];
    for (int i=0; i<nrms; ++i) {
      trace("i="+i);
      fandp = FakeData.seismicAndSlopes2d2014A(nsratio[i],F);
      f = fandp[0];
      p = fandp[1];
      lsf.findSlopes(f,pe);
      rms_error[i] = Util.rmsError(pe,p,false);
    }
    Util.writeBinary(rms_error,PATH+fileName);
  }

  public void testRmsErrorCurveLSF(Sampling sn, String fn2, String fn3, 
      float sigma1, float sigma2, float sigma3) {
    trace("Structure tensor:");
    int n1 = _s1.getCount();
    int n2 = _s2.getCount();
    int n3 = _s3.getCount();
    int nrms = sn.getCount();
    float[] nsratio = Util.f(sn.getValues());

    float[][][][] fandp = new float[3][n3][n2][n1];
    float[][][] f  = new float[n3][n2][n1];
    float[][][] p2 = new float[n3][n2][n1];
    float[][][] p3 = new float[n3][n2][n1];

    //Structure tensor
    LocalSlopeFinder lsf = new LocalSlopeFinder(sigma1,sigma2,sigma3,_pmax);
    float[][][] p2e = new float[n3][n2][n1];
    float[][][] p3e = new float[n3][n2][n1];
    float[] rms_error2 = new float[nrms];
    float[] rms_error3 = new float[nrms];
    for (int i=0; i<nrms; ++i) {
      trace("i="+i);
      fandp = FakeData.seismicAndSlopes3d2014A(nsratio[i],F);
      f  = fandp[0];
      p2 = fandp[1];
      p3 = fandp[2];
      lsf.findSlopes(f,p2e,p3e,null);
      rms_error2[i] = Util.rmsError(p2e,p2,false);
      rms_error3[i] = Util.rmsError(p3e,p3,false);
    }
    Util.writeBinary(rms_error2,PATH+fn2);
    Util.writeBinary(rms_error3,PATH+fn3);
  }

  public void testRmsErrorCurvePWD(Sampling sn, String fileName, 
      int rect1, int rect2) {
    trace("Plane-wave destructor:");
    int n1 = _s1.getCount();
    int n2 = _s2.getCount();
    int nrms = sn.getCount();
    float[] nsratio = Util.f(sn.getValues());

    float[][][] fandp = new float[2][n2][n1];
    float[][] f = new float[n2][n1];
    float[][] p = new float[n2][n1];

    //Plane-wave destruction filter
    Sfdip sd = new Sfdip(-_pmax,_pmax);
    sd.setRect(rect1,rect2);
    sd.setOrder(4);
    float[][] pe = new float[n2][n1];
    float[] rms_error = new float[nrms];
    for (int i=0; i<nrms; ++i) {
      trace("i="+i);
      fandp = FakeData.seismicAndSlopes2d2014A(nsratio[i],F);
      f = fandp[0];
      p = fandp[1];
      sd.findSlopes(_s1,_s2,f,pe);
      rms_error[i] = Util.rmsError(pe,p,false);
    }
    Util.writeBinary(rms_error,PATH+fileName);
  }

  public void testRmsErrorCurvePWD(Sampling sn, String fn2, String fn3, 
      int rect1, int rect2, int rect3) {
    trace("Plane-wave destructor:");
    int n1 = _s1.getCount();
    int n2 = _s2.getCount();
    int n3 = _s3.getCount();
    int nrms = sn.getCount();
    float[] nsratio = Util.f(sn.getValues());

    float[][][][] fandp = new float[3][n3][n2][n1];
    float[][][] f  = new float[n3][n2][n1];
    float[][][] p2 = new float[n3][n2][n1];
    float[][][] p3 = new float[n3][n2][n1];

    //Plane-wave destruction filter
    Sfdip sd = new Sfdip(-_pmax,_pmax);
    sd.setRect(rect1,rect2,rect3);
    sd.setOrder(4);
    sd.setN4(2);
    float[][][] p2e = new float[n3][n2][n1];
    float[][][] p3e = new float[n3][n2][n1];
    float[] rms_error2 = new float[nrms];
    float[] rms_error3 = new float[nrms];
    for (int i=0; i<nrms; ++i) {
      trace("i="+i);
      fandp = FakeData.seismicAndSlopes3d2014A(nsratio[i],F);
      f  = fandp[0];
      p2 = fandp[1];
      p3 = fandp[2];
      sd.findSlopes(_s1,_s2,_s3,f,p2e,p3e);
      rms_error2[i] = Util.rmsError(p2e,p2,false);
      rms_error3[i] = Util.rmsError(p3e,p3,false);
    }
    Util.writeBinary(rms_error2,PATH+fn2);
    Util.writeBinary(rms_error3,PATH+fn3);
  }

  public void testRmsErrorCurveSDW(Sampling sn, String fileName, 
      int k, double r1, double r2, double h1, double h2) {
    trace("Smooth dynamic warping:");
    int n1 = _s1.getCount();
    int n2 = _s2.getCount();
    int nrms = sn.getCount();
    float[] nsratio = Util.f(sn.getValues());

    float[][][] fandp = new float[2][n2][n1];
    float[][] f = new float[n2][n1];
    float[][] p = new float[n2][n1];
    Sampling ss1 = new Sampling(n1);
    Sampling ss2 = new Sampling(n2);

    //Smooth dynamic warping
    DynamicWarpingSlopes dws = new DynamicWarpingSlopes(k,_pmax,h1,h2,
                                       r1,r2,ss1,ss2);
    float[][] pe = new float[n2][n1];
    float[] rms_error = new float[nrms];
    for (int i=0; i<nrms; ++i) {
      trace("i="+i);
      fandp = FakeData.seismicAndSlopes2d2014A(nsratio[i],F);
      f = fandp[0];
      p = fandp[1];
      dws.findSmoothSlopes(_s1,f,pe);
      rms_error[i] = Util.rmsError(pe,p,false);
    }
    Util.writeBinary(rms_error,PATH+fileName);
  }

  public void testRmsErrorCurveSDW(Sampling sn, String fn2, String fn3, int k,
      double r1, double r2, double r3, double h1, double h2, double h3) {
    trace("Smooth dynamic warping:");
    int n1 = _s1.getCount();
    int n2 = _s2.getCount();
    int n3 = _s3.getCount();
    int nrms = sn.getCount();
    float[] nsratio = Util.f(sn.getValues());

    float[][][][] fandp = new float[3][n3][n2][n1];
    float[][][] f  = new float[n3][n2][n1];
    float[][][] p2 = new float[n3][n2][n1];
    float[][][] p3 = new float[n3][n2][n1];
    Sampling ss1 = new Sampling(n1);
    Sampling ss2 = new Sampling(n2);
    Sampling ss3 = new Sampling(n3);

    //Smooth dynamic warping
    DynamicWarpingSlopes dws = new DynamicWarpingSlopes(k,_pmax,h1,h2,h3,
                                       r1,r2,r3,ss1,ss2,ss3);
    float[][][] p2e = new float[n3][n2][n1];
    float[][][] p3e = new float[n3][n2][n1];
    float[] rms_error2 = new float[nrms];
    float[] rms_error3 = new float[nrms];
    for (int i=0; i<nrms; ++i) {
      trace("i="+i);
      fandp = FakeData.seismicAndSlopes3d2014A(nsratio[i],F);
      f  = fandp[0];
      p2 = fandp[1];
      p3 = fandp[2];
      dws.findSmoothSlopes(_s1,f,p2e,p3e);
      rms_error2[i] = Util.rmsError(p2e,p2,false);
      rms_error3[i] = Util.rmsError(p3e,p3,false);
    }
    Util.writeBinary(rms_error2,PATH+fn2);
    Util.writeBinary(rms_error3,PATH+fn3);
  }

  /**
   * nni = number noise images
   */
  public void testMeanCurveLSF(Sampling sn, String fileName, int nni,
      float sigma1, float sigma2) {
    boolean one = true;
    int n1 = _s1.getCount();
    int n2 = _s2.getCount();
    int n = sn.getCount();
    float[] nsratio = Util.f(sn.getValues());

    float[][] pbar = new float[n2][n1];
    float[][] p = FakeData.seismicAndSlopes2d2014A(_noise,F)[1];
    float[] mean_error = new float[n];
    for (int i=0; i<n; ++i) {
      trace("nsratio= "+nsratio[i]);
      pbar = testSampleMeanLSF(fileName,nni,nsratio[i],sigma1,sigma2);
      mean_error[i] = sum(sub(pbar,p))/(n1*n2);
    }
    Util.writeBinary(mean_error,PATH+fileName);
  }

  /**
   * nni = number noise images
   */
  public void testMeanCurvePWD(Sampling sn, String fileName, int nni,
      int rect1, int rect2) {
    boolean one = true;
    int n1 = _s1.getCount();
    int n2 = _s2.getCount();
    int n = sn.getCount();
    float[] nsratio = Util.f(sn.getValues());

    float[][] pbar = new float[n2][n1];
    float[][] p = FakeData.seismicAndSlopes2d2014A(_noise,F)[1];
    float[] mean_error = new float[n];
    for (int i=0; i<n; ++i) {
      trace("nsratio= "+nsratio[i]);
      pbar = testSampleMeanPWD(fileName,nni,nsratio[i],rect1,rect2);
      mean_error[i] = sum(sub(pbar,(p)))/(n1*n2);
    }
    Util.writeBinary(mean_error,PATH+fileName);
  }

  /**
   * nni = number noise images
   */
  public void testMeanCurveSDW(Sampling sn, String fileName, int nni,
      int k, double r1, double r2, double h1, double h2) {
    boolean one = true;
    int n1 = _s1.getCount();
    int n2 = _s2.getCount();
    int n = sn.getCount();
    float[] nsratio = Util.f(sn.getValues());

    float[][] pbar = new float[n2][n1];
    float[][] p = FakeData.seismicAndSlopes2d2014A(_noise,F)[1];
    float[] mean_error = new float[n];
    for (int i=0; i<n; ++i) {
      trace("nsratio= "+nsratio[i]);
      pbar = testSampleMeanSDW(fileName,nni,nsratio[i],k,r1,r2,h1,h2);
      mean_error[i] = sum(sub(pbar,(p)))/(n1*n2);
    }
    Util.writeBinary(mean_error,PATH+fileName);
  }

  public float[][] testSampleMeanLSF(String fileName, int n, float noise,
      float sigma1, float sigma2) {
    int n1 = _s1.getCount();
    int n2 = _s2.getCount();

    float[][][] fandp = new float[2][n2][n1];
    float[][] f = new float[n2][n1];
    float[][] pe = new float[n2][n1];
    float[][] psum = new float[n2][n1];
    LocalSlopeFinder lsf = new LocalSlopeFinder(sigma1,sigma2,_pmax);
    for (int i=0; i<n; ++i) {
      fandp = FakeData.seismicAndSlopes2d2014A(noise,T);
      f = fandp[0];
      trace("i= "+i);
      lsf.findSlopes((f),pe);
      psum = add(psum,pe);
    }
    float[][] mean = div(psum,n);
    Util.writeBinary(mean,PATH+fileName);
    return mean;
  }

  public float[][] testSampleMeanPWD(String fileName, int n, float noise,
      int rect1, int rect2) {
    int n1 = _s1.getCount();
    int n2 = _s2.getCount();

    float[][][] fandp = new float[2][n2][n1];
    float[][] f = new float[n2][n1];
    float[][] pe = new float[n2][n1];
    float[][] psum = new float[n2][n1];
    Sfdip sd = new Sfdip(-_pmax,_pmax);
    sd.setRect(rect1,rect2);
    sd.setOrder(4);
    for (int i=0; i<n; ++i) {
      fandp = FakeData.seismicAndSlopes2d2014A(noise,T);
      f = fandp[0];
      trace("i= "+i);
      sd.findSlopes(_s1,_s2,(f),pe);
      psum = add(psum,pe);
    }
    float[][] mean = div(psum,n);
    Util.writeBinary(mean,PATH+fileName);
    return mean;
  }

  public float[][] testSampleMeanSDW(String fileName, int n, float noise, 
      int k, double r1, double r2, double h1, double h2) {
    int n1 = _s1.getCount();
    int n2 = _s2.getCount();

    float[][][] fandp = new float[2][n2][n1];
    float[][] f = new float[n2][n1];
    float[][] pe = new float[n2][n1];
    float[][] psum = new float[n2][n1];
    Sampling ss1 = new Sampling(n1);
    Sampling ss2 = new Sampling(n2);
    DynamicWarpingSlopes dws = new DynamicWarpingSlopes(k,_pmax,h1,h2,
                                       r1,r2,ss1,ss2);
    for (int i=0; i<n; ++i) {
      fandp = FakeData.seismicAndSlopes2d2014A(noise,T);
      f = fandp[0];
      trace("i= "+i);
      dws.findSmoothSlopes(_s1,f,pe);
      psum = add(psum,pe);
    }
    float[][] mean = div(psum,n);
    Util.writeBinary(mean,PATH+fileName);
    return mean;
  }

  public void testSampleStdDevLSF(String fileName, int nni,
      float sigma1, float sigma2) {
    int n1 = _s1.getCount();
    int n2 = _s2.getCount();

    float[][][] fandp = new float[2][n2][n1];
    float[][] f = new float[n2][n1];
    float[][] p = new float[n2][n1];

    float[][] pe = new float[n2][n1];
    float[][] sumdiffsq = new float[n2][n1]; //sum of differences squared
    LocalSlopeFinder lsf = new LocalSlopeFinder(sigma1,sigma2,_pmax);
    for (int i=0; i<nni; ++i) {
      fandp = FakeData.seismicAndSlopes2d2014A(_noise,T);
      f = fandp[0];
      p = fandp[1];
      trace("i= "+i);
      lsf.findSlopes(f,pe);
      sumdiffsq = add(sumdiffsq,pow(sub(pe,p),2));
    }
    float[][] stddev = sqrt(div(sumdiffsq,nni));
    Util.writeBinary(stddev,PATH+fileName);
  }

  public void testSampleStdDevPWD(String fileName, int nni,
      int rect1, int rect2) {
    int n1 = _s1.getCount();
    int n2 = _s2.getCount();

    float[][][] fandp = new float[2][n2][n1];
    float[][] f = new float[n2][n1];
    float[][] p = new float[n2][n1];

    float[][] pe = new float[n2][n1];
    float[][] sumdiffsq = new float[n2][n1]; //sum of differences squared
    Sfdip sd = new Sfdip(-_pmax,_pmax);
    sd.setRect(rect1,rect2);
    sd.setOrder(4);
    for (int i=0; i<nni; ++i) {
      fandp = FakeData.seismicAndSlopes2d2014A(_noise,T);
      f = fandp[0];
      p = fandp[1];
      trace("i= "+i);
      sd.findSlopes(_s1,_s2,f,pe);
      sumdiffsq = add(sumdiffsq,pow(sub(pe,p),2));
    }
    float[][] stddev = sqrt(div(sumdiffsq,nni));
    Util.writeBinary(stddev,PATH+fileName);
  }

  public void testSampleStdDevSDW(String fileName, int nni,
      int k, double r1, double r2, double h1, double h2) {
    int n1 = _s1.getCount();
    int n2 = _s2.getCount();
    Sampling ss1 = new Sampling(n1);
    Sampling ss2 = new Sampling(n2);

    float[][][] fandp = new float[2][n2][n1];
    float[][] f = new float[n2][n1];
    float[][] p = new float[n2][n1];

    float[][] pe = new float[n2][n1];
    float[][] sumdiffsq = new float[n2][n1]; //sum of differences squared
    DynamicWarpingSlopes dws = new DynamicWarpingSlopes(k,_pmax,h1,h2,
                                       r1,r2,ss1,ss2);
    for (int i=0; i<nni; ++i) {
      fandp = FakeData.seismicAndSlopes2d2014A(_noise,T);
      f = fandp[0];
      p = fandp[1];
      trace("i= "+i);
      dws.findSmoothSlopes(_s1,f,pe);
      sumdiffsq = add(sumdiffsq,pow(sub(pe,p),2));
    }
    float[][] stddev = sqrt(div(sumdiffsq,nni));
    Util.writeBinary(stddev,PATH+fileName);
  }

  /****************************PLOTTING****************************/

  public void plotCurve(Sampling s, String fileName, String title, String hl, 
      String vl, float cmin, float cmax) {
    plotCurve(s,fileName,fileName,title,hl,vl,cmin,cmax);
  }

  public void plotCurve(Sampling s, String fileName1, String fileName2, 
      String title, String hl, String vl, float cmin, float cmax) {
    plotCurve(s,fileName1,fileName2,fileName1,title,hl,vl,cmin,cmax);
  }

  public void plotCurve(Sampling s, String fileName1, String fileName2,
      String fileName3, String title, String hl, String vl, 
      float cmin, float cmax) {
    boolean one = true;
    int n = s.getCount();
    float[] f1 = Util.readImage(n,PATH+fileName1);
    float[] f2 = Util.readImage(n,PATH+fileName2);
    float[] f3 = Util.readImage(n,PATH+fileName3);
    // paint 
    Plot.plot(s,f1,f2,f3,title,hl,vl,_fw,_fh,cmin,cmax,_slide,one,_paint);
  }

  public void plot2D(float[][] f, String title) {
    boolean one = true;
    int n1 = _s1.getCount();
    int n2 = _s2.getCount();
    float[][] p = Util.readImage(n1,n2,PATH+"data/"+title+"_p.dat");
    // clip, interp, title, paint, color
    String cbl = "slope (samples/trace)"; //colorbar label
    Plot.plot(_s1,_s2,f,p,title+num,cbl,_fw,_fh,-_cmax,_cmax,_clip,_title,
        _paint,_slide,one);
  }

  public void plot3DSub(float[][][] f, String title1, String title2) {
    boolean one = true;
    int n1 = _s1.getCount();
    int n2 = _s2.getCount();
    int n3 = _s3.getCount();
    float[][][] p2 = Util.readImage(n1,n2,n3,PATH+"data/"+title1+"_p2.dat");
    float[][][] p3 = Util.readImage(n1,n2,n3,PATH+"data/"+title1+"_p3.dat");
    float[][][] p2t = Util.readImage(n1,n2,n3,PATH+"data/"+title2+"_p2.dat");
    float[][][] p3t = Util.readImage(n1,n2,n3,PATH+"data/"+title2+"_p3.dat");
    //Plot.plot(_s1,_s2,_s3,f,sub(p2,p2t),title1+"_p2_slices",-_cmax,_cmax,
    //_paint);
    //Plot.plot(_s1,_s2,_s3,f,sub(p3,p3t),title1+"_p3_slices",-_cmax,_cmax,
    //_paint);
    Plot.plot(_s1,_s2,_s3,f,sub(p2,p2t),"slope (samples/trace)",title1+
        "_p2_panels",0.6f,0.8f,-_cmax,_cmax,_paint,_slide,one);
    Plot.plot(_s1,_s2,_s3,f,sub(p3,p3t),"slope (samples/trace)",title1+
       "_p3_panels",0.6f,0.8f,-_cmax,_cmax,_paint,_slide,one);
    trace("max and min values= "+max(sub(p3,p3t))+", "+min(sub(p2,p2t)));
  }

  public void plot3D(float[][][] f, String title) {
    boolean one = true;
    int n1 = _s1.getCount();
    int n2 = _s2.getCount();
    int n3 = _s3.getCount();
    float[][][] p2 = Util.readImage(n1,n2,n3,PATH+"data/"+title+"_p2.dat");
    float[][][] p3 = Util.readImage(n1,n2,n3,PATH+"data/"+title+"_p3.dat");
    Plot.plot(_s1,_s2,_s3,f,p2,title+"_p2_slices",-_cmax,_cmax,_paint);
    Plot.plot(_s1,_s2,_s3,f,p3,title+"_p3_slices",-_cmax,_cmax,_paint);
    Plot.plot(_s1,_s2,_s3,f,p2,"slope (samples/trace)",title+"_p2_panels",
        0.6f,0.8f,-_cmax,_cmax,_paint,_slide,one);
    Plot.plot(_s1,_s2,_s3,f,p3,"slope (samples/trace)",title+"_p3_panels",
        0.6f,0.8f,-_cmax,_cmax,_paint,_slide,one);
  }

  public void plotError(String t1, String t2, String t3) {
    int n  = 20;
    int n1 = _s1.getCount();
    int n2 = _s2.getCount();
    int n3 = _s3.getCount();
    float[][][] p2t1 = Util.readImage(n1,n2,n3,PATH+"data/"+t1+"_p2.dat");
    float[][][] p3t1 = Util.readImage(n1,n2,n3,PATH+"data/"+t1+"_p3.dat");
    float[][][] p2t2 = Util.readImage(n1,n2,n3,PATH+"data/"+t2+"_p2.dat");
    float[][][] p3t2 = Util.readImage(n1,n2,n3,PATH+"data/"+t2+"_p3.dat");
    float[][][] p2t3 = Util.readImage(n1,n2,n3,PATH+"data/"+t3+"_p2.dat");
    float[][][] p3t3 = Util.readImage(n1,n2,n3,PATH+"data/"+t3+"_p3.dat");
    Sampling s = new Sampling(n,1,1);
    float[] p   = new float[n];
    float[] pbar  = new float[n];
    float[] e1  = new float[n];
    float[] pe1 = new float[n];
    float[] e2  = new float[n];
    float[] pe2 = new float[n];
    float[] e3  = new float[n];
    float[] pe3 = new float[n];
    float[] e4  = new float[n];
    float med_error1 = Util.errorStatistic(p2t1,p3t1,e1,p,pe1,pbar);
    float med_error2 = Util.errorStatistic(p2t2,p3t2,e2,p,pe2,pbar);
    float med_error3 = Util.errorStatistic(p2t3,p3t3,e3,p,pe3,pbar);
    float med_error4 = Util.errorStatistic(p2t3,p3t3,e4,p,p,pbar);
    pbar = add(pe1,pe2);
    pbar = add(pbar,pe3);
    pbar = add(pbar,p);
    pbar = div(pbar,4.0f);
    med_error1 = Util.errorStatistic(p2t1,p3t1,e1,p,pe1,pbar);
    med_error2 = Util.errorStatistic(p2t2,p3t2,e2,p,pe2,pbar);
    med_error3 = Util.errorStatistic(p2t3,p3t3,e3,p,pe3,pbar);
    med_error4 = Util.errorStatistic(p2t3,p3t3,e4,p,p,pbar);
    trace(t1+" median error= "+med_error1);
    trace(t2+" median error= "+med_error2);
    trace(t3+" median error= "+med_error3);
    trace("elias median error= "+med_error4);
    //Plot.plot(s,pe1,pe2,pe3,p,"est_slope_plot","pick #",
    //    "slope (samples/trace)",_fw,_fh,-1f,1f,_slide,false,_paint);
    Plot.plot(s,e1,e2,e3,e4,"median_error_plot","pick #",
        "slope difference (samples/trace)",_fw,_fh,-0.2f,0.2f,_slide,false,
        _paint);
  }

  /**
   * Plots the mean and std dev images.
   */
  public void plotImage(String title, String hl, String vl, String cbl,
      String fileName, float cmin, float cmax) {
    boolean one = true;
    int n1 = _s1.getCount();
    int n2 = _s2.getCount();
    float[][] f = Util.readImage(n1,n2,PATH+fileName);
    // clip, interp, title, paint, color
    Plot.plot(_s1,_s2,f,title,cbl,_fw,_fh,cmin,cmax,T,F,_title,_paint,T,
        _slide,one);
  }

  public void plotTeaser(float[][] f, String title) {
    boolean paint = true;
    int n1 = _s1.getCount();
    int n2 = _s2.getCount();
    float[][] pe_lsf = Util.readImage(n1,n2,PATH+"data/lsf_"+title+".dat");
    float[][] pe_pwd = Util.readImage(n1,n2,PATH+"data/pwdm_"+title+".dat");
    float[][] pe_sdw = Util.readImage(n1,n2,PATH+"data/sdw_"+title+".dat");
    Plot.plot(_s1,_s2,f,pe_lsf,pe_pwd,pe_sdw,"teaser",_fw,_fh,_slide,paint);
  }

  public void plotOptimalParameters(String fileName, String hl, String vl, 
      Sampling sp1, Sampling sp2, float cmin, float cmax) {
    boolean one = true;
    int np1 = sp1.getCount();
    int np2 = sp2.getCount();
    float[] param1 = Util.f(sp1.getValues());
    float[] param2 = Util.f(sp2.getValues());
    int[] error_index = new int[2];    
    float[][] rmserror = Util.readImage(np1,np2,PATH+"data/"+fileName);
    float min_error = min(rmserror,error_index);
    trace(vl+"= "+param1[error_index[0]]+" "+hl+"= "+param2[error_index[1]]+
        " "+" Minimum Error Value= "+min_error);    
    // clip, interp, title, paint, color
    String cbl = "rms error (samples/trace)"; //colorbar label
    Plot.plot(sp1,sp2,rmserror,"rms_error_"+
        fileName.split("\\.(?=[^\\.]+$)")[0],cbl,_fw,_fh,cmin,cmax,
        T,T,_title,_paint,T,_slide,one);
  }

  /*
  private static void goTestSlopeVsError() {
    int n1 = 501;
    int n2 = 501;

    float noise = 0.5f;

    float[][][] synthAndSlope = FakeData.seismicAndSlopes2d2014A(noise,F);
    float[][] exact_slope = synthAndSlope[1];
    float[][] lsf_mean = Util.readImage(n1,n2,
        "/Users/earias/Home/git/ea/bench/src/slopes/data/lsf_mean.dat");
    float[][] mad_mean = Util.readImage(n1,n2,
        "/Users/earias/Home/git/ea/bench/src/slopes/data/mad_mean.dat");
    float[][] lsf_sd = Util.readImage(n1,n2,
        "/Users/earias/Home/git/ea/bench/src/slopes/data/lsf_sd.dat");
    float[][] mad_sd = Util.readImage(n1,n2,
        "/Users/earias/Home/git/ea/bench/src/slopes/data/mad_sd.dat");

    float[][] mad_abs = sub(mad_mean,exact_slope); 
    float[][] lsf_abs = sub(lsf_mean,exact_slope); 

    float[] mad_1Dslope = new float[n2*n1];
    float[] lsf_1Dslope = new float[n2*n1];
    float[] mad_1Dabs = new float[n2*n1];
    float[] lsf_1Dabs = new float[n2*n1];
    float[] exact_1Dslope= new float[n2*n1];
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        mad_1Dslope[i1+i2*n1] = abs(atan(mad_mean[i2][i1]));
        lsf_1Dslope[i1+i2*n1] = abs(atan(lsf_mean[i2][i1]));
        mad_1Dabs[i1+i2*n1] = abs(mad_sd[i2][i1]);
        lsf_1Dabs[i1+i2*n1] = abs(lsf_sd[i2][i1]);
        exact_1Dslope[i1+i2*n1] = exact_slope[i2][i1];
      }
    }
    Util.writeBinary(exact_1Dslope,
        "/Users/earias/Home/git/ea/bench/src/slopes/data/exact_slope1D.dat");
    Util.writeBinary(lsf_1Dslope,
        "/Users/earias/Home/git/ea/bench/src/slopes/data/struct_tens1D.dat");
    Util.writeBinary(mad_1Dslope,
        "/Users/earias/Home/git/ea/bench/src/slopes/data/plane_wave1D.dat");
    float fw = 0.75f; //fraction width for slide
    float fh = 0.8f; //fraction height for slide
    //Plot.plotHistogram(exact_1Dslope, lsf_1Dslope, mad_1Dslope,
    //  "Histogram",fw,fh,F);

    // paint
    //Plot.plot(mad_1Dslope,mad_1Dabs,"PWD Slope vs Absolute Error",
    //"Slope (samples/trace)","(samples/trace)",fw,fh,T);
    //Plot.plot(lsf_1Dslope,lsf_1Dabs,"LSF Slope vs Absolute Error",
    //"Slope (samples/trace)","(samples/trace)",fw,fh,T);
    Plot.plot(exact_1Dslope,mad_1Dabs,"PWD Slope vs Absolute Error",
    "Slope (samples/trace)","(samples/trace)",fw,fh,T);
    Plot.plot(exact_1Dslope,lsf_1Dabs,"LSF Slope vs Absolute Error",
    "Slope (samples/trace)","(samples/trace)",fw,fh,T);
  }*/

  private static void trace(String s) {
    System.out.println(s);
  }

  ///////////////////PRIVATE VARIABLES///////////////////////
  private static final String PATH = 
    "/users/elias.arias/Home/git/ea/bench/src/util/";
  private static final int _niter = 5;
  private static final float pi = FLT_PI;      
  private static final float _fw = 0.75f; //fraction width for slide
  private static final float _fh = 0.9f; //fraction height for slide
  private static final float _cmax = 2.0f;
  private static final boolean T = true;
  private static final boolean F = false;  
  private static final boolean _title = false;
  private static final boolean _paint = false;
  private static final boolean _clip = true;
  private static final boolean _slide = true;

  private Stopwatch _sw = new Stopwatch();
  private float _noise;
  private int num;
  private float _pmax;
  private Sampling _s1;
  private Sampling _s2;
  private Sampling _s3;
}
