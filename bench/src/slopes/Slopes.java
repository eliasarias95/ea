package slopes;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.Check;
import edu.mines.jtk.util.Stopwatch;

import static edu.mines.jtk.util.ArrayMath.*;

import util.*;
import warp.DynamicWarpingR;
import dnp.LocalSlopeFinder;
import dnp.PlaneWaveDestructor;


import javax.swing.*;

/**
 *  
 * @author Elias Arias, Colorado School of Mines, CWP
 * @version 27.1.2015
 */

/**
 * float[][][][] makeConstantSlope(float freq, float p2, float p3, 
 *                                 int n1, int n2, int n3)
 *
 * float[][][] superSample(int k, float[][][] f)
 */

public class Slopes{
  
  public Slopes(float noise, float pmax, Sampling s1, Sampling s2) {
    this(noise,pmax,s1,s2,null);
  }

  public Slopes(float noise, float pmax, Sampling s1, Sampling s2, Sampling s3){
    _noise = noise;
    _pmax = pmax;
    _s1 = s1;
    _s2 = s2;
    _s3 = s3;
  }

  /**
   * Creates a 2D synthetic seismic image with complex geologic structures and
   * known slope values.
   * @param noise RMS noise to RMS signal ratio for the generated image.
   * @param f array[501][501] for the seismic image.
   * @param p array[501][501] for the know slope values.
   */
  public static void makeSyntheticComplex(float noise, float[][] f, 
                                                       float[][] p) {
    int n2 = f.length;
    int n1 = f[0].length;
    float[][][] fandp = FakeData.seismicAndSlopes2d2014A(noise,false);
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        f[i2][i1] = fandp[0][i2][i1];
        p[i2][i1] = fandp[1][i2][i1];
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
  public static void makeSyntheticConstant(float freq, float p2, 
                                           float[][] f, float[][] p) {
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
          f[i3][i2][i1] = cos(2.0f*pi*(pc2*freq*i2+pc3*freq-freq*i1));
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
  public static void makeSyntheticComplex(float noise, 
                                          float[][][] f, float[][][] p) {
    float[][][][] fandp = FakeData.seismicAndSlopes3d2014A(noise);
    int n3 = f.length;
    int n2 = f[0].length;
    int n1 = f[0][0].length;
    for (int i3=0; i3<n3; ++i3){
      for (int i2=0; i2<n2; ++i2){
        for (int i1=0; i1<n1; ++i1){
          f[i3][i2][i1] = fandp[0][i3][i2][i1];
          p[i3][i2][i1] = fandp[1][i3][i2][i1];
        }
      }
    }
  }

  /**
   * Takes a real seismic image from the GOM and places it in the array f.
   * @param f array[920][301] that holds the data for the GOM seismic image.
   */
  public static void makeRealGOM(float[][] f) {
    int n1 = 301;
    int n2 = 920;
    float[][] temp = Util.readImage(n1,n2,PATH+"data/gom.dat");
    for (int i2=0; i2<n2; ++i2){
      for (int i1=0; i1<n1; ++i1){
        f[i2][i1] = temp[i2][i1];
      }
    }
    //float[][] fs = new float[n2][n1];
    //mul(f,.001f,fs);
    //Util.writeBinary(fs,PATH+"data/gom_scaled.dat");
  }

////////////////////////ESTIMATING SLOPES/////////////////////////////

  /**
   * Structure tensor: plots the estimated slopes.
   */
  public void plotLSF(float[][] f) {
    int n1 = _s1.getCount();
    int n2 = _s2.getCount();
    double d1 = _s1.getDelta();
    double d2 = _s2.getDelta();
    LocalSlopeFinder lsf = new LocalSlopeFinder(23.0f,1.0f,_pmax);
    float[][] pe = new float[n2][n1];
    _sw.restart();
    lsf.findSlopes(f,pe);
    pe = mul(pe,(float)(d1/d2));
    _sw.stop();
    trace("Structure tensor time = "+_sw.time());    

    // clip, interp, title, paint, colorbar, color
    String hl = "Traces"; //horizontal label
    String vl = "Samples"; //vertical label
    String cbl = "slope (samples/trace)"; //colorbar label
    Plot.plot(_s1,_s2,pe,"LSF noise= "+_noise,hl,vl,cbl,
        _fw,_fh,-_clipMax,_clipMax,
        _clip,F,_title,_paint,T,T);
  }

  /**
   * Structure tensor: plots the estimated slopes and RMS error.
   */
  public void plotLSF(float[][] f, float[][] p) {
    int n1 = _s1.getCount();
    int n2 = _s2.getCount();
    double d1 = _s1.getDelta();
    double d2 = _s2.getDelta();
    LocalSlopeFinder lsf = new LocalSlopeFinder(23.0f,1.0f,_pmax);
    float[][] pe = new float[n2][n1];
    _sw.restart();
    lsf.findSlopes(f,pe);
    pe = mul(pe,(float)(d1/d2));
    _sw.stop();
    trace("Structure tensor time = "+_sw.time());    
    trace("Structure tensor:");
    float error_lsf = Util.rmsError(pe,p,T);

    // clip, interp, title, paint, colorbar, color
    String hl = "Traces"; //horizontal label
    String vl = "Samples"; //vertical label
    String cbl = "slope (samples/trace)"; //colorbar label
    Plot.plot(_s1,_s2,pe,"LSF noise= "+_noise,hl,vl,cbl,
        _fw,_fh,-_clipMax,_clipMax,
        _clip,F,_title,_paint,T,T);
  }

  /**
   * PWD Madagascar: plots the estimated slopes.
   */
  public void plotPWDM(float[][] f) {
    int n1 = _s1.getCount();
    int n2 = _s2.getCount();
    double d1 = _s1.getDelta();
    double d2 = _s2.getDelta();
    Sfdip sd = new Sfdip(-_pmax,_pmax);
    sd.setRect(75,6);
    sd.setOrder(4);
    sd.setNiter(_niter);
    //sd.setBoth("y");
    float[][] pe = new float[n2][n1]; //pwd w/initial p=0 (Madagascar)
    _sw.restart();
    sd.findSlopes(_s1,_s2,f,pe);
    pe = mul(pe,(float)(d1/d2));
    _sw.stop();
    trace("Madagascar PWD time = "+_sw.time());    

    // interp, title, paint, colorbar, color
    String hl = "Traces"; //horizontal label
    String vl = "Samples"; //vertical label
    String cbl = "slope (samples/trace)"; //colorbar label
    Plot.plot(_s1,_s2,pe,"PWDM noise= "+_noise,hl,vl,cbl,
        _fw,_fh,-_clipMax,_clipMax,
        _clip,F,_title,_paint,T,T);
  }

  /**
   * PWD Madagascar: plots the estimated slopes and RMS error.
   */
  public void plotPWDM(float[][] f, float[][] p) {
    int n1 = _s1.getCount();
    int n2 = _s2.getCount();
    double d1 = _s1.getDelta();
    double d2 = _s2.getDelta();
    Sfdip sd = new Sfdip(-_pmax,_pmax);
    sd.setRect(75,6);
    sd.setOrder(4);
    sd.setNiter(_niter);
    //sd.setBoth("y");
    float[][] pe = new float[n2][n1]; //pwd w/initial p=0 (Madagascar)
    _sw.restart();
    sd.findSlopes(_s1,_s2,f,pe);
    pe = mul(pe,(float)(d1/d2));
    _sw.stop();
    trace("Madagascar PWD time = "+_sw.time());    
    trace("Madagascar PWD:");
    float error_pwdm = Util.rmsError(pe,p,T);

    // interp, title, paint, colorbar, color
    String hl = "Traces"; //horizontal label
    String vl = "Samples"; //vertical label
    String cbl = "slope (samples/trace)"; //colorbar label
    Plot.plot(_s1,_s2,pe,"PWDM noise= "+_noise,hl,vl,cbl,
        _fw,_fh,-_clipMax,_clipMax,
        _clip,F,_title,_paint,T,T);
  }

  /**
   * PWD Dave: plots the estimated slopes and RMS error.
   */
  public void plotPWDD(float[][] f) {
    int n1 = _s1.getCount();
    int n2 = _s2.getCount();
    double d1 = _s1.getDelta();
    double d2 = _s2.getDelta();
    PlaneWaveDestructor pwd = new PlaneWaveDestructor(-_pmax,_pmax);
    pwd.setSmoothness(6,1);
    pwd.setLateralBias(0.0);
    pwd.setOuterIterations(_niter); // default is 5
    pwd.setInnerIterations(20); // default is 20
    float[][] pe = new float[n2][n1]; //pwd w/initial p=0 (Dave)
    _sw.restart();
    pe = pwd.findSlopes(f);
    pwd.updateSlopes(f,pe);
    pe = mul(pe,(float)(d1/d2));
    _sw.stop();
    trace("Dave's PWD time = "+_sw.time());    

    // interp, title, paint, colorbar, color
    String hl = "Traces"; //horizontal label
    String vl = "Samples"; //vertical label
    String cbl = "slope (samples/trace)"; //colorbar label
    Plot.plot(_s1,_s2,pe,"PWDD noise= "+_noise,hl,vl,cbl,
        _fw,_fh,-_clipMax,_clipMax,
        _clip,F,_title,_paint,T,T);
  }

  /**
   * PWD Dave: plots the estimated slopes and RMS error.
   */
  public void plotPWDD(float[][] f, float[][] p) {
    int n1 = _s1.getCount();
    int n2 = _s2.getCount();
    double d1 = _s1.getDelta();
    double d2 = _s2.getDelta();
    PlaneWaveDestructor pwd = new PlaneWaveDestructor(-_pmax,_pmax);
    pwd.setSmoothness(6,1);
    pwd.setLateralBias(0.0);
    pwd.setOuterIterations(_niter); // default is 5
    pwd.setInnerIterations(20); // default is 20
    float[][] pe = new float[n2][n1]; //pwd w/initial p=0 (Dave)
    _sw.restart();
    pe = pwd.findSlopes(f);
    pwd.updateSlopes(f,pe);
    pe = mul(pe,(float)(d1/d2));
    _sw.stop();
    trace("Dave's PWD time = "+_sw.time());    
    trace("Dave's PWD:");
    float error_pwdd = Util.rmsError(pe,p,T);

    // interp, title, paint, colorbar, color
    String hl = "Traces"; //horizontal label
    String vl = "Samples"; //vertical label
    String cbl = "slope (samples/trace)"; //colorbar label
    Plot.plot(_s1,_s2,pe,"PWDD noise= "+_noise,hl,vl,cbl,
        _fw,_fh,-_clipMax,_clipMax,
        _clip,F,_title,_paint,T,T);
  }

  /**
   * Dynamic warping: plots the estimated slopes and RMS error.
   */
  public void plotDW(int k, float[][] f) {
    int n1 = _s1.getCount();
    int n2 = _s2.getCount();
    double d1 = _s1.getDelta();
    double d2 = _s2.getDelta();
    float strainMax1 = 0.2f;
    float strainMax2 = 1.0f;
    Sampling ss1 = new Sampling((n1-1)*k+1);
    DynamicWarpingSlopes dws = new DynamicWarpingSlopes((int)_pmax*k,ss1,_s2);
    dws.setK(k);
    _sw.restart();
    float[][] pe = dws.findSlopes(f);
    pe = mul(pe,(float)(d1/d2));
    _sw.stop();
    trace("Dynamic warping time = "+_sw.time());    

    // interp, title, paint, colorbar, color
    String hl = "Traces"; //horizontal label
    String vl = "Samples"; //vertical label
    String cbl = "slope (samples/trace)"; //colorbar label
    Plot.plot(_s1,_s2,pe,"DW noise= "+_noise,hl,vl,cbl,
        _fw,_fh,-_clipMax,_clipMax,
        _clip,F,_title,_paint,T,T);
  }

  /**
   * Dynamic warping: plots the estimated slopes and RMS error.
   */
  public void plotDW(int k, float[][] f, float[][] p) {
    int n1 = _s1.getCount();
    int n2 = _s2.getCount();
    double d1 = _s1.getDelta();
    double d2 = _s2.getDelta();
    float strainMax1 = 0.2f;
    float strainMax2 = 1.0f;
    Sampling ss1 = new Sampling((n1-1)*k+1);
    DynamicWarpingSlopes dws = new DynamicWarpingSlopes((int)_pmax*k,ss1,_s2);
    dws.setK(k);
    _sw.restart();
    float[][] pe = dws.findSlopes(f);
    pe = mul(pe,(float)(d1/d2));
    _sw.stop();
    trace("Dynamic warping time = "+_sw.time());    

    trace("Dynamic warping:");
    float error_dw = Util.rmsError(pe,p,T);

    // interp, title, paint, colorbar, color
    String hl = "Traces"; //horizontal label
    String vl = "Samples"; //vertical label
    String cbl = "slope (samples/trace)"; //colorbar label
    Plot.plot(_s1,_s2,pe,"DW noise= "+_noise,hl,vl,cbl,
        _fw,_fh,-_clipMax,_clipMax,
        _clip,F,_title,_paint,T,T);
  }

  /**
   * Smooth dynamic warping: plots the estimated slopes and RMS error.
   */
  public void plotSDW(int k, float[][] f) {
    int n1 = _s1.getCount();
    int n2 = _s2.getCount();
    double d1 = _s1.getDelta();
    double d2 = _s2.getDelta();
    double r1 = 0.1;
    double r2 = 0.4;
    double h1 = 72.0;
    double h2 = 12.0;
    DynamicWarpingSlopes dws = new DynamicWarpingSlopes(k,_pmax,h1,h2,
                                       -r1,r1,-r2,r2,_s1,_s2);
    _sw.restart();
    float[][] pe = dws.findSmoothSlopes(f);
    pe = mul(pe,(float)(d1/d2));
    _sw.stop();
    trace("Smooth dynamic warping time = "+_sw.time());    

    // interp, title, paint, colorbar, color
    String hl = "Traces"; //horizontal label
    String vl = "Samples"; //vertical label
    String cbl = "slope (samples/trace)"; //colorbar label
    Plot.plot(_s1,_s2,pe,"SDW noise= "+_noise,hl,vl,cbl,
        _fw,_fh,-_clipMax,_clipMax,
        _clip,F,_title,_paint,T,T);
  }

  /**
   * Smooth dynamic warping: plots the estimated slopes and RMS error.
   */
  public void plotSDW(int k, float[][] f, float[][] p) {
    int n1 = _s1.getCount();
    int n2 = _s2.getCount();
    double d1 = _s1.getDelta();
    double d2 = _s2.getDelta();
    double r1 = 0.1;
    double r2 = 0.4;
    double h1 = 72.0;
    double h2 = 12.0;
    DynamicWarpingSlopes dws = new DynamicWarpingSlopes(k,_pmax,h1,h2,
                                       -r1,r1,-r2,r2,_s1,_s2);
    _sw.restart();
    float[][] pe = dws.findSmoothSlopes(f);
    pe = mul(pe,(float)(d1/d2));
    _sw.stop();
    trace("Smooth dynamic warping time = "+_sw.time());    
    float error_sdw = Util.rmsError(pe,p,T);

    // interp, title, paint, colorbar, color
    String hl = "Traces"; //horizontal label
    String vl = "Samples"; //vertical label
    String cbl = "slope (samples/trace)"; //colorbar label
    Plot.plot(_s1,_s2,pe,"SDW noise= "+_noise,hl,vl,cbl,
        _fw,_fh,-_clipMax,_clipMax,
        _clip,F,_title,_paint,T,T);
  }

  /**
   * Plots the synthetic seismic image.
   */
  public void plotF(float[][] f) {
    // interp, title, paint, colorbar, color
    String hl = "Traces"; //horizontal label
    String vl = "Samples"; //vertical label
    String cbl = ""; //colorbar label
    Plot.plot(_s1,_s2,f,"Synthetic noise= "+_noise,hl,vl,cbl,
        _fw,_fh,0,0,
        F,F,_title,_paint,T,F);
  }

  /**
   * Plots the known slopes image.
   */
  public void plotP(float[][] p) {
    trace("max slope= "+max(p));
    // interp, title, paint, colorbar, color
    String hl = "Traces"; //horizontal label
    String vl = "Samples"; //vertical label
    String cbl = "slope (samples/trace)"; //colorbar label
    Plot.plot(_s1,_s2,p,"Known slopes",hl,vl,cbl,
        _fw,_fh,-_clipMax,_clipMax,
        _clip,F,_title,_paint,T,T);
  }

////////////////////TESTING/PLOTTING OPTIMAL PARAMETERS///////////////////

  /**
   * Tests for the parameters that produce the smallest RMS error.
   */
  public void testOptimalSmoothLSF(String method, float[][] f, float[][] p, 
                                   Sampling sp1, Sampling sp2) {
    int n1  = _s1.getCount();
    int n2  = _s2.getCount();
    int np1 = sp1.getCount();
    int np2 = sp2.getCount();
    double d1 = _s1.getDelta();
    double d2 = _s2.getDelta();
    float[] param1 = Util.f(sp1.getValues());
    float[] param2 = Util.f(sp2.getValues());

    float[][] pe = new float[n2][n1];
    float[][] rmserror = new float[np2][np1]; // RMS error
    LocalSlopeFinder lsf;
    for(int i2=0; i2<np2; ++i2) {
      for(int i1=0; i1<np1; ++i1) {
        lsf = new LocalSlopeFinder(param1[i1],param2[i2],_pmax);
        lsf.findSlopes(f,pe);
        pe = mul(pe,(float)(d1/d2));
        rmserror[i2][i1] = Util.rmsError(pe,p,F);
      }
    }
    Util.writeBinary(rmserror,PATH+method+"_errors.dat");
  }

  /**
   * Tests for the parameters that produce the smallest RMS error.
   */
  public void testOptimalSmoothPWDM(String method, float[][] f, float[][] p, 
                                    Sampling sp1, Sampling sp2) {
    int n1  = _s1.getCount();
    int n2  = _s2.getCount();
    int np1 = sp1.getCount();
    int np2 = sp2.getCount();
    double d1 = _s1.getDelta();
    double d2 = _s2.getDelta();
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
        pe = mul(pe,(float)(d1/d2));
        rmserror[i2][i1] = Util.rmsError(pe,p,F);
      }
    }
    Util.writeBinary(rmserror,PATH+method+"_errors.dat");
  }

  /**
   * Tests for the parameters that produce the smallest RMS error.
   */
  public void testOptimalSmoothPWDD(String method, float[][] f, float[][] p, 
                                    Sampling sp1, Sampling sp2) {
    int n1 = _s1.getCount();
    int n2 = _s2.getCount();
    int np1 = sp1.getCount();
    int np2 = sp2.getCount();
    double d1 = _s1.getDelta();
    double d2 = _s2.getDelta();
    float[] param1 = Util.f(sp1.getValues());
    float[] param2 = Util.f(sp2.getValues());

    float[][] pe = new float[n2][n1];
    float[][] rmserror = new float[np2][np1]; // RMS error
    PlaneWaveDestructor pwd = new PlaneWaveDestructor(-_pmax,_pmax);
    for(int i2=0; i2<np2; ++i2) {
      for(int i1=0; i1<np1; ++i1) {
        pwd.setSmoothness(param1[i1],param2[i2]);
        pe = pwd.findSlopes(f);
        pe = mul(pe,(float)(d1/d2));
        rmserror[i2][i1] = Util.rmsError(pe,p,F);
      }
    }
    Util.writeBinary(rmserror,PATH+method+"_errors.dat");
  }

  /**
   * Tests for the parameters that produce the smallest RMS error.
   */
  public void testOptimalSmoothSDW(String method, int k, 
              float[][] f, float[][] p, Sampling sp1, Sampling sp2) {
    int n1  = _s1.getCount();
    int n2  = _s2.getCount();
    int np1 = sp1.getCount();
    int np2 = sp2.getCount();
    double d1 = _s1.getDelta();
    double d2 = _s2.getDelta();
    double[] param1 = sp1.getValues();
    double[] param2 = sp2.getValues();
    double r1 = 0.1;
    double r2 = 0.4;

    float[][] pe = new float[n2][n1];
    float[][] rmserror = new float[np2][np1]; // RMS error
    DynamicWarpingSlopes dws;
    for(int i2=0; i2<np2; ++i2) {
      for(int i1=0; i1<np1; ++i1) {
        dws = new DynamicWarpingSlopes(k,_pmax,param1[i1],param2[i2],
                                       -r1,r1,-r2,r2,_s1,_s2);
        pe = dws.findSmoothSlopes(f);
        pe = mul(pe,(float)(d1/d2));
        rmserror[i2][i1] = Util.rmsError(pe,p,F);
      }
    }
    Util.writeBinary(rmserror,PATH+method+"_errors.dat");
  }

  /**
   * Tests for the parameters that produce the smallest RMS error.
   */
  public void testOptimalStrainSDW(String method, int k, 
              float[][] f, float[][] p, Sampling sp1, Sampling sp2) {
    int n1  = _s1.getCount();
    int n2  = _s2.getCount();
    int np1 = sp1.getCount();
    int np2 = sp2.getCount();
    double d1 = _s1.getDelta();
    double d2 = _s2.getDelta();
    float[] param1 = Util.f(sp1.getValues());
    float[] param2 = Util.f(sp2.getValues());
    double h1 = 72.0;
    double h2 = 12.0;

    float[][] pe = new float[n2][n1];
    float[][] rmserror = new float[np2][np1]; // RMS error
    DynamicWarpingSlopes dws;
    for(int i2=0; i2<np2; ++i2) {
      for(int i1=0; i1<np1; ++i1) {
        dws = new DynamicWarpingSlopes(k,_pmax,h1,h2,
                                       -param1[i1],param1[i1],
                                       -param2[i2],param2[i2],_s1,_s2);
        pe = dws.findSmoothSlopes(f);
        pe = mul(pe,(float)(d1/d2));
        rmserror[i2][i1] = Util.rmsError(pe,p,F);
      }
    }
    Util.writeBinary(rmserror,PATH+method+"_errors.dat");
  }

  public void plotOptimalParameters(String method, String hl, String vl,
                                    Sampling sp1, Sampling sp2) {
    int np1 = sp1.getCount();
    int np2 = sp2.getCount();
    float[] param1 = Util.f(sp1.getValues());
    float[] param2 = Util.f(sp2.getValues());
    int[] error_index = new int[2];    
    float[][] rmserror = Util.readImage(np1,np2,
        PATH+method+"_errors.dat");
    float min_error = min(rmserror,error_index);
    trace(vl+"= "+param1[error_index[0]]+" "+hl+"= "+param2[error_index[1]]+" "+
             " Minimum Error Value= "+min_error);    

    // clip, interp, title, paint, colorbar, color
    String cbl = "rms error (samples/trace)"; //colorbar label
    float cmin = 0.4f;
    float cmax = 0.75f;
    Plot.plot(sp1,sp2,rmserror,"RMS Error "+method,hl,vl,cbl,
        _fw,_fh,cmin,cmax,
        F,T,_title,_paint,T,T);
  }

  //////////////////OTHER METHODS FOR ALGORITHM EVALUATION///////////////

  public void testOrderVsTime(int norder, float[][] f, float[][] p) {
    int n1 = _s1.getCount();
    int n2 = _s2.getCount();
    double d1 = _s1.getDelta();
    double d2 = _s2.getDelta();
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
    Util.writeBinary(ovt,PATH+"orderVsTime.dat");
  }

  public void plotOrderVsTime(int norder) {
    float[] ovt = Util.readImage(norder,PATH+"orderVsTime.dat");
    // paint 
    Plot.plot(ovt,"Order vs time","Order","Time",_fw,_fh,F);
  }

  public void testRmsErrorCurveLSF(Sampling sn) {
    int n1 = _s1.getCount();
    int n2 = _s2.getCount();
    double d1 = _s1.getDelta();
    double d2 = _s2.getDelta();
    int nrms = sn.getCount();
    double drms = sn.getDelta();
    float[] noise = Util.f(sn.getValues());;

    float[][][] fandp = new float[2][n2][n1];
    float[][] f = new float[n2][n1];
    float[][] p = new float[n2][n1];

    //Structure tensor
    LocalSlopeFinder lsf = new LocalSlopeFinder(23.0f,1.0f,_pmax);
    float[][] p_lsf = new float[n2][n1];
    
    float[] rms_error_lsf = new float[nrms];
    for (int i=0; i<nrms; ++i) {
      fandp = FakeData.seismicAndSlopes2d2014A(noise[i],F);
      f = fandp[0];
      p = fandp[1];
      lsf.findSlopes(f,p_lsf);
      p_lsf = mul(p_lsf,(float)(d1/d2));
      rms_error_lsf[i] = Util.rmsError(p_lsf,p,false);
    }
    Util.writeBinary(rms_error_lsf,PATH+"rms_error_lsf.dat");
  }

  public void testRmsErrorCurvePWD(Sampling sn) {
    int n1 = _s1.getCount();
    int n2 = _s2.getCount();
    double d1 = _s1.getDelta();
    double d2 = _s2.getDelta();
    int nrms = sn.getCount();
    double drms = sn.getDelta();
    float[] noise = Util.f(sn.getValues());;

    float[][][] fandp = new float[2][n2][n1];
    float[][] f = new float[n2][n1];
    float[][] p = new float[n2][n1];

    //Plane-wave destruction filter
    Sfdip sd = new Sfdip(-_pmax,_pmax);
    sd.setRect(75,6);
    sd.setOrder(4);
    float[][] p_pwd = new float[n2][n1];

    float[] rms_error_pwd = new float[nrms];
    for (int i=0; i<nrms; ++i) {
      fandp = FakeData.seismicAndSlopes2d2014A(noise[i],F);
      f = fandp[0];
      p = fandp[1];
      sd.findSlopes(_s1,_s2,f,p_pwd);
      p_pwd = mul(p_pwd,(float)(d1/d2));
      rms_error_pwd[i] = Util.rmsError(p_pwd,p,false);
    }
    Util.writeBinary(rms_error_pwd,PATH+"rms_error_pwd.dat");
  }

  public void testRmsErrorCurveSDW(Sampling sn, int k) {
    int n1 = _s1.getCount();
    int n2 = _s2.getCount();
    double d1 = _s1.getDelta();
    double d2 = _s2.getDelta();
    int nrms = sn.getCount();
    double drms = sn.getDelta();
    float[] noise = Util.f(sn.getValues());;

    float[][][] fandp = new float[2][n2][n1];
    float[][] f = new float[n2][n1];
    float[][] p = new float[n2][n1];

    //Smooth dynamic warping
    double r1 = 0.1;
    double r2 = 0.4;
    double h1 = 40.0;
    double h2 =  9.0;
    DynamicWarpingSlopes dws = new DynamicWarpingSlopes(k,_pmax,h1,h2,
                                       -r1,r1,-r2,r2,_s1,_s2);
    float[][] p_sdw = new float[n2][n1];
    
    float[] rms_error_sdw = new float[nrms];
    for (int i=0; i<nrms; ++i) {
      fandp = FakeData.seismicAndSlopes2d2014A(noise[i],F);
      f = fandp[0];
      p = fandp[1];
      p_sdw = dws.findSmoothSlopes(f);
      p_sdw = mul(p_sdw,(float)(d1/d2));
      rms_error_sdw[i] = Util.rmsError(p_sdw,p,false);
    }
    Util.writeBinary(rms_error_sdw,PATH+"rms_error_sdw.dat");
  }

  public void plotRmsErrorCurves(Sampling sn) {
    float fw = 0.75f; //fraction width for slide
    float fh = 0.9f; //fraction height for slide
    int nrms = sn.getCount();
    float[] noise = Util.f(sn.getValues());;
    float[] rms_error_lsf = Util.readImage(nrms,PATH+"rms_error_lsf.dat");
    float[] rms_error_pwd = Util.readImage(nrms,PATH+"rms_error_pwd.dat");
    float[] rms_error_sdw = Util.readImage(nrms,PATH+"rms_error_sdw.dat");
    // paint
    Plot.plot(noise,rms_error_pwd,rms_error_lsf,rms_error_sdw,
        "PWD, LSF, and SDW RMS error vs NS ratio",fw,fh,F);
  }

  /*private static void goErrorLocation() {
    int i = 5;
    LocalSlopeFinder lsf = new LocalSlopeFinder(24.0f,1.0f,pmax);
    float[][] lsf_slope = new float[n2][n1];
    lsf.findSlopes(synth_data,lsf_slope);

    Sfdip sd = new Sfdip(-pmax,pmax);
    sd.setRect(76,6);
    sd.setOrder(2);
    float[][] mad_slope = new float[n2][n1];
    sd.findSlopes(s1,s2,synth_data,mad_slope);

    lsf_slope = mul(lsf_slope,d1/d2);
    mad_slope = mul(mad_slope,d1/d2);
    
    float[][] exact_grad = new float[n2][n1];
    float[][] lsf_grad = new float[n2][n1];
    float[][] mad_grad = new float[n2][n1];
    double rgf_sig = 2.0;
    RecursiveGaussianFilter rgf = new RecursiveGaussianFilter(rgf_sig);
    rgf.apply11(exact_slope,exact_grad);
    rgf.apply11(lsf_slope,lsf_grad);
    rgf.apply11(mad_slope,mad_grad);

    float exact_mag = sqrt(sum(pow(exact_grad,2))); 
    float lsf_mag = sqrt(sum(pow(lsf_grad,2))); 
    float mad_mag = sqrt(sum(pow(mad_grad,2))); 

    float[][] exact_rmsdiff = new float[n2][n1]; 
    float[][] lsf_rmsdiff = new float[n2][n1]; 
    float[][] mad_rmsdiff = new float[n2][n1]; 

    sub(lsf_slope,exact_slope,lsf_rmsdiff);
    sub(mad_slope,exact_slope,mad_rmsdiff);
    

    float fw = 0.9f; //fraction width for slide
    float fh = 0.7f; //fraction height for slide
    // title, paint, colorbar, color
    //Plot.plot(s1,s2,exact_rmsdiff,"RMS Diff Exact",fw,fh,F,F,T,T);
    Plot.plot(s1,s2,lsf_rmsdiff,"RMS Diff LOF w noise"+i,fw,fh,-2,2,F,F,T,T);
    Plot.plot(s1,s2,mad_rmsdiff,"RMS Diff MAD w noise"+i,fw,fh,-2,2,F,F,T,T);
  }*/

  /*private static void goTestSampleMeanSD() {
    float[][][] synthAndSlope = new float[2][n2][n1];
    float[][] synth_data = new float[n2][n1];
    float[][] exact_slope = new float[n2][n1];

    int N = 100;
    float[][] lsf_sum = new float[n2][n1];
    float[][] mad_sum = new float[n2][n1];
    float[][] dw_sum = new float[n2][n1];
    float[][][] lsf_p = new float[N][n2][n1];
    float[][][] mad_p = new float[N][n2][n1];
    float[][][] dw_p = new float[N][n2][n1];
    for (int i=0; i<N; ++i) {
      synthAndSlope = FakeData.seismicAndSlopes2d2014B(noise,T);
      synth_data = synthAndSlope[0];
      System.out.println("i= "+i);

      LocalSlopeFinder lsf = new LocalSlopeFinder(24.0f,1.0f,pmax);
      float[][] lsf_slope = new float[n2][n1];
      lsf.findSlopes(synth_data,lsf_slope);

      Sfdip sd = new Sfdip(-pmax,pmax);
      sd.setRect(76,6);
      sd.setOrder(2);
      float[][] mad_slope = new float[n2][n1];
      sd.findSlopes(s1,s2,synth_data,mad_slope);

      DynamicWarping dw = new DynamicWarping((int)-pmax,(int)pmax);
      dw.setShiftSmoothing(15.0,1.0);
      dw.setStrainMax(0.3,1.0); //only allow a max 50% stretch or squeeze
      float[][] dw_slope = new float[n2][n1];
      dw_slope = Slopes.DWSlopesAvg(dw,synth_data);

      lsf_slope = mul(lsf_slope,d1/d2);
      mad_slope = mul(mad_slope,d1/d2);
      dw_slope = mul(dw_slope,d1/d2);

      lsf_p[i] = lsf_slope;
      mad_p[i] = mad_slope;
      dw_p[i] = dw_slope;
      
      lsf_sum = add(lsf_sum,lsf_slope);
      mad_sum = add(mad_sum,mad_slope);
      dw_sum = add(dw_sum,dw_slope);
    }
    exact_slope = synthAndSlope[1];

    float[][] lsf_mean = div(lsf_sum,N);
    float[][] mad_mean = div(mad_sum,N);
    float[][] dw_mean = div(dw_sum,N);

    float[][] lsf_sumdiffsq = new float[n2][n1]; //sum of differences squared
    float[][] mad_sumdiffsq = new float[n2][n1]; //sum of differences squared
    float[][] dw_sumdiffsq = new float[n2][n1]; //sum of differences squared
    for (int i=0; i<N; ++i) {
      lsf_sumdiffsq = add(lsf_sumdiffsq,pow(sub(lsf_p[i],exact_slope),2));
      mad_sumdiffsq = add(mad_sumdiffsq,pow(sub(mad_p[i],exact_slope),2));
      dw_sumdiffsq = add(dw_sumdiffsq,pow(sub(dw_p[i],exact_slope),2));
      System.out.println("i= "+i);
    }
    float[][] lsf_sd= sqrt(div(lsf_sumdiffsq,N));
    float[][] mad_sd= sqrt(div(mad_sumdiffsq,N));
    float[][] dw_sd= sqrt(div(dw_sumdiffsq,N));
    Util.writeBinary(lsf_mean,
          "/Users/earias/Home/git/ea/bench/src/slopes/data/lsf_mean.dat");
    Util.writeBinary(mad_mean,
          "/Users/earias/Home/git/ea/bench/src/slopes/data/mad_mean.dat");
    Util.writeBinary(dw_mean,
          "/Users/earias/Home/git/ea/bench/src/slopes/data/dw_mean.dat");
    Util.writeBinary(lsf_sd,
          "/Users/earias/Home/git/ea/bench/src/slopes/data/lsf_sd.dat");
    Util.writeBinary(mad_sd,
          "/Users/earias/Home/git/ea/bench/src/slopes/data/mad_sd.dat");
    Util.writeBinary(dw_sd,
          "/Users/earias/Home/git/ea/bench/src/slopes/data/dw_sd.dat");
  }*/

  /*private static void goPlotSampleMeanSD() {
    float[][] lsf_mean = Util.readImage(n1,n2,
        "/Users/earias/Home/git/ea/bench/src/slopes/data/lsf_mean.dat");
    float[][] mad_mean = Util.readImage(n1,n2,
        "/Users/earias/Home/git/ea/bench/src/slopes/data/mad_mean.dat");
    float[][] dw_mean = Util.readImage(n1,n2,
        "/Users/earias/Home/git/ea/bench/src/slopes/data/dw_mean.dat");
    float[][] lsf_sd = Util.readImage(n1,n2,
        "/Users/earias/Home/git/ea/bench/src/slopes/data/lsf_sd.dat");
    float[][] mad_sd = Util.readImage(n1,n2,
        "/Users/earias/Home/git/ea/bench/src/slopes/data/mad_sd.dat");
    float[][] dw_sd = Util.readImage(n1,n2,
        "/Users/earias/Home/git/ea/bench/src/slopes/data/dw_sd.dat");
    // title, paint, colorbar, color
    //Plot.plot(s1,s2,lsf_mean,"Sample mean LOF",fw,fh,-4,4,F,T,T,T);
    //Plot.plot(s1,s2,mad_mean,"Sample mean PWD",fw,fh,-4,4,F,T,T,T);
    //Plot.plot(s1,s2,lsf_sd,"Sample variance LOF 03 max",fw,fh,0,0.3f,F,F,T,T);
    //Plot.plot(s1,s2,mad_sd,"Sample variance PWD",fw,fh,0,0.5f,F,F,T,T);
    Plot.plot(s1,s2,lsf_sd,"Sample standard deviation LOF",fw,fh,0,1,F,T,T,T);
    Plot.plot(s1,s2,mad_sd,"Sample standard deviation PWD",fw,fh,0,1,F,T,T,T);
    Plot.plot(s1,s2,dw_sd,"Sample standard deviation DW",fw,fh,0,1,F,T,T,T);
  }*/

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
    "/Users/earias/Home/git/ea/bench/src/util/";
  private static final int _niter = 5;
  private static final float pi = FLT_PI;      
  private static final float _fw = 0.75f; //fraction width for slide
  private static final float _fh = 0.9f; //fraction height for slide
  private static final float _clipMax = 4.0f;
  private static final boolean T = true;
  private static final boolean F = false;  
  private static final boolean _title = true;
  private static final boolean _paint = false;  
  private static final boolean _clip = true;  

  private Stopwatch _sw = new Stopwatch();
  private float _noise;
  private float _pmax;
  private Sampling _s1;
  private Sampling _s2;
  private Sampling _s3;
}
