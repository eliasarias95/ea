package slopes;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.Check;
import edu.mines.jtk.util.Stopwatch;

import static edu.mines.jtk.util.ArrayMath.*;

import util.*;
import dnp.LocalSlopeFinder;
import dnp.PlaneWaveDestructor;


import javax.swing.*;

/**
 *  
 * @author Elias Arias, Colorado School of Mines, CWP
 * @version 19.2.2015
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
  public static void makeSyntheticComplex(float noise, float[][] f, 
                                                       float[][] p, 
                                                       float[][] r) {
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
    int n1 = f[0].length;
    int n2 = f.length;
    float[][] temp = Util.sexp(Util.readImage(n1,n2,PATH+"data/gom.dat"));
    //mul(temp,.001f,temp);
    for (int i2=0; i2<n2; ++i2){
      for (int i1=0; i1<n1; ++i1){
        f[i2][i1] = temp[i2][i1];
      }
    }
    //Util.sexp(f);
    //float[][] fs = new float[n2][n1];
    //mul(f,.001f,f);
  }

  /**
   * Takes a real seismic image from the GOM and places it in the array f.
   * @param f array[920][301] that holds the data for the GOM seismic image.
   */
  public static void makeRealTp(float[][] f) {
    int n1 = f[0].length;
    int n2 = f.length;
    float[][] temp = (Util.readImage(n1,n2,PATH+"data/tp/tp73.dat"));
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
    float[][][] temp = Util.readImage(n1,n2,n3,PATH+
          //"data/tp/tpsz_subz_51_4_1400.dat");
          "data/tp/tpsz_subz_401_4_400.dat");
    for (int i3=0; i3<n3; ++i3){
      for (int i2=0; i2<n2; ++i2){
        for (int i1=0; i1<n1; ++i1){
          f[i3][i2][i1] = temp[i3][i2][i1];
        }
      }
    }
  }

////////////////////////ESTIMATING SLOPES/////////////////////////////

  /**
   * Structure tensor: plots the estimated slopes.
   */
  public void plotLSF(float[][] f, String title, String hl, String vl) {
    boolean one = true;
    int n1 = _s1.getCount();
    int n2 = _s2.getCount();
    LocalSlopeFinder lsf = new LocalSlopeFinder(23.0f,1.0f,_pmax);
    float[][] pe = new float[n2][n1];
    _sw.restart();
    lsf.findSlopes(f,pe);
    trace("lsf slopes max= "+max(pe));
    trace("lsf slopes min= "+min(pe));
    _sw.stop();
    trace("Structure tensor time = "+_sw.time());    

    Util.writeBinary(pe,PATH+"data/lsf_"+title+".dat");
    // clip, interp, title, paint, colorbar, color
    String cbl = "slope (samples/trace)"; //colorbar label
    Plot.plot(_s1,_s2,pe,"lsf_"+title,hl,vl,cbl,
        _fw,_fh,-_clipMax,_clipMax,
        _clip,F,_title,_paint,T,T,_slide,one);
  }

  /**
   * Structure tensor: plots the estimated slopes, RMS error, and time.
   */
  public void plotLSF(float[][] f, float[][] p) {
    boolean one = true;
    int n1 = _s1.getCount();
    int n2 = _s2.getCount();
    LocalSlopeFinder lsf = new LocalSlopeFinder(23.0f,1.0f,_pmax);
    float[][] pe = new float[n2][n1];
    _sw.restart();
    lsf.findSlopes(f,pe);
    _sw.stop();
    trace("Structure tensor time = "+_sw.time());    
    trace("Structure tensor:");
    float error_lsf = Util.rmsError(pe,p,T);

    // clip, interp, title, paint, colorbar, color
    String hl = "Traces"; //horizontal label
    String vl = "Samples"; //vertical label
    String cbl = "slope (samples/trace)"; //colorbar label
    Plot.plot(_s1,_s2,pe,"lsf_noise="+_noise,hl,vl,cbl,
        _fw,_fh,-_clipMax,_clipMax,
        _clip,F,_title,_paint,T,T,_slide,one);
  }

  /**
   * PWD Madagascar: plots the estimated slopes.
   */
  public void plotPWDM(float[][] f, String title, String hl, String vl) {
    boolean one = true;
    int n1 = _s1.getCount();
    int n2 = _s2.getCount();
    Sfdip sd = new Sfdip(-_pmax,_pmax);
    sd.setRect(75,6);
    sd.setOrder(4);
    sd.setNiter(_niter);
    sd.setNj(1);
    sd.setBoth("y");
    float[][] pe = new float[n2][n1]; //pwd w/initial p=0 (Madagascar)
    _sw.restart();
    sd.findSlopes(_s1,_s2,f,pe);
    trace("pwdm slopes max= "+max(pe));
    trace("pwdm slopes min= "+min(pe));
    _sw.stop();
    trace("Madagascar PWD time = "+_sw.time());    

    Util.writeBinary(pe,PATH+"data/pwdm_"+title+".dat");
    // interp, title, paint, colorbar, color
    String cbl = "slope (samples/trace)"; //colorbar label
    Plot.plot(_s1,_s2,pe,"pwdm_"+title,hl,vl,cbl,
        _fw,_fh,-_clipMax,_clipMax,
        _clip,F,_title,_paint,T,T,_slide,one);
  }

  /**
   * PWD Madagascar: plots the estimated slopes, RMS error, and time.
   */
  public void plotPWDM(float[][] f, float[][] p) {
    boolean one = true;
    int n1 = _s1.getCount();
    int n2 = _s2.getCount();
    Sfdip sd = new Sfdip(-_pmax,_pmax);
    sd.setRect(75,6);
    sd.setOrder(4);
    sd.setNiter(_niter);
    sd.setNj(1);
    sd.setBoth("y");
    float[][] pe = new float[n2][n1]; //pwd w/initial p=0 (Madagascar)
    _sw.restart();
    sd.findSlopes(_s1,_s2,f,pe);
    _sw.stop();
    trace("Madagascar PWD time = "+_sw.time());    
    trace("Madagascar PWD:");
    float error_pwdm = Util.rmsError(pe,p,T);

    // interp, title, paint, colorbar, color
    String hl = "Traces"; //horizontal label
    String vl = "Samples"; //vertical label
    String cbl = "slope (samples/trace)"; //colorbar label
    Plot.plot(_s1,_s2,pe,"pwdm_noise="+_noise,hl,vl,cbl,
        _fw,_fh,-_clipMax,_clipMax,
        _clip,F,_title,_paint,T,T,_slide,one);
  }

  /**
   * PWD Dave: plots the estimated slopes.
   */
  public void plotPWDD(float[][] f, String title, String hl, String vl) {
    boolean one = true;
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
    trace("pwdd slopes max= "+max(pe));
    trace("pwdd slopes min= "+min(pe));
    _sw.stop();
    trace("Dave's PWD time = "+_sw.time());    

    // interp, title, paint, colorbar, color
    String cbl = "slope (samples/trace)"; //colorbar label
    Plot.plot(_s1,_s2,pe,"pwdd_"+title,hl,vl,cbl,
        _fw,_fh,-_clipMax,_clipMax,
        _clip,F,_title,_paint,T,T,_slide,one);
  }

  /**
   * PWD Dave: plots the estimated slopes, RMS error, and time.
   */
  public void plotPWDD(float[][] f, float[][] p) {
    boolean one = true;
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
    Plot.plot(_s1,_s2,pe,"pwdd_noise="+_noise,hl,vl,cbl,
        _fw,_fh,-_clipMax,_clipMax,
        _clip,F,_title,_paint,T,T,_slide,one);
  }

  /**
   * Dynamic warping: plots the estimated slopes.
   */
  public void plotDW(int k, float[][] f, String title, String hl, String vl) {
    boolean one = true;
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
    trace("dw slopes max= "+max(pe));
    trace("dw slopes min= "+min(pe));
    _sw.stop();
    trace("Dynamic warping time = "+_sw.time());    

    // interp, title, paint, colorbar, color
    String cbl = "slope (samples/trace)"; //colorbar label
    Plot.plot(_s1,_s2,pe,"dw_"+title,hl,vl,cbl,
        _fw,_fh,-_clipMax,_clipMax,
        _clip,F,_title,_paint,T,T,_slide,one);
  }

  /**
   * Dynamic warping: plots the estimated slopes, RMS error, and time.
   */
  public void plotDW(int k, float[][] f, float[][] p) {
    boolean one = true;
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
    Plot.plot(_s1,_s2,pe,"dw_noise="+_noise,hl,vl,cbl,
        _fw,_fh,-_clipMax,_clipMax,
        _clip,F,_title,_paint,T,T,_slide,one);
  }

  /**
   * Smooth dynamic warping: plots the estimated slopes.
   */
  public void plotSDW(int k, float[][] f, String title, String hl, String vl) {
    boolean one = true;
    int n1 = _s1.getCount();
    int n2 = _s2.getCount();
    double r1 = 0.1;
    double r2 = 0.5;
    double h1 = 72.0;
    double h2 =  6.0;
    Sampling ss1 = new Sampling(n1);
    Sampling ss2 = new Sampling(n2);
    DynamicWarpingSlopes dws = new DynamicWarpingSlopes(k,_pmax,h1,h2,
                                       r1,r2,ss1,ss2);
    _sw.restart();
    float[][] pe = new float[n2][n1];
    dws.findSmoothSlopes(ss1,f,pe);
    trace("sdw slopes max= "+max(pe));
    trace("sdw slopes min= "+min(pe));
    _sw.stop();
    trace("Smooth dynamic warping time = "+_sw.time());    

    Util.writeBinary(pe,PATH+"data/sdw_"+title+".dat");
    // interp, title, paint, colorbar, color
    String cbl = "slope (samples/trace)"; //colorbar label
    Plot.plot(_s1,_s2,pe,"sdw_"+title,hl,vl,cbl,
        _fw,_fh,-_clipMax,_clipMax,
        _clip,F,_title,_paint,T,T,_slide,one);
  }

  /**
   * Smooth dynamic warping: plots the estimated slopes, RMS error, and time.
   */
  public void plotSDW(int k, float[][] f, float[][] p) {
    boolean one = true;
    int n1 = _s1.getCount();
    int n2 = _s2.getCount();
    double r1 = 0.1;
    double r2 = 0.3;
    double h1 = 72.0;
    double h2 = 12.0;
    Sampling ss1 = new Sampling(n1);
    Sampling ss2 = new Sampling(n2);
    DynamicWarpingSlopes dws = new DynamicWarpingSlopes(k,_pmax,h1,h2,
                                       r1,r2,ss1,ss2);
    _sw.restart();
    float[][] pe = new float[n2][n1];
    dws.findSmoothSlopes(ss1,f,pe);
    _sw.stop();
    trace("Smooth dynamic warping time = "+_sw.time());    
    float error_sdw = Util.rmsError(pe,p,T);

    // interp, title, paint, colorbar, color
    String hl = "Traces"; //horizontal label
    String vl = "Samples"; //vertical label
    String cbl = "slope (samples/trace)"; //colorbar label
    Plot.plot(_s1,_s2,pe,"sdw_noise="+_noise,hl,vl,cbl,
        _fw,_fh,-_clipMax,_clipMax,
        _clip,F,_title,_paint,T,T,_slide,one);
  }

   /**
   * Smooth dynamic warping: plots the estimated slopes, RMS error, and time.
   */
  public void plotSDW(int k, float[][][] f, String title) {
    boolean one = true;
    int n1 = _s1.getCount();
    int n2 = _s2.getCount();
    int n3 = _s3.getCount();
    double r1 = 0.1;
    double r2 = 0.3;
    double r3 = 0.3;
    double h1 = 72.0;
    double h2 = 12.0;
    double h3 = 12.0;
    Sampling ss1 = new Sampling(n1);
    Sampling ss2 = new Sampling(n2);
    Sampling ss3 = new Sampling(n3);
    DynamicWarpingSlopes dws = new DynamicWarpingSlopes(k,_pmax,h1,h2,h3,
                                       r1,r2,r3,ss1,ss2,ss3);
    _sw.restart();
    float[][][] p2 = new float[n3][n2][n1];
    float[][][] p3 = new float[n3][n2][n1];
    dws.findSmoothSlopes(_s1,f,p2,p3);
    _sw.stop();
    trace("Smooth dynamic warping time = "+_sw.time());    
    Util.writeBinary(p2,PATH+"data/"+title+"_p2.dat");
    Util.writeBinary(p3,PATH+"data/"+title+"_p3.dat");
  }

  /**
   * Plots the seismic image.
   */
  public void plotF(String title, String hl, String vl, float[][] f) {
    boolean one = true;
    // clip, interp, title, paint, colorbar, color
    Plot.plot(_s1,_s2,f,title,hl,vl,"",
        _fw,_fh,-2000,2000,
        F,F,_title,_paint,T,F,_slide,one);
  }

  public void plot3D(String title, float[][][] f) {
    int n1 = _s1.getCount();
    int n2 = _s2.getCount();
    int n3 = _s3.getCount();
    float[][][] p2 = Util.readImage(n1,n2,n3,PATH+"data/"+title+"_p2.dat");
    float[][][] p3 = Util.readImage(n1,n2,n3,PATH+"data/"+title+"_p3.dat");
    Plot.plot(_s1,_s2,_s3,f,p2,"slope (samples/trace)",title+"_p2_slices",
        -_clipMax,_clipMax,_paint);
    Plot.plot(_s1,_s2,_s3,f,p3,"slope (samples/trace)",title+"_p3_slices",
        -_clipMax,_clipMax,_paint);
    Plot.plotp(_s1,_s2,_s3,f,p2,"slope (samples/trace)",title+"_p2_panels",
        -_clipMax,_clipMax,_paint);
    Plot.plotp(_s1,_s2,_s3,f,p3,"slope (samples/trace)",title+"_p3_panels",
        -_clipMax,_clipMax,_paint);
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
    // clip, interp, title, paint, colorbar, color
    Plot.plot(_s1,_s2,f,title,hl,vl,cbl,
        _fw,_fh,cmin,cmax,
        T,F,_title,_paint,T,T,_slide,one);
  }

  /**
   * Plots the known slopes image.
   */
  public void plotP(float[][] p) {
    boolean one = true;
    trace("max slope= "+max(p));
    // clip, interp, title, paint, colorbar, color
    String hl = "Traces"; //horizontal label
    String vl = "Samples"; //vertical label
    String cbl = "slope (samples/trace)"; //colorbar label
    Plot.plot(_s1,_s2,p,"known_slopes",hl,vl,cbl,
        _fw,_fh,-_clipMax,_clipMax,
        _clip,F,_title,_paint,T,T,_slide,one);
  }

  public void plotTeaser(String title, float[][] f) {
    boolean paint = true;
    int n1 = _s1.getCount();
    int n2 = _s2.getCount();
    float[][] pe_lsf = Util.readImage(n1,n2,PATH+"data/lsf_"+title+".dat");
    float[][] pe_pwd = Util.readImage(n1,n2,PATH+"data/pwdm_"+title+".dat");
    float[][] pe_sdw = Util.readImage(n1,n2,PATH+"data/sdw_"+title+".dat");
    Plot.plot(_s1,_s2,f,pe_lsf,pe_pwd,pe_sdw,"teaser",_fw,_fh,_slide,paint);
  }

////////////////////TESTING/PLOTTING OPTIMAL PARAMETERS///////////////////

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
    double d1 = _s1.getDelta();
    double d2 = _s2.getDelta();
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
        pe = mul(pe,(float)(d1/d2));
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
    double d1 = _s1.getDelta();
    double d2 = _s2.getDelta();
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
        pe = mul(pe,(float)(d1/d2));
        rmserror[i2][i1] = Util.rmsError(pe,p,F);
      }
    }
    Util.writeBinary(rmserror,PATH+"data/"+fileName);
  }

  public void plotOptimalParameters(String fileName, String hl, String vl,
                                    Sampling sp1, Sampling sp2, 
                                    float cmin, float cmax) {
    boolean one = true;
    int np1 = sp1.getCount();
    int np2 = sp2.getCount();
    float[] param1 = Util.f(sp1.getValues());
    float[] param2 = Util.f(sp2.getValues());
    int[] error_index = new int[2];    
    float[][] rmserror = Util.readImage(np1,np2,PATH+"data/"+fileName);
    float min_error = min(rmserror,error_index);
    trace(vl+"= "+param1[error_index[0]]+" "+hl+"= "+param2[error_index[1]]+" "+
             " Minimum Error Value= "+min_error);    

    // clip, interp, title, paint, colorbar, color
    String cbl = "rms error (samples/trace)"; //colorbar label
    Plot.plot(sp1,sp2,rmserror,"rms_error_"+
        fileName.split("\\.(?=[^\\.]+$)")[0],hl,vl,cbl,
        _fw,_fh,cmin,cmax,
        T,T,_title,_paint,T,T,_slide,one);
  }

  //////////////////OTHER METHODS FOR ALGORITHM EVALUATION///////////////

  public void testOrderVsTime(Sampling s, String fileName, 
      float[][] f, float[][] p) {
    int norder = s.getCount();
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
    Util.writeBinary(ovt,PATH+fileName);
  }

  public void testRmsErrorCurveLSF(Sampling sn, String fileName, 
      float sigma1, float sigma2) {
    int n1 = _s1.getCount();
    int n2 = _s2.getCount();
    double d1 = _s1.getDelta();
    double d2 = _s2.getDelta();
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
      fandp = FakeData.seismicAndSlopes2d2014A(nsratio[i],F);
      f = fandp[0];
      p = fandp[1];
      lsf.findSlopes(f,pe);
      pe = mul(pe,(float)(d1/d2));
      rms_error[i] = Util.rmsError(pe,p,false);
    }
    Util.writeBinary(rms_error,PATH+fileName);
  }

  public void testRmsErrorCurvePWD(Sampling sn, String fileName, 
      int rect1, int rect2) {
    int n1 = _s1.getCount();
    int n2 = _s2.getCount();
    double d1 = _s1.getDelta();
    double d2 = _s2.getDelta();
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
      fandp = FakeData.seismicAndSlopes2d2014A(nsratio[i],F);
      f = fandp[0];
      p = fandp[1];
      sd.findSlopes(_s1,_s2,f,pe);
      pe = mul(pe,(float)(d1/d2));
      rms_error[i] = Util.rmsError(pe,p,false);
    }
    Util.writeBinary(rms_error,PATH+fileName);
  }

  public void testRmsErrorCurveSDW(Sampling sn, String fileName, 
      int k, double r1, double r2, double h1, double h2) {
    int n1 = _s1.getCount();
    int n2 = _s2.getCount();
    double d1 = _s1.getDelta();
    double d2 = _s2.getDelta();
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
      fandp = FakeData.seismicAndSlopes2d2014A(nsratio[i],F);
      f = fandp[0];
      p = fandp[1];
      dws.findSmoothSlopes(_s1,f,pe);
      pe = mul(pe,(float)(d1/d2));
      rms_error[i] = Util.rmsError(pe,p,false);
    }
    Util.writeBinary(rms_error,PATH+fileName);
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
      //Plot.plot(_s1,_s2,sub(pbar,Util.flip2(p)),"mean error","Traces","Samples",
       //"slope",_fw,_fh,-_clipMax,_clipMax,false,F,T,F,T,T,_slide,one);
      mean_error[i] = sum(sub(pbar,(p)))/(n1*n2);
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
      //Plot.plot(_s1,_s2,sub(pbar,Util.flip2(p)),"mean error","Traces","Samples",
       //"slope",_fw,_fh,-_clipMax,_clipMax,false,F,T,F,T,T,_slide,one);
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
      //Plot.plot(_s1,_s2,sub(pbar,Util.flip2(p)),"mean error","Traces","Samples",
       //"slope",_fw,_fh,-_clipMax,_clipMax,false,F,T,F,T,T,_slide,one);
      mean_error[i] = sum(sub(pbar,(p)))/(n1*n2);
    }
    Util.writeBinary(mean_error,PATH+fileName);
  }

  public float[][] testSampleMeanLSF(String fileName, int n, float noise,
      float sigma1, float sigma2) {
    int n1 = _s1.getCount();
    int n2 = _s2.getCount();
    double d1 = _s1.getDelta();
    double d2 = _s2.getDelta();

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
      pe = mul(pe,(float)(d1/d2));
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
    double d1 = _s1.getDelta();
    double d2 = _s2.getDelta();

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
      pe = mul(pe,(float)(d1/d2));
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
    double d1 = _s1.getDelta();
    double d2 = _s2.getDelta();

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
      pe = mul(pe,(float)(d1/d2));
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
    double d1 = _s1.getDelta();
    double d2 = _s2.getDelta();

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
      pe = mul(pe,(float)(d1/d2));
      sumdiffsq = add(sumdiffsq,pow(sub(pe,p),2));
    }
    float[][] stddev = sqrt(div(sumdiffsq,nni));
    Util.writeBinary(stddev,PATH+fileName);
  }

  public void testSampleStdDevPWD(String fileName, int nni,
      int rect1, int rect2) {
    int n1 = _s1.getCount();
    int n2 = _s2.getCount();
    double d1 = _s1.getDelta();
    double d2 = _s2.getDelta();

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
      pe = mul(pe,(float)(d1/d2));
      sumdiffsq = add(sumdiffsq,pow(sub(pe,p),2));
    }
    float[][] stddev = sqrt(div(sumdiffsq,nni));
    Util.writeBinary(stddev,PATH+fileName);
  }

  public void testSampleStdDevSDW(String fileName, int nni,
      int k, double r1, double r2, double h1, double h2) {
    int n1 = _s1.getCount();
    int n2 = _s2.getCount();
    double d1 = _s1.getDelta();
    double d2 = _s2.getDelta();
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
      pe = mul(pe,(float)(d1/d2));
      sumdiffsq = add(sumdiffsq,pow(sub(pe,p),2));
    }
    float[][] stddev = sqrt(div(sumdiffsq,nni));
    Util.writeBinary(stddev,PATH+fileName);
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

////////////////////PLOTTING///////////////////////////

  public void plotCurve(Sampling s, String fileName, String title, 
                                    String hl, String vl, 
                                    float cmin, float cmax) {
    plotCurve(s,fileName,fileName,title,hl,vl,cmin,cmax);
  }

  public void plotCurve(Sampling s, String fileName1, String fileName2, 
                                    String title, String hl, String vl,
                                    float cmin, float cmax) {
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
  private static final float _clipMax = 1.5f;
  private static final boolean T = true;
  private static final boolean F = false;  
  private static final boolean _title = false;
  private static final boolean _paint = false;  
  private static final boolean _clip = true;
  private static final boolean _slide = false;

  private Stopwatch _sw = new Stopwatch();
  private float _noise;
  private int num;
  private float _pmax;
  private Sampling _s1;
  private Sampling _s2;
  private Sampling _s3;
}
