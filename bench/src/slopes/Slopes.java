package slopes;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.Check;
import edu.mines.jtk.util.Stopwatch;

import static edu.mines.jtk.util.ArrayMath.*;

import dnp.LocalSlopeFinder;
import dnp.PlaneWaveDestructor;
import util.FakeData;
import util.Util;
import util.Plot;
import warp.DynamicWarpingR;


import javax.swing.*;

/**
 *  
 * @author Elias Arias, Colorado School of Mines, CWP
 * @version 27.1.2015
 */

/**
 * float[][][] makeConstantSlope(float freq, float p2, int n1, int n2)
 * float[][][][] makeConstantSlope(float freq, float p2, float p3, 
 *                                 int n1, int n2, int n3)
 *
 * float[][][] superSample(int k, float[][][] f)
 */

public class Slopes{

  public static void setChickenTestParameters(float pk) {
    float number = 01.0f;
    _ng = (int)(501*number);
    _n1 = 501;
    _n2 = 501;
    _dg = 1.0f/number;
    _d1 = 1.0f;
    _d2 = 1.0f;
    _fg = 0.0f;
    _f1 = 0.0f;
    _f2 = 0.0f;
    _s1 = new Sampling(_n1,_d1,_f1);
    _s2 = new Sampling(_n2,_d2,_f2);
    _sg = new Sampling(_ng,_dg,_fg);
    float w = 0.1f; //frequency (cycles/sample)
    _f = new float[_n2][_n1];
    _pk = new float[_n2][_n1];
    for (int i2=0; i2<_n2; ++i2){
      for (int i1=0; i1<_n1; ++i1){
        _f[i2][i1] = cos(2.0f*pi*(pk*w*i2-w*i1));
        _pk[i2][i1] = pk;
      }
    }
    for (int i=0; i<_n2; ++i) {
      _f[i] = Util.reSampleSinc(_s1,_f[i],_sg);
      _pk[i] = Util.reSampleSinc(_s1,_pk[i],_sg);
    }
  }

  public static void setSynthParameters(float noise) {
    float number = 1.0f;
    _ng = (int)(501*number);
    _n1 = 501;
    _n2 = 501;
    _dg = 1.0f/number;
    _d1 = 1.0f;
    _d2 = 1.0f;
    _fg = 0.0f;
    _f1 = 0.0f;
    _f2 = 0.0f;
    _s1 = new Sampling(_n1,_d1,_f1);
    _s2 = new Sampling(_n2,_d2,_f2);
    _sg = new Sampling(_ng,_dg,_fg);
    _noise = noise;
    float[][][] fandpk = FakeData.seismicAndSlopes2d2014B(_noise,F);
    _f = fandpk[0];  //synthetic seismic data
    _pk = fandpk[1]; //exact slope values
    for (int i=0; i<_n2; ++i) {
      _f[i] = Util.reSampleSinc(_s1,_f[i],_sg);
      _pk[i] = Util.reSampleSinc(_s1,_pk[i],_sg);
    }
  }

  public static void setGOMParameters() {
    float number = 1.0f;
    _ng = (int)(501*number);
    _n1 = 301;
    _n2 = 920;
    _dg = 1.0f/number;
    _d1 = 0.004f;
    _d2 = .02667f;
    _fg = 0.0f;
    _f1 = 1.6f;
    _f2 = 0;
    _s1 = new Sampling(_n1,_d1,_f1);
    _s2 = new Sampling(_n2,_d2,_f2);
    _sg = new Sampling(_ng,_dg,_fg);
    _f = Util.readImage(_n1,_n2,"data/gom.dat");
    float[][] fs = new float[_n2][_n1];
    mul(_f,.001f,fs);
    Util.writeBinary(fs,"data/gom_scaled.dat");
  }

  public static void setSmoothingParameters(float l_param) {
    _dparam= 1.0f;
    _fparam = 1.0f;
    _n1param = (int)((l_param-_fparam+1)/_dparam);
    _n2param = 10;
    _s1param = new Sampling(_n1param,_dparam,_fparam);
    _s2param = new Sampling(_n2param,_dparam,_fparam);
    _param1 = new float[_n1param];
    _param2 = new float[_n2param];
    for(int i=0; i<_n1param; ++i) {
      _param1[i] = i*_dparam+_fparam;
    }
    for(int i=0; i<_n2param; ++i) {
      _param2[i] = i*_dparam+_fparam;
    }
  }

  public static void setStrainParameters() {
    _dparam = 0.1f; // strain sampling rate
    _fparam = 0.1f;
    _n1param = (int)((1.0f-_fparam+0.1f)/_dparam); //# strain vals to test
    _n2param = _n1param;
    _s1param = new Sampling(_n1param,_dparam,_fparam);
    _s2param = new Sampling(_n2param,_dparam,_fparam);
    _param1 = new float[_n1param];
    _param2 = new float[_n2param];
    for(int i=0; i<_n1param; ++i) {
      _param1[i] = i*_dparam+_fparam;
      _param2[i] = i*_dparam+_fparam;
    }
  }

////////////////////////ESTIMATING SLOPES/////////////////////////////

  /**
   * Structure tensor: plots the estimated slopes and RMS error.
   */
  public static void plotLSF(boolean error) {
    LocalSlopeFinder lsf = new LocalSlopeFinder(24.0f,1.0f,_pmax);
    float[][] p_lsf = new float[_n2][_n1];
    lsf.findSlopes(_f,p_lsf);
    p_lsf = mul(p_lsf,_dg/_d2);

    if (error) {
      System.out.println("Structure tensor:");
      float error_lsf = Util.rmsError(p_lsf,_pk,_dg,_d2,T);
    }

    // clip, interp, title, paint, colorbar, color
    String hl = "Traces"; //horizontal label
    String vl = "Samples"; //vertical label
    String cbl = "slope (samples/trace)"; //colorbar label
    Plot.plot(_sg,_s2,p_lsf,"LSF noise= "+_noise,hl,vl,cbl,
        _fw,_fh,-_clipMax,_clipMax,
        _clip,F,_title,_paint,T,T);
  }

  /**
   * PWD Madagascar: plots the estimated slopes and RMS error.
   */
  public static void plotPWDM(boolean error) {
    Sfdip sd = new Sfdip(-_pmax,_pmax);
    sd.setRect(76,6);
    sd.setOrder(5);
    sd.setNiter(_niter);
    //sd.setBoth("y");
    float[][] p_pwdm = new float[_n2][_ng]; //pwd w/initial p=0 (Madagascar)
    sd.findSlopes(_sg,_s2,_f,p_pwdm);
    p_pwdm = mul(p_pwdm,_dg/_d2);

    if (error) {
      System.out.println("Madagascar PWD:");
      float error_pwdm = Util.rmsError(p_pwdm,_pk,_dg,_d2,T);
    }

    // interp, title, paint, colorbar, color
    String hl = "Traces"; //horizontal label
    String vl = "Samples"; //vertical label
    String cbl = "slope (samples/trace)"; //colorbar label
    Plot.plot(_sg,_s2,p_pwdm,"PWDM noise= "+_noise,hl,vl,cbl,
        _fw,_fh,-_clipMax,_clipMax,
        _clip,F,_title,_paint,T,T);
  }

  /**
   * Dynamic warping: plots the estimated slopes and RMS error.
   */
  /*public static void plotDW(boolean error) {
    DynamicWarping dw = new DynamicWarping((int)-_pmax,(int)_pmax);
    dw.setShiftSmoothing(18.0,1.0);
    dw.setStrainMax(0.1,0.1);
    dw.setErrorSmoothing(2);
    float[][] p_dw = DWSlopesAvg(dw);
    p_dw = mul(p_dw,_dg/_d2);

    if (error) {
      System.out.println("Dynamic warping:");
      float error_dw = Util.rmsError(p_dw,_pk,_dg,_d2,T);
    }
    System.out.println("slopes: "+p_dw[230][225]);

    // interp, title, paint, colorbar, color
    String hl = "Traces"; //horizontal label
    String vl = "Samples"; //vertical label
    String cbl = "slope (samples/trace)"; //colorbar label
    Plot.plot(_sg,_s2,p_dw,"DW noise= "+_noise,hl,vl,cbl,
        _fw,_fh,-_clipMax,_clipMax,
        _clip,F,_title,_paint,T,T);
  }*/

  /**
   * Dynamic warping: plots the estimated slopes and RMS error.
   */
  public static void plotDW(boolean error, int k) {
    float strainMax1 = 0.2f;
    float strainMax2 = 1.0f;
    Sampling s1 = new Sampling(_n1);
    Sampling s2 = new Sampling(_n2);
    DynamicWarpingSlopes dws = new DynamicWarpingSlopes((int)_pmax*k,s1,s2);
    dws.setK(k);
    float[][] p_dw = dws.findSlopes(_f);

    if (error) {
      System.out.println("Dynamic warping:");
      float error_dw = Util.rmsError(p_dw,_pk,_dg,_d2,T);
    }
    System.out.println("slopes: "+p_dw[230][225]);

    // interp, title, paint, colorbar, color
    String hl = "Traces"; //horizontal label
    String vl = "Samples"; //vertical label
    String cbl = "slope (samples/trace)"; //colorbar label
    Plot.plot(_sg,_s2,p_dw,"DW noise= "+_noise,hl,vl,cbl,
        _fw,_fh,-_clipMax,_clipMax,
        _clip,F,_title,_paint,T,T);
  }

  /**
   * Smooth dynamic warping: plots the estimated slopes and RMS error.
   */
  public static void plotSDW(boolean error, int k) {
    double r1min = -0.2;
    double r1max = 0.2;
    double r2min = -0.4;
    double r2max = 0.4;
    double h1 = 20.0;
    double h2 =  9.0;
    Sampling s1 = new Sampling(_n1); //shift sampling in 1st dimension
    Sampling s2 = new Sampling(_n2); //shift sampling in 2nd dimension
    DynamicWarpingSlopes dws = new DynamicWarpingSlopes(k,_pmax,h1,h2,
                                       r1min,r1max,r2min,r2max,s1,s2);
    float[][] p_sdw = dws.findSmoothSlopes(_f);
    //DynamicWarpingK dwk = new DynamicWarpingK(k,-_pmax,_pmax,s1,s2);
    //dwk.setSmoothness(h1,h2);
    //dwk.setStrainLimits(r1min,r1max,r2min,r2max);
    
    //float[][] p_dwk = SDWSlopesAvg(dwk);
    //p_dwk = mul(p_dwk,_dg/_d2);

    if (error) {
      System.out.println("Smooth dynamic warping:");
      float error_sdw = Util.rmsError(p_sdw,_pk,_dg,_d2,T);
    }
    System.out.println("slopes: "+p_sdw[230][225]);

    // interp, title, paint, colorbar, color
    String hl = "Traces"; //horizontal label
    String vl = "Samples"; //vertical label
    String cbl = "slope (samples/trace)"; //colorbar label
    Plot.plot(_sg,_s2,p_sdw,"SDW noise= "+_noise,hl,vl,cbl,
        _fw,_fh,-_clipMax,_clipMax,
        _clip,F,_title,_paint,T,T);
  }

  /**
   * PWD Dave: plots the estimated slopes and RMS error.
   */
  public static void plotPWDD(boolean error) {
    PlaneWaveDestructor pwd = new PlaneWaveDestructor(-_pmax,_pmax);
    pwd.setSmoothness(6,1);
    pwd.setLateralBias(0.0);
    pwd.setOuterIterations(_niter); // default is 5
    pwd.setInnerIterations(20); // default is 20
    float[][] p_pwdd = new float[_n2][_ng]; //pwd w/initial p=0 (Dave)
    p_pwdd= pwd.findSlopes(_f);
    pwd.updateSlopes(_f,p_pwdd);
    p_pwdd = mul(p_pwdd,_dg/_d2);

    if (error) {
      System.out.println("Dave's PWD:");
      float error_pwdd = Util.rmsError(p_pwdd,_pk,_dg,_d2,T);
    }

    // interp, title, paint, colorbar, color
    String hl = "Traces"; //horizontal label
    String vl = "Samples"; //vertical label
    String cbl = "slope (samples/trace)"; //colorbar label
    Plot.plot(_sg,_s2,p_pwdd,"PWDD noise= "+_noise,hl,vl,cbl,
        _fw,_fh,-_clipMax,_clipMax,
        _clip,F,_title,_paint,T,T);
  }

  /**
   * Plots the synthetic seismic image and the corresponding slopes image.
   */
  public static void plotFandPk() {
    System.out.println("max slope= "+max(_pk));
    // interp, title, paint, colorbar, color
    String hl = "Traces"; //horizontal label
    String vl = "Samples"; //vertical label
    String cbl = ""; //colorbar label
    Plot.plot(_sg,_s2,_f,"Synthetic noise= "+_noise,hl,vl,cbl,
        _fw,_fh,-_clipMax,_clipMax,
        F,F,_title,_paint,T,F);

    cbl = "slope (samples/trace)"; //colorbar label
    Plot.plot(_sg,_s2,_pk,"Known slopes",hl,vl,cbl,
        _fw,_fh,-_clipMax,_clipMax,
        _clip,F,_title,_paint,T,T);
  }

////////////////////TESTING/PLOTTING OPTIMAL PARAMETERS///////////////////

  /**
   * Tests for the parameters that produce the smallest RMS error.
   * The m variable requires a number 1,2,3 or 4.
   * For optimal structure tensor parameters, method = 1
   * For optimal Madagascar PWD parameters,   method = 2
   * For optimal dynamic warping parameters,  method = 3
   * For optimal Dave's PWD parameters,       method = 4
   */
  public static void testOptimalParameters(int m, String method) {
    Check.argument(m==1 || m==2 || m==3 || m==4 || m==5,"valid method");

    _n1param = _s1param.getCount();
    _n2param = _s2param.getCount();
    float[][] pe = new float[_n2][_n1];
    float[][] rmserror = new float[_n2param][_n1param]; // RMS error
    if (m==1) {
      LocalSlopeFinder lsf;
      for(int i2=0; i2<_n2param; ++i2) {
        for(int i1=0; i1<_n1param; ++i1) {
          lsf = new LocalSlopeFinder(_param1[i1],_param2[i2],_pmax);
          lsf.findSlopes(_f,pe);
          rmserror[i2][i1] = Util.rmsError(pe,_pk,_d1,_d2,false);
        }
      }
    }

    else if (m==2) {
      Sfdip sd = new Sfdip(-_pmax,_pmax);
      sd.setOrder(2);
      for(int i2=0; i2<_n2param; ++i2) {
        for(int i1=0; i1<_n1param; ++i1) {
          sd.setRect((int)_param1[i1],(int)_param2[i2]);
          sd.findSlopes(_s1,_s2,_f,pe);
          rmserror[i2][i1] = Util.rmsError(pe,_pk,_d1,_d2,false);
        }
      }
    }

    else if (m==3) {
      DynamicWarping dw = new DynamicWarping((int)-_pmax,(int)_pmax);
      dw.setStrainMax(1.0,1.0);
      dw.setErrorSmoothing(2);
      for(int i2=0; i2<_n2param; ++i2) {
        for(int i1=0; i1<_n1param; ++i1) {
          dw.setShiftSmoothing(_param1[i1],_param2[i2]);
          //pe = DWSlopesAvg(dw);
          //rmserror[i2][i1] = Util.rmsError(pe,_pk,_d1,_d2,false);
        }
      }
    }

    else if (m==4) {
      DynamicWarping dw = new DynamicWarping((int)-_pmax,(int)_pmax);
      dw.setErrorSmoothing(2);
      //dw.setShiftSmoothing(18,1);
      for(int i2=0; i2<_n2param; ++i2) {
        for(int i1=0; i1<_n1param; ++i1) {
          dw.setStrainMax(_param1[i1],_param2[i2]);
          //pe = DWSlopesAvg(dw);
          //rmserror[i2][i1] = Util.rmsError(pe,_pk,_d1,_d2,false);
        }
      }
    }

    else if (m==5) {
      PlaneWaveDestructor pwd = new PlaneWaveDestructor(-_pmax,_pmax);
      for(int i2=0; i2<_n2param; ++i2) {
        for(int i1=0; i1<_n1param; ++i1) {
          pwd.setSmoothness(_param1[i1],_param2[i2]);
          pe = pwd.findSlopes(_f);
          rmserror[i2][i1] = Util.rmsError(pe,_pk,_d1,_d2,false);
        }
      }
    }
    Util.writeBinary(rmserror,_path+method+"_errors.dat");
  }

  public static void plotOptimalParameters(String method, String hl, String vl){
    _n1param = _s1param.getCount();
    _n2param = _s2param.getCount();

    int[] error_index = new int[2];    
    float[][] rmserror = Util.readImage(_n1param,_n2param,
        _path+method+"_errors.dat");
    float min_error = min(rmserror,error_index);
    System.out.println(vl+"= "+_param1[error_index[0]]+" "+
                       hl+"= "+_param2[error_index[1]]+" "+
                      " Minimum Error Value= "+min_error);    

    // interp, title, paint, colorbar, color
    String cbl = "rms error (samples/trace)"; //colorbar label
    float cmin = 0.4f;
    float cmax = 0.75f;
    Plot.plot(_s1param,_s2param,rmserror,"RMS Error "+method,hl,vl,cbl,
        _fw,_fh,cmin,cmax,
        F,T,F,_paint,T,T);
  }

  //////////////////OTHER METHODS FOR ALGORITHM EVALUATION///////////////

  public static void testOrderVsTime(int norder) {
    Stopwatch sw = new Stopwatch();
    float[] ovt = new float[norder];
    Sfdip sd = new Sfdip(-_pmax,_pmax);
    sd.setRect(76,6);
    sd.setNiter(_niter);
    float[][] p_pwdm = new float[_n2][_n1]; //pwd w/initial p=0 (Madagascar)
    for (int i=1; i<norder; ++i) {
        sw.restart();
        sd.setOrder(i);
        sd.findSlopes(_s1,_s2,_f,p_pwdm);
        ovt[i] = (float)sw.time();
        System.out.println("ovt["+i+"]= "+ovt[i]);
    }
    Util.writeBinary(ovt,_path+"orderVsTime.dat");
  }

  public static void plotOrderVsTime(int norder) {
    float[] ovt = Util.readImage(norder,_path+"orderVsTime.dat");
    // paint 
    Plot.plot(ovt,"Order vs time","Order","Time",_fw,_fh,F);
  }

  /*public static void rmsErrorCurve(int nrms) {
    for (int i=0; i<nrms; ++i) {
      _noise = i;
    }

  }
  */
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

  ////////////////ESTIMATING SLOPES USING (SMOOTH) DYNAMIC WARPING/////////////

  public static float[][] SDWSlopesAvg(DynamicWarpingK dwk) {
    float[][] pp = new float[_n2][_ng];
    float[][] pm = new float[_n2][_ng];
    float[][] pa = new float[_n2][_ng];
    float[][] fp = new float[_n2][_ng];
    float[][] fm = new float[_n2][_ng];

    fp[0]    = _f[1];
    fm[0]    = _f[0];
    fp[_n2-1] = _f[_n2-1];
    fm[_n2-1] = _f[_n2-2];

    for (int i2=1; i2<_n2-1; ++i2) {
      fp[i2] = _f[i2+1];
      fm[i2] = _f[i2-1];
    }

    pp = dwk.findShifts(_sg,_f,_sg,fp);
    pm = dwk.findShifts(_sg,_f,_sg,fm);
    pa = sub(pp,pm);
    pa = mul(pa,0.5f);
    //pm = mul(pm,-1.0f);
    return pa;
  }

  /*
  public static float[][] SDWSlopesAvg(DynamicWarpingR dwr, int nl) {
    float[][] g = copy(_ng,_n2,_f);
    g = copy(_ng+nl,_n2,g);
    Sampling sf = new Sampling(_ng);
    Sampling sg = new Sampling(_ng+nl);
    float[][] p = dwr.findShifts(_sg,_f,_sg,g);
    return p;
  }*/

  /*
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
  */

  ///////////////////PRIVATE VARIABLES///////////////////////
  private static final String _path = 
    "/Users/earias/Home/git/ea/bench/src/slopes/data/";
  private static final int _niter = 5;
  private static final float _fw = 0.75f; //fraction width for slide
  private static final float _fh = 0.9f; //fraction height for slide
  private static final float _pmax = 5.0f;
  private static final float _clipMax = 4.0f;
  private static final float pi = FLT_PI;
  private static final boolean T = true;
  private static final boolean F = false;  
  private static final boolean _title = false;
  private static final boolean _paint = true;  
  private static final boolean _clip = true;  
  private static Sampling _s1 = new Sampling(1,1,1);  
  private static Sampling _s2 = new Sampling(1,1,1);
  private static Sampling _sg = new Sampling(1,1,1);
  private static Sampling _s1param = new Sampling(1,1,1);
  private static Sampling _s2param = new Sampling(1,1,1);

  private static int _n1,_n2,_ng,_n1param,_n2param;
  private static float _d1,_d2,_dg,_f1,_f2,_fg,_dparam,_fparam,_noise;
  private static float[] _param1 = new float[_n1param];
  private static float[] _param2 = new float[_n2param];
  private static float[][] _f = new float[_n2][_n1];
  private static float[][] _pk = new float[_n2][_n1];

}
