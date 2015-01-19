package slopes;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.Stopwatch;

import static edu.mines.jtk.util.ArrayMath.*;

import dnp.LocalSlopeFinder;
import dnp.PlaneWaveDestructor;
import util.FakeData;
import util.*;

import javax.swing.*;

/**
 * This software acts as a machine for testing slope estimation algorithms.
 * Using other pieces of software, I generate and save a synthetic image 
 * containing structural features that may be encountered in real seismic 
 * images. Since I generated the synthetic image, I know the exact slopes
 * at every sample in that image. This allows me to test how well different
 * slope estimating algorithms work. The three methods tested in this software
 * are the plane-wave destruction filter method, structure tensor method, and
 * dynamic warping method. 
 * @author Elias Arias, Colorado School of Mines, CWP
 * @version 15.1.2015
 */

public class SlopeAlgorithmEval {

  private static void timeAndPlotLSF() {
    sw.start();
    Slopes.goLSF(s1,s2,f,pk,error);
    sw.stop();
    System.out.println("Structure tensor time = "+sw.time());
  }

  private static void timeAndPlotPWDM() {
    sw.restart();
    Slopes.goPWDM(s1,s2,f,pk,error);
    sw.stop();
    System.out.println("Madagascar PWD time = "+sw.time());
  }

  private static void timeAndPlotDW() {
    sw.restart();
    Slopes.goDW(s1,s2,f,pk,error);
    sw.stop();
    System.out.println("Dynamic warping time = "+sw.time());
  }

  private static void timeAndPlotPWDD() {
    sw.restart();
    Slopes.goPWDD(s1,s2,f,pk,error);
    sw.stop();
    System.out.println("Dave's PWD time = "+sw.time());
  }

  private static void setSynthParameters() {
    n1 = 501;
    n2 = 501;

    d1 = 1.0f;
    d2 = 1.0f;
    f1 = 0.0f;
    f2 = 0.0f;

    s1 = new Sampling(n1,d1,f1);
    s2 = new Sampling(n2,d2,f2);

    //Synth data
    float[][][] fandpk = FakeData.seismicAndSlopes2d2014B(noise,F);
    f = fandpk[0];  //synthetic seismic data
    pk = fandpk[1]; //exact slope values
  }

  private static void setGOMParameters() {
    n1 = 301;
    n2 = 920;

    d1 = 0.004f;
    d2 = .02667f;
    f1 = 1.6f;
    f2 = 0;

    s1 = new Sampling(n1,d1,f1);
    s2 = new Sampling(n2,d2,f2);

    f = Util.readImage(n1,n2,"data/gom.dat");
    float[][] fs = new float[n2][n1];
    mul(f,.001f,fs);
    Util.writeBinary(fs,"data/gom_scaled.dat");
  }
/////////////////////////////NOT REFACTORED////////////////////////////  
  private static void goTestStrainValuesDW(Sampling s1, Sampling s2, 
      float[][] f, float[][] pk) {

    int n1 = s1.getCount();
    int n2 = s2.getCount();
    float d1 = (float)s1.getDelta();
    float d2 = (float)s2.getDelta();

    float dstrain = 0.1f; // strain sampling rate
    float lstrain = 1.0f; // last strain 
    float fstrain = 0.1f;
    int nstrain1 = (int)((lstrain-fstrain)/dstrain); //# strain vals to test
    int nstrain2 = nstrain1;
    Sampling sstrain1 = new Sampling(nstrain1,dstrain,fstrain);
    Sampling sstrain2 = new Sampling(nstrain2,dstrain,fstrain);

    float[][] rmserror = new float[nstrain2][nstrain1]; // RMS error

    float[] strain1 = new float[nstrain1];
    float[] strain2 = new float[nstrain2];
    for(int i=0; i<nstrain1; ++i) {
      strain1[i] = i*dstrain+fstrain;
      strain2[i] = i*dstrain+fstrain;
    }
    
    float[][] pe = new float[n2][n1];

    int pmax = 5;
    DynamicWarping dw = new DynamicWarping(-pmax,pmax);
    dw.setShiftSmoothing(15.0,1.0);
    dw.setErrorSmoothing(2);
    
    for(int i2=0; i2<nstrain2; ++i2) {
      for(int i1=0; i1<nstrain1; ++i1) {
        dw.setStrainMax(strain1[i1],strain2[i2]);
        pe = Util.DWSlopesAvg(dw,f);
        rmserror[i2][i1] = Util.rmsError(pe,pk,d1,d2,true);
      }
    }
    Util.writeBinary(rmserror,
      "/Users/earias/Home/git/ea/bench/src/slopes/data/strainErrValsDWavg.dat");
  }

  private static void goPlotStrainValuesDW() {
    float dstrain= 0.1f; // strain sampling rate
    float lstrain= 1.0f; // last strain
    float fstrain= 0.1f;
    int nstrain1 = (int)((lstrain-fstrain)/dstrain); // # strain vals to test
    int nstrain2 = nstrain1;
    int[] err_index = new int[2];

    float[] strain1= new float[nstrain1];
    float[] strain2= new float[nstrain2];
    for(int i=0; i<nstrain1; ++i) {
      strain1[i] = i*dstrain+fstrain;
    }

    for(int i=0; i<nstrain2; ++i) {
      strain2[i] = i*dstrain+fstrain;
    }

    Sampling sstrain1 = new Sampling(nstrain1,dstrain,fstrain);
    Sampling sstrain2 = new Sampling(nstrain2,dstrain,fstrain);

    float[][] err = Util.readImage(nstrain1,nstrain2,
     "/Users/earias/Home/git/ea/bench/src/slopes/data/strainErrValsDWavg.dat");
    float min_err = min(err,err_index);
    System.out.println("Strain1= "+strain1[err_index[0]]+
                      " Strain2= "+strain2[err_index[1]]+
                      " Minimum Error Value= "+min_err);

    float fw = 0.75f; //fraction width for slide
    float fh = 0.9f; //fraction height for slide
    // strain, interp, title, paint, colorbar, color
    Plot.plot(sstrain1,sstrain2,err,"RMS Error Strain DW",fw,fh,T,T,F,T,T,T);
  }

  private static void goRmsErrorCurves() {
    int n1 = 501;
    int n2 = 501;

    float d1 = 1.0f;
    float d2 = 1.0f;

    float f1 = 0;
    float f2 = 0;

    Sampling s1 = new Sampling(n1,d1,f1);
    Sampling s2 = new Sampling(n2,d2,f2);

    float pmax = 5.0f;

    float[][][] synthAndSlope = new float[2][n2][n1];
    float[][] synth_data = new float[n2][n1];
    float[][] exact_slope = new float[n2][n1];

    float[][] mad_diff = new float[n2][n1];
    float[][] temp_mad = new float[n2][n1];

    float[][] lsf_diff = new float[n2][n1];
    float[][] temp_lsf = new float[n2][n1];

    float[][] dw_diff = new float[n2][n1];
    float[][] temp_dw = new float[n2][n1];

    //Structure tensor
    LocalSlopeFinder lsf = new LocalSlopeFinder(24.0f,1.0f,pmax);
    float[][] lsf_slope = new float[n2][n1];

    //Plane-wave destruction filter
    Sfdip sd = new Sfdip(-pmax,pmax);
    sd.setRect(76,6);
    sd.setOrder(2);
    float[][] mad_slope = new float[n2][n1];

    //Dynamic warping
    DynamicWarping dw = new DynamicWarping((int)-pmax,(int)pmax);
    dw.setShiftSmoothing(15.0,1.0);
    dw.setStrainMax(0.1,0.5);
    dw.setErrorSmoothing(2);
    float[][] dw_slope = new float[n2][n1];

    int nerr = 21;
    float[] noise = new float[nerr];
    float[] rmserr_lsf = new float[nerr];
    float[] rmserr_mad = new float[nerr];
    float[] rmserr_dw = new float[nerr];
    for (int i=0; i<nerr; ++i) {
      noise[i] = 0.05f*i;
      synthAndSlope = FakeData.seismicAndSlopes2d2014B(noise[i],F);
      synth_data = synthAndSlope[0];
      exact_slope = synthAndSlope[1];
      lsf.findSlopes(synth_data,lsf_slope);
      sd.findSlopes(s1,s2,synth_data,mad_slope);
      dw_slope = Util.DWSlopesAvg(dw,synth_data);
      lsf_slope = mul(lsf_slope,d1/d2);
      mad_slope = mul(mad_slope,d1/d2);
      dw_slope = mul(dw_slope,d1/d2);
      lsf_diff = sub(lsf_slope,exact_slope);
      mad_diff = sub(mad_slope,exact_slope);
      dw_diff = sub(dw_slope,exact_slope);
      temp_lsf = pow(lsf_diff,2);
      temp_mad = pow(mad_diff,2);
      temp_dw = pow(dw_diff,2);
      rmserr_lsf[i] = sqrt(sum(temp_lsf)/(n2*n1));
      rmserr_mad[i] = sqrt(sum(temp_mad)/(n2*n1));
      rmserr_dw[i] = sqrt(sum(temp_dw)/(n2*n1));
    }

    float fw = 0.75f; //fraction width for slide
    float fh = 0.9f; //fraction height for slide
    // paint
    Plot.plot(noise,rmserr_mad,rmserr_lsf,rmserr_dw,
        "PWD, LOF, and DW RMS error vs NS ratio",fw,fh,T);
  }

  /**
   * Dp = delta p = change in slope
   */
  /**
  private static void goDpImages() {
    int n1 = 501;
    int n2 = 501;

    float d1 = 1.0f;
    float d2 = 1.0f;

    float f1 = 0;
    float f2 = 0;
    float pmax = 5.0f;
    float noise = 0.5f;
    float w = 10.0f; // fpeak=0.1   1/fpeak=10
    float c = 1.0f/FLT_PI;
    float sigma1 = 13.0f;
    float sigma2 = 1.0f;

    Sampling s1 = new Sampling(n1,d1,f1);
    Sampling s2 = new Sampling(n2,d2,f2);

    float[][][] synthAndSlope = FakeData.seismicAndSlopes2d2014B(noise,T);

    float[][] synth_data = synthAndSlope[0];
    float[][] p = synthAndSlope[1];
    System.out.println("max slope= "+max(p));

    float[][] dp = new float[n2][n1];
    dp = mul(w,add(c/sigma2,mul(p,c/sigma1))); 

    float fw = 0.9f; //fraction width for slide
    float fh = 0.8f; //fraction height for slide
    // title, paint, colorbar, color
    Plot.plot(s1,s2,dp,"DP image",fw,fh,-4,4,F,T,T,T);
    Plot.plot(s1,s2,p,"Exact Slopes",fw,fh,-4,4,F,T,T,T);
  }*/

  private static void goTestSlopeVsError() {
    int n1 = 501;
    int n2 = 501;

    float noise = 0.5f;

    float[][][] synthAndSlope = FakeData.seismicAndSlopes2d2014B(noise,F);
    float[][] exact_slope = synthAndSlope[1];
    float[][] lsf_mean = Util.readImage(n1,n2,
        "/Users/earias/Home/git/ea/bench/src/pwd/data/lsf_mean.dat");
    float[][] mad_mean = Util.readImage(n1,n2,
        "/Users/earias/Home/git/ea/bench/src/pwd/data/mad_mean.dat");
    float[][] lsf_sd = Util.readImage(n1,n2,
        "/Users/earias/Home/git/ea/bench/src/pwd/data/lsf_sd.dat");
    float[][] mad_sd = Util.readImage(n1,n2,
        "/Users/earias/Home/git/ea/bench/src/pwd/data/mad_sd.dat");

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
        "/Users/earias/Home/git/ea/bench/src/pwd/data/exact_slope1D.dat");
    Util.writeBinary(lsf_1Dslope,
        "/Users/earias/Home/git/ea/bench/src/pwd/data/struct_tens1D.dat");
    Util.writeBinary(mad_1Dslope,
        "/Users/earias/Home/git/ea/bench/src/pwd/data/plane_wave1D.dat");
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
  }

  private static void goTestSampleMeanSD() {
    int n1 = 501;
    int n2 = 501;

    float d1 = 1.0f;
    float d2 = 1.0f;

    float f1 = 0;
    float f2 = 0;
    float pmax = 5.0f;
    float noise = 0.5f;

    Sampling s1 = new Sampling(n1,d1,f1);
    Sampling s2 = new Sampling(n2,d2,f2);

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
      dw_slope = Util.DWSlopesAvg(dw,synth_data);

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
  }

  private static void goPlotSampleMeanSD() {
    int n1 = 501;
    int n2 = 501;

    float d1 = 1.0f;
    float d2 = 1.0f;

    float f1 = 0;
    float f2 = 0;

    Sampling s1 = new Sampling(n1,d1,f1);
    Sampling s2 = new Sampling(n2,d2,f2);

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
    float fw = 0.75f; //fraction width for slide
    float fh = 0.9f; //fraction height for slide
    // title, paint, colorbar, color
    //Plot.plot(s1,s2,lsf_mean,"Sample mean LOF",fw,fh,-4,4,F,T,T,T);
    //Plot.plot(s1,s2,mad_mean,"Sample mean PWD",fw,fh,-4,4,F,T,T,T);
    //Plot.plot(s1,s2,lsf_sd,"Sample variance LOF 03 max",fw,fh,0,0.3f,F,F,T,T);
    //Plot.plot(s1,s2,mad_sd,"Sample variance PWD",fw,fh,0,0.5f,F,F,T,T);
    Plot.plot(s1,s2,lsf_sd,"Sample standard deviation LOF",fw,fh,0,1,F,T,T,T);
    Plot.plot(s1,s2,mad_sd,"Sample standard deviation PWD",fw,fh,0,1,F,T,T,T);
    Plot.plot(s1,s2,dw_sd,"Sample standard deviation DW",fw,fh,0,1,F,T,T,T);
  }

  private static void goErrorLocation() {
    int n1 = 501;
    int n2 = 501;

    float d1 = 1.0f;
    float d2 = 1.0f;

    float f1 = 0;
    float f2 = 0;
    float pmax = 5.0f;
    float noise = 0.5f;

    Sampling s1 = new Sampling(n1,d1,f1);
    Sampling s2 = new Sampling(n2,d2,f2);

    int i = 5;
    float[][][] synthAndSlope = FakeData.seismicAndSlopes2d2014B(noise,T);

    float[][] synth_data = synthAndSlope[0];
    float[][] exact_slope = synthAndSlope[1];
    System.out.println("max slope= "+max(exact_slope));

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
  }

///////////////////VARIABLES///////////////////////
  private static final boolean T = true;
  private static final boolean F = false;
  private static final boolean title = true;
  private static final boolean paint = true;
  private static final boolean error = false;
  private static final float pi = FLT_PI;
  private static final float noise = 0.5f;
  private static Stopwatch sw = new Stopwatch();

  private static int n1,n2;
  private static float d1,d2,f1,f2;
  private static Sampling s1 = new Sampling(1,1,1);
  private static Sampling s2 = new Sampling(1,1,1);
  private static float[][] f = new float[n2][n1];  //synthetic seismic data
  private static float[][] pk = new float[n2][n1]; //exact slope values

  public static void main(String[] args) {
    SwingUtilities.invokeLater(new Runnable() {
      public void run() {
        float d_param= 1.0f;
        float l_param = 20.0f;
        float f_param = 1.0f;
        int n1_param = (int)((l_param-f_param+1)/d_param);
        int n2_param = 10;
        Sampling s1_param = new Sampling(n1_param,d_param,f_param);
        Sampling s2_param = new Sampling(n2_param,d_param,f_param);

        float[] param1 = new float[n1_param];
        float[] param2 = new float[n2_param];
        for(int i=0; i<n1_param; ++i) {
          param1[i] = i*d_param+f_param;
        }
        for(int i=0; i<n2_param; ++i) {
          param2[i] = i*d_param+f_param;
        }

////////////////////////////////PLOTTING////////////////////////////////////
        /*
        //Testing for optimal parameters
        Util.testOptimalParameters(s1_param,s2_param,s1,s2,param1,param2,1,
            "Structure tensor");
        Util.testOptimalParameters(s1_param,s2_param,s1,s2,param1,param2,2,
            "Madagascar PWD");
        Util.testOptimalParameters(s1_param,s2_param,s1,s2,param1,param2,3,
            "Dynamic warping");
        Util.testOptimalParameters(s1_param,s2_param,s1,s2,param1,param2,4,
            "Dave's PWD");
        */
        //setSynthParameters();
        setGOMParameters();

        timeAndPlotLSF();
        timeAndPlotPWDM();
        timeAndPlotDW();
        timeAndPlotPWDD();

        setSynthParameters();
        Slopes.goFandPk(s1,s2);
        
        /*
        //Plotting for optimal parameters
        System.out.println("Optimal sigma");
        Util.plotOptimalParameters(s1_param,s2_param,s1,s2,
            param1,param2,"Structure tensor");
        System.out.println("Optimal rect");
        Util.plotOptimalParameters(s1_param,s2_param,s1,s2,
            param1,param2,"Madagascar PWD");
        System.out.println("Optimal shift smoothing");
        Util.plotOptimalParameters(s1_param,s2_param,s1,s2,
            param1,param2,"Dynamic warping");
        System.out.println("Optimal smoothing");
        Util.plotOptimalParameters(s1_param,s2_param,s1,s2,
            param1,param2,"Dave's PWD");
        */

        //goDpImages();
        //goTestSlopeVsError();
        //goRmsErrorCurves();
        //goTestSampleMeanSD();
        //goPlotSampleMeanSD();
        //goErrorLocation();
        //goTestStrainValuesDW();
        //goPlotStrainValuesDW();
        //goTestPmaxValues();
      }
    });
  }
}
