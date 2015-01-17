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

  private static void goLSF() {
    int n1 = 501;
    int n2 = 501;

    float d1 = 1.0f;
    float d2 = 1.0f;

    float f1 = 0;
    float f2 = 0;
    float pmax = 5.0f;

    Sampling s1 = new Sampling(n1,d1,f1);
    Sampling s2 = new Sampling(n2,d2,f2);

    float[][][] fandpk = FakeData.seismicAndSlopes2d2014B(noise,F);

    float[][] f = fandpk[0];  //synthetic seismic data
    float[][] pk = fandpk[1]; //exact slope values
    //System.out.println("max slope= "+max(exact_slope));

    //Structure tensor method
    LocalSlopeFinder lsf = new LocalSlopeFinder(24.0f,1.0f,pmax);
    float[][] p_lsf = new float[n2][n1];
    lsf.findSlopes(f,p_lsf);

    System.out.println("Structure tensor:");
    float error_lsf = Util.rmsError(p_lsf,pk,d1,d2,T);

    float fw = 0.75f; //fraction width for slide
    float fh = 0.9f; //fraction height for slide
    // title, paint, colorbar, color
    Plot.plot(s1,s2,p_lsf,"LSF noise= "+noise,fw,fh,-4,4,T,F,T,T);
  }

  private static void goPWDM() {
    int n1 = 501;
    int n2 = 501;

    float d1 = 1.0f;
    float d2 = 1.0f;

    float f1 = 0;
    float f2 = 0;
    float pmax = 5.0f;
    int niter = 5;

    Sampling s1 = new Sampling(n1,d1,f1);
    Sampling s2 = new Sampling(n2,d2,f2);

    float[][][] fandpk = FakeData.seismicAndSlopes2d2014B(noise,F);

    float[][] f = fandpk[0];  //synthetic seismic data
    float[][] pk = fandpk[1]; //exact slope values
    //System.out.println("max slope= "+max(exact_slope));

    //Plane-wave destruction filter
    Sfdip sd = new Sfdip(-pmax,pmax);
    sd.setRect(76,6);
    sd.setOrder(5);
    sd.setNiter(niter);
    //sd.setBoth("y");
    float[][] p_pwdm = new float[n2][n1]; //pwd w/initial p=0 (Madagascar)
    sd.findSlopes(s1,s2,f,p_pwdm);

    System.out.println("Madagascar PWD:");
    float error_pwdm = Util.rmsError(p_pwdm,pk,d1,d2,T);

    float fw = 0.75f; //fraction width for slide
    float fh = 0.9f; //fraction height for slide
    // title, paint, colorbar, color
    Plot.plot(s1,s2,p_pwdm,"PWDM noise= "+noise,fw,fh,-4,4,T,F,T,T);
  }

  private static void goDW() {
    int n1 = 501;
    int n2 = 501;

    float d1 = 1.0f;
    float d2 = 1.0f;

    float f1 = 0;
    float f2 = 0;
    float pmax = 5.0f;

    Sampling s1 = new Sampling(n1,d1,f1);
    Sampling s2 = new Sampling(n2,d2,f2);

    float[][][] fandpk = FakeData.seismicAndSlopes2d2014B(noise,F);

    float[][] f = fandpk[0];  //synthetic seismic data
    float[][] pk = fandpk[1]; //exact slope values
    //System.out.println("max slope= "+max(exact_slope));

    //Dynamic warping
    DynamicWarping dw = new DynamicWarping((int)-pmax,(int)pmax);
    dw.setShiftSmoothing(40.0,1.0);
    dw.setStrainMax(1.0,1.0);
    dw.setErrorSmoothing(2);
    float[][] p_dw= Util.DWSlopesAvg(dw,f);

    System.out.println("Dynamic warping:");
    float error_dw = Util.rmsError(p_dw,pk,d1,d2,T);

    float fw = 0.75f; //fraction width for slide
    float fh = 0.9f; //fraction height for slide
    // title, paint, colorbar, color
    Plot.plot(s1,s2,p_dw,"DW noise= "+noise,fw,fh,-4,4,T,F,T,T);
  }

  private static void goPWDD() {
    int n1 = 501;
    int n2 = 501;

    float d1 = 1.0f;
    float d2 = 1.0f;

    float f1 = 0;
    float f2 = 0;
    float pmax = 5.0f;
    int niter = 5;

    Sampling s1 = new Sampling(n1,d1,f1);
    Sampling s2 = new Sampling(n2,d2,f2);

    float[][][] fandpk = FakeData.seismicAndSlopes2d2014B(noise,F);

    float[][] f = fandpk[0];  //synthetic seismic data
    float[][] pk = fandpk[1]; //exact slope values
    //System.out.println("max slope= "+max(exact_slope));

    //Dave's plane-wave destruction filter
    PlaneWaveDestructor pwd = new PlaneWaveDestructor(-pmax,pmax);
    pwd.setSmoothness(6,1);
    pwd.setLateralBias(0.0);
    pwd.setOuterIterations(niter); // default is 5
    pwd.setInnerIterations(20); // default is 20
    float[][] p_pwdd = new float[n2][n1]; //pwd w/initial p=0 (Dave)
    p_pwdd= pwd.findSlopes(f);
    pwd.updateSlopes(f,p_pwdd);

    System.out.println("Dave's PWD:");
    float error_pwdd = Util.rmsError(p_pwdd,pk,d1,d2,T);

    float fw = 0.75f; //fraction width for slide
    float fh = 0.9f; //fraction height for slide
    // title, paint, colorbar, color
    Plot.plot(s1,s2,p_pwdd,"PWDD noise= "+noise,fw,fh,-4,4,T,F,T,T);
  }

  private static void goFandPk() {
    int n1 = 501;
    int n2 = 501;

    float d1 = 1.0f;
    float d2 = 1.0f;

    float f1 = 0;
    float f2 = 0;

    Sampling s1 = new Sampling(n1,d1,f1);
    Sampling s2 = new Sampling(n2,d2,f2);

    float[][][] fandpk = FakeData.seismicAndSlopes2d2014B(noise,F);

    float[][] f = fandpk[0];  //synthetic seismic data
    float[][] pk = fandpk[1]; //exact slope values
    //System.out.println("max slope= "+max(exact_slope));

    float fw = 0.75f; //fraction width for slide
    float fh = 0.9f; //fraction height for slide
    // title, paint, colorbar, color
    Plot.plot(s1,s2,f,"Synthetic seismic noise= "+noise,fw,fh,F,T,T,F);
    Plot.plot(s1,s2,pk,"Known Slopes",fw,fh,-4.0f,4.0f,F,F,T,T);
  }
  
  /*
  private static void goGom() {
    int n1 = 301;
    int n2 = 920;
    int n  = n1*n2;

    float d1 = 0.004f;
    float d2 = .02667f;

    float f1 = 1.6f;
    float f2 = 0;
    float pmax = 5.0f;

    Sampling s1 = new Sampling(n1,d1,f1);
    Sampling s2 = new Sampling(n2,d2,f2);

    float[][] orig_data = Util.readImage(n1,n2,"data/gom.dat");
    float[][] scaled_data = new float[n2][n1];
    mul(orig_data,.001f,scaled_data);
    Util.writeBinary(scaled_data,"data/gom_scaled.dat");

    //Structure tensor
    LocalSlopeFinder lsf = new LocalSlopeFinder(8.0f,4.0f,pmax);
    float[][] lsf_slope = new float[n2][n1];
    lsf.findSlopes(orig_data,lsf_slope);

    //Plane-wave destruction filter
    Sfdip sd = new Sfdip(-pmax,pmax);
    sd.setRect(12,6);
    float[][] mad_slope = new float[n2][n1];
    sd.findSlopes(s1,s2,orig_data,mad_slope);

    //Dynamic warping
    DynamicWarping dw = new DynamicWarping((int)-pmax,(int)pmax);
    dw.setShiftSmoothing(15.0,1.0);
    dw.setStrainMax(0.1,0.5);
    dw.setErrorSmoothing(2);
    float[][] dw_slope = DWSlopesAvg(dw,orig_data);
    
    lsf_slope = mul(lsf_slope,d1/d2);
    mad_slope = mul(mad_slope,d1/d2);
    dw_slope = mul(dw_slope,d1/d2);

    System.out.println("lsf slope= "+lsf_slope[234][54]);
    System.out.println("mad slope= "+mad_slope[234][54]);
    System.out.println("dw slope= "+dw_slope[234][54]);

    float fw = 0.75f; //fraction width for slide
    float fh = 0.9f; //fraction height for slide
    // title, paint, colorbar, color
    Plot.plot(s1,s2,Util.sexp(scaled_data),
        "GOM Near Offset Data (Gained)",fw,fh,T,F,F,F);
    Plot.plot(s1,s2,lsf_slope,"LSF Slopes",fw,fh,T,F,T,T);
    Plot.plot(s1,s2,mad_slope,"Madagascar Slopes",fw,fh,T,F,T,T);
    Plot.plot(s1,s2,dw_slope,"DW Slopes",fw,fh,T,F,T,T);
    //Plot.plot(s1,s2,sexp(scaled_data),lsf_slope,mad_slope,"teaser",F);
  }*/

  private static void goPlotDavePWD() {
    float dsigma = 1.0f; // sigma sampling rate
    float lsigma = 40.0f; // last sigma
    float fsigma = 1.0f;
    int nsigma1 = (int)((lsigma-fsigma+1)/dsigma); // # of sigma vals to test
    int nsigma2 = 10;
    int[] err_index = new int[2];

    float[] sigma1 = new float[nsigma1];
    float[] sigma2 = new float[nsigma2];
    for(int i=0; i<nsigma1; ++i) {
      sigma1[i] = i*dsigma+fsigma;
    }
    for(int i=0; i<nsigma2; ++i) {
      sigma2[i] = i*dsigma+fsigma;
    }

    Sampling ssigma1 = new Sampling(nsigma1,dsigma,fsigma);
    Sampling ssigma2 = new Sampling(nsigma2,dsigma,fsigma);

    float[][] err_vals = Util.readImage(nsigma1,nsigma2,
          "/Users/earias/Home/git/ea/bench/src/slopes/data/DavePWDtest.dat");
    float min_err = min(err_vals,err_index);
    System.out.println("Sigma1= "+sigma1[err_index[0]]+
                      " Sigma2= "+sigma2[err_index[1]]+
                      " Minimum Error Value= "+min_err);

    float fw = 0.75f; //fraction width for slide
    float fh = 0.9f; //fraction height for slide
    // sigma, interp, title, paint, colorbar, color
    Plot.plot(ssigma1,ssigma2,err_vals,"RMS Error Sigma",fw,fh,T,T,F,T,T,T);
  }

  private static void goTestStrainValuesDW() {
    int n1 = 501;
    int n2 = 501;

    float dstrain = 0.1f; // strain sampling rate
    float lstrain = 1.0f; // last strain 
    float fstrain = 0.1f;
    int nstrain1 = (int)((lstrain-fstrain)/dstrain); //# strain vals to test
    int nstrain2 = nstrain1;

    Sampling sstrain1 = new Sampling(nstrain1,dstrain,fstrain);
    Sampling sstrain2 = new Sampling(nstrain2,dstrain,fstrain);

    float d1 = 1;
    float d2 = 1;

    float f1 = 0;
    float f2 = 0;

    int[] err_index = new int[2];
    float[][] err = new float[nstrain2][nstrain1]; // RMS error

    float[] strain1 = new float[nstrain1];
    float[] strain2 = new float[nstrain2];
    for(int i=0; i<nstrain1; ++i) {
      strain1[i] = i*dstrain+fstrain;
      //System.out.println("strain1 vals= "+strain1[i]);
    }
    
    for(int i=0; i<nstrain2; ++i) {
      strain2[i] = i*dstrain+fstrain;
      //System.out.println("strain2 vals= "+strain2[i]);
    }
    
    Sampling s1 = new Sampling(n1,d1,f1);
    Sampling s2 = new Sampling(n2,d2,f2);

    float noise = 0.5f;
    float[][][] synthAndSlope = FakeData.seismicAndSlopes2d2014B(noise,F);
    float[][] synth_data = synthAndSlope[0];
    float[][] exact_slope = synthAndSlope[1];
    float[][] dw_slope = new float[n2][n1];

    int pmax = 5;
    DynamicWarping dw = new DynamicWarping(-pmax,pmax);
    dw.setShiftSmoothing(15.0,1.0);
    dw.setErrorSmoothing(4);
    float[][] dw_diff;
    
    for(int i2=0; i2<nstrain2; ++i2) {
      for(int i1=0; i1<nstrain1; ++i1) {
        dw.setStrainMax(strain1[i1],strain2[i2]);
        dw_slope = Util.DWSlopesAvg(dw,synth_data);
        dw_slope = mul(dw_slope,d1/d2);
        dw_diff = sub(dw_slope,exact_slope);
        float[][] temp = pow(dw_diff,2);
        err[i2][i1] = sqrt(sum(temp)/(n2*n1));
      }
    }
    Util.writeBinary(err,
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

  private static void goPlotSigmaValuesDW() {
    float dsigma = 1.0f; // sigma sampling rate
    float lsigma = 80.0f; // last sigma
    float fsigma = 1.0f;
    int nsigma1 = (int)((lsigma-fsigma+1)/dsigma); // # of sigma vals to test
    int nsigma2 = 10;
    int[] err_index = new int[2];

    float[] sigma1 = new float[nsigma1];
    float[] sigma2 = new float[nsigma2];
    for(int i=0; i<nsigma1; ++i) {
      sigma1[i] = i*dsigma+fsigma;
    }
    for(int i=0; i<nsigma2; ++i) {
      sigma2[i] = i*dsigma+fsigma;
    }

    Sampling ssigma1 = new Sampling(nsigma1,dsigma,fsigma);
    Sampling ssigma2 = new Sampling(nsigma2,dsigma,fsigma);

    float[][] err = Util.readImage(nsigma1,nsigma2,
      "/Users/earias/Home/git/ea/bench/src/slopes/data/sigmaErrValsDWavg.dat");
    float min_err = min(err,err_index);
    System.out.println("Sigma1= "+sigma1[err_index[0]]+
                      " Sigma2= "+sigma2[err_index[1]]+
                      " Minimum Error Value= "+min_err);

    float fw = 0.75f; //fraction width for slide
    float fh = 0.9f; //fraction height for slide
    // sigma, interp, title, paint, colorbar, color
    Plot.plot(ssigma1,ssigma2,err,"RMS Error Sigma DW avg",fw,fh,T,T,F,T,T,T);
  }

  private static void goPlotSigmaValues() {
    float dsigma = 1.0f; // sigma sampling rate
    float lsigma = 80.0f; // last sigma
    float fsigma = 1.0f;
    int nsigma1 = (int)((lsigma-fsigma+1)/dsigma); // # of sigma vals to test
    int nsigma2 = 10;
    int[] err_index = new int[2];

    float[] sigma1 = new float[nsigma1];
    float[] sigma2 = new float[nsigma2];
    for(int i=0; i<nsigma1; ++i) {
      sigma1[i] = i*dsigma+fsigma;
    }
    for(int i=0; i<nsigma2; ++i) {
      sigma2[i] = i*dsigma+fsigma;
    }

    Sampling ssigma1 = new Sampling(nsigma1,dsigma,fsigma);
    Sampling ssigma2 = new Sampling(nsigma2,dsigma,fsigma);

    float[][] err_vals = Util.readImage(nsigma1,nsigma2,
          "/Users/earias/Home/git/ea/bench/src/slopes/data/sigmaErrVals.dat");
    float min_err = min(err_vals,err_index);
    System.out.println("Sigma1= "+sigma1[err_index[0]]+
                      " Sigma2= "+sigma2[err_index[1]]+
                      " Minimum Error Value= "+min_err);

    float fw = 0.75f; //fraction width for slide
    float fh = 0.9f; //fraction height for slide
    // sigma, interp, title, paint, colorbar, color
    Plot.plot(ssigma1,ssigma2,err_vals,"RMS Error Sigma",fw,fh,T,T,F,T,T,T);
  }

  private static void goPlotRectValues() {
    int frect = 1;
    int lrect = 80;
    int drect = 1;
    int nrect1 = (lrect-frect+1)/drect;
    int nrect2 = 10;

    int[] rect1 = new int[nrect1];
    int[] rect2 = new int[nrect2];

    for(int i=0; i<nrect1; ++i) {
      rect1[i] = drect*i+frect;
    }
    for(int i=0; i<nrect2; ++i) {
      rect2[i] = drect*i+frect;
    }

    int[] err_index = new int[2];

    Sampling srect1 = new Sampling(nrect1,drect,frect);
    Sampling srect2 = new Sampling(nrect2,drect,frect);

    float[][] err_vals = Util.readImage(nrect1,nrect2,
        "/Users/earias/Home/git/ea/bench/src/slopes/data/rectErrVals.dat");
    float min_err = min(err_vals,err_index);
    System.out.println("Rect1= "+rect1[err_index[0]]+
                      " Rect2= "+rect2[err_index[1]]+
                      " Minimum Error Value= "+min_err);

    float fw = 0.75f; //fraction width for slide
    float fh = 0.9f; //fraction height for slide
    // sigma, interp, title, paint, colorbar, color
    Plot.plot(srect1,srect2,err_vals,"RMS Error Rect",fw,fh,F,T,F,T,T,T);
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

  public static final String DATA_NAME = "gom_scaled.dat";
  public static final boolean T = true;
  public static final boolean F = false;
  public static final float pi = FLT_PI;
  private static final float noise = 0.5f;
  public static void main(String[] args) {
    SwingUtilities.invokeLater(new Runnable() {
      public void run() {
        // have user input argument to show which plot they want to see
        // e.g. 1 for synthetic 
        //      2 for std dev....
        int n1 = 501;
        int n2 = 501;

        float d1 = 1;
        float d2 = 1;
        float f1 = 0;
        float f2 = 0;

        float d_param= 1.0f; // parameter sampling rate
        float l_param = 20.0f; // last parameter
        float f_param = 1.0f;
        int n1_param = (int)((l_param-f_param+1)/d_param); // # of parameters
        int n2_param = 10;

        Sampling s1_param = new Sampling(n1_param,d_param,f_param);
        Sampling s2_param = new Sampling(n2_param,d_param,f_param);
        Sampling s1 = new Sampling(n1,d1,f1);
        Sampling s2 = new Sampling(n2,d2,f2);

        Stopwatch sw = new Stopwatch();
        sw.start();
        goLSF();
        sw.stop();
        System.out.println("Structure tensor time = "+sw.time());
        sw.restart();
        goPWDM();
        sw.stop();
        System.out.println("Madagascar PWD time = "+sw.time());
        sw.restart();
        goDW();
        sw.stop();
        System.out.println("Dynamic warping time = "+sw.time());
        sw.restart();
        goPWDD();
        sw.stop();
        System.out.println("Dave's PWD time = "+sw.time());

        //Util.testOptimalParameters(s1_param,s2_param,s1,s2,1,
        //    "Structure tensor");
        //Util.testOptimalParameters(s1_param,s2_param,s1,s2,2,
        //    "Madagascar PWD");
        //Util.testOptimalParameters(s1_param,s2_param,s1,s2,3,
        //    "Dynamic warping");
        //Util.testOptimalParameters(s1_param,s2_param,s1,s2,4,
        //    "Dave's PWD");
        goFandPk();
        //goDynamicWarpingSlopes();
        //goDpImages();
        //goTestSlopeVsError();
        //goRmsErrorCurves();
        //goTestSampleMeanSD();
        //goPlotSampleMeanSD();
        //goErrorLocation();
        //goTestSigmaValues();
        //goPlotSigmaValues();
        //goTestRectValues();
        //goPlotRectValues();
        //goTestSigmaValuesDW();
        //goPlotSigmaValuesDW();
        //goTestStrainValuesDW();
        //goPlotStrainValuesDW();
        //goTestDavePWD();
        //goPlotDavePWD();
        //goTestPmaxValues();
        //goGom();
      }
    });
  }
}
