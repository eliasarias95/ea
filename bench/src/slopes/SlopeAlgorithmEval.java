package slopes;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.mosaic.*;
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
 * @version 19.1.2015
 */

public class SlopeAlgorithmEval {

  private static void timeAndPlotLSF(boolean error) {
    _sw.start();
    Slopes.plotLSF(error);
    _sw.stop();
    System.out.println("Structure tensor time = "+_sw.time());
  }

  private static void timeAndPlotPWDM(boolean error) {
    _sw.restart();
    Slopes.plotPWDM(error);
    _sw.stop();
    System.out.println("Madagascar PWD time = "+_sw.time());
  }

  private static void timeAndPlotDW(boolean error, int k) {
    _sw.restart();
    Slopes.plotDW(error,k);
    _sw.stop();
    System.out.println("Dynamic warping time = "+_sw.time());
  }

  private static void timeAndPlotSDW(boolean error, int k) {
    _sw.restart();
    Slopes.plotSDW(error,k);
    _sw.stop();
    System.out.println("Smooth dynamic warping time = "+_sw.time());
  }

  private static void timeAndPlotPWDD(boolean error) {
    _sw.restart();
    Slopes.plotPWDD(error);
    _sw.stop();
    System.out.println("Dave's PWD time = "+_sw.time());
  }
/////////////////////////////NOT REFACTORED////////////////////////////  

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
      //dw_slope = Slopes.DWSlopesAvg(dw);
      lsf_slope = mul(lsf_slope,d1/d2);
      mad_slope = mul(mad_slope,d1/d2);
      //dw_slope = mul(dw_slope,d1/d2);
      lsf_diff = sub(lsf_slope,exact_slope);
      mad_diff = sub(mad_slope,exact_slope);
      //dw_diff = sub(dw_slope,exact_slope);
      temp_lsf = pow(lsf_diff,2);
      temp_mad = pow(mad_diff,2);
      //temp_dw = pow(dw_diff,2);
      rmserr_lsf[i] = sqrt(sum(temp_lsf)/(n2*n1));
      rmserr_mad[i] = sqrt(sum(temp_mad)/(n2*n1));
      //rmserr_dw[i] = sqrt(sum(temp_dw)/(n2*n1));
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
  }

  private static Stopwatch _sw = new Stopwatch();
  private static final boolean T = true;
  private static final boolean F = false;    
  public static void main(String[] args) {
    SwingUtilities.invokeLater(new Runnable() {
      public void run() {     
        float pk = -0.7f;
        float noise = 0.5f;
        float l_param = 10.0f;
        boolean error = true;
        int norder = 10;
        int k = 5;

        //set parameters for testing
        Slopes.setSynthParameters(noise);
        //Slopes.setChickenTestParameters(pk);
        //Slopes.setGOMParameters();
        
        //time methods and plot estimated slopes
        //timeAndPlotLSF(error);
        //timeAndPlotPWDM(error);
        timeAndPlotDW(error,k);
        //timeAndPlotSDW(error,k);
        //timeAndPlotPWDD(error);
        //Slopes.testOrderVsTime(norder);
        //Slopes.plotOrderVsTime(norder);
        //Slopes.plotFandPk(); //plot synthetic seismic and exact slopes

        //Testing and plot optimal parameters
        //Slopes.setSmoothingParameters(l_param);
        //Slopes.testOptimalParameters(1,"Structure tensor");
        //Slopes.plotOptimalParameters("Structure tensor","Sigma2","Sigma1");

        //Slopes.testOptimalParameters(2,"Madagascar PWD");
        //Slopes.plotOptimalParameters("Madagascar PWD","Rect2","Rect1");

        //Slopes.testOptimalParameters(3,"Dynamic warping smooth");
        //Slopes.plotOptimalParameters("Dynamic warping smooth","Usmooth2",
        //    "Usmooth1");

        //Slopes.testOptimalParameters(5,"Dave's PWD");
        //Slopes.plotOptimalParameters("Dave's PWD","Smooth2","Smooth1");

        //Slopes.setStrainParameters();
        //Slopes.testOptimalParameters(4,"Dynamic warping strain");
        //Slopes.plotOptimalParameters("Dynamic warping strain","Strain2",
        //    "Strain1");

        //goDpImages();
        //goTestSlopeVsError();
        //goRmsErrorCurves();
        //goTestPmaxValues();
      }
    });
  }
}
