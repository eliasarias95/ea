package slopes;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.mosaic.*;
import edu.mines.jtk.util.Stopwatch;

import static edu.mines.jtk.util.ArrayMath.*;

import dnp.LocalSlopeFinder;
import dnp.PlaneWaveDestructor;
import util.Plot;
import util.Util;
import util.FakeData;

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
      synthAndSlope = FakeData.seismicAndSlopes2d2014A(noise[i],F);
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
  }

  private static void trace(String s) {
    System.out.println(s);
  }

  private static final int norder = 15;
  private static final int k = 10;
  private static final float pmax = 5.0f;
  private static final float noise = 0.5f;
  private static final float freq = 0.1f;
  private static final float pc2   = -0.7f;//constant slope
  private static final float pc3   = 1.3f;//constant slope
  private static final String PATH = 
    "/Users/earias/Home/git/ea/bench/src/util/";
  private static final boolean T = true;
  private static final boolean F = false;    
  public static void main(String[] args) {
    SwingUtilities.invokeLater(new Runnable() {
      public void run() {     
        int n1 = 501;
        int n2 = 501;
        float[][] f2D  = new float[n2][n1];//2D seismic image
        float[][] p2D = new float[n2][n1]; //2D slope values
        Sampling s1 = new Sampling(n1);
        Sampling s2 = new Sampling(n2);
        Slopes s2D = new Slopes(noise,pmax,s1,s2);

        n1 = 101;
        n2 = 102;
        int n3 = 103;
        float[][][] f3  = new float[n3][n2][n1];//3D seismic image
        float[][][] p3  = new float[n3][n2][n1];//3D slope values 
        s1 = new Sampling(n1);
        s2 = new Sampling(n2);
        Sampling s3 = new Sampling(n3);
        Slopes s3D = new Slopes(noise,pmax,s1,s2,s3);
        //Slopes.makeSyntheticConstant(freq,pc2,pc3,f3,p3);
        //Slopes.makeSyntheticComplex(noise,f3,p3);

        /*
        n1 = 301;
        n2 = 920;
        float d1 = 0.004f;
        float d2 = .02667f;
        float f1 = 1.6f;
        float f2 = 0;
        f2D  = new float[n2][n1];//for GOM image
        s1 = new Sampling(n1,d1,f1);
        s2 = new Sampling(n2,d2,f2);
        Slopes.makeRealGOM(f2D);
        */

        //Slopes.makeSyntheticConstant(freq,pc2,f2D,p2D);
        Slopes.makeSyntheticComplex(noise,f2D,p2D);

        double dp = 1.0f;
        double fp = 1.0f;
        double lp = 80.0f;
        int np1 = (int)((lp)/dp);
        int np2 = 20;
        Sampling sp1 = new Sampling(np1,dp,fp);
        Sampling sp2 = new Sampling(np2,dp,fp);

        //Testing and plot optimal parameters
        //s2D.testOptimalSmoothLSF("Structure tensor",f2D,p2D,sp1,sp2);
        //s2D.plotOptimalParameters("Structure tensor","Sigma2","Sigma1",sp1,sp2);

        //s2D.testOptimalSmoothPWDM("Madagascar PWD",f2D,p2D,sp1,sp2);
        //s2D.plotOptimalParameters("Madagascar PWD","Rect2","Rect1",sp1,sp2);

        //s2D.testOptimalSmoothPWDD("Dave's PWD",f2D,p2D,sp1,sp2);
        //s2D.plotOptimalParameters("Dave's PWD","Smooth2","Smooth1",sp1,sp2);

        //s2D.testOptimalSmoothSDW("Smooth dynamic warping h",k,f2D,p2D,sp1,sp2);
        //s2D.plotOptimalParameters("Smooth dynamic warping h","h2","h1",sp1,sp2);

        dp = 0.1f; // strain sampling rate
        fp = 0.1f;
        lp = 2.0f;
        np1 = (int)((lp)/dp); //# strain vals to test
        np2 = np1;
        sp1 = new Sampling(np1,dp,fp);
        sp2 = new Sampling(np2,dp,fp);
        //s2D.testOptimalStrainSDW("Smooth dynamic warping r",k,f2D,p2D,sp1,sp2);
        s2D.plotOptimalParameters("Smooth dynamic warping r","r2","r1",sp1,sp2);

        //s2D.plotLSF(f2D,p2D);
        //s2D.plotLSF(f2D); //exact slopes unknown
        //s2D.plotPWDM(f2D,p2D);
        //s2D.plotPWDM(f2D); //exact slopes unknown
        //s2D.plotPWDD(f2D,p2D);
        //s2D.plotPWDD(f2D); //exact slopes unknown
        //s2D.plotDW(k,f2D,p2D);
        //s2D.plotDW(k,f2D); //exact slopes unknown
        //s2D.plotSDW(k,f2D,p2D);
        //s2D.plotSDW(k,f2D); //exact slopes unknown
        //s2D.testOrderVsTime(norder,f2D,p2D);
        //s2D.plotOrderVsTime(norder);
        //s2D.plotF(f2D); //plot synthetic seismic
        //s2D.plotP(p2D); //plot exact slopes

        //goDpImages();
        //goTestSlopeVsError();
        //goRmsErrorCurves();
        //goTestPmaxValues();
      }
    });
  }
}
