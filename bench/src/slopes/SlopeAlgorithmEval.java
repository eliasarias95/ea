package slopes;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.*;
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
 * @version 19.2.2015
 */

public class SlopeAlgorithmEval {

  /**
   * Runs the method of choice for the complex 2D synthetic where:
   * method == 1, for structure tensor
   * method == 2, for plane-wave destructor
   * method == 3, for smooth dynamic warping
   * method == 4, for all methods.
   */
  private static void complex2D(int method) {
    Check.argument(method==1 || method==2 || method==3 || method==4,
        "not valid input");
    int n1 = 501;
    int n2 = 501;
    float[][] f2D = new float[n2][n1];//2D seismic image
    float[][] p2D = new float[n2][n1]; //2D slope values
    float[][] r2D = new float[n2][n1]; //2D reflectivity 
    Sampling s1 = new Sampling(n1);
    Sampling s2 = new Sampling(n2);
    Slopes s2D = new Slopes(noise,pmax,s1,s2);
    Slopes.makeSyntheticComplex(noise,f2D,p2D,r2D);
    if (method==1) s2D.plotLSF(f2D,p2D);
    if (method==2) s2D.plotPWDM(f2D,p2D);
    if (method==3) s2D.plotSDW(k,f2D,p2D);
    if (method==4) {
      s2D.plotLSF(f2D,p2D);
      s2D.plotPWDM(f2D,p2D);
      s2D.plotSDW(k,f2D,p2D);
    }
    s2D.plotF("synthetic_noise="+noise,f2D); //plot synthetic seismic
    s2D.plotF("reflectivity_noise="+noise,r2D); //plot reflectivity 
    s2D.plotP(p2D); //plot exact slopes
  }

  /**
   * Runs the method of choice for the complex 3D synthetic where:
   * method == 1, for structure tensor
   * method == 2, for plane-wave destructor
   * method == 3, for smooth dynamic warping
   * method == 4, for all methods.
   */
  private static void complex3D(int method) {
    Check.argument(method==1 || method==2 || method==3 || method==4,
        "not valid input");
    int n1 = 101;
    int n2 = 102;
    int n3 = 103;
    float[][][] f3  = new float[n3][n2][n1];//3D seismic image
    float[][][] p3  = new float[n3][n2][n1];//3D slope values 
    Sampling s1 = new Sampling(n1);
    Sampling s2 = new Sampling(n2);
    Sampling s3 = new Sampling(n3);
    Slopes s3D = new Slopes(noise,pmax,s1,s2,s3);
    //Slopes.makeSyntheticConstant(freq,pc2,pc3,f3,p3);
    Slopes.makeSyntheticComplex(noise,f3,p3);
  }

  /**
   * Runs the method of choice for the constant 2D synthetic where:
   * method == 1, for structure tensor
   * method == 2, for plane-wave destructor
   * method == 3, for smooth dynamic warping
   * method == 4, for all methods.
   */
  private static void constant2D(int method) {
    Check.argument(method==1 || method==2 || method==3 || method==4,
        "not valid input");
    int n1 = 501;
    int n2 = 501;
    float[][] f2D  = new float[n2][n1];//2D seismic image
    float[][] p2D = new float[n2][n1]; //2D slope values
    Sampling s1 = new Sampling(n1);
    Sampling s2 = new Sampling(n2);
    Slopes s2D = new Slopes(noise,pmax,s1,s2);
    Slopes.makeSyntheticConstant(freq,pc2,f2D,p2D);
    if (method==1) s2D.plotLSF(f2D,p2D);
    if (method==2) s2D.plotPWDM(f2D,p2D);
    if (method==3) s2D.plotSDW(k,f2D,p2D);
    if (method==4) {
      s2D.plotLSF(f2D,p2D);
      s2D.plotPWDM(f2D,p2D);
      s2D.plotSDW(k,f2D,p2D);
    }
    s2D.plotF("constant2D_p2="+pc2,f2D); //plot synthetic seismic
    s2D.plotP(p2D); //plot exact slopes
  }

  /**
   * Runs the method of choice for the GOM real data where:
   * method == 1, for structure tensor
   * method == 2, for plane-wave destructor
   * method == 3, for smooth dynamic warping
   * method == 4, for all methods.
   */
  private static void GOM(int method) {
    Check.argument(method==1 || method==2 || method==3 || method==4,
        "not valid input");
    int n1 = 301;
    int n2 = 920;
    float d1 = 0.004f;
    float d2 = .02667f;
    float f1 = 1.6f;
    float f2 = 0;
    float[][] f2D  = new float[n2][n1];//for GOM image
    Sampling s1 = new Sampling(n1,d1,f1);
    Sampling s2 = new Sampling(n2,d2,f2);
    //Sampling s1 = new Sampling(n1);
    //Sampling s2 = new Sampling(n2);
    Slopes s2D = new Slopes(noise,pmax,s1,s2);
    Slopes.makeRealGOM(f2D);
    if (method==1) s2D.plotLSF(f2D);
    if (method==2) s2D.plotPWDM(f2D);
    if (method==3) s2D.plotSDW(k,f2D);
    if (method==4) {
      s2D.plotLSF(f2D);
      s2D.plotPWDM(f2D);
      s2D.plotSDW(k,f2D);
    }
    s2D.plotF("GOM2D",f2D); //plot seismic
  }

  private static void optimal(int method, boolean test, double lp, int np2) {
    Check.argument(method==1 || method==2 || method==3 || method==4,
        "not valid input");
    int n1 = 501;
    int n2 = 501;
    float[][] f2D = new float[n2][n1]; //2D seismic image
    float[][] p2D = new float[n2][n1]; //2D slope values
    float[][] r2D = new float[n2][n1]; //2D reflectivity 
    Sampling s1 = new Sampling(n1);
    Sampling s2 = new Sampling(n2);
    Slopes s2D = new Slopes(noise,pmax,s1,s2);
    Slopes.makeSyntheticComplex(noise,f2D,p2D,r2D);

    double dp = 1.0f;
    double fp = 1.0f;
    int np1 = (int)((lp)/dp);
    Sampling sp1 = new Sampling(np1,dp,fp);
    Sampling sp2 = new Sampling(np2,dp,fp);
    String ttl;
    if (test) {
      if (method==1) {
        ttl = "Structure tensor noise="+noise;
        s2D.testOptimalSmoothLSF(ttl,f2D,p2D,sp1,sp2);
        s2D.plotOptimalParameters(ttl,"sigma2","sigma1",sp1,sp2);
      }
      if (method==2) {
        ttl = "Madagascar PWD noise="+noise;
        s2D.testOptimalSmoothPWDM(ttl,f2D,p2D,sp1,sp2);
        s2D.plotOptimalParameters(ttl,"rect2","rect1",sp1,sp2);
      }
      if (method==3) {
        ttl = "Smooth dynamic warping noise="+noise;
        s2D.testOptimalSmoothSDW(ttl,k,f2D,p2D,sp1,sp2);
        s2D.plotOptimalParameters(ttl,"h2","h1",sp1,sp2);
        dp = 0.1f; // strain sampling rate
        fp = 0.1f;
        np1 = (int)((2.0)/dp); //# strain vals to test
        sp1 = new Sampling(np1,dp,fp);
        sp2 = new Sampling(np1,dp,fp);
        ttl = "Smooth dynamic warping r noise="+noise;
        s2D.testOptimalStrainSDW(ttl,k,f2D,p2D,sp1,sp2);
        s2D.plotOptimalParameters(ttl,"r2","r1",sp1,sp2);
      }
      if (method==4) {
        ttl = "Structure tensor noise="+noise;
        s2D.testOptimalSmoothLSF(ttl,f2D,p2D,sp1,sp2);
        s2D.plotOptimalParameters(ttl,"sigma2","sigma1",sp1,sp2);
        ttl = "Madagascar PWD noise="+noise;
        s2D.testOptimalSmoothPWDM(ttl,f2D,p2D,sp1,sp2);
        s2D.plotOptimalParameters(ttl,"rect2","rect1",sp1,sp2);
        ttl = "Smooth dynamic warping noise="+noise;
        s2D.testOptimalSmoothSDW(ttl,k,f2D,p2D,sp1,sp2);
        s2D.plotOptimalParameters(ttl,"h2","h1",sp1,sp2);
        dp = 0.1f; // strain sampling rate
        fp = 0.1f;
        np1 = (int)((2.0)/dp); //# strain vals to test
        sp1 = new Sampling(np1,dp,fp);
        sp2 = new Sampling(np1,dp,fp);
        ttl = "Smooth dynamic warping r noise="+noise;
        s2D.testOptimalStrainSDW(ttl,k,f2D,p2D,sp1,sp2);
        s2D.plotOptimalParameters(ttl,"r2","r1",sp1,sp2);
      }
    }
    else {
      if (method==1) {
        ttl = "Structure tensor noise="+noise;
        s2D.plotOptimalParameters(ttl,"sigma2","sigma1",sp1,sp2);
      }
      if (method==2) {
        ttl = "Madagascar PWD noise="+noise;
        s2D.plotOptimalParameters(ttl,"rect2","rect1",sp1,sp2);
      }
      if (method==3) {
        ttl = "Smooth dynamic warping noise="+noise;
        s2D.plotOptimalParameters(ttl,"h2","h1",sp1,sp2);
        dp = 0.1f; // strain sampling rate
        fp = 0.1f;
        np1 = (int)((2.0)/dp); //# strain vals to test
        sp1 = new Sampling(np1,dp,fp);
        sp2 = new Sampling(np1,dp,fp);
        ttl = "Smooth dynamic warping r noise="+noise;
        s2D.plotOptimalParameters(ttl,"r2","r1",sp1,sp2);
      }
      if (method==4) {
        ttl = "Structure tensor noise="+noise;
        s2D.plotOptimalParameters(ttl,"sigma2","sigma1",sp1,sp2);
        ttl = "Madagascar PWD noise="+noise;
        s2D.plotOptimalParameters(ttl,"rect2","rect1",sp1,sp2);
        ttl = "Smooth dynamic warping noise="+noise;
        s2D.plotOptimalParameters(ttl,"h2","h1",sp1,sp2);
        dp = 0.1f; // strain sampling rate
        fp = 0.1f;
        np1 = (int)((2.0)/dp); //# strain vals to test
        sp1 = new Sampling(np1,dp,fp);
        sp2 = new Sampling(np1,dp,fp);
        ttl = "Smooth dynamic warping r noise="+noise;
        s2D.plotOptimalParameters(ttl,"r2","r1",sp1,sp2);
      }
    }
  }

  private static void errorCurves(int method, boolean test) {
    Check.argument(method==1 || method==2 || method==3 || method==4,
        "not valid input");
    int n1 = 501;
    int n2 = 501;
    float[][] f2D = new float[n2][n1];//2D seismic image
    float[][] p2D = new float[n2][n1]; //2D slope values
    float[][] r2D = new float[n2][n1]; //2D reflectivity 
    Sampling s1 = new Sampling(n1);
    Sampling s2 = new Sampling(n2);
    Slopes s2D = new Slopes(noise,pmax,s1,s2);
    Slopes.makeSyntheticComplex(noise,f2D,p2D,r2D);
    double dn = 0.05;
    int nn = 21;
    Sampling sn = new Sampling(nn,dn,0.0);
    float sigma1,sigma2;
    int rect1,rect2;
    double r1,r2,h1,h2;

    if (noise==0.0) {
      sigma1=14.0f; sigma2=1.0f;
      rect1=34; rect2=2;
      h1=30.0; h2=9.0;
      r1=0.1; r2=0.6;
    }

    else {
      sigma1=23.0f; sigma2=1.0f;
      rect1=75; rect2=6;
      h1=72.0; h2=12.0;
      r1=0.1; r2=0.4;
    }

    if (test) {
      if (method==1) s2D.testRmsErrorCurveLSF(sn,sigma1,sigma2);
      if (method==2) s2D.testRmsErrorCurvePWD(sn,rect1,rect2);
      if (method==3) s2D.testRmsErrorCurveSDW(sn,k,r1,r2,h1,h2);
      if (method==4) {
        s2D.testRmsErrorCurveLSF(sn,sigma1,sigma2);
        s2D.testRmsErrorCurvePWD(sn,rect1,rect2);
        s2D.testRmsErrorCurveSDW(sn,k,r1,r2,h1,h2);
      }
    }
    s2D.plotRmsErrorCurves(sn,1);
  }

  private static void trace(String s) {
    System.out.println(s);
  }

  private static final int norder = 15;
  private static final int k = 10;
  private static final float pmax = 9.0f;
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
        float[][] f2D = new float[n2][n1];//2D seismic image
        float[][] p2D = new float[n2][n1]; //2D slope values
        float[][] r2D = new float[n2][n1]; //2D reflectivity 
        Sampling s1 = new Sampling(n1);
        Sampling s2 = new Sampling(n2);
        Slopes s2D = new Slopes(noise,pmax,s1,s2);
        Slopes.makeSyntheticComplex(noise,f2D,p2D,r2D);

        float sigma1 = 23.0f;
        float sigma2 = 1.0f;
        int N = 100;
        double dn = 0.05;
        int n = 21;
        Sampling sn = new Sampling(n,dn,0.0);
        //s2D.testMeanErrorCurveLSF(sn,N,sigma1,sigma2);
        //s2D.plotRmsErrorCurves(sn,N);
        //s2D.testSampleMeanLSF(nn);
        //s2D.testSampleMeanPWD(nn);
        //s2D.testSampleMeanSDW(nn,k);
        //s2D.testSampleStdDevLSF(nn);
        //s2D.meanErrorCurveLSF();
        //s2D.meanErrorCurvePWD();
        //s2D.meanErrorCurveSDW();
      
        //1=lsf  2=pwd  3=sdw  4=all
        complex2D(1);
        //1=lsf  2=pwd  3=sdw  4=all
        //constant2D(4);
        //1=lsf  2=pwd  3=sdw  4=all
        //GOM(1);

        //1=lsf  2=pwd  3=sdw  4=all, test?
        //errorCurves(4,true);
        //1=lsf  2=pwd  3=sdw  4=all, test?, n1, n2
        //optimal(4,true,40,15);

        //s2D.testOrderVsTime(norder,f2D,p2D);
        //s2D.plotOrderVsTime(norder);

        //goDpImages();
        //goTestSlopeVsError();
        //goRmsErrorCurves();
        //goTestPmaxValues();
      }
    });
  }
}
