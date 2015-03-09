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
    String hl = "Traces";
    String vl = "Samples";
    s2D.plotF("synthetic_noise="+noise,hl,vl,f2D); //plot synthetic seismic
    s2D.plotF("reflectivity_noise="+noise,hl,vl,r2D); //plot reflectivity 
    s2D.plotP(p2D); //plot exact slopes
    //s2D.plotF("flipped_synthetic_noise="+noise,hl,vl,flip2(f2D)); //plot synthetic seismic
    //s2D.plotP(flip2(p2D)); //plot exact slopes
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
    float[][][] p2  = new float[n3][n2][n1];//3D slope values inline 
    float[][][] p3  = new float[n3][n2][n1];//3D slope values xline
    Sampling s1 = new Sampling(n1);
    Sampling s2 = new Sampling(n2);
    Sampling s3 = new Sampling(n3);
    Slopes s3D = new Slopes(noise,pmax,s1,s2,s3);
    Slopes.makeSyntheticConstant(freq,pc2,pc3,f3,p2,p3);
    s3D.plot3D(f3);
    //Slopes.makeSyntheticComplex(noise,f3,p3);
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
    double d1 = 0.5;
    double d2 = 0.25;
    float[][] f2D  = new float[n2][n1];//2D seismic image
    float[][] p2D = new float[n2][n1]; //2D slope values
    Sampling s1 = new Sampling(n1,d1,0.0);
    Sampling s2 = new Sampling(n2,d2,0.0);
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
    String hl = "Traces";
    String vl = "Samples";
    s2D.plotF("constant2D_p2="+pc2,hl,vl,f2D); //plot synthetic seismic
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
    s1 = new Sampling(n1);
    s2 = new Sampling(n2);
    Slopes s2D = new Slopes(noise,pmax,s1,s2);
    Slopes.makeRealGOM(f2D);
    String title = "GOM2D";
    String hl = "Distance (km)";
    String vl = "Time (s)";
    if (method==1) s2D.plotLSF(f2D,title,hl,vl);
    if (method==2) s2D.plotPWDM(f2D,title,hl,vl);
    if (method==3) s2D.plotSDW(k,f2D,title,hl,vl);
    if (method==4) {
      s2D.plotLSF(f2D,title,hl,vl);
      s2D.plotPWDM(f2D,title,hl,vl);
      s2D.plotSDW(k,f2D,title,hl,vl);
    }
    s2D.plotF(title,hl,vl,f2D); //plot seismic
  }

  /**
   * Runs the method of choice for the GOM real data where:
   * method == 1, for structure tensor
   * method == 2, for plane-wave destructor
   * method == 3, for smooth dynamic warping
   * method == 4, for all methods.
   */
  private static void TP(int method) {
    Check.argument(method==1 || method==2 || method==3 || method==4,
        "not valid input");
    int n1 = 251;
    int n2 = 357;
    float d1 = 0.004f;
    float d2 = .025f;
    float f1 = 0.5f;
    float f2 = 0;
    float[][] f2D  = new float[n2][n1];//for GOM image
    Sampling s1 = new Sampling(n1,d1,f1);
    Sampling s2 = new Sampling(n2,d2,f2);
    //s1 = new Sampling(n1);
    //s2 = new Sampling(n2);
    Slopes s2D = new Slopes(noise,pmax,s1,s2);
    Slopes.makeRealTp(f2D);
    String title = "TP2D";
    String hl = "Distance (km)";
    String vl = "Time (s)";
    if (method==1) s2D.plotLSF(f2D,title,hl,vl);
    if (method==2) s2D.plotPWDM(f2D,title,hl,vl);
    if (method==3) s2D.plotSDW(k,f2D,title,hl,vl);
    if (method==4) {
      s2D.plotLSF(f2D,title,hl,vl);
      s2D.plotPWDM(f2D,title,hl,vl);
      s2D.plotSDW(k,f2D,title,hl,vl);
    }
    s2D.plotF(title,hl,vl,f2D); //plot seismic
  }

  private static void optimal(int method, boolean test, double lp, int np2) {
    Check.argument(method==1 || method==2 || method==3 || method==4,
        "not valid input");
    if (noise==0.0f)      num = 1;
    else if (noise==0.5f) num = 2;
    else if (noise==1.0f) num = 3;
    else                  num = 4;    
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
    float cmin_lsf = 0.4f;
    float cmax_lsf = 0.7f;
    float cmin_pwd = 0.63f;
    float cmax_pwd = 0.7f;
    float cmin_sdw = 0.58f;
    float cmax_sdw = 0.7f;
    float cmin_sdwr = 0.55f;
    float cmax_sdwr = 0.6f;
    if (test) {
      if (method==1) {
        ttl = "structure_tensor_errors_"+num+".dat";
        s2D.testOptimalSmoothLSF(ttl,f2D,p2D,sp1,sp2);
        s2D.plotOptimalParameters(ttl,"sigma2","sigma1",sp1,sp2,cmin_lsf,cmax_lsf);
      }
      if (method==2) {
        ttl = "madagascar_pwd_errors_"+num+".dat";
        s2D.testOptimalSmoothPWDM(ttl,f2D,p2D,sp1,sp2);
        s2D.plotOptimalParameters(ttl,"rect2","rect1",sp1,sp2,cmin_pwd,cmax_pwd);
      }
      if (method==3) {
        ttl = "smooth_dynamic_warping_h_errors_"+num+".dat";
        s2D.testOptimalSmoothSDW(ttl,k,f2D,p2D,sp1,sp2);
        s2D.plotOptimalParameters(ttl,"h2","h1",sp1,sp2,cmin_sdw,cmax_sdw);
        dp = 0.1f; // strain sampling rate
        fp = 0.1f;
        np1 = (int)((2.0)/dp); //# strain vals to test
        sp1 = new Sampling(np1,dp,fp);
        sp2 = new Sampling(np1,dp,fp);
        ttl = "smooth_dynamic_warping_r_errors_"+num+".dat";
        s2D.testOptimalStrainSDW(ttl,k,f2D,p2D,sp1,sp2);
        s2D.plotOptimalParameters(ttl,"r2","r1",sp1,sp2,cmin_sdwr,cmax_sdwr);
      }
      if (method==4) {
        ttl = "structure_tensor_errors_"+num+".dat";
        s2D.testOptimalSmoothLSF(ttl,f2D,p2D,sp1,sp2);
        s2D.plotOptimalParameters(ttl,"sigma2","sigma1",sp1,sp2,cmin_lsf,cmax_lsf);
        ttl = "madagascar_pwd_errors_"+num+".dat";
        s2D.testOptimalSmoothPWDM(ttl,f2D,p2D,sp1,sp2);
        s2D.plotOptimalParameters(ttl,"rect2","rect1",sp1,sp2,cmin_pwd,cmax_pwd);
        ttl = "smooth_dynamic_warping_h_errors_"+num+".dat";
        s2D.testOptimalSmoothSDW(ttl,k,f2D,p2D,sp1,sp2);
        s2D.plotOptimalParameters(ttl,"h2","h1",sp1,sp2,cmin_sdw,cmax_sdw);
        dp = 0.1f; // strain sampling rate
        fp = 0.1f;
        np1 = (int)((2.0)/dp); //# strain vals to test
        sp1 = new Sampling(np1,dp,fp);
        sp2 = new Sampling(np1,dp,fp);
        ttl = "smooth_dynamic_warping_r_errors_"+num+".dat";
        s2D.testOptimalStrainSDW(ttl,k,f2D,p2D,sp1,sp2);
        s2D.plotOptimalParameters(ttl,"r2","r1",sp1,sp2,cmin_sdwr,cmax_sdwr);
      }
    }
    else {
      if (method==1) {
        ttl = "structure_tensor_errors_"+num+".dat";
        s2D.plotOptimalParameters(ttl,"sigma2","sigma1",sp1,sp2,cmin_lsf,cmax_lsf);
      }
      if (method==2) {
        ttl = "madagascar_pwd_errors_"+num+".dat";
        s2D.plotOptimalParameters(ttl,"rect2","rect1",sp1,sp2,cmin_pwd,cmax_pwd);
      }
      if (method==3) {
        ttl = "smooth_dynamic_warping_h_errors_"+num+".dat";
        s2D.plotOptimalParameters(ttl,"h2","h1",sp1,sp2,cmin_sdw,cmax_sdw);
        dp = 0.1f; // strain sampling rate
        fp = 0.1f;
        np1 = (int)((2.0)/dp); //# strain vals to test
        sp1 = new Sampling(np1,dp,fp);
        sp2 = new Sampling(np1,dp,fp);
        ttl = "smooth_dynamic_warping_r_errors_"+num+".dat";
        s2D.plotOptimalParameters(ttl,"r2","r1",sp1,sp2,cmin_sdwr,cmax_sdwr);
      }
      if (method==4) {
        ttl = "structure_tensor_errors_"+num+".dat";
        s2D.plotOptimalParameters(ttl,"sigma2","sigma1",sp1,sp2,cmin_lsf,cmax_lsf);
        ttl = "madagascar_pwd_errors_"+num+".dat";
        s2D.plotOptimalParameters(ttl,"rect2","rect1",sp1,sp2,cmin_pwd,cmax_pwd);
        ttl = "smooth_dynamic_warping_h_errors_"+num+".dat";
        s2D.plotOptimalParameters(ttl,"h2","h1",sp1,sp2,cmin_sdw,cmax_sdw);
        dp = 0.1f; // strain sampling rate
        fp = 0.1f;
        np1 = (int)((2.0)/dp); //# strain vals to test
        sp1 = new Sampling(np1,dp,fp);
        sp2 = new Sampling(np1,dp,fp);
        ttl = "smooth_dynamic_warping_r_errors_"+num+".dat";
        s2D.plotOptimalParameters(ttl,"r2","r1",sp1,sp2,cmin_sdwr,cmax_sdwr);
      }
    }
  }

  private static void rmsErrorCurves(int method, boolean test) {
    Check.argument(method==1 || method==2 || method==3 || method==4,
        "not valid input");
    if (noise==0.0f)      num = 1;
    else if (noise==0.5f) num = 2;
    else if (noise==1.0f) num = 3;
    else                  num = 4;    
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
    String fileName1 = "data/rms_error_lsf_"+num+".dat";
    String fileName2 = "data/rms_error_pwd_"+num+".dat";
    String fileName3 = "data/rms_error_sdw_"+num+".dat";
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
      if (method==1) s2D.testRmsErrorCurveLSF(sn,fileName1,sigma1,sigma2);
      if (method==2) s2D.testRmsErrorCurvePWD(sn,fileName2,rect1,rect2);
      if (method==3) s2D.testRmsErrorCurveSDW(sn,fileName3,k,r1,r2,h1,h2);
      if (method==4) {
        s2D.testRmsErrorCurveLSF(sn,fileName1,sigma1,sigma2);
        s2D.testRmsErrorCurvePWD(sn,fileName2,rect1,rect2);
        s2D.testRmsErrorCurveSDW(sn,fileName3,k,r1,r2,h1,h2);
      }
    }
    String title = "pwd_lsf_sdw_rmserror_vs_nsratio"+num;
    String hl = "Noise/signal";
    String vl = "RMS error";
    s2D.plotCurve(sn,fileName1,fileName2,fileName3,title,hl,vl,0,1);
  }

  private static void meanErrorCurves(int method, boolean test, int nni) {
    Check.argument(method==1 || method==2 || method==3 || method==4,
        "not valid input");
    if (noise==0.0f)      num = 1;
    else if (noise==0.5f) num = 2;
    else if (noise==1.0f) num = 3;
    else                  num = 4;    
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
    String fileName1 = "data/mean_error_lsf_"+nni+"_"+num+".dat";
    String fileName2 = "data/mean_error_pwd_"+nni+"_"+num+".dat";
    String fileName3 = "data/mean_error_sdw_"+nni+"_"+num+".dat";
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
      if (method==1) s2D.testMeanCurveLSF(sn,fileName1,nni,sigma1,sigma2);
      if (method==2) s2D.testMeanCurvePWD(sn,fileName2,nni,rect1,rect2);
      if (method==3) s2D.testMeanCurveSDW(sn,fileName3,nni,k,r1,r2,h1,h2);
      if (method==4) {
        s2D.testMeanCurveLSF(sn,fileName1,nni,sigma1,sigma2);
        s2D.testMeanCurvePWD(sn,fileName2,nni,rect1,rect2);
        s2D.testMeanCurveSDW(sn,fileName3,nni,k,r1,r2,h1,h2);
      }
    }
    String title = "pwd_lsf_sdw_mean_error_vs_nsratio"+nni;
    String hl = "Noise/signal";
    String vl = "Mean error";
    s2D.plotCurve(sn,fileName1,fileName2,fileName3,title,hl,vl,0,0.02f);
  }

  private static void mean() {

  }

  private static void stdDev() {
    
  }

  private static void trace(String s) {
    System.out.println(s);
  }

  private static int num;
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

        n1 = 2762;
        n2 = 357;
        int n3 = 161;
        double d1 = .002;
        double d2 = .025;
        double d3 = .025;
        s1 = new Sampling(n1,d1,0.0);
        s2 = new Sampling(n2,d2,0.0);
        Sampling s3 = new Sampling(n3,d3,0.0);
        Slopes s3D = new Slopes(noise,pmax,s1,s2,s3);
        //float[][][] teapot = Util.readImage(n1,n2,n3,PATH+"data/tp/tpsz.dat");
        //float[][][] teapot = Util.readImage(n1,n2,n3,PATH+
        //    "data/tp/tpsz_subz_401_4_400.dat");
        //s3D.plot3D(teapot);

        float sigma1 = 23.0f; float sigma2 = 1.0f;
        int rect1 = 34;       int rect2 = 2;
        double h1 = 30.0;     double h2 = 9.0;
        double r1 = 0.1;      double r2 = 0.6;

        int nni = 100;
        double dn = 0.05;
        int n = 21;
        Sampling sn = new Sampling(n,dn,0.0);
        String fn_lsf = "data/mean_lsf.dat";
        String fn_pwd = "data/mean_pwd.dat";
        String fn_sdw = "data/mean_sdw.dat";

        String title = "lsf_mean"+nni;
        String hl = "Traces";
        String vl = "Samples";
        String cbl = "std (samples/trace)";
        //s2D.testSampleMeanLSF(fn_lsf,nni,noise,sigma1,sigma2);
        //s2D.testSampleMeanPWD(fn_pwd,nni,noise,rect1,rect2);
        //s2D.testSampleMeanSDW(fn_sdw,nni,noise,k,r1,r2,h1,h2);
        //s2D.plotImage(title,hl,vl,cbl,fn_lsf);
        title = "pwd_mean"+nni;
        //s2D.plotImage(title,hl,vl,cbl,fn_pwd);
        title = "sdw_mean"+nni;
        //s2D.plotImage(title,hl,vl,cbl,fn_sdw);

        fn_lsf = "data/stddev_lsf.dat";
        fn_pwd = "data/stddev_pwd.dat";
        fn_sdw = "data/stddev_sdw.dat";
        //s2D.testSampleStdDevLSF(fn_lsf,nni,sigma1,sigma2);
        //s2D.testSampleStdDevPWD(fn_pwd,nni,rect1,rect2);
        //s2D.testSampleStdDevSDW(fn_sdw,nni,k,r1,r2,h1,h2);
        title = "lsf_stddev"+nni;
        s2D.plotImage(title,hl,vl,cbl,fn_lsf);
        title = "pwd_stddev"+nni;
        s2D.plotImage(title,hl,vl,cbl,fn_pwd);
        title = "sdw_stddev"+nni;
        s2D.plotImage(title,hl,vl,cbl,fn_sdw);

        //complex3D(1);
      
        //1=lsf  2=pwd  3=sdw  4=all
        //complex2D(4);
        //1=lsf  2=pwd  3=sdw  4=all
        //constant2D(1);
        //1=lsf  2=pwd  3=sdw  4=all
        //GOM(4);
        //1=lsf  2=pwd  3=sdw  4=all
        //TP(4);

        //1=lsf  2=pwd  3=sdw  4=all, test?
        //rmsErrorCurves(4,false);
        //1=lsf  2=pwd  3=sdw  4=all, test?
        //meanErrorCurves(4,false,nni);
        //1=lsf  2=pwd  3=sdw  4=all, test?, n1, n2
        //optimal(4,false,80,15);

        Sampling sovt = new Sampling(norder);
        String fn = "data/order_vs_time_"+norder+".dat";
        title = "order_vs_time";
        hl = "Order";
        vl = "Time";
        //s2D.testOrderVsTime(sovt,fn,f2D,p2D);
        //s2D.plotCurve(sovt,fn,title,hl,vl,0,78);

        //goDpImages();
        //goTestSlopeVsError();
        //goRmsErrorCurves();
        //goTestPmaxValues();
      }
    });
  }
}
