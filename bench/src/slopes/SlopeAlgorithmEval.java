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
    float[][] f = new float[n2][n1];//2D seismic image
    float[][] p = new float[n2][n1]; //2D slope values
    float[][] r = new float[n2][n1]; //2D reflectivity 
    Sampling s1 = new Sampling(n1);
    Sampling s2 = new Sampling(n2);
    Slopes s2D = new Slopes(noise,pmax,s1,s2);
    Slopes.makeSyntheticComplex(noise,f,p,r);
    String lsf_title = "lsf_complex2D";
    String pwd_title = "pwd_complex2D";
    String sdw_title = "sdw_complex2D";
    if (method==1) {
      s2D.estimateLSF(f,p,lsf_title);
      s2D.plot2D(f,lsf_title);
    }
    if (method==2) {
      s2D.estimatePWDM(f,p,pwd_title);
      s2D.plot2D(f,pwd_title);
    }
    if (method==3) {
      s2D.estimateSDW(k,f,p,sdw_title);
      s2D.plot2D(f,sdw_title);
    }
    if (method==4) {
      //s2D.estimateLSF(f,p,lsf_title);
      s2D.plot2D(f,lsf_title);
      //s2D.estimatePWDM(f,p,pwd_title);
      s2D.plot2D(f,pwd_title);
      //s2D.estimateSDW(k,f,p,sdw_title);
      s2D.plot2D(f,sdw_title);
    }
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
    float[][][] f   = new float[n3][n2][n1];//3D seismic image
    float[][][] p2  = new float[n3][n2][n1];//3D slope values inline 
    float[][][] p3  = new float[n3][n2][n1];//3D slope values xline
    Sampling s1 = new Sampling(n1);
    Sampling s2 = new Sampling(n2);
    Sampling s3 = new Sampling(n3);
    Slopes s3D = new Slopes(noise,pmax,s1,s2,s3);
    String lsf_title = "lsf_complex3D";
    String pwd_title = "pwd_complex3D";
    String sdw_title = "sdw_complex3D";
    Slopes.makeSyntheticComplex(noise,f,p2,p3);
    if (method==1) {
      s3D.estimateLSF(f,p2,p3,lsf_title);
      s3D.plot3D(f,lsf_title);
    }
    if (method==2) {
      s3D.estimatePWDM(f,p2,p3,pwd_title);
      s3D.plot3D(f,pwd_title);
    }
    if (method==3) {
      s3D.estimateSDW(k,f,p2,p3,sdw_title);
      s3D.plot3D(f,sdw_title);
    }
    if (method==4) {
      s3D.estimateLSF(f,p2,p3,lsf_title);
      s3D.plot3D(f,lsf_title);
      s3D.estimatePWDM(f,p2,p3,pwd_title);
      s3D.plot3D(f,pwd_title);
      s3D.estimateSDW(k,f,p2,p3,sdw_title);
      s3D.plot3D(f,sdw_title);
    }
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
    float[][] f  = new float[n2][n1];//2D seismic image
    float[][] p = new float[n2][n1]; //2D slope values
    Sampling s1 = new Sampling(n1);
    Sampling s2 = new Sampling(n2);
    Slopes s2D = new Slopes(noise,pmax,s1,s2);
    Slopes.makeSyntheticConstant(freq,pc2,f,p);
    String lsf_title = "lsf_constant2D";
    String pwd_title = "pwd_constant2D";
    String sdw_title = "sdw_constant2D";
    if (method==1) {
      s2D.estimateLSF(f,p,lsf_title);
      s2D.plot2D(f,lsf_title);
    }
    if (method==2) {
      s2D.estimatePWDM(f,p,pwd_title);
      s2D.plot2D(f,pwd_title);
    }
    if (method==3) {
      s2D.estimateSDW(k,f,p,sdw_title);
      s2D.plot2D(f,sdw_title);
    }
    if (method==4) {
      //s2D.estimateLSF(f,p,lsf_title);
      s2D.plot2D(f,lsf_title);
      //s2D.estimatePWDM(f,p,pwd_title);
      s2D.plot2D(f,pwd_title);
      //s2D.estimateSDW(k,f,p,sdw_title);
      s2D.plot2D(f,sdw_title);
    }
  }

  /**
   * Runs the method of choice for the constant 3D synthetic where:
   * method == 1, for structure tensor
   * method == 2, for plane-wave destructor
   * method == 3, for smooth dynamic warping
   * method == 4, for all methods.
   */
  private static void constant3D(int method) {
    Check.argument(method==1 || method==2 || method==3 || method==4,
        "not valid input");
    int n1 = 101;
    int n2 = 102;
    int n3 = 103;
    float[][][] f   = new float[n3][n2][n1];//2D seismic image
    float[][][] p2  = new float[n3][n2][n1]; //2D slope values
    float[][][] p3  = new float[n3][n2][n1]; //2D slope values
    Sampling s1 = new Sampling(n1);
    Sampling s2 = new Sampling(n2);
    Sampling s3 = new Sampling(n3);
    Slopes s3D = new Slopes(noise,pmax,s1,s2,s3);
    Slopes.makeSyntheticConstant(freq,pc2,pc3,f,p2,p3);
    String lsf_title = "lsf_constant3D";
    String pwd_title = "pwd_constant3D";
    String sdw_title = "sdw_constant3D";
    if (method==1) {
      //s3D.estimateLSF(f,p2,p3,lsf_title);
      s3D.plot3D(f,lsf_title);
    }
    if (method==2) {
      //s3D.estimatePWDM(f,p2,p3,pwd_title);
      s3D.plot3D(f,pwd_title);
    }
    if (method==3) { 
      //s3D.estimateSDW(k,f,p2,p3,sdw_title);
      s3D.plot3D(f,sdw_title);
    }
    if (method==4) {
      //s3D.estimateLSF(f,p2,p3,lsf_title);
      s3D.plot3D(f,lsf_title);
      //s3D.estimatePWDM(f,p2,p3,pwd_title);
      s3D.plot3D(f,pwd_title);
      //s3D.estimateSDW(k,f,p2,p3,sdw_title);
      s3D.plot3D(f,sdw_title);
    }
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
    float[][] f  = new float[n2][n1];//for GOM image
    Sampling s1 = new Sampling(n1,d1,f1);
    Sampling s2 = new Sampling(n2,d2,f2);
    s1 = new Sampling(n1);
    s2 = new Sampling(n2);
    Slopes s2D = new Slopes(noise,pmax,s1,s2);
    Slopes.makeRealGOM(f);
    String lsf_title = "lsf_GOM2D";
    String pwd_title = "pwd_GOM2D";
    String sdw_title = "sdw_GOM2D";
    if (method==1) {
      s2D.estimateLSF(f,null,lsf_title);
      s2D.plot2D(f,lsf_title);
    }
    if (method==2) {
      s2D.estimatePWDM(f,null,pwd_title);
      s2D.plot2D(f,pwd_title);
    }
    if (method==3) {
      s2D.estimateSDW(k,f,null,sdw_title);
      s2D.plot2D(f,sdw_title);
    }
    if (method==4) {
      //s2D.estimateLSF(f,null,lsf_title);
      s2D.plot2D(f,lsf_title);
      //s2D.estimatePWDM(f,null,pwd_title);
      s2D.plot2D(f,pwd_title);
      //s2D.estimateSDW(k,f,null,sdw_title);
      s2D.plot2D(f,sdw_title);
    }
  }

  /**
   * Runs the method of choice for the GOM real data where:
   * method == 1, for structure tensor
   * method == 2, for plane-wave destructor
   * method == 3, for smooth dynamic warping
   * method == 4, for all methods.
   */
  private static void TP2(int method) {
    Check.argument(method==1 || method==2 || method==3 || method==4,
        "not valid input");
    int n1 = 251;
    int n2 = 357;
    //int n1 = 401;
    //int n2 = 161;
    float d1 = 0.004f;
    float d2 = .025f;
    float f1 = 0.5f;
    float f2 = 0;
    float[][] f  = new float[n2][n1];//for TP image
    Sampling s1 = new Sampling(n1);
    Sampling s2 = new Sampling(n2);
    Slopes s2D = new Slopes(noise,pmax,s1,s2);
    Slopes.makeRealTp(f);
    String lsf_title = "lsf_TP2D";
    String pwd_title = "pwd_TP2D";
    String sdw_title = "sdw_TP2D";
    if (method==1) {
      s2D.estimateLSF(f,null,lsf_title);
      s2D.plot2D(f,lsf_title);
    }
    if (method==2) {
      s2D.estimatePWDM(f,null,pwd_title);
      s2D.plot2D(f,pwd_title);
    }
    if (method==3) {
      s2D.estimateSDW(k,f,null,sdw_title);
      s2D.plot2D(f,sdw_title);
    }
    if (method==4) {
      s2D.estimateLSF(f,null,lsf_title);
      s2D.plot2D(f,lsf_title);
      s2D.estimatePWDM(f,null,pwd_title);
      s2D.plot2D(f,pwd_title);
      s2D.estimateSDW(k,f,null,sdw_title);
      s2D.plot2D(f,sdw_title);
    }
    //s2D.plotTeaser(f,"teaser");
  }

  /**
   * Runs the method of choice for the GOM real data where:
   * method == 1, for structure tensor
   * method == 2, for plane-wave destructor
   * method == 3, for smooth dynamic warping
   * method == 4, for all methods.
   */
  private static void TP3(int method) {
    Check.argument(method==1 || method==2 || method==3 || method==4,
        "not valid input");
    int n1 = 251;
    int n2 = 357;
    int n3 = 161;
    float d1 = 0.002f;
    float d2 = .025f;
    float d3 = .025f;
    float f1 = 0;
    float f2 = 0;
    float f3 = 0;
    float[][][] f  = new float[n3][n2][n1]; //for TP image
    f = Util.readImage(n1,n2,n3,PATH+"data/tp/tpst_subt_251_4_500.dat");
    float[][][] g  = new float[n2][n3][n1]; //for TP image
    Sampling s1 = new Sampling(n1);
    Sampling s2 = new Sampling(n2);
    Sampling s3 = new Sampling(n3);
    Slopes s3D = new Slopes(noise,pmax,s1,s2,s3);
    Util.transpose23(f,g);
    //String lsf_title = "tp/lsf";
    //String pwd_title = "tp/pwd";
    //String sdw_title = "tp/sdw";
    String lsf_title = "lsf_TP3D";
    String pwd_title = "pwd_TP3D";
    String sdw_title = "sdw_TP3D";
    //String lsf_title = "lsf_temp";
    //String lsf_title = "lsf_gom3D";
    if (method==1) { 
      //s3D.estimateLSF(f,null,null,lsf_title);
      //s3D.estimateTransLSF(g,lsf_title+"_trans");
      s3D.plot3D(f,lsf_title);
      //s3D.retranspose(lsf_title+"_trans");
      //s3D.plot3D(f,lsf_title+"_trans");
      //s3D.plot3DSub(f,lsf_title,lsf_title+"_trans");
      //float[][][] p3 = Util.readImage(n1,n2,n3,PATH+"data/"+lsf_title+"_p3.dat");
      //s2D.plot2D(Util.slice13(317,p3),"lsf_TP2D");
    }
    if (method==2) { 
      //s3D.estimatePWDM(f,null,null,pwd_title);
      //s3D.estimateTransPWDM(g,pwd_title+"_trans");
      s3D.plot3D(f,pwd_title);
      //s3D.retranspose(pwd_title+"_trans");
      //s3D.plot3D(f,pwd_title+"_trans");
      //s3D.plot3DSub(f,pwd_title,pwd_title+"_trans");
      //float[][][] p3 = Util.readImage(n1,n2,n3,PATH+"data/"+pwd_title+"_p3.dat");
      //s2D.plot2D(Util.slice13(317,p3),"pwd_TP2D");
    }
    if (method==3) {
      //s3D.estimateSDW(k,f,null,null,sdw_title);
      //s3D.estimateTransSDW(k,g,sdw_title+"_trans");
      s3D.plot3D(f,sdw_title);
      //s3D.retranspose(sdw_title+"_trans");
      //s3D.plot3D(f,sdw_title+"_trans");
      //s3D.plot3DSub(f,sdw_title,sdw_title+"_trans");
      //float[][][] p3 = Util.readImage(n1,n2,n3,PATH+"data/"+sdw_title+"_p3.dat");
      //s2D.plot2D(Util.slice13(317,p3),"sdw_TP2D");
    }
    if (method==4) {
      s3D.estimateLSF(f,null,null,lsf_title);
      //s3D.estimateTransLSF(g,lsf_title+"_trans");
      s3D.plot3D(f,lsf_title);
      s3D.estimatePWDM(f,null,null,pwd_title);
      //s3D.estimateTransPWDM(g,pwd_title+"_trans");
      s3D.plot3D(f,pwd_title);
      s3D.estimateSDW(k,f,null,null,sdw_title);
      //s3D.estimateTransSDW(k,g,sdw_title+"_trans");
      s3D.plot3D(f,sdw_title);
      s3D.plotError(lsf_title,pwd_title,sdw_title);
    }
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
    float[][] f = new float[n2][n1]; //2D seismic image
    float[][] p = new float[n2][n1]; //2D slope values
    float[][] r = new float[n2][n1]; //2D reflectivity 
    Sampling s1 = new Sampling(n1);
    Sampling s2 = new Sampling(n2);
    Slopes s2D = new Slopes(noise,pmax,s1,s2);
    Slopes.makeSyntheticComplex(noise,f,p,r);

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
        s2D.testOptimalSmoothLSF(ttl,f,p,sp1,sp2);
        s2D.plotOptimalParameters(ttl,"sigma2","sigma1",sp1,sp2,cmin_lsf,cmax_lsf);
      }
      if (method==2) {
        ttl = "madagascar_pwd_errors_"+num+".dat";
        s2D.testOptimalSmoothPWDM(ttl,f,p,sp1,sp2);
        s2D.plotOptimalParameters(ttl,"rect2","rect1",sp1,sp2,cmin_pwd,cmax_pwd);
      }
      if (method==3) {
        ttl = "smooth_dynamic_warping_h_errors_"+num+".dat";
        s2D.testOptimalSmoothSDW(ttl,k,f,p,sp1,sp2);
        s2D.plotOptimalParameters(ttl,"h2","h1",sp1,sp2,cmin_sdw,cmax_sdw);
        dp = 0.1f; // strain sampling rate
        fp = 0.1f;
        np1 = (int)((2.0)/dp); //# strain vals to test
        sp1 = new Sampling(np1,dp,fp);
        sp2 = new Sampling(np1,dp,fp);
        ttl = "smooth_dynamic_warping_r_errors_"+num+".dat";
        s2D.testOptimalStrainSDW(ttl,k,f,p,sp1,sp2);
        s2D.plotOptimalParameters(ttl,"r2","r1",sp1,sp2,cmin_sdwr,cmax_sdwr);
      }
      if (method==4) {
        ttl = "structure_tensor_errors_"+num+".dat";
        s2D.testOptimalSmoothLSF(ttl,f,p,sp1,sp2);
        s2D.plotOptimalParameters(ttl,"sigma2","sigma1",sp1,sp2,cmin_lsf,cmax_lsf);
        ttl = "madagascar_pwd_errors_"+num+".dat";
        s2D.testOptimalSmoothPWDM(ttl,f,p,sp1,sp2);
        s2D.plotOptimalParameters(ttl,"rect2","rect1",sp1,sp2,cmin_pwd,cmax_pwd);
        ttl = "smooth_dynamic_warping_h_errors_"+num+".dat";
        s2D.testOptimalSmoothSDW(ttl,k,f,p,sp1,sp2);
        s2D.plotOptimalParameters(ttl,"h2","h1",sp1,sp2,cmin_sdw,cmax_sdw);
        dp = 0.1f; // strain sampling rate
        fp = 0.1f;
        np1 = (int)((2.0)/dp); //# strain vals to test
        sp1 = new Sampling(np1,dp,fp);
        sp2 = new Sampling(np1,dp,fp);
        ttl = "smooth_dynamic_warping_r_errors_"+num+".dat";
        s2D.testOptimalStrainSDW(ttl,k,f,p,sp1,sp2);
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

  private static void rmsErrorCurves2(int method, boolean test) {
    Check.argument(method==1 || method==2 || method==3 || method==4,
        "not valid input");
    if (noise==0.0f)      num = 1;
    else if (noise==0.5f) num = 2;
    else if (noise==1.0f) num = 3;
    else                  num = 4;    
    int n1 = 501;
    int n2 = 501;
    float[][] f = new float[n2][n1];//2D seismic image
    float[][] p = new float[n2][n1]; //2D slope values
    float[][] r = new float[n2][n1]; //2D reflectivity 
    Sampling s1 = new Sampling(n1);
    Sampling s2 = new Sampling(n2);
    Slopes s2D = new Slopes(noise,pmax,s1,s2);
    Slopes.makeSyntheticComplex(noise,f,p,r);
    double dn = 0.05;
    int nn = 21;
    Sampling sn = new Sampling(nn,dn,0.0);
    float sigma1,sigma2;
    int rect1,rect2;
    double r1,r2,h1,h2;
    String fileName1 = "data/rms_error2_lsf_"+num+".dat";
    String fileName2 = "data/rms_error2_pwd_"+num+".dat";
    String fileName3 = "data/rms_error2_sdw_"+num+".dat";
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
    String vl = "slope error (samples/trace)";
    s2D.plotCurve(sn,fileName1,fileName2,fileName3,title,hl,vl,0,1);
  }

  private static void rmsErrorCurves3(int method, boolean test) {
    Check.argument(method==1 || method==2 || method==3 || method==4,
        "not valid input");
    if (noise==0.0f)      num = 1;
    else if (noise==0.5f) num = 2;
    else if (noise==1.0f) num = 3;
    else                  num = 4;    
    int n1 = 101;
    int n2 = 102;
    int n3 = 103;
    float[][][] f  = new float[n3][n2][n1];//3D seismic image
    float[][][] p2 = new float[n3][n2][n1];//3D slope values
    float[][][] p3 = new float[n3][n2][n1];//3D slope values
    Sampling s1 = new Sampling(n1);
    Sampling s2 = new Sampling(n2);
    Sampling s3 = new Sampling(n3);
    Slopes s3D = new Slopes(noise,pmax,s1,s2,s3);
    Slopes.makeSyntheticComplex(noise,f,p2,p3);
    double dn = 0.05;
    int nn = 21;
    Sampling sn = new Sampling(nn,dn,0.0);
    float sigma1,sigma2,sigma3;
    int rect1,rect2,rect3;
    double r1,r2,r3,h1,h2,h3;
    String fn_lsf2 = "data/rms_error_p2_lsf_"+num+".dat";
    String fn_lsf3 = "data/rms_error_p3_lsf_"+num+".dat";
    String fn_pwd2 = "data/rms_error_p2_pwd_"+num+".dat";
    String fn_pwd3 = "data/rms_error_p3_pwd_"+num+".dat";
    String fn_sdw2 = "data/rms_error_p2_sdw_"+num+".dat";
    String fn_sdw3 = "data/rms_error_p3_sdw_"+num+".dat";
    if (noise==0.0) {
      sigma1=14.0f; sigma2=1.0f; sigma3=sigma2;
      rect1=34; rect2=2; rect3=rect2;
      h1=30.0; h2=9.0; h3=h2;
      r1=0.1;  r2=0.6; r3=r2;
    }

    else {
      sigma1=23.0f; sigma2=1.0f; sigma3=sigma2;
      rect1=20; rect2=9; rect3=rect2;
      h1=20.0; h2=9.0; h3=h2;
      r1=0.1; r2=0.6; r3=r2;
    }

    if (test) {
      if (method==1) 
        s3D.testRmsErrorCurveLSF(sn,fn_lsf2,fn_lsf3,sigma1,sigma2,sigma3);
      if (method==2) 
        s3D.testRmsErrorCurvePWD(sn,fn_pwd2,fn_pwd3,rect1,rect2,rect3);
      if (method==3) 
        s3D.testRmsErrorCurveSDW(sn,fn_sdw2,fn_sdw3,k,r1,r2,r3,h1,h2,h3);
      if (method==4) {
        //s3D.testRmsErrorCurveLSF(sn,fn_lsf2,fn_lsf3,sigma1,sigma2,sigma3);
        s3D.testRmsErrorCurvePWD(sn,fn_pwd2,fn_pwd3,rect1,rect2,rect3);
        s3D.testRmsErrorCurveSDW(sn,fn_sdw2,fn_sdw3,k,r1,r2,r3,h1,h2,h3);
      }
    }
    String title = "pwd_lsf_sdw_rmserror_p2_vs_nsratio"+num;
    String hl = "Noise/signal";
    String vl = "slope error (samples/trace)";
    s3D.plotCurve(sn,fn_lsf2,fn_pwd2,fn_sdw2,title,
        hl,vl,0,0.3f);
    title = "pwd_lsf_sdw_rmserror_p3_vs_nsratio"+num;
    s3D.plotCurve(sn,fn_lsf3,fn_pwd3,fn_sdw3,title,
        hl,vl,0,0.3f);
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
    float[][] f = new float[n2][n1];//2D seismic image
    float[][] p = new float[n2][n1]; //2D slope values
    float[][] r = new float[n2][n1]; //2D reflectivity 
    Sampling s1 = new Sampling(n1);
    Sampling s2 = new Sampling(n2);
    Slopes s2D = new Slopes(noise,pmax,s1,s2);
    Slopes.makeSyntheticComplex(noise,f,p,r);
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
    String vl = "slope error (samples/trace)";
    s2D.plotCurve(sn,fileName1,fileName2,fileName3,title,hl,vl,0,0.02f);
  }

  private static void mean(int method, boolean test, int nni) {
    Check.argument(method==1 || method==2 || method==3 || method==4,
        "not valid input");
    int n1 = 501;
    int n2 = 501;
    float[][] f = new float[n2][n1];//2D seismic image
    float[][] p = new float[n2][n1]; //2D slope values
    float[][] r = new float[n2][n1]; //2D reflectivity 
    Sampling s1 = new Sampling(n1);
    Sampling s2 = new Sampling(n2);
    Slopes s2D = new Slopes(noise,pmax,s1,s2);
    Slopes.makeSyntheticComplex(noise,f,p,r);

    float sigma1 = 23.0f; float sigma2 = 1.0f;
    int rect1 = 75;       int rect2 = 6;
    double h1 = 72.0;     double h2 = 12.0;
    double r1 = 0.1;      double r2 = 0.4;
    //fn = filename
    String fn_lsf = "data/mean_lsf.dat";
    String fn_pwd = "data/mean_pwd.dat";
    String fn_sdw = "data/mean_sdw.dat";
    String title = "lsf_mean"+nni;
    String hl = "Traces";
    String vl = "Samples";
    String cbl = "slope (samples/trace)";
    float cmax = 4.0f;
    if (method==1) {
      if (test) s2D.testSampleMeanLSF(fn_lsf,nni,noise,sigma1,sigma2);
      s2D.plotImage(title,hl,vl,cbl,fn_lsf,-cmax,cmax);
    }
    if (method==2) {
      if (test) s2D.testSampleMeanPWD(fn_pwd,nni,noise,rect1,rect2);
      title = "pwd_mean"+nni;
      s2D.plotImage(title,hl,vl,cbl,fn_pwd,-cmax,cmax);
    }
    if (method==3) {
      if (test) s2D.testSampleMeanSDW(fn_sdw,nni,noise,k,r1,r2,h1,h2);
      title = "sdw_mean"+nni;
      s2D.plotImage(title,hl,vl,cbl,fn_sdw,-cmax,cmax);
    }
    if (method==4) {
      if (test) {
        s2D.testSampleMeanLSF(fn_lsf,nni,noise,sigma1,sigma2);
        s2D.testSampleMeanPWD(fn_pwd,nni,noise,rect1,rect2);
        s2D.testSampleMeanSDW(fn_sdw,nni,noise,k,r1,r2,h1,h2);
      }
      s2D.plotImage(title,hl,vl,cbl,fn_lsf,-cmax,cmax);
      title = "pwd_mean"+nni;
      s2D.plotImage(title,hl,vl,cbl,fn_pwd,-cmax,cmax);
      title = "sdw_mean"+nni;
      s2D.plotImage(title,hl,vl,cbl,fn_sdw,-cmax,cmax);
    }
  }

  private static void stdDev(int method, boolean test, int nni) {
    Check.argument(method==1 || method==2 || method==3 || method==4,
        "not valid input");
    int n1 = 501;
    int n2 = 501;
    float[][] f = new float[n2][n1];//2D seismic image
    float[][] p = new float[n2][n1]; //2D slope values
    float[][] r = new float[n2][n1]; //2D reflectivity 
    Sampling s1 = new Sampling(n1);
    Sampling s2 = new Sampling(n2);
    Slopes s2D = new Slopes(noise,pmax,s1,s2);
    Slopes.makeSyntheticComplex(noise,f,p,r);

    float sigma1 = 23.0f; float sigma2 = 1.0f;
    int rect1 = 75;       int rect2 = 6;
    double h1 = 72.0;     double h2 = 12.0;
    double r1 = 0.1;      double r2 = 0.4;
    //fn = filename
    String fn_lsf = "data/stddev_lsf.dat";
    String fn_pwd = "data/stddev_pwd.dat";
    String fn_sdw = "data/stddev_sdw.dat";
    String title = "lsf_stddev"+nni;
    String hl = "Traces";
    String vl = "Samples";
    String cbl = "std (samples/trace)";
    float cmin = 0.0f;
    float cmax = 1.0f;
    if (method==1) {
      if (test) s2D.testSampleStdDevLSF(fn_lsf,nni,sigma1,sigma2);
      s2D.plotImage(title,hl,vl,cbl,fn_lsf,cmin,cmax);
    }
    if (method==2) {
      if (test) s2D.testSampleStdDevPWD(fn_pwd,nni,rect1,rect2);
      title = "pwd_stddev"+nni;
      s2D.plotImage(title,hl,vl,cbl,fn_pwd,cmin,cmax);
    }
    if (method==3) {
      if (test) s2D.testSampleStdDevSDW(fn_sdw,nni,k,r1,r2,h1,h2);
      title = "sdw_stddev"+nni;
      s2D.plotImage(title,hl,vl,cbl,fn_sdw,cmin,cmax);
    }
    if (method==4) {
      if (test) {
        s2D.testSampleStdDevLSF(fn_lsf,nni,sigma1,sigma2);
        s2D.testSampleStdDevPWD(fn_pwd,nni,rect1,rect2);
        s2D.testSampleStdDevSDW(fn_sdw,nni,k,r1,r2,h1,h2);
      }
      s2D.plotImage(title,hl,vl,cbl,fn_lsf,cmin,cmax);
      title = "pwd_stddev"+nni;
      s2D.plotImage(title,hl,vl,cbl,fn_pwd,cmin,cmax);
      title = "sdw_stddev"+nni;
      s2D.plotImage(title,hl,vl,cbl,fn_sdw,cmin,cmax);
    }
  }

  private static void trace(String s) {
    System.out.println(s);
  }

  private static int num;
  private static final int norder = 15;
  private static final int k = 10;
  private static final float pmax = 9.0f;
  private static final float noise = 0.0f;
  private static final float freq = 0.1f;
  private static final float pc2   = -0.7f;//constant slope
  private static final float pc3   = 1.3f;//constant slope
  private static final String PATH = 
    "/users/elias.arias/Home/git/ea/bench/src/util/";
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

        //1=lsf  2=pwd  3=sdw  4=all
        //complex2D(1);
        //1=lsf  2=pwd  3=sdw  4=all
        complex3D(1);
        //1=lsf  2=pwd  3=sdw  4=all
        //constant2D(4);
        //1=lsf  2=pwd  3=sdw  4=all
        //constant3D(4);
        //1=lsf  2=pwd  3=sdw  4=all
        //GOM(4);
        //1=lsf  2=pwd  3=sdw  4=all
        //TP2(4);
        //1=lsf  2=pwd  3=sdw  4=all
        //TP3(4);

        int nni = 100;
        //1=lsf  2=pwd  3=sdw  4=all, test?
        //rmsErrorCurves2(4,true);
        //1=lsf  2=pwd  3=sdw  4=all, test?
        //rmsErrorCurves3(4,false);
        //1=lsf  2=pwd  3=sdw  4=all, test?
        //meanErrorCurves(4,false,nni);
        //1=lsf  2=pwd  3=sdw  4=all, test?, n1, n2
        //optimal(1,true,15, 3);
        //1=lsf  2=pwd  3=sdw  4=all, test?
        //mean(4,false,nni);
        //1=lsf  2=pwd  3=sdw  4=all, test?
        //stdDev(4,false,nni);

        Sampling sovt = new Sampling(norder);
        String fn = "data/order_vs_time_"+norder+".dat";
        String title = "order_vs_time";
        String hl = "Order";
        String vl = "Time";
        //s2D.testOrderVsTime(sovt,fn,f2D,p2D);
        //s2D.plotCurve(sovt,fn,title,hl,vl,0,78);

        //goTestSlopeVsError();
      }
    });
  }
}
