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
        //Slopes.makeSyntheticConstant(freq,pc2,f2D,p2D);
        Slopes.makeSyntheticComplex(noise,f2D,p2D);

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

        double dp = 1.0f;
        double fp = 1.0f;
        double lp = 80.0f;
        int np1 = (int)((lp)/dp);
        int np2 = 20;
        Sampling sp1 = new Sampling(np1,dp,fp);
        Sampling sp2 = new Sampling(np2,dp,fp);

        double dn = 0.05;
        int nn = 21;
        Sampling sn = new Sampling(nn,dn,0.0);
        //s2D.testRmsErrorCurveLSF(sn);
        //s2D.testRmsErrorCurvePWD(sn);
        //s2D.testRmsErrorCurveSDW(sn,k);
        s2D.plotRmsErrorCurves(sn);

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
        //s2D.plotOptimalParameters("Smooth dynamic warping r","r2","r1",sp1,sp2);

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
