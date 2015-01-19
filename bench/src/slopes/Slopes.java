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

public class Slopes{

  public static void goLSF(Sampling s1, Sampling s2, 
      float[][] f, float[][] pk, boolean error) {
    float d1 = (float)s1.getDelta();
    float d2 = (float)s2.getDelta();
    int n1 = s1.getCount();
    int n2 = s2.getCount();

    float pmax = 5.0f;
    //Structure tensor method
    LocalSlopeFinder lsf = new LocalSlopeFinder(24.0f,1.0f,pmax);
    float[][] p_lsf = new float[n2][n1];
    lsf.findSlopes(f,p_lsf);

    if (error) {
      System.out.println("Structure tensor:");
      float error_lsf = Util.rmsError(p_lsf,pk,d1,d2,T);
    }

    float fw = 0.75f; //fraction width for slide
    float fh = 0.9f; //fraction height for slide
    // title, paint, colorbar, color
    Plot.plot(s1,s2,p_lsf,"LSF noise= "+noise,fw,fh,-4,4,title,paint,T,T);
  }

  public static void goPWDM(Sampling s1, Sampling s2, 
      float[][] f, float[][] pk, boolean error) {
    float d1 = (float)s1.getDelta();
    float d2 = (float)s2.getDelta();
    int n1 = s1.getCount();
    int n2 = s2.getCount();

    float pmax = 5.0f;
    int niter = 5;

    //Plane-wave destruction filter
    Sfdip sd = new Sfdip(-pmax,pmax);
    sd.setRect(76,6);
    sd.setOrder(5);
    sd.setNiter(niter);
    //sd.setBoth("y");
    float[][] p_pwdm = new float[n2][n1]; //pwd w/initial p=0 (Madagascar)
    sd.findSlopes(s1,s2,f,p_pwdm);

    if (error) {
      System.out.println("Madagascar PWD:");
      float error_pwdm = Util.rmsError(p_pwdm,pk,d1,d2,T);
    }

    float fw = 0.75f; //fraction width for slide
    float fh = 0.9f; //fraction height for slide
    // title, paint, colorbar, color
    Plot.plot(s1,s2,p_pwdm,"PWDM noise= "+noise,fw,fh,-4,4,title,paint,T,T);
  }

  public static void goDW(Sampling s1, Sampling s2, 
      float[][] f, float[][] pk, boolean error) {
    float d1 = (float)s1.getDelta();
    float d2 = (float)s2.getDelta();
    int n1 = s1.getCount();
    int n2 = s2.getCount();
    float pmax = 5.0f;

    //Dynamic warping
    DynamicWarping dw = new DynamicWarping((int)-pmax,(int)pmax);
    dw.setShiftSmoothing(40.0,1.0);
    dw.setStrainMax(1.0,1.0);
    dw.setErrorSmoothing(2);
    float[][] p_dw= Util.DWSlopesAvg(dw,f);

    if (error) {
      System.out.println("Dynamic warping:");
      float error_dw = Util.rmsError(p_dw,pk,d1,d2,T);
    }

    float fw = 0.75f; //fraction width for slide
    float fh = 0.9f; //fraction height for slide
    // title, paint, colorbar, color
    Plot.plot(s1,s2,p_dw,"DW noise= "+noise,fw,fh,-4,4,title,paint,T,T);
  }

  public static void goPWDD(Sampling s1, Sampling s2, 
      float[][] f, float[][] pk, boolean error) {
    float d1 = (float)s1.getDelta();
    float d2 = (float)s2.getDelta();
    int n1 = s1.getCount();
    int n2 = s2.getCount();
    float pmax = 5.0f;
    int niter = 5;

    //Dave's plane-wave destruction filter
    PlaneWaveDestructor pwd = new PlaneWaveDestructor(-pmax,pmax);
    pwd.setSmoothness(6,1);
    pwd.setLateralBias(0.0);
    pwd.setOuterIterations(niter); // default is 5
    pwd.setInnerIterations(20); // default is 20
    float[][] p_pwdd = new float[n2][n1]; //pwd w/initial p=0 (Dave)
    p_pwdd= pwd.findSlopes(f);
    pwd.updateSlopes(f,p_pwdd);

    if (error) {
      System.out.println("Dave's PWD:");
      float error_pwdd = Util.rmsError(p_pwdd,pk,d1,d2,T);
    }

    float fw = 0.75f; //fraction width for slide
    float fh = 0.9f; //fraction height for slide
    // title, paint, colorbar, color
    Plot.plot(s1,s2,p_pwdd,"PWDD noise= "+noise,fw,fh,-4,4,title,paint,T,T);
  }

  public static void goFandPk(Sampling s1, Sampling s2) {
    float d1 = (float)s1.getDelta();
    float d2 = (float)s2.getDelta();
    int n1 = s1.getCount();
    int n2 = s2.getCount();

    float[][][] fandpk = FakeData.seismicAndSlopes2d2014B(noise,F);
    float[][] f = fandpk[0];  //synthetic seismic data
    float[][] pk = fandpk[1]; //exact slope values    
    System.out.println("max slope= "+max(pk));

    float fw = 0.75f; //fraction width for slide
    float fh = 0.9f; //fraction height for slide
    // title, paint, colorbar, color
    Plot.plot(s1,s2,f,"Synthetic seismic noise= "+noise,fw,fh,F,T,T,F);
    Plot.plot(s1,s2,pk,"Known Slopes",fw,fh,-4.0f,4.0f,title,paint,T,T);
  }

  private static final float noise = 0.5f;  
  private static final boolean T = true;
  private static final boolean F = false;  
  private static final boolean title = true;
  private static final boolean paint = true;  
}
