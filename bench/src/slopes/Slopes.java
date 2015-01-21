package slopes;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.Check;
import edu.mines.jtk.util.Stopwatch;

import static edu.mines.jtk.util.ArrayMath.*;

import dnp.LocalSlopeFinder;
import dnp.PlaneWaveDestructor;
import util.FakeData;
import util.*;

import javax.swing.*;

/**
 *  
 * @author Elias Arias, Colorado School of Mines, CWP
 * @version 19.1.2015
 */

public class Slopes{

  public static void setChickenTestParameters(float pk) {
    _n1 = 501;
    _n2 = 501;
    _d1 = 1.0f;
    _d2 = 1.0f;
    _f1 = 0.0f;
    _f2 = 0.0f;
    _s1 = new Sampling(_n1,_d1,_f1);
    _s2 = new Sampling(_n2,_d2,_f2);
    float w = 0.1f; //frequency
    float k = pk*w; //wavenumber
    _f = new float[_n2][_n1];
    _pk = new float[_n2][_n1];
    for (int i2=0; i2<_n2; ++i2){
      for (int i1=0; i1<_n1; ++i1){
        _f[i2][i1] = cos(k*i2-w*i1);
        _pk[i2][i1] = pk;
      }
    }
  }

  public static void setSynthParameters(float noise) {
    _n1 = 501;
    _n2 = 501;
    _d1 = 1.0f;
    _d2 = 1.0f;
    _f1 = 0.0f;
    _f2 = 0.0f;
    _s1 = new Sampling(_n1,_d1,_f1);
    _s2 = new Sampling(_n2,_d2,_f2);
    _noise = noise;
    float[][][] fandpk = FakeData.seismicAndSlopes2d2014B(_noise,F);
    _f = fandpk[0];  //synthetic seismic data
    _pk = fandpk[1]; //exact slope values
  }

  public static void setGOMParameters() {
    _n1 = 301;
    _n2 = 920;
    _d1 = 0.004f;
    _d2 = .02667f;
    _f1 = 1.6f;
    _f2 = 0;
    _s1 = new Sampling(_n1,_d1,_f1);
    _s2 = new Sampling(_n2,_d2,_f2);
    _f = Util.readImage(_n1,_n2,"data/gom.dat");
    float[][] fs = new float[_n2][_n1];
    mul(_f,.001f,fs);
    Util.writeBinary(fs,"data/gom_scaled.dat");
  }

  public static void setSmoothingParameters(float l_param) {
    d_param= 1.0f;
    f_param = 1.0f;
    n1_param = (int)((l_param-f_param+1)/d_param);
    n2_param = 10;
    s1_param = new Sampling(n1_param,d_param,f_param);
    s2_param = new Sampling(n2_param,d_param,f_param);
    param1 = new float[n1_param];
    param2 = new float[n2_param];
    for(int i=0; i<n1_param; ++i) {
      param1[i] = i*d_param+f_param;
    }
    for(int i=0; i<n2_param; ++i) {
      param2[i] = i*d_param+f_param;
    }
  }

  public static void setStrainParameters() {
    d_param = 0.1f; // strain sampling rate
    f_param = 0.1f;
    n1_param = (int)((1.0f-f_param+0.1f)/d_param); //# strain vals to test
    n2_param = n1_param;
    s1_param = new Sampling(n1_param,d_param,f_param);
    s2_param = new Sampling(n2_param,d_param,f_param);
    param1 = new float[n1_param];
    param2 = new float[n2_param];
    for(int i=0; i<n1_param; ++i) {
      param1[i] = i*d_param+f_param;
      param2[i] = i*d_param+f_param;
    }
  }

  /**
   * Structure tensor: plots the estimated slopes and RMS error.
   */
  public static void goLSF(boolean error) {
    LocalSlopeFinder lsf = new LocalSlopeFinder(10.0f,1.0f,pmax);
    float[][] p_lsf = new float[_n2][_n1];
    lsf.findSlopes(_f,p_lsf);

    if (error) {
      System.out.println("Structure tensor:");
      float error_lsf = Util.rmsError(p_lsf,_pk,_d1,_d2,T);
    }

    // title, paint, colorbar, color
    Plot.plot(_s1,_s2,p_lsf,"LSF noise= "+_noise,fw,fh,-clipMax,clipMax,
        title,paint,T,T);
  }

  /**
   * PWD Madagascar: plots the estimated slopes and RMS error.
   */
  public static void goPWDM(boolean error) {
    Sfdip sd = new Sfdip(-pmax,pmax);
    sd.setRect(76,6);
    sd.setOrder(5);
    sd.setNiter(niter);
    //sd.setBoth("y");
    float[][] p_pwdm = new float[_n2][_n1]; //pwd w/initial p=0 (Madagascar)
    sd.findSlopes(_s1,_s2,_f,p_pwdm);

    if (error) {
      System.out.println("Madagascar PWD:");
      float error_pwdm = Util.rmsError(p_pwdm,_pk,_d1,_d2,T);
    }

    // title, paint, colorbar, color
    Plot.plot(_s1,_s2,p_pwdm,"PWDM noise= "+_noise,fw,fh,-clipMax,clipMax,
        title,paint,T,T);
  }

  /**
   * Dynamic warping: plots the estimated slopes and RMS error.
   */
  public static void goDW(boolean error) {
    DynamicWarping dw = new DynamicWarping((int)-pmax,(int)pmax);
    dw.setShiftSmoothing(10.0,1.0);
    dw.setStrainMax(0.1,0.1);
    dw.setErrorSmoothing(2);
    float[][] p_dw= Util.DWSlopesAvg(dw,_f);

    if (error) {
      System.out.println("Dynamic warping:");
      float error_dw = Util.rmsError(p_dw,_pk,_d1,_d2,T);
    }

    // title, paint, colorbar, color
    Plot.plot(_s1,_s2,p_dw,"DW noise= "+_noise,fw,fh,-clipMax,clipMax,
        title,paint,T,T);
  }

  /**
   * PWD Dave: plots the estimated slopes and RMS error.
   */
  public static void goPWDD(boolean error) {
    PlaneWaveDestructor pwd = new PlaneWaveDestructor(-pmax,pmax);
    pwd.setSmoothness(6,1);
    pwd.setLateralBias(0.0);
    pwd.setOuterIterations(niter); // default is 5
    pwd.setInnerIterations(20); // default is 20
    float[][] p_pwdd = new float[_n2][_n1]; //pwd w/initial p=0 (Dave)
    p_pwdd= pwd.findSlopes(_f);
    pwd.updateSlopes(_f,p_pwdd);

    if (error) {
      System.out.println("Dave's PWD:");
      float error_pwdd = Util.rmsError(p_pwdd,_pk,_d1,_d2,T);
    }

    // title, paint, colorbar, color
    Plot.plot(_s1,_s2,p_pwdd,"PWDD noise= "+_noise,fw,fh,-clipMax,clipMax,
        title,paint,T,T);
  }

  /**
   * Plots the synthetic seismic image and the corresponding slopes image.
   */
  public static void plotFandPk() {
    System.out.println("max slope= "+max(_pk));
    // title, paint, colorbar, color
    Plot.plot(_s1,_s2,_f,"Synthetic noise= "+_noise,fw,fh,title,paint,T,F);
    Plot.plot(_s1,_s2,_pk,"Known Slopes",fw,fh,-clipMax,clipMax,title,paint,T,T);
  }

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

    n1_param = s1_param.getCount();
    n2_param = s2_param.getCount();
    float[][] pe = new float[_n2][_n1];
    float[][] rmserror = new float[n2_param][n1_param]; // RMS error
    if (m==1) {
      LocalSlopeFinder lsf;
      for(int i2=0; i2<n2_param; ++i2) {
        for(int i1=0; i1<n1_param; ++i1) {
          lsf = new LocalSlopeFinder(param1[i1],param2[i2],pmax);
          lsf.findSlopes(_f,pe);
          rmserror[i2][i1] = Util.rmsError(pe,_pk,_d1,_d2,false);
        }
      }
    }

    else if (m==2) {
      Sfdip sd = new Sfdip(-pmax,pmax);
      sd.setOrder(2);
      for(int i2=0; i2<n2_param; ++i2) {
        for(int i1=0; i1<n1_param; ++i1) {
          sd.setRect((int)param1[i1],(int)param2[i2]);
          sd.findSlopes(_s1,_s2,_f,pe);
          rmserror[i2][i1] = Util.rmsError(pe,_pk,_d1,_d2,false);
        }
      }
    }

    else if (m==3) {
      DynamicWarping dw = new DynamicWarping((int)-pmax,(int)pmax);
      dw.setStrainMax(1.0,1.0);
      dw.setErrorSmoothing(2);
      for(int i2=0; i2<n2_param; ++i2) {
        for(int i1=0; i1<n1_param; ++i1) {
          dw.setShiftSmoothing(param1[i1],param2[i2]);
          pe = Util.DWSlopesAvg(dw,_f);
          rmserror[i2][i1] = Util.rmsError(pe,_pk,_d1,_d2,false);
        }
      }
    }

    else if (m==4) {
      DynamicWarping dw = new DynamicWarping((int)-pmax,(int)pmax);
      dw.setErrorSmoothing(2);
      for(int i2=0; i2<n2_param; ++i2) {
        for(int i1=0; i1<n1_param; ++i1) {
          dw.setStrainMax(1.0,1.0);
          dw.setShiftSmoothing(param1[i1],param2[i2]);
          pe = Util.DWSlopesAvg(dw,_f);
          rmserror[i2][i1] = Util.rmsError(pe,_pk,_d1,_d2,false);
        }
      }
    }

    else if (m==5) {
      PlaneWaveDestructor pwd = new PlaneWaveDestructor(-pmax,pmax);
      for(int i2=0; i2<n2_param; ++i2) {
        for(int i1=0; i1<n1_param; ++i1) {
          pwd.setSmoothness(param1[i1],param2[i2]);
          pe = pwd.findSlopes(_f);
          rmserror[i2][i1] = Util.rmsError(pe,_pk,_d1,_d2,false);
        }
      }
    }
    Util.writeBinary(rmserror,path+method+"_errors.dat");
  }

  public static void plotOptimalParameters(String method, String hl, String vl) {
    n1_param = s1_param.getCount();
    n2_param = s2_param.getCount();

    int[] error_index = new int[2];    
    float[][] rmserror = Util.readImage(n1_param,n2_param,
        path+method+"_errors.dat");
    float min_error = min(rmserror,error_index);
    System.out.println(vl+"= "+param1[error_index[0]]+" "+
                       hl+"= "+param2[error_index[1]]+" "+
                      " Minimum Error Value= "+min_error);    

    // interp, title, paint, colorbar, color
    Plot.plot(s1_param,s2_param,rmserror,"RMS Error "+method,hl,vl,fw,fh,
        T,F,T,T,T);
  }

  ///////////////////VARIABLES///////////////////////
  private static final String path = 
    "/Users/earias/Home/git/ea/bench/src/slopes/data/";
  private static final int niter = 5;
  private static final float fw = 0.75f; //fraction width for slide
  private static final float fh = 0.9f; //fraction height for slide
  private static final float pmax = 5.0f;
  private static final float clipMax = 4.0f;
  private static final boolean T = true;
  private static final boolean F = false;  
  public static final boolean title = true;
  public static final boolean paint = false;  
  public static Sampling _s1 = new Sampling(1,1,1);  
  private static Sampling _s2 = new Sampling(1,1,1);
  private static Sampling s1_param = new Sampling(1,1,1);
  private static Sampling s2_param = new Sampling(1,1,1);

  private static int _n1,_n2,n1_param,n2_param;
  private static float _d1,_d2,_f1,_f2,d_param,f_param,_noise;
  private static float[] param1 = new float[n1_param];
  private static float[] param2 = new float[n2_param];
  private static float[][] _f = new float[_n2][_n1];
  private static float[][] _pk = new float[_n2][_n1];

}
