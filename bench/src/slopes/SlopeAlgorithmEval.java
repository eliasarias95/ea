package slopes;

import edu.mines.jtk.io.*;
import edu.mines.jtk.dsp.*;
import edu.mines.jtk.awt.ColorMap;
import edu.mines.jtk.mosaic.SimplePlot;

import static edu.mines.jtk.util.ArrayMath.*;

import dnp.LocalSlopeFinder;
import util.FakeData;
import utils.Plot;

import java.io.*;
import java.nio.*;
import java.util.Random;

import javax.swing.*;

/**
 * This software acts as a machine for testing slope estimation algorithms.
 * Using other pieces of software, I generate and save a synthetic image 
 * containing structural features that may be encountered in real seismic 
 * images. 
 * @author Elias Arias, Colorado School of Mines, CWP
 * @version 09.11.2014
 */

public class SlopeAlgorithmEval {

  private static float[] addNoise(double nrms, float[] f) {
    int n1 = f.length;
    Random r = new Random(1);
    float[] g = mul(2.0f,sub(randfloat(r,n1),0.5f));
    RecursiveGaussianFilter rgf = new RecursiveGaussianFilter(2.0);
    rgf.apply1(g,g); // 1st derivative enhances high-frequencies
    g = mul(g,(float)nrms*rms(f)/rms(g));
    return add(f,g);
  }

  private static float rms(float[] f) {
    int n1 = f.length;
    double sum = 0.0;
    for (int i1=0; i1<n1; ++i1) {
      float fi = f[i1];
      sum += fi*fi;
    }
    return (float)sqrt(sum/n1);
  }

  private static float[][] DWSlopesCluster(DynamicWarping dw, 
      float[][] synth_data, int ntraces) {

    int n1 = synth_data[0].length;
    int n2 = synth_data.length;

    int[] index = new int[ntraces];

    float[][] dw_slope = new float[n2][n1];
    float[][] errors = new float[n2][n1];

    float[][] intermed_shifts = new float[ntraces][n1];
    float[][] cluster = new float[ntraces][n1];
    float[][] cluster_shifts = new float[ntraces][n1];

    float[] error_sum = new float[ntraces];

    //shift n traces to the same position which is the nth traces position
    for (int it=ntraces-1; it>0; --it) {
      intermed_shifts[it] = dw.findShifts(synth_data[it],
                                         synth_data[it-1]);
      cluster[it] = dw.applyShifts(intermed_shifts[it],synth_data[it-1]);
    }

    //now have n traces at the same level and attempting to shift each trace
    //to the (n+1)th trace 
    for (int i2=ntraces; i2<n2; ++i2) {
      for (int it=0; it<ntraces; ++it) {
        errors = dw.computeErrors(cluster[it],synth_data[i2]);
        cluster_shifts[it] = dw.findShifts(cluster[it],synth_data[i2]);
        error_sum[it] = dw.sumErrors(errors,cluster_shifts[it]);
      }

      min(error_sum,index);
      dw_slope[i2] = cluster_shifts[index[0]];

      for (int it=ntraces-1; it>0; --it) {
        intermed_shifts[it] = dw.findShifts(synth_data[i2],
                                        synth_data[i2+it-ntraces]);
        cluster[it] = dw.applyShifts(intermed_shifts[it],
                                         synth_data[i2+it-ntraces]);
      }
    }

    for (int i=0; i<ntraces; ++i)
      dw_slope[i] = dw.findShifts(synth_data[i],synth_data[i+1]);

    return dw_slope;
  }

  /**
   * 
   */
  private static void goSynth() {
    int n1 = 501;
    int n2 = 501;

    //float d1 = 0.002f;
    //float d2 = 0.0016f;
    float d1 = 1.0f;
    float d2 = 1.0f;

    float f1 = 0;
    float f2 = 0;
    float pmax = 5.0f;
    float noise = 1.0f;
    int niter = 5;

    Sampling s1 = new Sampling(n1,d1,f1);
    Sampling s2 = new Sampling(n2,d2,f2);

    float[][][] synthAndSlope = FakeData.seismicAndSlopes2d2014B(noise,F);

    float[][] synth_data = synthAndSlope[0];
    float[][] exact_slope = synthAndSlope[1];
    //System.out.println("max slope= "+max(exact_slope));

    LocalSlopeFinder lsf = new LocalSlopeFinder(24.0f,1.0f,pmax);
    float[][] lsf_slope = new float[n2][n1];
    lsf.findSlopes(synth_data,lsf_slope);

    Sfdip sd = new Sfdip(-pmax,pmax);
    sd.setRect(76,6);
    sd.setOrder(2);
    sd.setNiter(niter);
    float[][] mad_slope_init = new float[n2][n1];
    float[][] mad_slope = new float[n2][n1];
    sd.findSlopes(s1,s2,synth_data,exact_slope,mad_slope_init);
    sd.findSlopes(s1,s2,synth_data,mad_slope);

    int shiftmax = 8;
    int ntraces = 10;
    DynamicWarping dw = new DynamicWarping(-shiftmax,shiftmax);
    dw.setShiftSmoothing(15.0,1.0);
    dw.setStrainMax(0.5,0.5); //only allow a max 50% stretch or squeeze
    float[][] dw_slope = DWSlopesCluster(dw,synth_data,ntraces);

    lsf_slope = mul(lsf_slope,d1/d2);
    mad_slope = mul(mad_slope,d1/d2);
    dw_slope = mul(dw_slope,d1/d2);

    float[][] lsf_diff = sub(lsf_slope,exact_slope);
    float[][] lsf_temp = pow(lsf_diff,2);
    float lsf_err = sqrt(sum(lsf_temp)/(n2*n1));

    float[][] mad_diff = sub(mad_slope,exact_slope);
    float[][] mad_temp = pow(mad_diff,2);
    float mad_err = sqrt(sum(mad_temp)/(n2*n1));

    float[][] dw_diff = sub(dw_slope,exact_slope);
    float[][] dw_temp = pow(dw_diff,2);
    float dw_err = sqrt(sum(dw_temp)/(n2*n1));

    System.out.println("Error lsf="+lsf_err);
    System.out.println("Error mad="+mad_err);
    System.out.println("Error dw="+dw_err);

    float fw = 0.75f; //fraction width for slide
    float fh = 0.9f; //fraction height for slide
    // title, paint, colorbar, color
    //Plot.plot(s1,s2,synth_data,
    //    "Synthetic seismic noise= "+noise,fw,fh,F,T,T,F);
    //Plot.plot(s1,s2,lsf_diff,dw_diff,
    //"Absolute differences with optimized "+ 
    //    "rect,sigma pmax= "+pmax+" and noise= "+noise,fw,fh,T);
    //Plot.plot(s1,s2,exact_slope,"Exact Slopes",fw,fh,-4.0f,4.0f,F,T,T,T);
    //Plot.plot(s1,s2,lsf_slope,"LOF Slopes noise= "+noise,fw,fh,-4,4,F,F,T,T);
    //Plot.plot(s1,s2,mad_slope_init,"PWD Slopes with initial values noise= "
    //    +noise+"iter= "+niter,fw,fh,-4,4,T,F,T,T);
    //Plot.plot(s1,s2,mad_slope,"PWD Slopes noise= "+noise,fw,fh,-4,4,F,F,T,T);
    //Plot.plot(s1,s2,mad_slope,"test PWD Slopes noise= "+noise,fw,fh,-4,4,F,T,T,T);
    //Plot.plot(s1,s2,exact_slope,lsf_slope,mad_slope,"Slopes noise= "+noise,fw,fh,T);
  }

 /**
  * Perhaps the shifts are too generic. For example, when noise is introduced,
  * we get the shift value to every single point in the image, including noise
  * points. We should instead do the entire DW process to get only the optimal
  * shift points to things that are not noise. The optimal shift path. Also,
  * why are the two ways I estimate slopes different? They should yield the
  * same results but don't.
  */ 
  private static void goDynamicWarpingSlopes() {
    int n1 = 501;
    int n2 = 501;

    float d1 = 1.0f;
    float d2 = 1.0f;

    float f1 = 0;
    float f2 = 0;
    float noise = 1.0f;

    Sampling s1 = new Sampling(n1,d1,f1);
    Sampling s2 = new Sampling(n2,d2,f2);

    float[][][] synthAndSlope = FakeData.seismicAndSlopes2d2014B(noise,F);

    float[][] synth_data = synthAndSlope[0];
    float[][] exact_slope = synthAndSlope[1];

    int shiftmax = 10;
    DynamicWarping dw = new DynamicWarping(-shiftmax,shiftmax);
    dw.setShiftSmoothing(10.0,1.0);
    dw.setStrainMax(0.5,0.5);

    float[][] dw_slope1 = new float[n2][n1];
    float[][] dw_slope = new float[n2][n1];
    float[][] synth_data_shifted = new float[n2][n1];
    synth_data_shifted[n2-1] = synth_data[n2-1];
    for (int i2=0; i2<n2-1; ++i2) {
      dw_slope[i2] = dw.findShifts(synth_data[i2],synth_data[i2+1]);
      synth_data_shifted[i2] = synth_data[i2+1];
    }
    dw_slope1 = dw.findShifts(synth_data,synth_data_shifted);
    System.out.println("real slope value= "+exact_slope[470][470]);
    System.out.println("dw slope value= "+dw_slope[470][470]);

    float fw = 0.75f; //fraction width for slide
    float fh = 0.9f; //fraction height for slide
    // title, paint, colorbar, color
    //Plot.plot(s1,s2,synth_data,"Synthetic seismic noise= "+noise,fw,fh,F,T,T,F);
    Plot.plot(s1,s2,exact_slope,"Exact Slopes",fw,fh,-4f,4f,T,F,T,T);
    Plot.plot(s1,s2,dw_slope,"1D DW Slopes noise= "+noise,fw,fh,-4f,4f,T,F,T,T);
    Plot.plot(s1,s2,dw_slope1,"2D DW Slopes noise= "+noise,fw,fh,-4f,4f,T,F,T,T);
  }

   private static void goDynamicWarpingSlopesTesting() {
    int n1 = 301;
    int shift = 10;
    float[] a = SlopeAlgorithmEval.readImage(n1,
        "/Users/earias/Home/git/ea/bench/src/slopes/data/chris.dat");
    float[] b = new float[n1];
    for (int i=0; i<n1-shift; ++i)
      b[i] = a[i+shift];

    int count = 0;
    for (int i=n1-shift; i<n1; ++i) {
      b[i] = a[count];
      ++count;
    }

    float[] c = new float[n1];
    c = addNoise(0.1,a);

    int shiftmax = 50;
    DynamicWarping dw = new DynamicWarping(-shiftmax,shiftmax);
    dw.setShiftSmoothing(1.0,1.0);
    //dw.setStrainMax(0.9,0.9);

    float[] dw_slope = new float[n1];
    dw_slope = dw.findShifts(a,c);

    for (int i=0; i<n1-3; ++i)
      System.out.println("dw_slope= "+dw_slope[i]);

    float fw = 0.75f; //fraction width for slide
    float fh = 0.9f; //fraction height for slide
    Plot.plot(a,c,"Curves","Index","Value",fw,fh,false);
    //Plot.plot(dw_slope,"Slope","Index","Value",fw,fh,false);
  }

  /**
   * need to compute errors compared to exact answer (std dev and rms).
   * need to compute optimal smoothing (and possible other) parameters 
   * similar to how I did for PWD and ST.
   * Need to make sure I can produce correct units for different values for 
   * different sampling values (try the GOM example)
   */
  private static void goDynamicWarpingSlopesCluster() {
    int n1 = 501;
    int n2 = 501;

    float d1 = 1.0f;
    float d2 = 1.0f;

    float f1 = 0;
    float f2 = 0;
    float noise = 0.0f;

    Sampling s1 = new Sampling(n1,d1,f1);
    Sampling s2 = new Sampling(n2,d2,f2);

    float[][][] synthAndSlope = FakeData.seismicAndSlopes2d2014B(noise,F);

    float[][] synth_data = synthAndSlope[0];
    float[][] exact_slope = synthAndSlope[1];

    int shiftmax = 8;
    int ntraces = 10;
    DynamicWarping dw = new DynamicWarping(-shiftmax,shiftmax);
    dw.setShiftSmoothing(15.0,1.0);
    dw.setStrainMax(0.5,0.5); //only allow a max 50% stretch or squeeze
    float[][] dw_slope = DWSlopesCluster(dw,synth_data,ntraces);

    //dw_slope = mul(dw_slope,d1/d2);
    float fw = 0.75f; //fraction width for slide
    float fh = 0.9f; //fraction height for slide
    // title, paint, colorbar, color
    Plot.plot(s1,s2,dw_slope,"DW Slopes noise= "+noise,fw,fh,-4f,4f,F,F,T,T);
    Plot.plot(s1,s2,exact_slope,"Exact Slopes",fw,fh,-4f,4f,F,F,T,T);
  }

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

    float[][] orig_data = readImage(n1,n2,"data/gom.dat");
    float[][] scaled_data = new float[n2][n1];
    mul(orig_data,.001f,scaled_data);
    writeBinary(scaled_data,"data/gom_scaled.dat");

    LocalSlopeFinder lsf = new LocalSlopeFinder(8.0f,4.0f,pmax);
    float[][] lsf_slope = new float[n2][n1];
    lsf.findSlopes(orig_data,lsf_slope);

    Sfdip sd = new Sfdip(-pmax,pmax);
    sd.setRect(12,6);
    float[][] mad_slope = new float[n2][n1];
    sd.findSlopes(s1,s2,orig_data,mad_slope);

    int shiftmax = 8;
    int ntraces = 10;
    DynamicWarping dw = new DynamicWarping(-shiftmax,shiftmax);
    dw.setShiftSmoothing(15.0,1.0);
    dw.setStrainMax(0.5,0.5); //only allow a max 50% stretch or squeeze
    float[][] dw_slope = DWSlopesCluster(dw,orig_data,ntraces);
    
    lsf_slope = mul(lsf_slope,d1/d2);
    mad_slope = mul(mad_slope,d1/d2);
    dw_slope = mul(dw_slope,d1/d2);

    // title, paint, colorbar, color
    //Plot.plot(s1,s2,sexp(scaled_data),"GOM Near Offset Data (Gained)",T,F,F,F);
    //Plot.plot(s1,s2,lsf_slope,"LSF Slopes",T,F,F,F);
    //Plot.plot(s1,s2,mad_slope,"Madagascar Slopes",T,F,F,F);
    //Plot.plot(s1,s2,dw_slope,"DW Slopes",T,F,F,F);
    //Plot.plot(s1,s2,sexp(scaled_data),lsf_slope,mad_slope,"teaser",F);
  }

  private static void goTestSigmaValuesDW() {
    int n1 = 501;
    int n2 = 501;

    float dsigma = 1.0f; // sigma sampling rate
    float lsigma = 80.0f; // last sigma
    float fsigma = 1.0f;
    int nsigma1 = (int)((lsigma-fsigma+1)/dsigma); // # of sigma vals to test
    int nsigma2 = 10;

    Sampling ssigma1 = new Sampling(nsigma1,dsigma,fsigma);
    Sampling ssigma2 = new Sampling(nsigma2,dsigma,fsigma);

    float d1 = 1;
    float d2 = 1;

    float f1 = 0;
    float f2 = 0;
    float pmax = 5.0f;

    int[] err_index = new int[2];
    float[][] err = new float[nsigma2][nsigma1]; // RMS error

    float[] sigma1 = new float[nsigma1];
    float[] sigma2 = new float[nsigma2];
    for(int i=0; i<nsigma1; ++i) {
      sigma1[i] = i*dsigma+fsigma;
    }
    for(int i=0; i<nsigma2; ++i) {
      sigma2[i] = i*dsigma+fsigma;
    }
    
    Sampling s1 = new Sampling(n1,d1,f1);
    Sampling s2 = new Sampling(n2,d2,f2);

    float noise = 0.5f;
    float[][][] synthAndSlope = FakeData.seismicAndSlopes2d2014B(noise,F);
    float[][] synth_data = synthAndSlope[0];
    float[][] exact_slope = synthAndSlope[1];
    float[][] dw_slope = new float[n2][n1];

    int shiftmax = 8;
    int ntraces = 10;
    DynamicWarping dw = new DynamicWarping(-shiftmax,shiftmax);
    dw.setStrainMax(0.5,0.5); //only allow a max 50% stretch or squeeze
    float[][] dw_diff;
    
    for(int i2=0; i2<nsigma2; ++i2) {
      for(int i1=0; i1<nsigma1; ++i1) {
        dw.setShiftSmoothing(sigma1[i1],sigma2[i2]);
        dw_slope = DWSlopesCluster(dw,synth_data,ntraces);
        dw_slope = mul(dw_slope,d1/d2);
        dw_diff = sub(dw_slope,exact_slope);
        float[][] temp = pow(dw_diff,2);
        err[i2][i1] = sqrt(sum(temp)/(n2*n1));
      }
    }
    writeBinary(err,
          "/Users/earias/Home/git/ea/bench/src/slopes/data/sigmaErrValsDW.dat");
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

    float[][] err = readImage(nsigma1,nsigma2,
          "/Users/earias/Home/git/ea/bench/src/slopes/data/sigmaErrValsDW.dat");
    float min_err = min(err,err_index);
    System.out.println("Sigma1= "+sigma1[err_index[0]]+
                      " Sigma2= "+sigma2[err_index[1]]+
                      " Minimum Error Value= "+min_err);

    float fw = 0.75f; //fraction width for slide
    float fh = 0.9f; //fraction height for slide
    // sigma, interp, title, paint, colorbar, color
    Plot.plot(ssigma1,ssigma2,err,"RMS Error Sigma",fw,fh,T,T,F,T,T,T);
  }

  private static void goTestSigmaValues() {
    int n1 = 501;
    int n2 = 501;

    float dsigma = 1.0f; // sigma sampling rate
    float lsigma = 80.0f; // last sigma
    float fsigma = 1.0f;
    int nsigma1 = (int)((lsigma-fsigma+1)/dsigma); // # of sigma vals to test
    int nsigma2 = 10;

    Sampling ssigma1 = new Sampling(nsigma1,dsigma,fsigma);
    Sampling ssigma2 = new Sampling(nsigma2,dsigma,fsigma);

    float d1 = 1;
    float d2 = 1;

    float f1 = 0;
    float f2 = 0;
    float pmax = 5.0f;

    int[] err_index = new int[2];
    float[][] err_vals = new float[nsigma2][nsigma1]; // RMS error

    float[] sigma1 = new float[nsigma1];
    float[] sigma2 = new float[nsigma2];
    for(int i=0; i<nsigma1; ++i) {
      sigma1[i] = i*dsigma+fsigma;
    }
    for(int i=0; i<nsigma2; ++i) {
      sigma2[i] = i*dsigma+fsigma;
    }
    
    Sampling s1 = new Sampling(n1,d1,f1);
    Sampling s2 = new Sampling(n2,d2,f2);

    float noise = 0.5f;
    float[][][] synthAndSlope = FakeData.seismicAndSlopes2d2014B(noise,F);
    float[][] synth_data = synthAndSlope[0];
    float[][] exact_slope = synthAndSlope[1];
    float[][] lsf_slope = new float[n2][n1];
    LocalSlopeFinder lsf;
    float[][] lsf_diff;
    
    for(int i2=0; i2<nsigma2; ++i2) {
      for(int i1=0; i1<nsigma1; ++i1) {
        lsf = new LocalSlopeFinder(sigma1[i1],sigma2[i2],pmax);
        lsf.findSlopes(synth_data,lsf_slope);
        lsf_slope = mul(lsf_slope,d1/d2);
        lsf_diff = sub(lsf_slope,exact_slope);
        float[][] temp = pow(lsf_diff,2);
        err_vals[i2][i1] = sqrt(sum(temp)/(n2*n1));
      }
    }
    writeBinary(err_vals,
          "/Users/earias/Home/git/ea/bench/src/slopes/data/sigmaErrVals.dat");
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

    float[][] err_vals = readImage(nsigma1,nsigma2,
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

  private static void goTestRectValues() {
    int n1 = 501;
    int n2 = 501;

    int frect = 1;
    int lrect = 80;
    int drect = 1;
    int nrect1 = (lrect-frect+1)/drect;
    int nrect2 = 10;

    Sampling srect1 = new Sampling(nrect1,drect,frect);
    Sampling srect2 = new Sampling(nrect2,drect,frect);

    float d1 = 1;
    float d2 = 1;

    float f1 = 0;
    float f2 = 0;
    float pmax = 5.0f;

    int[] err_index = new int[2];
    float[][] err_vals = new float[nrect2][nrect1]; // RMS error

    int[] rect1 = new int[nrect1];
    int[] rect2 = new int[nrect2];

    for(int i=0; i<nrect1; ++i) {
      rect1[i] = drect*i+frect;
    }
    for(int i=0; i<nrect2; ++i) {
      rect2[i] = drect*i+frect;
    }
    
    Sampling s1 = new Sampling(n1,d1,f1);
    Sampling s2 = new Sampling(n2,d2,f2);

    float noise = 0.5f;
    float[][][] synthAndSlope = FakeData.seismicAndSlopes2d2014B(noise,F);
    float[][] synth_data = synthAndSlope[0];
    float[][] exact_slope = synthAndSlope[1];

    float[][] mad_slope = new float[n2][n1];
    Sfdip sd = new Sfdip(-pmax,pmax);
    float[][] mad_diff;
    
    for(int i2=0; i2<nrect2; ++i2) {
      for(int i1=0; i1<nrect1; ++i1) {
        sd.setRect(rect1[i1],rect2[i2]);
        sd.findSlopes(s1,s2,synth_data,mad_slope);
        mad_slope = mul(mad_slope,d1/d2);
        mad_diff = sub(mad_slope,exact_slope);
        float[][] temp = pow(mad_diff,2);
        err_vals[i2][i1] = sqrt(sum(temp)/(n2*n1));
      }
    }
    writeBinary(err_vals,
          "/Users/earias/Home/git/ea/bench/src/slopes/data/rectErrVals.dat");
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

    float[][] err_vals = readImage(nrect1,nrect2,
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

    LocalSlopeFinder lsf = new LocalSlopeFinder(24.0f,1.0f,pmax);
    float[][] lsf_slope = new float[n2][n1];

    Sfdip sd = new Sfdip(-pmax,pmax);
    sd.setRect(76,6);
    sd.setOrder(2);
    float[][] mad_slope = new float[n2][n1];

    int shiftmax = 8;
    int ntraces = 10;
    DynamicWarping dw = new DynamicWarping(-shiftmax,shiftmax);
    dw.setShiftSmoothing(15.0,1.0);
    dw.setStrainMax(0.5,0.5); //only allow a max 50% stretch or squeeze
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
      dw_slope = DWSlopesCluster(dw,synth_data,ntraces);
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
        "PWD, LOF, and DW RMS error vs NS ratio",fw,fh,F);
  }

  /*
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
    float[][] lsf_mean = readImage(n1,n2,
        "/Users/earias/Home/git/ea/bench/src/pwd/data/lsf_mean.dat");
    float[][] mad_mean = readImage(n1,n2,
        "/Users/earias/Home/git/ea/bench/src/pwd/data/mad_mean.dat");
    float[][] lsf_sd = readImage(n1,n2,
        "/Users/earias/Home/git/ea/bench/src/pwd/data/lsf_sd.dat");
    float[][] mad_sd = readImage(n1,n2,
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
    writeBinary(exact_1Dslope,
        "/Users/earias/Home/git/ea/bench/src/pwd/data/exact_slope1D.dat");
    writeBinary(lsf_1Dslope,
        "/Users/earias/Home/git/ea/bench/src/pwd/data/struct_tens1D.dat");
    writeBinary(mad_1Dslope,
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
    float[][][] lsf_p = new float[N][n2][n1];
    float[][][] mad_p = new float[N][n2][n1];
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

      lsf_slope = mul(lsf_slope,d1/d2);
      mad_slope = mul(mad_slope,d1/d2);

      lsf_p[i] = lsf_slope;
      mad_p[i] = mad_slope;
      
      lsf_sum = add(lsf_sum,lsf_slope);
      mad_sum = add(mad_sum,mad_slope);
    }
    exact_slope = synthAndSlope[1];

    float[][] lsf_mean = div(lsf_sum,N);
    float[][] mad_mean = div(mad_sum,N);

    float[][] lsf_sumdiffsq = new float[n2][n1]; //sum of differences squared
    float[][] mad_sumdiffsq = new float[n2][n1]; //sum of differences squared
    for (int i=0; i<N; ++i) {
      lsf_sumdiffsq = add(lsf_sumdiffsq,pow(sub(lsf_p[i],exact_slope),2));
      mad_sumdiffsq = add(mad_sumdiffsq,pow(sub(mad_p[i],exact_slope),2));
      System.out.println("i= "+i);
    }
    float[][] lsf_sd= sqrt(div(lsf_sumdiffsq,N));
    float[][] mad_sd= sqrt(div(mad_sumdiffsq,N));
    writeBinary(lsf_mean,
          "/Users/earias/Home/git/ea/bench/src/pwd/data/lsf_mean.dat");
    writeBinary(mad_mean,
          "/Users/earias/Home/git/ea/bench/src/pwd/data/mad_mean.dat");
    writeBinary(lsf_sd,
          "/Users/earias/Home/git/ea/bench/src/pwd/data/lsf_sd.dat");
    writeBinary(mad_sd,
          "/Users/earias/Home/git/ea/bench/src/pwd/data/mad_sd.dat");
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

    float[][] lsf_mean = readImage(n1,n2,
        "/Users/earias/Home/git/ea/bench/src/pwd/data/lsf_mean.dat");
    float[][] mad_mean = readImage(n1,n2,
        "/Users/earias/Home/git/ea/bench/src/pwd/data/mad_mean.dat");
    float[][] lsf_sd = readImage(n1,n2,
        "/Users/earias/Home/git/ea/bench/src/pwd/data/lsf_sd.dat");
    float[][] mad_sd = readImage(n1,n2,
        "/Users/earias/Home/git/ea/bench/src/pwd/data/mad_sd.dat");
    float fw = 0.75f; //fraction width for slide
    float fh = 0.9f; //fraction height for slide
    // title, paint, colorbar, color
    //Plot.plot(s1,s2,lsf_mean,"Sample mean LOF",fw,fh,-4,4,F,T,T,T);
    //Plot.plot(s1,s2,mad_mean,"Sample mean PWD",fw,fh,-4,4,F,T,T,T);
    //Plot.plot(s1,s2,lsf_sd,"Sample variance LOF 03 max",fw,fh,0,0.3f,F,F,T,T);
    //Plot.plot(s1,s2,mad_sd,"Sample variance PWD",fw,fh,0,0.5f,F,F,T,T);
    Plot.plot(s1,s2,lsf_sd,"Sample standard deviation LOF",fw,fh,0,1,F,T,T,T);
    Plot.plot(s1,s2,mad_sd,"Sample standard deviation PWD",fw,fh,0,1,F,T,T,T);
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

  /**
   * Ungains the data input.
   * @param x array[n2][n1] of data to be ungained
   * @return array[n2][n1] of ungained data
   */
  private static float[][] slog(float[][] x) {
    return mul(sgn(x),sub(exp(abs(x)),1.0f));
  }

  /**
   * Gains the data input.
   * @param x array[n2][n1] of data to be gained
   * @return array[n2][n1] of gained data
   */
  private static float[][] sexp(float[][] x) {
    return mul(sgn(x),log(add(abs(x),1.0f)));
  }

  /**
   * Reads a binary file.
   * @param n1 the length of floats in the 1st-dimension
   * @param n2 the length of floats in the 2nd-dimension
   * @param fileName the name of the file to be read
   * @return array[n2][n1] of floats read from file
   */
  public static float[][] readImage(int n1, int n2, String fileName) {
    try {
      ArrayInputStream ais = new ArrayInputStream(fileName);
      float[][] x = new float[n2][n1];
      ais.readFloats(x);
      ais.close();
      return x;
    } catch (IOException e) {
      throw new RuntimeException(e);
    }
  }

  public static float[] readImage(int n1, String fileName) {
    try {
      ArrayInputStream ais = new ArrayInputStream(fileName);
      float[] x = new float[n1];
      ais.readFloats(x);
      ais.close();
      return x;
    } catch (IOException e) {
      throw new RuntimeException(e);
    }
  }

  /**
   * Writes seismic data to binary file.
   * @param x array[n2][n1] of data to write to the binary file
   * @param fileName name of output binary file
   */
  private static void writeBinary(float[][] x, String fileName) {
    ByteOrder byteOrder = ByteOrder.BIG_ENDIAN;
    try {
      ArrayOutputStream aos = new ArrayOutputStream(fileName,byteOrder);
      aos.writeFloats(x);
      aos.close();
    } catch (IOException e) {
      throw new RuntimeException(e);
    }
  }

  /**
   * Writes seismic data to binary file.
   * @param x array[n2][n1] of data to write to the binary file
   * @param fileName name of output binary file
   */
  public static void writeBinary(float[] x, String fileName) {
    ByteOrder byteOrder = ByteOrder.BIG_ENDIAN;
    try {
      ArrayOutputStream aos = new ArrayOutputStream(fileName,byteOrder);
      aos.writeFloats(x);
      aos.close();
    } catch (IOException e) {
      throw new RuntimeException(e);
    }
  }

  public static final String DATA_NAME = "gom_scaled.dat";
  public static final boolean T = true;
  public static final boolean F = false;
  public static final float pi = FLT_PI;
  public static void main(String[] args) {
    SwingUtilities.invokeLater(new Runnable() {
      public void run() {
        // have user input argument to show which plot they want to see
        // e.g. 1 for synthetic 
        //      2 for std dev....
        //goSynth();
        //goDynamicWarpingSlopes();
        //goDynamicWarpingSlopesTesting();
        //goDynamicWarpingSlopesCluster();
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
        goTestSigmaValuesDW();
        goPlotSigmaValuesDW();
        //goTestPmaxValues();
        //goGom();
        //go3D();
      }
    });
  }
}
