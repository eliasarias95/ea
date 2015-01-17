package util;

import java.util.Random;
import java.io.IOException;

import javax.swing.*;

import dnp.*;
import edu.mines.jtk.dsp.*;
import edu.mines.jtk.awt.*;
import edu.mines.jtk.mosaic.*;
import edu.mines.jtk.io.ArrayInputStream;
import edu.mines.jtk.io.ArrayOutputStream;

import static edu.mines.jtk.util.ArrayMath.*;

/**
 * @author Stefan Compton
 * @edited by Elias Arias
 */
public class Synthetic {

  /**
   * Writes seismic data to binary file.
   * @param x array[n2][n1] of data to write to the binary file
   * @param fileName name of output binary file
   */
  public static void writeBinary(float[][] x, String fileName) {
    try {
      ArrayOutputStream aos = new ArrayOutputStream(fileName);
      aos.writeFloats(x);
      aos.close();
    } catch (IOException e) {
      throw new RuntimeException(e);
    }
  }

  /**
   * Get a synthetic trace with the given Sampling and peak frequency.
   * @param s Sampling of synthetic trace.
   * @param fpeak peak frequency in Hz.
   * @param seed for random events.
   * @return a synthetic trace.
   */
  public static float[] getTrace(Sampling s, float fpeak, long seed) {
    int n = s.getCount();
    double d = s.getDelta();
    fpeak = fpeak*(float)d; // cycles/s * s/sample
    float[] f = makeRandomEvents(n,seed);
    return addRickerWavelet(fpeak,f);
  }

  /**
   * Creates 2D synthetic seismic data that follows a sine wave shape,
   * (time down, distance right).
   * @param st the time sampling
   * @param trc the array[3*nt] of input floats to represent a trace 
   * @param A the amplitude (0 to peak/trough) of the sine wave you want
   * @param nx the number of samples in the x direction
   */
  public static float[][] sineData(Sampling st, Sampling sx, float ft,
      float[] trc, float A, float[][] slope) {
    int nx = sx.getCount();
    int nt = st.getCount();
    int ntrc = trc.length;
    double dx = sx.getDelta();
    double dt = st.getDelta();
    //float interval = (float)();
    float[][] synth = new float[nx][nt];
    SincInterp si = new SincInterp();
    Sampling sst;
    float in;
    float b;
    float f1 = 3.0f*(float)dt;
    float f2 = 12.0f*(float)dt;
    for (int ix=0; ix<nx; ++ix) {
      b = 2.0f*pi*(f2-f1)/nx;
      in = -b*ix*ix/10.0f;
      sst = new Sampling(ntrc,dt,A*sin(in)+ft);
      si.interpolate(sst,trc,st,synth[ix]);
      for (int it=0; it<nt; ++it) {
        slope[ix][it] = -2.0f*b*ix*A*cos(in)/10.0f;
      }
    }
    return synth;
  }

  /**
   * Creates a synthetic 2D seismic image with a specified shift in each
   * following trace using a supplied trace.
   * @param trc an array of floats representing a single seismic trace.
   * @param n_traces the number of traces you would like for your image.
   * @param shift the amount of samples you would like each trace
   * shifted from the trace before it.
   * @return a 2D synthetic seismic image.
   */
  public static float[][] make2D(Sampling st, Sampling sx, float ft, 
      float[] trc, float shift, float[][] slope) {
    int nx = sx.getCount();
    int nt = st.getCount();
    int ntrc = trc.length;
    double dx = sx.getDelta();
    double dt = st.getDelta();
    float[][] synth = new float[nx][nt];
    SincInterp si = new SincInterp();
    Sampling sst;
    for (int ix=0; ix<nx; ++ix) {
      sst = new Sampling(ntrc,dt,shift*(float)ix/(nx-1)+ft);
      si.interpolate(sst,trc,st,synth[ix]);
      for (int it=0; it<nt; ++it) {
        slope[ix][it] = shift/(nx-1);
      }
    }
    return synth;
  }

  /**
   * Get random events for a synthetic trace. 
   * These are just impulses without a wavelet.
   * @param n length of the output array.
   * @param seed for random events.
   * @return an array of random impulse responses.
   */
  private static float[] makeRandomEvents(int n, long seed) {
    Random r = new Random(seed);
    return pow(mul(2.0f,sub(randfloat(r,n),0.5f)),15.0f);
  }

  /**
   * Adds a ricker wavelet with specified peak frequency to array 
   * <code>f</code>.
   * @param fpeak peak frequency in Cycles/Sample.
   * @param f array of impulse responses.
   * @return an array of random impulse responses convolved with a ricker
   *   wavelet.
   */
  private static float[] addRickerWavelet(float fpeak, float[] f) {
    int n = f.length;
    int ih = (int)(3.0f/fpeak);
    int nh = 1+2*ih;
    float[] h = new float[nh];
    for (int jh=0; jh<nh; jh++)
      h[jh] = ricker(fpeak,jh-ih);
    float[] g = new float[n];
    Conv.conv(nh,-ih,h,n,0,f,n,0,g);
    return g;
  }

  /**
   * Plots a 2D array of floats with specified title.
   * @param f the 2D array of floats to be plotted.
   * @param title the title of the image generated.
   */
  public static void plot(float[][] f, String title) {
    PixelsView pv = new PixelsView(f);
    pv.setOrientation(PixelsView.Orientation.X1DOWN_X2RIGHT);
    pv.setColorModel(ColorMap.JET);
    //pv.setClips(1.99f,2.01f);
    PlotPanel pp = new PlotPanel();
    pp.addTiledView(pv);
    pp.setTitle(title);
    pp.setHLabel("# Traces");
    pp.setVLabel("# Samples");
    pp.addColorBar();
    pp.setColorBarWidthMinimum(70);
    PlotFrame pf = new PlotFrame(pp);
    pf.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
    pf.setVisible(true);
  }

  /**
   * Plots a 2D array of floats with specified title.
   * @param s1 the sampling in the 1st-dimension
   * @param s2 the sampling in the 2nd-dimension
   * @param f the 2D array of floats to be plotted.
   * @param title the title of the image generated.
   */
  public static void plot(float[][] f, String title,
      boolean ttl, boolean paint, boolean cb, boolean color) {
    PlotPanel pp = new PlotPanel(1,1,PlotPanel.Orientation.X1DOWN_X2RIGHT);
    PixelsView pv = pp.addPixels(f);
    pv.setOrientation(PixelsView.Orientation.X1DOWN_X2RIGHT);
    //pv.setClips(-0.1f,0.1f);
    pp.addTiledView(pv);
    pp.setHLabel("Distance (km)");
    pp.setVLabel("Time (s)");
    pp.setColorBarWidthMinimum(70);
    PlotFrame pf = new PlotFrame(pp);
    pf.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
    pf.setVisible(true);
    if (color) pv.setColorModel(ColorMap.JET);
    if (ttl) pp.setTitle(title);
    if (cb) pp.addColorBar();
    if (paint) { 
      pf.paintToPng(720,6.51,
          "/Users/earias/Documents/latex/figures/"+title+".png");
    }
  }

  /**
   * Plots a 2D array of floats with specified title.
   * @param s1 the sampling in the 1st-dimension
   * @param s2 the sampling in the 2nd-dimension
   * @param f the 2D array of floats to be plotted.
   * @param title the title of the image generated.
   */
  public static void plot(Sampling s1, Sampling s2, float[][] f, String title,
      boolean ttl, boolean paint, boolean cb, boolean color) {
    PlotPanel pp = new PlotPanel(1,1,PlotPanel.Orientation.X1DOWN_X2RIGHT);
    PixelsView pv = pp.addPixels(s1,s2,f);
    pv.setOrientation(PixelsView.Orientation.X1DOWN_X2RIGHT);
    //pv.setClips(-0.1f,0.1f);
    pp.addTiledView(pv);
    pp.setHLabel("Distance (km)");
    pp.setVLabel("Time (s)");
    pp.setColorBarWidthMinimum(70);
    PlotFrame pf = new PlotFrame(pp);
    pf.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
    pf.setVisible(true);
    if (color) pv.setColorModel(ColorMap.JET);
    if (ttl) pp.setTitle(title);
    if (cb) pp.addColorBar();
    if (paint) { 
      pf.paintToPng(720,6.51,
          "/Users/earias/Documents/latex/figures/"+title+".png");
    }
  }
  
  /**
   * Plots a 2D array of floats with specified title.
   * @param s1 the sampling in the 1st-dimension
   * @param s2 the sampling in the 2nd-dimension
   * @param f array[n2][n1] of floats to be plotted
   * @param title the title of the image generated
   */
  private static void plot(Sampling s1, Sampling s2, 
      float[][] f1, float[][] f2, float[][] f3, String title) {
    PlotPanel pp = new PlotPanel(1,3,PlotPanel.Orientation.X1DOWN_X2RIGHT);
    PixelsView pv1 = pp.addPixels(0,0,s1,s2,f1);
    PixelsView pv2 = pp.addPixels(0,1,s1,s2,f2);
    PixelsView pv3 = pp.addPixels(0,2,s1,s2,f3);
    pv2.setColorModel(ColorMap.JET);
    pv3.setColorModel(ColorMap.JET);
    pp.addColorBar("Slope (samples/trace)");
    pp.setColorBarWidthMinimum(70);
    pp.setHLabel(0,"Distance (km)");
    pp.setVLabel(0,"Time (s)");
    pp.setHLabel(1,"Distance (km)");
    pp.setVLabel(0,"Time (s)");
    pp.setHLabel(2,"Distance (km)");
    pp.setVLabel(0,"Time (s)");
    pp.getMosaic().getTileAxisTop(0).setFormat("%1.6G");
    PlotFrame pf = new PlotFrame(pp);
    pf.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE);
    pf.setFontSizeForPrint(8.0,469.0);
    pf.setSize(1200,400);
    pf.setVisible(true);
    //pf.paintToPng(720,6.51,"/Users/earias/Documents/latex/figures/"
    //    +title+".png");
  }

  public static void plot(float[] f1, float[] f2, String title) {
    PointsView pv = new PointsView(f1,f2);
    PlotPanel pp = new PlotPanel();
    pp.addTiledView(pv);
    pp.setTitle(title);
    pp.setHLabel("# Traces");
    pp.setVLabel("# Samples");
    pp.setColorBarWidthMinimum(70);
    PlotFrame pf = new PlotFrame(pp);
    pf.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
    pf.setFontSizeForPrint(8.0,469.0);
    pf.setVisible(true);
  }

  //////////////////////////////////////////////////////////////////////////
  // Private

  private static final float pi = FLT_PI;

  private static float ricker(float fpeak, int time) {
    float x = FLT_PI*fpeak*time;
    return (1.0f-2.0f*x*x)*exp(-x*x);
  }

}
