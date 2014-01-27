package utils;

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
   * Creates a synthetic 2D seismic image with a specified shift in each
   * following trace using a supplied trace.
   * @param trc an array of floats representing a single seismic trace.
   * @param n_traces the number of traces you would like for your image.
   * @param shift the amount of samples you would like each trace
   * shifted from the trace before it.
   * @return a 2D synthetic seismic image.
   */
  public static float[][] make2D(float[] trc, int n_traces, int shift) {
    int n1 = trc.length;
    float[][] synth = new float[n_traces][n1];
    synth[0] = trc;
    for (int i1=1; i1<n_traces; ++i1)
      synth[i1] = shift(synth[i1-1],shift);
    return synth;
  }

  /**
   * Shifts a given trace by a number of samples specified by the user.
   * @param trc an array of floats representing the trace to be shifted.
   * @param shift the amount of samples you would like each trace
   * shifted from the trace before it.
   * @return the same trace shifted by the # of samples specified.
   */
  public static float[] shift(float[] trc, int shift) {
    int n1 = trc.length;
    int j1 = 0;
    float[] strc = new float[n1];
    for (int i1=0; i1<n1; ++i1){
      j1 = i1-shift;
      if (j1<0) j1 += n1;
      if (j1<n1)
        strc[i1] = trc[j1];
    }
    return strc;
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
    PlotPanel pp = new PlotPanel();
    pp.addTiledView(pv);
    pp.setTitle(title);
    pp.setHLabel("# Traces");
    pp.setVLabel("# Samples");
    pp.addColorBar();
    pp.setColorBarWidthMinimum(100);
    PlotFrame pf = new PlotFrame(pp);
    pf.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
    pf.setVisible(true);
  }

  //////////////////////////////////////////////////////////////////////////
  // Private

  private static float ricker(float fpeak, int time) {
    float x = FLT_PI*fpeak*time;
    return (1.0f-2.0f*x*x)*exp(-x*x);
  }

}
