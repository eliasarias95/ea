package util;

import edu.mines.jtk.mosaic.*;
import edu.mines.jtk.awt.ColorMap;
import edu.mines.jtk.dsp.*;

import static edu.mines.jtk.util.ArrayMath.*;

import javax.swing.*;

/**
 * Plotting class to hold various plot methods.
 * @author Elias Arias, Colorado School of Mines, CWP
 * @version 19.01.2015
 */
public class Plot {

  private static final double ratio = 16.0/9.0;

  public static void plot(float[] f1, String title, String hl,
      String vl, float fw, float fh, boolean print) {
    int fwi = round(1920*fw/2+1);
    int fhi = round(1080*fh/2+1);
    PointsView pv = new PointsView(f1);
    pv.setLineColor(java.awt.Color.RED);
    pv.setLineWidth(5);
    pv.setOrientation(PointsView.Orientation.X1DOWN_X2RIGHT);
    PlotPanel pp = new PlotPanel();
    pp.addTiledView(pv);
    //pp.setTitle(title);
    pp.setHLabel(hl);
    pp.setVLabel(vl);
    pp.setColorBarWidthMinimum(70);
    //pp.setHInterval(0.5);
    //pp.setVInterval(0.2);
    //pp.setVLimits(0,1);
    PlotFrame pf = new PlotFrame(pp);
    pf.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
    pf.setFontSizeForSlide(fw,fh,ratio);
    pf.setVisible(true);
    if (print) {
      int dpi = 720;
      pf.paintToPng(dpi,(1920f*fw-1)/dpi,
          "/Users/earias/Documents/research/figures/presentation/"+title+".png");
    }
  }

  public static void plot(float[] f1, float[] f2, String title, String hl,
      String vl, float fw, float fh, boolean print) {
    int fwi = round(1920*fw/2+1);
    int fhi = round(1080*fh/2+1);
    PointsView pv1 = new PointsView(f1);
    PointsView pv2 = new PointsView(f2);
    pv1.setLineColor(java.awt.Color.RED);
    pv1.setLineWidth(5);
    pv2.setLineWidth(5);
    pv1.setOrientation(PointsView.Orientation.X1DOWN_X2RIGHT);
    pv2.setOrientation(PointsView.Orientation.X1DOWN_X2RIGHT);
    PlotPanel pp1 = new PlotPanel();
    PlotPanel pp2 = new PlotPanel();
    pp1.addTiledView(pv1);
    pp2.addTiledView(pv2);
    //pp.setTitle(title);
    pp1.setHLabel(hl);
    pp1.setVLabel(vl);
    pp1.setColorBarWidthMinimum(70);
    pp2.setHLabel(hl);
    pp2.setVLabel(vl);
    pp2.setColorBarWidthMinimum(70);
    //pp.setHInterval(0.5);
    //pp.setVInterval(0.2);
    //pp.setVLimits(0,1);
    PlotFrame pf = new PlotFrame(pp1,pp2,PlotFrame.Split.HORIZONTAL);
    pf.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
    pf.setFontSizeForSlide(fw,fh,ratio);
    pf.setVisible(true);
    if (print) {
      int dpi = 720;
      pf.paintToPng(dpi,(1920f*fw-1)/dpi,
          "/Users/earias/Documents/research/figures/presentation/"+title+".png");
    }
  }

  public static void plotHistogram(float[] f1, float[] f2, float[] f3,
      String title, float fw, float fh, boolean print) {
    int n = f1.length;
    int fwi = round(1920*fw/2+1);
    int fhi = round(1080*fh/2+1);
    Histogram f1_hist = new Histogram(f1);
    float[] f1_dens = mul(n,f1_hist.getDensities());
    Histogram f2_hist = new Histogram(f2);
    float[] f2_dens = mul(n,f2_hist.getDensities());
    Histogram f3_hist = new Histogram(f3);
    float[] f3_dens = mul(n,f3_hist.getDensities());

    SequenceView sv1 = new SequenceView(f1_dens);
    sv1.setColor(java.awt.Color.BLACK);
    SequenceView sv2 = new SequenceView(f2_dens);
    sv2.setColor(java.awt.Color.RED);
    SequenceView sv3 = new SequenceView(f3_dens);
    sv3.setColor(java.awt.Color.BLUE);

    PlotPanel pp = new PlotPanel();
    pp.addTiledView(sv1);
    pp.addTiledView(sv2);
    pp.addTiledView(sv3);
    pp.setHLabel("Slope (deg)");
    pp.setVLabel("Frequency");
    PlotFrame pf = new PlotFrame(pp);
    pf.setVisible(true);
    if (print) {
      int dpi = 720;
      pf.paintToPng(dpi,(1920f*fw-1)/dpi,
          "/Users/earias/Documents/research/figures/presentation/"+title+".png");
    }
  }

  public static void plot(float[] f1, float[] f2, float[] f3, float[] f4, 
      String title, float fw, float fh, boolean print) {
    int fwi = round(1920*fw/2+1);
    int fhi = round(1080*fh/2+1);
    PointsView pv1 = new PointsView(f1,f2);
    PointsView pv2 = new PointsView(f1,f3);
    PointsView pv3 = new PointsView(f1,f4);
    pv1.setLineColor(java.awt.Color.RED);
    pv1.setLineWidth(5);
    pv2.setLineColor(java.awt.Color.BLUE);
    pv2.setLineWidth(5);
    pv3.setLineColor(java.awt.Color.GREEN);
    pv3.setLineWidth(5);
    PlotPanel pp = new PlotPanel();
    pp.addTiledView(pv1);
    pp.addTiledView(pv2);
    pp.addTiledView(pv3);
    //pp.setTitle(title);
    pp.setHLabel("Noise/signal");
    pp.setVLabel("RMS error");
    pp.setColorBarWidthMinimum(70);
    pp.setHInterval(0.5);
    pp.setVInterval(0.2);
    pp.setVLimits(0,1);
    PlotFrame pf = new PlotFrame(pp);
    pf.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
    pf.setFontSizeForSlide(fw,fh,ratio);
    pf.setVisible(true);
    if (print) {
      int dpi = 720;
      pf.paintToPng(dpi,(1920f*fw-1)/dpi,
          "/Users/earias/Documents/research/figures/presentation/"+title+".png");
    }
  }

  /**
   * Plots a 2D array of floats with specified title.
   * @param s1 the sampling in the 1st-dimension
   * @param s2 the sampling in the 2nd-dimension
   * @param f the 2D array of floats to be plotted.
   * @param title the title of the image generated.
   */
  public static void plot(float[][] f, String title, float fw, float fh,
      boolean ttl, boolean paint, boolean cb, boolean color) {
    int fwi = round(1920*fw/2+1);
    int fhi = round(1080*fh/2+1);
    PlotPanel pp = new PlotPanel(1,1,PlotPanel.Orientation.X1DOWN_X2RIGHT);
    PixelsView pv = pp.addPixels(f);
    pv.setOrientation(PixelsView.Orientation.X1DOWN_X2RIGHT);
    //pv.setClips(-6.0f,6.0f);
    pp.addTiledView(pv);
    pp.setVLabel("Samples");
    pp.setHLabel("Traces");
    pp.setColorBarWidthMinimum(70);
    PlotFrame pf = new PlotFrame(pp);
    pf.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
    pf.setSize(fwi,fhi);
    pf.setFontSizeForSlide(fw,fh,ratio);
    pf.setVisible(true);
    if (color) pv.setColorModel(ColorMap.JET);
    if (ttl) pp.setTitle(title);
    if (cb) pp.addColorBar();
    if (paint) { 
      int dpi = 720;
      pf.paintToPng(dpi,(1920f*fw-1)/dpi,
          "/Users/earias/Documents/research/figures/presentation/"+title+".png");
    }
  }

  /**
   * Plots a 2D array of floats with specified title.
   * @param s1 the sampling in the 1st-dimension
   * @param s2 the sampling in the 2nd-dimension
   * @param f the 2D array of floats to be plotted.
   * @param title the title of the image generated.
   */
  public static void plot(Sampling s1, Sampling s2, String title, float[][] f,
      float fw, float fh, boolean paint) {
    int fwi = round(1920*fw/2+1);
    int fhi = round(1080*fh/2+1);
    PlotPanel pp = new PlotPanel(1,1,PlotPanel.Orientation.X1DOWN_X2RIGHT);
    PixelsView pv = pp.addPixels(s1,s2,f);
    pv.setOrientation(PixelsView.Orientation.X1RIGHT_X2UP);
    pp.addTiledView(pv);
    pp.setColorBarWidthMinimum(100);
    pp.setHLabel("a");
    pp.setVLabel("b");
    PlotFrame pf = new PlotFrame(pp);
    pf.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
    pf.setSize(fwi,fhi);
    pf.setFontSizeForSlide(fw,fh,ratio);
    pf.setVisible(true);
    pv.setInterpolation(PixelsView.Interpolation.NEAREST);
    pv.setColorModel(ColorMap.JET);
    pv.setClips(0f,7000f);
    pp.setTitle(title);
    pp.addColorBar("Distance (samples)");
    if (paint) { 
      int dpi = 720;
      pf.paintToPng(dpi,(1920f*fw-1)/dpi,
          "/Users/earias/Documents/research/figures/presentation/"+title+".png");
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
      float fw, float fh,
      boolean ttl, boolean paint, boolean cb, boolean color) {
    int fwi = round(1920*fw/2+1);
    int fhi = round(1080*fh/2+1);
    PlotPanel pp = new PlotPanel(1,1,PlotPanel.Orientation.X1DOWN_X2RIGHT);
    PixelsView pv = pp.addPixels(s1,s2,f);
    pv.setOrientation(PixelsView.Orientation.X1DOWN_X2RIGHT);
    pp.addTiledView(pv);
    pp.setVLabel("Samples");
    pp.setHLabel("Traces");
    pp.setColorBarWidthMinimum(90);
    PlotFrame pf = new PlotFrame(pp);
    pf.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
    pf.setSize(fwi,fhi);
    pf.setFontSizeForSlide(fw,fh,ratio);
    pf.setVisible(true);
    if (color) pv.setColorModel(ColorMap.JET);
    if (ttl) pp.setTitle(title);
    if (cb) pp.addColorBar("slope (samples/trace)");
    if (paint) { 
      int dpi = 720;
      pf.paintToPng(dpi,(1920f*fw-1)/dpi,
          "/Users/earias/Documents/research/figures/presentation/"+title+".png");
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
      String hl, String vl, float fw, float fh,
      boolean interp, boolean ttl, boolean paint, boolean cb, boolean color) {
    int fwi = round(1920*fw/2+1);
    int fhi = round(1080*fh/2+1);
    PlotPanel pp = new PlotPanel(1,1,PlotPanel.Orientation.X1DOWN_X2RIGHT);
    PixelsView pv = pp.addPixels(s1,s2,f);
    pv.setOrientation(PixelsView.Orientation.X1DOWN_X2RIGHT);
    pp.addTiledView(pv);
    pp.setHLabel(hl);
    pp.setVLabel(vl);
    //pv.setClips(0.665f,0.70f);
    pp.setColorBarWidthMinimum(120);
    PlotFrame pf = new PlotFrame(pp);
    pf.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
    pf.setSize(fwi,fhi);
    pf.setFontSizeForSlide(fw,fh,ratio);
    pf.setVisible(true);
    if (interp) pv.setInterpolation(PixelsView.Interpolation.NEAREST);
    if (color) pv.setColorModel(ColorMap.JET);
    if (ttl) pp.setTitle(title);
    if (cb) pp.addColorBar("rms error (samples/trace)");
    if (paint) { 
      int dpi = 720;
      pf.paintToPng(dpi,(1920f*fw-1)/dpi,
          "/Users/earias/Documents/research/figures/presentation/"+title+".png");
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
      float fw, float fh, float clipMin, float clipMax,
      boolean ttl, boolean paint, boolean cb, boolean color) {
    int fwi = round(1920*fw/2+1);
    int fhi = round(1080*fh/2+1);
    PlotPanel pp = new PlotPanel(1,1,PlotPanel.Orientation.X1DOWN_X2RIGHT);
    PixelsView pv = pp.addPixels(s1,s2,f);
    pv.setOrientation(PixelsView.Orientation.X1DOWN_X2RIGHT);
    //pv.setClips(clipMin,clipMax);
    pp.addTiledView(pv);
    pp.setVLabel("Samples");
    pp.setHLabel("Traces");
    pp.setColorBarWidthMinimum(90);
    PlotFrame pf = new PlotFrame(pp);
    pf.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
    pf.setSize(fwi,fhi);
    pf.setFontSizeForSlide(fw,fh,ratio);
    pf.setVisible(true);
    if (color) pv.setColorModel(ColorMap.JET);
    if (ttl) pp.setTitle(title);
    if (cb) pp.addColorBar("slope (samples/trace)");
    if (paint) { 
      int dpi = 720;
      pf.paintToPng(dpi,(1920f*fw-1)/dpi,
          "/Users/earias/Documents/research/figures/presentation/"+title+".png");
    }
  }

  /**
   * Plots a 2D array of floats with specified title.
   * @param s1 the sampling in the 1st-dimension
   * @param s2 the sampling in the 2nd-dimension
   * @param f the 2D array of floats to be plotted.
   * @param title the title of the image generated.
   */
  public static void plot(Sampling s1, Sampling s2, float[][] f, float[] f1, 
      float[] f2, String title, float fw, float fh, boolean paint) {
    int fwi = round(1920*fw/2+1);
    int fhi = round(1080*fh/2+1);
    PlotPanel pp = new PlotPanel(1,1,PlotPanel.Orientation.X1DOWN_X2RIGHT);
    PixelsView pv = pp.addPixels(s1,s2,f);
    PointsView pvl = new PointsView(f1,f2);
    pvl.setLineWidth(3);
    pv.setOrientation(PixelsView.Orientation.X1RIGHT_X2UP);
    pp.addTiledView(pv);
    pp.addTiledView(pvl);
    pp.setColorBarWidthMinimum(100);
    pp.setHLabel("a");
    pp.setVLabel("b");
    PlotFrame pf = new PlotFrame(pp);
    pf.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
    pf.setSize(fwi,fhi);
    pf.setFontSizeForSlide(fw,fh,ratio);
    pf.setVisible(true);
    pv.setInterpolation(PixelsView.Interpolation.NEAREST);
    pv.setColorModel(ColorMap.JET);
    //pv.setClips(0f,50f);
    pp.setTitle(title);
    pp.addColorBar("Distance (samples)");
    if (paint) { 
      int dpi = 720;
      pf.paintToPng(dpi,(1920f*fw-1)/dpi,
          "/Users/earias/Documents/research/figures/presentation/"+title+".png");
    }
  }

  /**
   * Plots a 2D array of floats with specified title.
   * @param s1 the sampling in the 1st-dimension
   * @param s2 the sampling in the 2nd-dimension
   * @param f array[n2][n1] of floats to be plotted
   * @param title the title of the image generated
   */
  public static void plot(Sampling s1, Sampling s2, 
      float[][] f1, float[][] f2, String title, float fw, float fh, 
      boolean print) {
    int fwi = round(1920*fw/2+1);
    int fhi = round(1080*fh/2+1);
    PlotPanel pp = new PlotPanel(1,2,PlotPanel.Orientation.X1DOWN_X2RIGHT);
    PixelsView pv1 = pp.addPixels(0,0,s1,s2,f1);
    PixelsView pv2 = pp.addPixels(0,1,s1,s2,f2);
    pv1.setColorModel(ColorMap.RED_WHITE_BLUE);
    pv2.setColorModel(ColorMap.RED_WHITE_BLUE);
    pv1.setClips(-1.5f,1.5f);
    pv2.setClips(-1.5f,1.5f);
    pp.addColorBar("Slope error (samples/trace)");
    pp.setColorBarWidthMinimum(90);
    pp.setHLabel(0,"Traces");
    pp.setVLabel(0,"Samples");
    pp.setHLabel(1,"Traces");
    pp.setVLabel(0,"Samples");
    pp.getMosaic().getTileAxisTop(0).setFormat("%1.6G");
    PlotFrame pf = new PlotFrame(pp);
    pf.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE);
    pf.setSize(fwi,fhi);
    pf.setFontSizeForSlide(fw,fh,ratio);
    pf.setVisible(true);
    if (print) {
      int dpi = 720;
      pf.paintToPng(dpi,(1920f*fw-1)/dpi,
          "/Users/earias/Documents/research/figures/presentation/"+title+".png");
    }
  }

  /**
   * Plots a 2D array of floats with specified title.
   * @param s1 the sampling in the 1st-dimension
   * @param s2 the sampling in the 2nd-dimension
   * @param f array[n2][n1] of floats to be plotted
   * @param title the title of the image generated
   */
  public static void plot(Sampling s1, Sampling s2, 
      float[][] f1, float[][] f2, float[][] f3, String title, 
      float fw, float fh, boolean print) {
    int fwi = round(1920*fw/2+1);
    int fhi = round(1080*fh/2+1);
    PlotPanel pp = new PlotPanel(1,3,PlotPanel.Orientation.X1DOWN_X2RIGHT);
    PixelsView pv1 = pp.addPixels(0,0,s1,s2,f1);
    PixelsView pv2 = pp.addPixels(0,1,s1,s2,f2);
    PixelsView pv3 = pp.addPixels(0,2,s1,s2,f3);
    pv1.setColorModel(ColorMap.JET);
    pv2.setColorModel(ColorMap.JET);
    pv3.setColorModel(ColorMap.JET);
    //pv1.setClips(-3.0f,3.0f);
    pv2.setClips(-3.0f,3.0f);
    pv3.setClips(-3.0f,3.0f);
    pp.addColorBar("Slope (samples/traces)");
    pp.setColorBarWidthMinimum(100);
    pp.setHLabel(0,"Traces");
    pp.setVLabel(0,"Samples");
    pp.setHLabel(1,"Traces");
    pp.setVLabel(0,"Samples");
    pp.setHLabel(2,"Traces");
    pp.setVLabel(0,"Samples");
    pp.getMosaic().getTileAxisTop(0).setFormat("%1.6G");
    PlotFrame pf = new PlotFrame(pp);
    pf.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE);
    pf.setSize(fwi,fhi);
    pf.setFontSizeForSlide(fw,fh,ratio);
    pf.setVisible(true);
    if (print) {
      int dpi = 720;
      pf.paintToPng(dpi,(1920f*fw-1)/dpi,
          "/Users/earias/Documents/research/figures/presentation/"+title+".png");
    }
  }
}
