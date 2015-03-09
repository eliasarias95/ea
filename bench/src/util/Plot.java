package util;

import edu.mines.jtk.mosaic.*;
import edu.mines.jtk.awt.ColorMap;
import edu.mines.jtk.dsp.*;
import edu.mines.jtk.sgl.*;

import static edu.mines.jtk.util.ArrayMath.*;

import java.text.DecimalFormat;
import javax.swing.*;

/**
 * Plotting class to hold various plot methods.
 * @author Elias Arias, Colorado School of Mines, CWP
 * @version 26.01.2015
 */
public class Plot {

  private static final double _ratio = 16.0/9.0;
  private static final String _paths = 
    "/Users/earias/Documents/research/figures/slides/";
  private static final String _pathp = 
    "/Users/earias/Documents/research/figures/paper/";

  public static void plot(Sampling s, float[] f1, float[] f2, float[] f3, 
      String title, String hl, String vl, 
      float fw, float fh, float cmin, float cmax,
      boolean slide, boolean one, boolean paint) {
    PointsView pv1 = new PointsView(s,f1);
    PointsView pv2 = new PointsView(s,f2);
    PointsView pv3 = new PointsView(s,f3);
    pv1.setLineColor(java.awt.Color.RED);
    pv1.setLineWidth(5);
    pv2.setLineColor(java.awt.Color.BLUE);
    pv2.setLineWidth(5);
    pv3.setLineColor(java.awt.Color.BLACK);
    pv3.setLineWidth(5);
    PlotPanel pp = new PlotPanel();
    pp.addTiledView(pv1);
    pp.addTiledView(pv2);
    pp.addTiledView(pv3);
    pp.setHLabel(hl);
    pp.setVLabel(vl);
    pp.setColorBarWidthMinimum(70);
    //pp.setHInterval(0.5);
    //pp.setVInterval(0.01);
    pp.setVLimits(cmin,cmax);
    PlotFrame pf = new PlotFrame(pp);
    pf.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
    pf.setVisible(true);
    int dpi = 720;
    if (slide) {
      pf.setFontSizeForSlide(fw,fh,_ratio);
      if (paint) {
        pf.paintToPng(dpi,(1920f*fw-1)/dpi,_paths+title+".png");
      }
    }
    else {
      if (one) {
        pf.setFontSizeForPrint(8.0,222.0);
        if (paint) {
          pf.paintToPng(dpi,3.08,_pathp+title+".png");
        }
      }
      else {
        pf.setFontSizeForPrint(8.0,469.0);
        if (paint) {
          pf.paintToPng(dpi,6.51,_pathp+title+".png");
        }
      }
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
      String hl, String vl, String cbl, float fw, float fh, 
      float clipMin, float clipMax, boolean clip, boolean interp, boolean ttl, 
      boolean paint, boolean cb, boolean color, boolean slide, boolean one) {
    int fwi = round(1920*fw/2+1);
    int fhi = round(1080*fh/2+1);
    PlotPanel pp = new PlotPanel(1,1,PlotPanel.Orientation.X1DOWN_X2RIGHT);
    PixelsView pv = pp.addPixels(s1,s2,f);
    //PixelsView pv = pp.addPixels(f);
    pv.setOrientation(PixelsView.Orientation.X1DOWN_X2RIGHT);
    pp.addTiledView(pv);
    pp.setHLabel(hl);
    pp.setVLabel(vl);
    //pp.setVInterval(20.0);
    pp.setColorBarWidthMinimum(100);
    PlotFrame pf = new PlotFrame(pp);
    pf.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
    //pf.setSize(fwi,fhi);
    pf.setFontSizeForSlide(fw,fh,_ratio);
    pf.setVisible(true);
    if (clip) pv.setClips(clipMin,clipMax);
    if (interp) pv.setInterpolation(PixelsView.Interpolation.NEAREST);
    if (color) pv.setColorModel(ColorMap.JET);
    if (ttl) pp.setTitle(title);
    if (cb) pp.addColorBar(cbl);
    int dpi = 720;
    if (slide) {
      pf.setFontSizeForSlide(fw,fh,_ratio);
      if (paint) {
        pf.paintToPng(dpi,(1920f*fw-1)/dpi,_paths+title+".png");
      }
    }
    else {
      if (one) {
        pf.setFontSizeForPrint(8.0,222.0);
        if (paint) {
          pf.paintToPng(dpi,3.08,_pathp+title+".png");
        }
      }
      else {
        pf.setFontSizeForPrint(8.0,469.0);
        if (paint) {
          pf.paintToPng(dpi,6.51,_pathp+title+".png");
        }
      }
    }
  }

  /**
   * For absolute difference images.
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
    pf.setFontSizeForSlide(fw,fh,_ratio);
    pf.setVisible(true);
    if (print) {
      int dpi = 720;
      pf.paintToPng(dpi,(1920f*fw-1)/dpi,_paths+title+".png");
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
    pf.setFontSizeForSlide(fw,fh,_ratio);
    pf.setVisible(true);
    if (print) {
      int dpi = 720;
      pf.paintToPng(dpi,(1920f*fw-1)/dpi,_paths+title+".png");
    }
  }

  public static void plot(Sampling s1, Sampling s2, Sampling s3, float[][][] f) {
    ImagePanelGroup img = new ImagePanelGroup(s1,s2,s3,f);
		SimpleFrame sf = new SimpleFrame();
		sf.addImagePanels(img);
    PlotPanelPixels3 ppp31 = new PlotPanelPixels3(
				PlotPanelPixels3.Orientation.X1DOWN_X2RIGHT,
				PlotPanelPixels3.AxesPlacement.LEFT_BOTTOM,s1,s2,s3,f);
		ppp31.setLineColor(java.awt.Color.GREEN);
		ppp31.addColorBar();
		PlotFrame pf = new PlotFrame(ppp31);
		pf.setVisible(true);
  }
}
