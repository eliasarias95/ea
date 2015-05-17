package util;

import edu.mines.jtk.mosaic.*;
import edu.mines.jtk.awt.ColorMap;
import edu.mines.jtk.dsp.*;
import edu.mines.jtk.sgl.*;

import static edu.mines.jtk.util.ArrayMath.*;

import java.awt.Color;
import java.text.DecimalFormat;
import java.io.IOException;
import javax.swing.*;

/**
 * Plotting class to hold various plot methods.
 * @author Elias Arias, Colorado School of Mines, CWP
 * @version 26.01.2015
 */
public class Plot {

  private static int k1 = 185;
  private static int k2 = 0;
  private static int k3 = 83;
  //private static int k1 = 99;
  //private static int k2 = 90;
  //private static int k3 = 90;
  private static final double _ratio = 16.0/9.0;
  private static final String _paths = 
    "/Users/earias/Documents/research/figures/testing/";
  private static final String _pathp = 
    "/Users/earias/Documents/research/figures/thesis/";

  public static void plot(Sampling s, float[] f1, float[] f2, float[] f3, 
      String title, String hl, String vl, 
      float fw, float fh, float cmin, float cmax,
      boolean slide, boolean one, boolean paint) {
    float[] zero = new float[s.getCount()];
    fill(0.00000001f,zero);
    PointsView pv1 = new PointsView(s,f1);
    PointsView pv2 = new PointsView(s,f2);
    PointsView pv3 = new PointsView(s,f3);
    PointsView pv4 = new PointsView(s,zero);
    //pv1.setLineWidth(5);
    //pv1.setLineColor(java.awt.Color.green);
    pv1.setMarkColor(java.awt.Color.GREEN);
    pv1.setMarkStyle(PointsView.Mark.FILLED_CIRCLE);
    pv1.setLineStyle(PointsView.Line.NONE);
    pv1.setMarkSize(12);
    //pv2.setLineWidth(5);
    //pv2.setLineColor(java.awt.Color.BLUE);
    pv2.setMarkColor(java.awt.Color.BLUE);
    pv2.setMarkStyle(PointsView.Mark.FILLED_CIRCLE);
    pv2.setLineStyle(PointsView.Line.NONE);
    pv2.setMarkSize(12);
    //pv3.setLineWidth(5);
    //pv3.setLineColor(java.awt.Color.RED);
    pv3.setMarkColor(java.awt.Color.RED);
    pv3.setMarkStyle(PointsView.Mark.FILLED_CIRCLE);
    pv3.setLineStyle(PointsView.Line.NONE);
    pv3.setMarkSize(12);
    pv4.setLineColor(java.awt.Color.BLACK);
    pv4.setLineWidth(2);
    PlotPanel pp = new PlotPanel();
    pp.addTiledView(pv4);
    pp.addTiledView(pv1);
    pp.addTiledView(pv2);
    pp.addTiledView(pv3);
    pp.setHLabel(hl);
    pp.setVLabel(vl);
    pp.setColorBarWidthMinimum(70);
    //pp.setHInterval(0.5);
    pp.setVInterval(0.1);
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

  public static void plot(
      Sampling s, float[] f1, float[] f2, float[] f3, float[] f4, 
      String title, String hl, String vl, float fw, float fh, 
      float cmin, float cmax, boolean slide, boolean one, boolean paint) {
    float[] zero = new float[s.getCount()];
    fill(0.00000001f,zero);
    PointsView pv5 = new PointsView(s,zero);
    PointsView pv1 = new PointsView(s,f1);
    PointsView pv2 = new PointsView(s,f2);
    PointsView pv3 = new PointsView(s,f3);
    PointsView pv4 = new PointsView(s,f4);
    pv1.setMarkColor(java.awt.Color.GREEN);
    pv1.setMarkStyle(PointsView.Mark.FILLED_CIRCLE);
    pv1.setLineStyle(PointsView.Line.NONE);
    pv1.setMarkSize(12);
    pv2.setMarkColor(java.awt.Color.BLUE);
    pv2.setMarkStyle(PointsView.Mark.FILLED_CIRCLE);
    pv2.setLineStyle(PointsView.Line.NONE);
    pv2.setMarkSize(12);
    pv3.setMarkColor(java.awt.Color.RED);
    pv3.setMarkStyle(PointsView.Mark.FILLED_CIRCLE);
    pv3.setLineStyle(PointsView.Line.NONE);
    pv3.setMarkSize(12);
    pv4.setMarkColor(java.awt.Color.BLACK);
    pv4.setMarkStyle(PointsView.Mark.FILLED_CIRCLE);
    pv4.setLineStyle(PointsView.Line.NONE);
    pv4.setMarkSize(12);
    pv5.setLineColor(java.awt.Color.BLACK);
    pv5.setLineWidth(2);
    PlotPanel pp = new PlotPanel();
    pp.addTiledView(pv5);
    pp.addTiledView(pv4);
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
      String cbl, float fw, float fh, float cmin, float cmax, 
      boolean clip, boolean interp, boolean ttl, boolean paint, boolean color, 
      boolean slide, boolean one) {
    int fwi = round(1920*fw/2+1);
    int fhi = round(1080*fh/2+1);
    PlotPanel pp = new PlotPanel(1,1,PlotPanel.Orientation.X1DOWN_X2RIGHT);
    PixelsView pv = pp.addPixels(s1,s2,f);
    pv.setOrientation(PixelsView.Orientation.X1DOWN_X2RIGHT);
    pp.addTiledView(pv);
    pp.setHLabel("Traces");
    pp.setVLabel("Samples");
    pp.addColorBar(cbl);
    pp.setColorBarWidthMinimum(100);
    PlotFrame pf = new PlotFrame(pp);
    pf.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
    pf.setSize(fwi,fhi);
    pf.setFontSizeForSlide(fw,fh,_ratio);
    pf.setVisible(true);
    if (clip) pv.setClips(cmin,cmax);
    if (interp) pv.setInterpolation(PixelsView.Interpolation.NEAREST);
    if (color) pv.setColorModel(ColorMap.JET);
    if (ttl) pp.setTitle(title);
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
  public static void plot(Sampling s1, Sampling s2, float[][] f, float[][] g, 
      String title, String cbl, float fw, float fh, float cmin, float cmax, 
      boolean clip, boolean ttl, boolean paint, 
      boolean slide, boolean one) {
    int fwi = round(1920*fw/2+1);
    int fhi = round(1080*fh/2+1);
    PlotPanel pp = new PlotPanel(1,1,PlotPanel.Orientation.X1DOWN_X2RIGHT);
    PixelsView pv1 = pp.addPixels(s1,s2,f);
    PixelsView pv2 = pp.addPixels(s1,s2,g);
    pv1.setOrientation(PixelsView.Orientation.X1DOWN_X2RIGHT);
    //pv1.setClips(-2000,2000);
    pv2.setOrientation(PixelsView.Orientation.X1DOWN_X2RIGHT);
    pv2.setColorModel(ColorMap.setAlpha(ColorMap.JET,0.4));
    pv1.setInterpolation(PixelsView.Interpolation.NEAREST);
    pv2.setInterpolation(PixelsView.Interpolation.NEAREST);
    pp.addTiledView(pv1);
    pp.addTiledView(pv2);
    pp.setHLabel("Traces");
    pp.setVLabel("Samples");
    pp.addColorBar(cbl);
    pp.setColorBarWidthMinimum(110);
    pp.setHInterval(100);
    pp.setVInterval(100);
    PlotFrame pf = new PlotFrame(pp);
    pf.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
    //pf.setSize(fwi,fhi);
    pf.setSize(850,600);
    pf.setVisible(true);
    if (clip) pv2.setClips(cmin,cmax);
    if (ttl) pp.setTitle(title);
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
   * @param f array[n2][n1] of floats to be plotted
   * @param title the title of the image generated
   */
  public static void plot(Sampling s1, Sampling s2, 
      float[][] f1, float[][] f2, float[][] f3, float[][] f4, String title, 
      float fw, float fh, boolean slide, boolean paint) {
    int fwi = round(1920*fw/2+1);
    int fhi = round(1080*fh/2+1);
    PlotPanel pp = new PlotPanel(1,4,PlotPanel.Orientation.X1DOWN_X2RIGHT);
    PixelsView pv1 = pp.addPixels(0,0,s1,s2,f1);
    PixelsView pv2 = pp.addPixels(0,1,s1,s2,f2);
    PixelsView pv3 = pp.addPixels(0,2,s1,s2,f3);
    PixelsView pv4 = pp.addPixels(0,3,s1,s2,f4);
    pv2.setColorModel(ColorMap.JET);
    pv3.setColorModel(ColorMap.JET);
    pv4.setColorModel(ColorMap.JET);
    pv2.setClips(-1.0f,1.0f);
    pv3.setClips(-1.0f,1.0f);
    pv4.setClips(-1.0f,1.0f);
    pv1.setInterpolation(PixelsView.Interpolation.NEAREST);
    pv2.setInterpolation(PixelsView.Interpolation.NEAREST);
    pv3.setInterpolation(PixelsView.Interpolation.NEAREST);
    pv4.setInterpolation(PixelsView.Interpolation.NEAREST);
    pp.addColorBar("Slope (samples/traces)");
    pp.setColorBarWidthMinimum(100);
    pp.setVLabel(0,"Samples");
    pp.setHLabel(0,"Traces");
    pp.setHLabel(1,"Traces");
    pp.setHLabel(2,"Traces");
    pp.setHLabel(3,"Traces");
    PlotFrame pf = new PlotFrame(pp);
    pf.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE);
    pf.setSize(fwi*4,fhi);
    pf.setVisible(true);
    int dpi = 720;
    if (slide) {
      pf.setFontSizeForSlide(fw,fh,_ratio);
      if (paint) {
        pf.paintToPng(dpi,(1920f*fw-1)/dpi,_paths+title+".png");
      }
    }
    else {
      pf.setFontSizeForPrint(8.0,469.0);
      if (paint) {
        pf.paintToPng(dpi,6.51,_pathp+title+".png");
      }
    }
  }

  public static void plot(Sampling s1, Sampling s2, Sampling s3, 
      float[][][] f, float[][][] g, String title, 
      float cmin, float cmax, boolean paint) {
    int n1 = s1.getCount();
    int n2 = s2.getCount();
    int n3 = s3.getCount();
    SimpleFrame sf = new SimpleFrame(AxesOrientation.XRIGHT_YIN_ZDOWN);
    ImagePanelGroup2 ipg = new ImagePanelGroup2(s1,s2,s3,f,g);
    ipg.setClips1(-4,4);
    //ipg.setClips1(-2000,2000);
    ipg.setClips2(cmin,cmax);
    ipg.setColorModel2(ColorMap.setAlpha(ColorMap.JET,0.4));
    sf.getWorld().addChild(ipg);
    ipg.setSlices(k1,k2,k3);
    sf.setSize(985,700);   // for sch data
    ViewCanvas vc = sf.getViewCanvas();
    vc.setBackground(Color.WHITE);
    double radius = 0.5*sqrt(n1*n1+n2*n2+n3*n3);
    OrbitView ov = sf.getOrbitView();
    ov.setEyeToScreenDistance(3018.87); // for consistency with brooks
    ov.setWorldSphere(new BoundingSphere(0.5*n1,0.5*n2,0.5*n3,radius));
    ov.setAzimuthAndElevation(50.0,15.0);
    ov.setScale(1.0);
    ov.setTranslate(new Vector3(0.090,0.238,0.012));
    sf.setVisible(true);
    if (paint) {
      sf.paintToFile(_paths+title+".png");
    }
  }

  public static void plot(Sampling s1, Sampling s2, Sampling s3, 
      float[][][] f, float[][][] g, String cbl, String title, 
      float fw, float fh, float cmin, float cmax, 
      boolean paint, boolean slide, boolean one) {
    int fwi = round(1920*fw/2+1);
    int fhi = round(1080*fh/2+1);
    Color background = Color.WHITE;
    int n1 = s1.getCount();
    int n2 = s2.getCount();
    int n3 = s3.getCount();
    float[][][] ft = new float[n2][n3][n1];
    float[][][] gt = new float[n2][n3][n1];
    Util.transpose23(f,ft);
    Util.transpose23(g,gt);
    PlotPanelPixels3 pp = new PlotPanelPixels3(
      PlotPanelPixels3.Orientation.X1DOWN_X2RIGHT,
      PlotPanelPixels3.AxesPlacement.LEFT_BOTTOM,
      s1,s3,s2,ft);
    pp.setClips(-4,4);
    pp.setSlices(k1,k3,k2);
    pp.setLabel1("Samples");
    pp.setLabel3("Crossline");
    pp.setLabel2("Inline");
    //pp.getMosaic().setHeightElastic(1, 85);
    pp.setLineColor(Color.YELLOW);
    ColorBar cb = pp.addColorBar(cbl);
    cb.setInterval(1.0);
    pp.setInterval1(100);
    pp.setInterval2(100);
    pp.setInterval3(100);
    PixelsView pv12 = new PixelsView(s1,s3,Util.slice12(k2,gt));
    pv12.setOrientation(PixelsView.Orientation.X1DOWN_X2RIGHT);
    pv12.setInterpolation(PixelsView.Interpolation.NEAREST);
    PixelsView pv13 = new PixelsView(s1,s2,Util.slice13(k3,gt));
    pv13.setOrientation(PixelsView.Orientation.X1DOWN_X2RIGHT);
    pv13.setInterpolation(PixelsView.Interpolation.NEAREST);
    PixelsView pv23 = new PixelsView(s3,s2,Util.slice23(k1,gt));
    pv23.setOrientation(PixelsView.Orientation.X1RIGHT_X2UP);
    pv23.setInterpolation(PixelsView.Interpolation.NEAREST);
    pv12.setColorModel(ColorMap.setAlpha(ColorMap.JET,0.4));
    pv13.setColorModel(ColorMap.setAlpha(ColorMap.JET,0.4));
    pv23.setColorModel(ColorMap.setAlpha(ColorMap.JET,0.4));
    float ccmin = -0.0001f;
    float ccmax = 0.0001f;
    pv12.setClips(cmin,cmax);
    pv13.setClips(cmin,cmax);
    pv23.setClips(cmin,cmax);
    pp.getPixelsView12().getTile().addTiledView(pv12);
    pp.getPixelsView13().getTile().addTiledView(pv13);
    pp.getPixelsView23().getTile().addTiledView(pv23);
    PlotFrame pf = new PlotFrame(pp);
    pf.setBackground(background);
    pp.setColorBarWidthMinimum(110);
    //pf.setFontSize(18) //for print
    //pf.setFontSize(30); //for slices
    //pf.setFontSizeForPrint(1.0,0.8)
    pf.setSize(1150,800);
    pp.getMosaic().setWidthElastic(1,90);
    pf.setVisible(true);
    int dpi = 720;
    if (slide) {
      pf.setFontSizeForSlide(fw,fh,_ratio);
      if (paint) {
        pf.paintToPng(dpi,(1920f*fh-1)/dpi,_paths+title+".png");
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
}
