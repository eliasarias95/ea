package util;

import edu.mines.jtk.mosaic.*;
import edu.mines.jtk.awt.ColorMap;
import edu.mines.jtk.dsp.*;
import edu.mines.jtk.sgl.*;
import edu.mines.jtk.util.SimpleFloat3;

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
      String cbl, float fw, float fh, float cmin, float cmax, 
      boolean clip, boolean interp, boolean ttl, boolean paint, boolean color, 
      boolean slide, boolean one) {
    int fwi = round(1920*fw/2+1);
    int fhi = round(1080*fh/2+1);
    PlotPanel pp = new PlotPanel(1,1,PlotPanel.Orientation.X1DOWN_X2RIGHT);
    PixelsView pv = pp.addPixels(s1,s2,f);
    //PixelsView pv = pp.addPixels(f);
    pv.setOrientation(PixelsView.Orientation.X1DOWN_X2RIGHT);
    pp.addTiledView(pv);
    pp.setHLabel("Traces");
    pp.setVLabel("Samples");
    //pp.setVInterval(20.0);
    pp.addColorBar(cbl);
    pp.setColorBarWidthMinimum(100);
    PlotFrame pf = new PlotFrame(pp);
    pf.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
    //pf.setSize(1200,900);
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
    ipg.setClips2(cmin,cmax);
    ipg.setColorModel2(ColorMap.setAlpha(ColorMap.JET,0.4));
    sf.getWorld().addChild(ipg);
    int k1 = 220;
    int k2 = 317;
    int k3 = 41;
    ipg.setSlices(k1,k2,k3);
    sf.setSize(985,700);   // for sch data
    //sf.setSize(837,700); // for fake data
    ViewCanvas vc = sf.getViewCanvas();
    vc.setBackground(Color.WHITE);
    double radius = 0.5*sqrt(n1*n1+n2*n2+n3*n3);
    OrbitView ov = sf.getOrbitView();
    ov.setEyeToScreenDistance(3018.87); // for consistency with brooks
    ov.setWorldSphere(new BoundingSphere(0.5*n1,0.5*n2,0.5*n3,radius));
    ov.setAzimuthAndElevation(35.0,20.0);
    ov.setScale(1.2);
    ov.setTranslate(new Vector3(0.090,0.238,0.012));
    sf.setVisible(true);
    if (paint) {
      sf.paintToFile(_paths+title+".png");
    }
  }

  public static void plot(Sampling s1, Sampling s2, Sampling s3, 
      float[][][] f, float[][][] g, String cbl, String title, 
      float cmin, float cmax, boolean paint) {
    Color background = Color.WHITE;
    int n1 = s1.getCount();
    int n2 = s2.getCount();
    int n3 = s3.getCount();
    PlotPanelPixels3 pp = new PlotPanelPixels3(
      PlotPanelPixels3.Orientation.X1DOWN_X2RIGHT,
      PlotPanelPixels3.AxesPlacement.LEFT_BOTTOM,
      s1,s2,s3,f);
    int k1 = 220;
    int k2 = 317;
    int k3 = 41;
    pp.setSlices(k1,k2,k3);
    pp.setLabel1("Samples");
    pp.setLabel2("Inline (traces)");
    pp.setLabel3("Crossline (traces)");
    //pp.mosaic.setHeightElastic(0,100);
    pp.getMosaic().setHeightElastic(1, 85);
    pp.setLineColor(Color.YELLOW);
    ColorBar cb = pp.addColorBar(cbl);
    pp.setInterval1(50);
    pp.setInterval2(50);
    pp.setInterval3(50);
    PixelsView pv12 = new PixelsView(s1,s2,slice12(k3,g));
    pv12.setOrientation(PixelsView.Orientation.X1DOWN_X2RIGHT);
    pv12.setInterpolation(PixelsView.Interpolation.NEAREST);
    PixelsView pv13 = new PixelsView(s1,s3,slice13(k2,g));
    pv13.setOrientation(PixelsView.Orientation.X1DOWN_X2RIGHT);
    pv13.setInterpolation(PixelsView.Interpolation.NEAREST);
    PixelsView pv23 = new PixelsView(s2,s3,slice23(k1,g));
    pv23.setOrientation(PixelsView.Orientation.X1RIGHT_X2UP);
    pv23.setInterpolation(PixelsView.Interpolation.NEAREST);
    pv12.setColorModel(ColorMap.setAlpha(ColorMap.JET,0.4));
    pv13.setColorModel(ColorMap.setAlpha(ColorMap.JET,0.4));
    pv23.setColorModel(ColorMap.setAlpha(ColorMap.JET,0.4));
    pv12.setClips(cmin,cmax);
    pv13.setClips(cmin,cmax);
    pv23.setClips(cmin,cmax);
    pp.getPixelsView12().getTile().addTiledView(pv12);
    pp.getPixelsView13().getTile().addTiledView(pv13);
    pp.getPixelsView23().getTile().addTiledView(pv23);
    PlotFrame pf = new PlotFrame(pp);
    pf.setBackground(background);
    pp.setColorBarWidthMinimum(120);
    //pf.setFontSize(18) //for print
    pf.setFontSize(27); //for slices
    //pf.setFontSizeForPrint(1.0,0.8)
    pf.setSize(1150,800);
    pf.setVisible(true);
    if (paint) {
      pf.paintToPng(360,7.0,_paths+title+".png");
    }
  }
  
  private static float[][] slice12(int k3, float[][][] f) {
    int n1 = f[0][0].length;
    int n2 = f[0].length;
    float[][] s = new float[n2][n1];
    new SimpleFloat3(f).get12(n1,n2,0,0,k3,s);
    return s;
  }

  private static float[][] slice13(int k2, float[][][] f) {
    int n1 = f[0][0].length;
    int n3 = f.length;
    float[][] s = new float[n3][n1];
    new SimpleFloat3(f).get13(n1,n3,0,k2,0,s);
    return s;
  }

  private static float[][] slice23(int k1, float[][][] f) {
    int n2 = f[0].length;
    int n3 = f.length;
    float[][] s = new float[n3][n2];
    new SimpleFloat3(f).get23(n2,n3,k1,0,0,s);
    return s;
  }
}
