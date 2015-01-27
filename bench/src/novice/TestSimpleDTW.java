package novice;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.mosaic.*;
import edu.mines.jtk.awt.ColorMap;
import edu.mines.jtk.util.Stopwatch;
import edu.mines.jtk.util.RandomFloat;

import static edu.mines.jtk.util.ArrayMath.*;

import util.*;

import java.util.Random;
import javax.swing.*;

public class TestSimpleDTW {

  private static void daveDW() {
    int n1 = 101;
    int shift = 10;
    float[] a = Util.readImage(n1,
        "/Users/earias/Home/git/ea/bench/src/novice/data/chris.dat");
    float[] b = new float[n1];
    for (int i=0; i<n1-shift; ++i)
      b[i] = a[i+shift];

    int count = 0;
    for (int i=n1-shift; i<n1; ++i) {
      b[i] = a[count];
      ++count;
    }

    float[] c = new float[n1];
    c = Util.addNoise(0.1,a);

    int shiftmax = 50;
    DynamicWarping dw = new DynamicWarping(-shiftmax,shiftmax);
    dw.setShiftSmoothing(1.0,1.0);
    //dw.setStrainMax(0.9,0.9);

    float[] dw_slope = new float[n1];
    dw_slope = dw.findShifts(a,c);

    float fw = 0.75f; //fraction width for slide
    float fh = 0.9f; //fraction height for slide
    //Plot.plot(a,c,"Curves","Index","Value",fw,fh,false);
    //Plot.plot(dw_slope,"Slope","Index","Value",fw,fh,false);
  }

  private static void goDTW(boolean bounded) {
    int n1 = 101;
    int shift = 10;
    float[] a = Util.readImage(n1,
        "/Users/earias/Home/git/ea/bench/src/novice/data/chris.dat");
    float[] b = new float[n1];
    for (int i=0; i<n1-shift; ++i)
      b[i] = a[i+shift];

    int count = 0;
    for (int i=n1-shift; i<n1; ++i) {
      b[i] = a[count];
      ++count;
    }

    float[] c = new float[n1];
    c = Util.addNoise(0.1,a);

    int na = a.length;
    int nb = b.length; 
    float[][] dist = new float[nb][na];
    float[][] accum = new float[nb][na];

    int bound = (int)(na*0.05f);

    if (bounded) {
      //SimpleDTW.bounded(a,b,dist,accum,bound);
      SimpleDTW.bounded(a,c,dist,accum,bound);
    }
    else {
      //SimpleDTW.unbounded(a,b,dist,accum);
      SimpleDTW.unbounded(a,c,dist,accum);
    }

    /**
     * For displaying the path of least error on accumulation matrix
    float[] xf = new float[nb+na];
    float[] yf = new float[nb+na];
    int[][] path = SimpleDTW.getPath(dist,accum,xf,yf);
    
    int[] map_x = new int[nb+na];
    int[] map_y = new int[nb+na];
    for (int i=0; i<nb+na; ++i) {
      map_x[i] = path[i][0];
      map_y[i] = path[i][0];
    }
    */

    int fd = 0;
    int dd = 1;
    int nd1 = na; 
    int nd2 = nb;
    Sampling sd1 = new Sampling(nd1,dd,fd);
    Sampling sd2 = new Sampling(nd2,dd,fd);

    float fw = 0.70f; //fraction width for slide
    float fh = 0.9f; //fraction height for slide
    //Plot.plot(a,c,"Curves","Index","Value",fw,fh,false);

    if (bounded) {
      plot(sd1,sd2,"distances (bounded)",dist,fw,fh,false);
      plot(sd1,sd2,"accumulated cost (bounded)",accum,fw,fh,false);
      //plot(sd1,sd2,"accumulated cost (bounded)",accum,xf,yf,fw,fh,
      //    false); // this plot to show path
    }
    else {
      plot(sd1,sd2,"distances (unbounded)",dist,fw,fh,false);
      plot(sd1,sd2,"accumulated cost (unbounded)",accum,fw,fh,false);
      //plot(sd1,sd2,"accumulated cost (unbounded)",accum,xf,yf,fw,fh,
      //    false); // this plot to show path

    }

    System.out.println("cost[4][6]"+accum[4][6]);
    //System.out.println("cost[3][6]"+accum[300][n-1]);
    //System.out.println("cost[3][6]"+accum[350][n-1]);
  }

  /**
   * Plots a 2D array of floats with specified title.
   * @param s1 the sampling in the 1st-dimension
   * @param s2 the sampling in the 2nd-dimension
   * @param title the title of the image generated.
   * @param f the 2D array of floats to be plotted.
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
      pf.paintToPng(dpi,(1920f*fw-1)/dpi,path+title+".png");
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
      pf.paintToPng(dpi,(1920f*fw-1)/dpi,path+title+".png");
    }
  }

  private static final double ratio = 16.0/9.0;  
  private static final String path = 
    "/Users/earias/Documents/junk_figures/";  
  public static void main(String[] args) {
    SwingUtilities.invokeLater(new Runnable() {
      public void run() {
        Stopwatch sw = new Stopwatch();
        sw.start();
        goDTW(false);
        sw.stop();
        double ubtime = sw.time();
        System.out.println("Unbounded time= "+ubtime);
        sw.restart();
        goDTW(true);
        sw.stop();
        double btime = sw.time();
        System.out.println("Bounded time= "+btime);
        System.out.println("Bounded is "+ubtime/btime+" times faster.");
        daveDW();
      }
    });
  }
}
