package novice;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.RandomFloat;
import edu.mines.jtk.util.Stopwatch;

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
      Plot.plot(sd1,sd2,"distances (bounded)",dist,fw,fh,false);
      Plot.plot(sd1,sd2,"accumulated cost (bounded)",accum,fw,fh,false);
      //Plot.plot(sd1,sd2,"accumulated cost (bounded)",accum,xf,yf,fw,fh,
      //    false); // this plot to show path
    }
    else {
      Plot.plot(sd1,sd2,"distances (unbounded)",dist,fw,fh,false);
      Plot.plot(sd1,sd2,"accumulated cost (unbounded)",accum,fw,fh,false);
      //Plot.plot(sd1,sd2,"accumulated cost (unbounded)",accum,xf,yf,fw,fh,
      //    false); // this plot to show path

    }

    System.out.println("cost[4][6]"+accum[4][6]);
    //System.out.println("cost[3][6]"+accum[300][n-1]);
    //System.out.println("cost[3][6]"+accum[350][n-1]);
  }

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
