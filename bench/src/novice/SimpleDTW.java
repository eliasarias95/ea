package novice;

import edu.mines.jtk.io.*;
import edu.mines.jtk.mosaic.*;
import edu.mines.jtk.awt.ColorMap;

import static edu.mines.jtk.util.ArrayMath.*;

import utils.Plot;

import java.io.*;
import java.nio.*;

/**
 * This software very simply shows the uses and applicability of dynamic time 
 * warping in trying to match two signals to one another.
 * @author Elias Arias, Colorado School of Mines, CWP
 * @version 06.10.2014
 */
public class SimpleDTW {

  public static void unbounded(float[] a, float[] b, float[][] dist, 
      float[][] accum) {
    int na = a.length;
    int nb = b.length;

    for (int i2=0; i2<nb; ++i2)
      for (int i1=0; i1<na; ++i1)
        dist[i2][i1] = pow(a[i1]-b[i2],2);

    accum[0][0] = dist[0][0];
    for (int i1=1; i1<na; ++i1)
      accum[0][i1] = dist[0][i1] + accum[0][i1-1];
    for (int i2=1; i2<nb; ++i2)
      accum[i2][0] = dist[i2][0] + accum[i2-1][0];
    for (int i2=1; i2<nb; ++i2) {
      for (int i1=1; i1<na; ++i1) {
        accum[i2][i1] = min(accum[i2-1][i1-1],
                                       accum[i2-1][i1  ],
                                       accum[i2  ][i1-1])
                                      + dist[i2][i1];
      }
    }
  }

  /**
   * For now, a and b must be of equal length.
   */
  public static void bounded(float[] a, float[] b, float[][] dist, 
      float[][] accum, int bound) {
    int na = a.length;

    //building the distance matrix
    int j = 0;
    for (int i2=0; i2<bound; ++i2) {
      for (int i1=0; i1<bound+j; ++i1) {
        dist[i2][i1] = pow(a[i1]-b[i2],2);
      }
      ++j;
    }

    j = 1;
    for (int i2=bound; i2<na-bound; ++i2) {
      for (int i1=j; i1<2*bound+j; ++i1) {
        dist[i2][i1] = pow(a[i1]-b[i2],2);
      }
      ++j;
    }

    j = 0;
    for (int i2=na-bound; i2<na; ++i2) {
      for (int i1=na-2*bound+j; i1<na; ++i1) {
        dist[i2][i1] = pow(a[i1]-b[i2],2);
      }
      ++j;
    }

    //building the accumulation matrix
    accum[0][0] = dist[0][0];
    for (int i1=1; i1<na; ++i1)
      accum[0][i1] = 99999f;
    for (int i2=1; i2<na; ++i2)
      accum[i2][0] = 99999f;

    j = 0;
    for (int i2=1; i2<bound; ++i2) {
      for (int i1=1; i1<bound+j; ++i1) {
        accum[i2][i1] = min(accum[i2-1][i1-1],
                                       accum[i2-1][i1  ],
                                       accum[i2  ][i1-1])
                                      + dist[i2][i1];
      }
      ++j;
    }

    j = 1;
    for (int i2=bound; i2<na-bound; ++i2) {
      for (int i1=j; i1<2*bound+j; ++i1) {
        accum[i2][i1] = min(accum[i2-1][i1-1],
                                       accum[i2-1][i1  ],
                                       accum[i2  ][i1-1])
                                      + dist[i2][i1];
      }
      ++j;
    }

    j = 0;
    for (int i2=na-bound; i2<na; ++i2) {
      for (int i1=na-2*bound+j; i1<na; ++i1) {
        accum[i2][i1] = min(accum[i2-1][i1-1],
                                       accum[i2-1][i1  ],
                                       accum[i2  ][i1-1])
                                      + dist[i2][i1];
      }
      ++j;
    }
  }

  public static int[][] getPath(float[][] dist, float[][] accum, 
      float[] xf, float[] yf) {
    int na = accum[0].length;
    int nb = accum.length;

    int[][] path = new int[na+nb][2];
    int i1 = na-1;
    int i2 = nb-1;
    int i3 = 0;
    float cost = 0;
    int[] x = new int[nb+na];
    int[] y = new int[nb+na];
    while (i2>0 && i1>0) {
      if (accum[i2-1][i1] == min(accum[i2-1][i1-1], 
                                 accum[i2-1][i1  ], 
                                 accum[i2  ][i1-1])) {
        i2 = i2-1;
    }
      else if (accum[i2][i1-1] == min(accum[i2-1][i1-1], 
                                      accum[i2-1][i1  ], 
                                      accum[i2  ][i1-1]))
      {
        i1 = i1-1;
      }
      else {
        i2 = i2-1;
        i1 = i1-1;
      }
    
      path[i3][0]=i1;
      path[i3][1]=i2;
      i3++;
      for (int i=0; i<nb+na; ++i) {
        x[i] = path[i][0];
        y[i] = path[i][1];
        xf[i] = (float)path[i][0];
        yf[i] = (float)path[i][1];
        cost += dist[x[i]][y[i]];
      }
    }
    //path[na-1][0]=0;
    //path[na-1][1]=0;
    return path;
  }
}
