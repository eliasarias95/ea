import pwd.*;

import edu.mines.jtk.io.*;
import edu.mines.jtk.dsp.*;
import edu.mines.jtk.mosaic.*;
import edu.mines.jtk.awt.ColorMap;

import static edu.mines.jtk.util.ArrayMath.*;

import dnp.LocalSlopeFinder;

import java.io.*;
import java.nio.*;

import javax.swing.*;

/**
 * 3D plane-wave destruction.
 * @author Elias Arias, Colorado School of Mines, CWP
 * @version 18.12.2013
 */

public class PWD {
  private static void go() {
    int n1 = 301;
    int n2 = 920;
    int n  = n1*n2;

    float d1 = 0.004f;
    float d2 = .02667f;
    float f1 = 1.6f;
    float f2 = 0;

    Sampling s1 = new Sampling(n1,d1,f1);
    Sampling s2 = new Sampling(n2,d2,f2);

    float[][] x1 = readImage(n1,n2,"src/pwd/gom.dat");
    mul(x1,.001f,x1);
    writeBinary(x1,"src/pwd/"+DATA_NAME);
    float[][] x2 = readImage(n1,n2,"src/pwd/"+DATA_NAME);

    LocalSlopeFinder lsf = new LocalSlopeFinder(8.0f,4.0f,2);
    float[][] lsf_output2 = new float[n2][n1];
    float[][] filter_output2 = new float[n2][n1];
    lsf.findSlopes(x2,lsf_output2);

    float[] x = new float[n];
    float[] lsf_output31 = new float[n];
    float[] filter_output31 = new float[n];
    float[][] filter_output32 = new float[n2][n1];

    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        lsf_output31[i1+i2*n1] = lsf_output2[i2][i1];
        x[i1+i2*n1] = x2[i2][i1];
      }
    }

    Filter2 fil2 = new Filter2(n1,n2,2,1,lsf_output2);
    fil2.destructor(false,x2,filter_output2);
    //mul(filter_output2,0.1f,filter_output2);


    Filter3 fil3 = new Filter3(n1,n2,0,2,1,lsf_output31);
    fil3.inlineDestructor(true,false,x,filter_output31);
    mul(filter_output31,0.1f,filter_output31);

    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        filter_output32[i2][i1] = filter_output31[i1+i2*n1];
      }
    }
    System.out.println("values= "+filter_output32[2][1]);

    plot(s1,s2,sexp(x2),"GOM Near Offset Data");
    plot(s1,s2,lsf_output2,"Slopes");

    plot(s1,s2,(filter_output2),"Plane-wave Destruction 2D Filter");
    plot(s1,s2,(sub(x2,filter_output2)),"Cleaned Data");
    //plot(s1,s2,filter_output32,"Plane-wave Destruction 3D Filter");
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
  private static float[][] readImage(int n1, int n2, String fileName) {
    ByteOrder byteOrder = ByteOrder.LITTLE_ENDIAN;
    try {
      ArrayInputStream ais = new ArrayInputStream(fileName,byteOrder);
      float[][] x = new float[n2][n1];
      ais.readFloats(x);
      ais.close();
      return x;
    } catch (IOException e) {
      throw new RuntimeException(e);
    }
  }

  /**
   * Writes seismic data to binary file.
   * @param x1 array[n2][n1] of data to write to the binary file
   * @param fileName name of output binary file
   */
  private static void writeBinary(float[][] x, String fileName) {
    ByteOrder byteOrder = ByteOrder.LITTLE_ENDIAN;
    try {
      ArrayOutputStream aos = new ArrayOutputStream(fileName,byteOrder);
      aos.writeFloats(x);
      aos.close();
    } catch (IOException e) {
      throw new RuntimeException(e);
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
      float[][] f, String title) {
    PlotPanel pp = new PlotPanel(1,1,PlotPanel.Orientation.X1DOWN_X2RIGHT);
    PixelsView pv = pp.addPixels(s1,s2,f);
    //pv.setColorModel(ColorMap.JET);
    pp.addColorBar();
    pp.setColorBarWidthMinimum(70);
    pp.setTitle(title);
    pp.setHLabel("Distance (km)");
    pp.setVLabel("Time (s)");
    pp.getMosaic().getTileAxisTop(0).setFormat("%1.6G");
    PlotFrame pf = new PlotFrame(pp);
    pf.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE);
    pf.setFontSize(20);
    pf.setSize(800,600);
    pf.setVisible(true);
    //pf.paintToPng(720,6.51,"C:/pro/ea/bench/src/pwd/"+title+".png");
  }

  public static final String DATA_NAME = "gom_scaled.dat";
  public static void main(String[] args) {
    SwingUtilities.invokeLater(new Runnable() {
      public void run() {
        go();
      }
    });
  }
}
