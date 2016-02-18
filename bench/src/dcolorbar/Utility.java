package dcolorbar;

import java.io.*;
import java.nio.*;
import java.awt.*;
import java.awt.image.*;
import java.util.Random;
import javax.imageio.*;
import javax.swing.*;

import edu.mines.jtk.io.*;
import edu.mines.jtk.awt.*;
import edu.mines.jtk.dsp.*;
import edu.mines.jtk.sgl.*;
import edu.mines.jtk.util.*;
import edu.mines.jtk.interp.*;
import edu.mines.jtk.mosaic.*;
import edu.mines.jtk.ogl.Gl.*;

import static edu.mines.jtk.util.ArrayMath.*;
/*
 *  
 * @author Elias Arias, Colorado School of Mines, CWP
 * @version 19.1.2015
 */

public class Utility {

  private static String path = 
    "/Users/earias/Home/git/ea/bench/src/dcolorbar/data/";

  /**
   * Converts a 1D array of doubles into an array of floats.
   * @param d array[n1] of doubles to be converted.
   * @return array[n1] of floats.
   */
  public static float[] f(double[] x) {
    int n1 = x.length;
    float[] f = new float[n1];
    for (int i1=0; i1<n1; ++i1)
      f[i1] = (float)x[i1];
    return f;
  }

  /**
   * Reads a binary file.
   * @param n1 the length of floats in the 1st-dimension
   * @param fileName the name of the file to be read
   * @return array[n1] of floats read from file
   */
  public static float[] read(int n1, String fileName) {
    ByteOrder byteOrder = ByteOrder.BIG_ENDIAN;
    try {
      ArrayInputStream ais = new ArrayInputStream(fileName,byteOrder);
      float[] x = new float[n1];
      ais.readFloats(x);
      ais.close();
      return x;
    } catch (IOException e) {
      throw new RuntimeException(e);
    }
  }

  /**
   * Reads a binary file.
   * @param n1 the length of floats in the 1st-dimension
   * @param fileName the name of the file to be read
   * @return array[n1] of floats read from file
   */
  public static float[] readL(int n1, String fileName) {
    ByteOrder byteOrder = ByteOrder.LITTLE_ENDIAN;
    try {
      ArrayInputStream ais = new ArrayInputStream(fileName,byteOrder);
      float[] x = new float[n1];
      ais.readFloats(x);
      ais.close();
      return x;
    } catch (IOException e) {
      throw new RuntimeException(e);
    }
  }


  /**
   * Reads a binary file.
   * @param n1 the length of floats in the 1st-dimension
   * @param n2 the length of floats in the 2nd-dimension
   * @param fileName the name of the file to be read
   * @return array[n2][n1] of floats read from file
   */
  public static float[][] read(int n1, int n2, String fileName) {
    ByteOrder byteOrder = ByteOrder.BIG_ENDIAN;
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
   * Reads a binary file.
   * @param n1 the length of floats in the 1st-dimension
   * @param n2 the length of floats in the 2nd-dimension
   * @param fileName the name of the file to be read
   * @return array[n2][n1] of floats read from file
   */
  public static float[][] readL(int n1, int n2, String fileName) {
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
   * Reads a binary file.
   * @param n1 the length of floats in the 1st-dimension
   * @param n2 the length of floats in the 2nd-dimension
   * @param n3 the length of floats in the 3rd-dimension
   * @param fileName the name of the file to be read
   * @return array[n3][n2][n1] of floats read from file
   */
  public static float[][][] read(int n1, int n2, int n3, String fileName) {
    ByteOrder byteOrder = ByteOrder.BIG_ENDIAN;
    try {
      ArrayInputStream ais = new ArrayInputStream(fileName,byteOrder);
      float[][][] x = new float[n3][n2][n1];
      ais.readFloats(x);
      ais.close();
      return x;
    } catch (IOException e) {
      throw new RuntimeException(e);
    }
  }

  /**
   * Reads a binary file.
   * @param n1 the length of floats in the 1st-dimension
   * @param n2 the length of floats in the 2nd-dimension
   * @param n3 the length of floats in the 3rd-dimension
   * @param fileName the name of the file to be read
   * @return array[n3][n2][n1] of floats read from file
   */
  public static float[][][] readL(int n1, int n2, int n3, String fileName) {
    ByteOrder byteOrder = ByteOrder.LITTLE_ENDIAN;
    try {
      ArrayInputStream ais = new ArrayInputStream(fileName,byteOrder);
      float[][][] x = new float[n3][n2][n1];
      ais.readFloats(x);
      ais.close();
      return x;
    } catch (IOException e) {
      throw new RuntimeException(e);
    }
  }

  /**
   * Writes seismic data to binary file.
   * @param x array[n1] of data to write to the binary file
   * @param fileName name of output binary file
   */
  public static void write(float[] x, String fileName) {
    ByteOrder byteOrder = ByteOrder.BIG_ENDIAN;
    try {
      ArrayOutputStream aos = new ArrayOutputStream(fileName,byteOrder);
      aos.writeFloats(x);
      aos.close();
    } catch (IOException e) {
      throw new RuntimeException(e);
    }
  }

  /**
   * Writes seismic data to binary file.
   * @param x array[n2][n1] of data to write to the binary file
   * @param fileName name of output binary file
   */
  public static void writeL(float[] x, String fileName) {
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
   * Writes seismic data to binary file.
   * @param x array[n2][n1] of data to write to the binary file
   * @param fileName name of output binary file
   */
  public static void write(float[][] x, String fileName) {
    ByteOrder byteOrder = ByteOrder.BIG_ENDIAN;
    try {
      ArrayOutputStream aos = new ArrayOutputStream(fileName,byteOrder);
      aos.writeFloats(x);
      aos.close();
    } catch (IOException e) {
      throw new RuntimeException(e);
    }
  }

  /**
   * Writes seismic data to binary file.
   * @param x array[n2][n1] of data to write to the binary file
   * @param fileName name of output binary file
   */
  public static void writeL(float[][] x, String fileName) {
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
   * Writes seismic data to binary file.
   * Java default byteorder is BIG_ENDIAN.
   * @param x array[n3][n2][n1] of data to write to the binary file
   * @param fileName name of output binary file
   */
  public static void write(float[][][] x, String fileName) {
    ByteOrder byteOrder = ByteOrder.BIG_ENDIAN;
    try {
      ArrayOutputStream aos = new ArrayOutputStream(fileName,byteOrder);
      aos.writeFloats(x);
      aos.close();
    } catch (IOException e) {
      throw new RuntimeException(e);
    }
  }

  /**
   * Writes seismic data to binary file.
   * Java default byteorder is BIG_ENDIAN.
   * @param x array[n3][n2][n1] of data to write to the binary file
   * @param fileName name of output binary file
   */
  public static void writeL(float[][][] x, String fileName) {
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
   * Makes a 1D array into a 2D array assuming the 1D array is column major.
   */
  public static float[][] make2D(int n1, int n2, float[] x) {
    Check.argument(n1*n2 == x.length, "n1*n2 != x.length dimension error");
    float[][] y = new float[n2][n1];
    for (int i2=0; i2<n2; ++i2)
      for (int i1=0; i1<n1; ++i1)
        y[i2][i1] = x[i1+i2*n1];
    return y;
  }

  /**
	 * Returns RGB (red, green, blue) components of an image read from a file.
	 * @param fileName the name of the file containing a JPEG image.
	 * @return array of RGB components; the components (red,green,blue) are 
	 * stored in {rgb[0],rgb[1],rgb[2]}.
	 */
	public static float[][][] rgbFromFile(String fileName) {
		BufferedImage bi = null;
		try {
			File file = new File(fileName);
			bi = ImageIO.read(file);
		} catch (IOException ioe) {
			throw new RuntimeException(ioe);
		}
		int w = bi.getWidth();
		int h = bi.getHeight();
		int[] pixels = bi.getRGB(0,0,w,h,null,0,w);
		float[][][] rgb = new float[3][h][w];
		for (int i=0,k=0; i<h; ++i) {
			for (int j=0; j<w; ++j,++k) {
				int p = pixels[k];
				int r = (p>>16)&0xff;
				int g = (p>> 8)&0xff;
				int b = (p>> 0)&0xff;
				rgb[0][i][j] = (float)r/256.0f;
				rgb[1][i][j] = (float)g/256.0f;
				rgb[2][i][j] = (float)b/256.0f;
			}
		}
		return rgb;
	}

  public static int minIndex(float[] array) {
    int narray = array.length;
    float min = FLT_MAX;
    int index = 0;
    for (int i=0; i<narray; ++i) {
      if (array[i] < min) {
        min = array[i];
        index = i;
      }
    }
    return index;
  }

  public static float[][] imageToData(
      float[] optcb, float[][] rI2, float[][] gI2, float[][] bI2) {
    int nd1 = rI2[0].length;
    int nd2 = rI2.length;
    float[] dist = new float[256];
    float[][] data = new float[nd1][nd2];
    float rdist = 0.0f;
    float gdist = 0.0f;
    float bdist = 0.0f;
    Sampling cbsamp = new Sampling(256,2.0/256.0,-1.0);
    for (int i2=0; i2<nd2; ++i2) {
      for (int i1=0; i1<nd1; ++i1) {
        int j = 0;
        for (int i=0; i<256; ++i) {
          rdist = optcb[j  ] - rI2[i2][i1];
          gdist = optcb[j+1] - gI2[i2][i1];
          bdist = optcb[j+2] - bI2[i2][i1];
          dist[i] = rdist*rdist+gdist*gdist+bdist*bdist;
          j += 3;
        }
        int k = minIndex(dist);
        data[i1][i2] = (float)cbsamp.getValue(k);
      }
    }
    return data;
  }

  public static PointGroup[] buildAxis() {
    int nPoints = 1000;
    int threeNPoints = 3*nPoints;
    float[] x = zerofloat(threeNPoints);
    float[] x1 = zerofloat(threeNPoints);
    float[] x2 = zerofloat(threeNPoints);
    float[] x3 = zerofloat(threeNPoints);
    float[] x4 = zerofloat(threeNPoints);
    float[] y = zerofloat(threeNPoints);
    float[] y1 = zerofloat(threeNPoints);
    float[] y2 = zerofloat(threeNPoints);
    float[] y3 = zerofloat(threeNPoints);
    float[] y4 = zerofloat(threeNPoints);
    float[] z = zerofloat(threeNPoints);
    float[] z1 = zerofloat(threeNPoints);
    float[] z2 = zerofloat(threeNPoints);
    float[] z3 = zerofloat(threeNPoints);
    float[] z4 = zerofloat(threeNPoints);
    float[] red = zerofloat(threeNPoints);
    float[] green = zerofloat(threeNPoints);
    float[] blue = zerofloat(threeNPoints);
    float[] white = fillfloat(1.0f,threeNPoints);
    float[] zeroToOneRGBSmallAxis = rampfloat(0.0f,0.0002f,nPoints);
    float[] zeroToOneWhite = rampfloat(0.0f,0.001f,nPoints);
    float valueRGB = 0.0f;
    float valueWhite = 0.0f;
    int j = 0;
    for (int i=0; i<nPoints; ++i) {
      valueRGB = zeroToOneRGBSmallAxis[i];
      valueWhite = zeroToOneWhite[i];
      //x
      x[j] = valueRGB;
      x[j+1] = 0.0f;
      x[j+2] = 0.0f; 
      //x1
      x1[j] = valueWhite;
      x1[j+1] = 0.0f;
      x1[j+2] = 0.0f; 
      //x2
      x2[j] = valueWhite;
      x2[j+1] = 1.0f;
      x2[j+2] = 0.0f; 
      //x3
      x3[j] = valueWhite;
      x3[j+1] = 0.0f;
      x3[j+2] = 1.0f; 
      //x4
      x4[j] = valueWhite;
      x4[j+1] = 1.0f;
      x4[j+2] = 1.0f; 
      //y
      y[j] = 0.0f;
      y[j+1] = valueRGB;
      y[j+2] = 0.0f; 
      //y1
      y1[j] = 0.0f;
      y1[j+1] = valueWhite;
      y1[j+2] = 0.0f; 
      //y2
      y2[j] = 1.0f;
      y2[j+1] = valueWhite;
      y2[j+2] = 0.0f; 
      //y3
      y3[j] = 0.0f;
      y3[j+1] = valueWhite;
      y3[j+2] = 1.0f; 
      //y4
      y4[j] = 1.0f;
      y4[j+1] = valueWhite;
      y4[j+2] = 1.0f; 
      //z
      z[j] = 0.0f; 
      z[j+1] = 0.0f;
      z[j+2] = valueRGB;
      //z1
      z1[j] = 0.0f; 
      z1[j+1] = 0.0f;
      z1[j+2] = valueWhite;
      //z2
      z2[j] = 1.0f; 
      z2[j+1] = 0.0f;
      z2[j+2] = valueWhite;
      //z3
      z3[j] = 0.0f; 
      z3[j+1] = 1.0f;
      z3[j+2] = valueWhite;
      //z4
      z4[j] = 1.0f; 
      z4[j+1] = 1.0f;
      z4[j+2] = valueWhite;
      red[j] = 1.0f;
      red[j+1] = 0.0f;
      red[j+2] = 0.0f; 
      green[j] = 0.0f;
      green[j+1] = 1.0f;
      green[j+2] = 0.0f; 
      blue[j] = 0.0f; 
      blue[j+1] = 0.0f;
      blue[j+2] = 1.0f;
      j += 3;
    }
    PointGroup pgx = new PointGroup(0.011f,x,red);
    PointGroup pgy = new PointGroup(0.011f,y,green);
    PointGroup pgz = new PointGroup(0.011f,z,blue);
    PointGroup pgx1 = new PointGroup(0.01f,x1,white);
    PointGroup pgx2 = new PointGroup(0.01f,x2,white);
    PointGroup pgx3 = new PointGroup(0.01f,x3,white);
    PointGroup pgx4 = new PointGroup(0.01f,x4,white);
    PointGroup pgy1 = new PointGroup(0.01f,y1,white);
    PointGroup pgy2 = new PointGroup(0.01f,y2,white);
    PointGroup pgy3 = new PointGroup(0.01f,y3,white);
    PointGroup pgy4 = new PointGroup(0.01f,y4,white);
    PointGroup pgz1 = new PointGroup(0.01f,z1,white);
    PointGroup pgz2 = new PointGroup(0.01f,z2,white);
    PointGroup pgz3 = new PointGroup(0.01f,z3,white);
    PointGroup pgz4 = new PointGroup(0.01f,z4,white);
    return new PointGroup[]{pgx,pgx1,pgx2,pgx3,pgx4,
                            pgy,pgy1,pgy2,pgy3,pgy4,
                            pgz,pgz1,pgz2,pgz3,pgz4};
  }


  public static void plot(float[] rgb1, float[] rgb2) {
    SimpleFrame sf = new SimpleFrame(AxesOrientation.XRIGHT_YIN_ZDOWN);
    PointGroup pg1 = new PointGroup(0.03f,rgb1,rgb1);
    PointGroup pg2 = new PointGroup(0.01f,rgb2,rgb2);
    PointGroup[] pgxyz = buildAxis();

    sf.getWorld().addChild(pg1);
    sf.getWorld().addChild(pg2);
    for (int i=0; i<15; ++i) {
      sf.getWorld().addChild(pgxyz[i]);
    }
    ViewCanvas vc = sf.getViewCanvas();
    vc.setBackground(Color.BLACK);
    OrbitView ov = sf.getOrbitView();
    //ov.setEyeToScreenDistance(3018.87); // for consistency with brooks
    //ov.setWorldSphere(new BoundingSphere(0.5*n1,0.5*n1,0.5*n1,radius));
    ov.setAzimuthAndElevation(220.0,45.0);
    ov.setScale(1.0);
    sf.setVisible(true);
    //sf.paintToFile(path+title+".png");
  }

  public static void plot(Sampling s1, Sampling s2, float[][] f) {
    PlotPanel pp = new PlotPanel(1,1,PlotPanel.Orientation.X1DOWN_X2RIGHT);
    PixelsView pv = pp.addPixels(s1,s2,f);
    //pv.setClips(-1.0f,1.0f);

    java.awt.image.IndexColorModel icm = ColorMap.GRAY;
    String fn = "gray";
    //java.awt.image.IndexColorModel icm = ColorMap.RED_WHITE_BLUE;
    //String fn = "red_white_blue";
    //java.awt.image.IndexColorModel icm = ColorMap.GRAY_YELLOW_RED;
    //String fn = "gray_yellow_red";
    //java.awt.image.IndexColorModel icm = ColorMap.PRISM;
    //String fn = "prism";
    ColorMap cm = new ColorMap(icm);
    pv.setColorModel(ColorMap.setAlpha(icm,1.0));
    int ns = 256;
    Sampling s = new Sampling(ns,1.0/ns,0.0);
    float[] vals = f(s.getValues());
    float[] rgb = cm.getRgbFloats(vals);
    //writeL(rgb,path+fn+".bin");

    //pv.setInterpolation(PixelsView.Interpolation.NEAREST);
    pp.addTiledView(pv);
    pp.setHLabel("Traces");
    pp.setVLabel("Samples");
    pp.addColorBar("amplitude");
    pp.setColorBarWidthMinimum(120);
    PlotFrame pf = new PlotFrame(pp);
    pf.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
    pf.setSize(850,600);
    pf.setFontSizeForSlide(0.8f,0.8f,16.0/9.0);
    //pf.paintToPng(720,(1920f-1)/720,path+fn+".png");
    pf.setVisible(true);
  }

  /**
   * Makes a 2D array into a 1D array such that the 1D array is column major.
   */
  public static float[] make1D(float[][] x) {
    int n2 = x.length;
    int n1 = x[0].length;
    float[] y = new float[n1*n2];
    for (int i2=0; i2<n2; ++i2)
      for (int i1=0; i1<n1; ++i1)
        y[i1+i2*n1] = x[i2][i1];
    return y;
  }

  public static void trace(String s) {
    System.out.println(s);
  }
}
