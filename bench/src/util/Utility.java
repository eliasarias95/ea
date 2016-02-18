import edu.mines.jtk.io.*;

import java.io.*;
import java.nio.*;
import java.util.Random;

/**
 *  
 * @author Elias Arias, Colorado School of Mines, CWP
 * @version 19.1.2015
 */

public class Utility {

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

  public static void trace(String s) {
    System.out.println(s);
  }
}
