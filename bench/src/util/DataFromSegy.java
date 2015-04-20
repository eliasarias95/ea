package util;

import edu.mines.jtk.io.ArrayFile;
import java.io.IOException;
import java.nio.ByteOrder;

/**
 * A class to turn seismic .segy files into Java usable float arrays.
 * @author Elias Arias, Colorado School of Mines, CWP
 * @version 27.03.2015
 */
public class DataFromSegy {

  /**
   * Gets 2D data from a file and puts it into a 2D array of floats.
   * @param fileName name of the file including the path to the file.
   * @param n1 samples in the 1st dimension.
   * @param n2 samples in the 2nd dimension.
   * @return array[n2][n1] of floats.
   */
	public static float[][] grab2DFromFile(String fileName, int n1, int n2) {
		float[][] data = new float[n2][n1];
		try {
			ArrayFile af = new ArrayFile(
          fileName,"r",ByteOrder.BIG_ENDIAN,ByteOrder.LITTLE_ENDIAN);
			af.readFloats(data);
			System.out.println("Data extracted from "+fileName);
			af.close();
		} catch (IOException e) {
			System.out.println("File Cannot be found!");
			throw new RuntimeException(e);
		}
		return data;
	}

  /**
   * Gets 3D data from a file and puts it into a 3D array of floats.
   * @param fileName name of the file including the path to the file.
   * @param n1 samples in the 1st dimension.
   * @param n2 samples in the 2nd dimension.
   * @param n3 samples in the 3rd dimension.
   * @return array[n3][n2][n1] of floats.
   */
  public static float[][][] grab3DFromFile(
      String fileName, int n1, int n2, int n3) {
    float[][][] data = new float[n3][n2][n1];
    try {
      ArrayFile af = new ArrayFile(
          fileName,"r",ByteOrder.BIG_ENDIAN,ByteOrder.LITTLE_ENDIAN);
      af.readFloats(data);
      System.out.println("Data extracted from "+fileName);
      af.close();
    } catch (IOException e) {
      System.out.println("File Cannot be found!");
      throw new RuntimeException(e);
    }
    return data;
  }
}
