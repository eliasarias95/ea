package pwd;

import edu.mines.jtk.io.ArrayInputStream;
import edu.mines.jtk.io.ArrayOutputStream;

import java.io.*;
import java.nio.ByteOrder;

/**
 * Java class to allow for more java-like use of the sfdip code in Madagascar.
 * @author Elias Arias, Colorado School of Mines, CWP
 * @version 20.1.2014
 */
public class Sfdip {

  /**
   * Default constructor that initializes all fields to the default values 
   * used in Madagascar.
   */
  public Sfdip(float pmin, float pmax) {
    _pmin = pmin;
    _pmax = pmax;
  }

  /**
   * Set the in-line antialiasing parameter value.
   * @param nj1 in-line antialiasing value
   */
  public void setNj(int nj1) {
    setNj(nj1,1);
  }

  /**
   * Set the in-line and cross-line antialiasing parameter values.
   * @param nj1 in-line antialiasing value
   * @param nj2 cross-line antialiasing value
   */
  public void setNj(int nj1, int nj2) {
    _nj1 = nj1;
    _nj2 = nj2;
  }

  /**
   * Set the dip smoothness in the 1st dimension.
   * @param rect1 dip smoothness in the 1st dimension
   */
  public void setRect(int rect1) {
    setRect(rect1,1,1);
  }

  /**
   * Set the dip smoothness in the 1st and 2nd dimensions.
   * @param rect1 dip smoothness in the 1st dimension
   * @param rect2 dip smoothness in the 2nd dimension
   */
  public void setRect(int rect1, int rect2) {
    setRect(rect1,rect2,1);
  }

  /**
   * Set the dip smoothness in the 1st, 2nd, and 3rd dimensions.
   * @param rect1 dip smoothness in the 1st dimension
   * @param rect2 dip smoothness in the 2nd dimension
   * @param rect3 dip smoothness in the 3rd dimension
   */
  public void setRect(int rect1, int rect2, int rect3) {
    _rect1 = rect1;
    _rect2 = rect2;
    _rect3 = rect3;
  }

  /**
   * Set the order of accuracy parameter value.
   * @param order accuracy order
   */
  public void setOrder(int order) {
    _order = order;
  }

  /**
   * Set the number of outer iterations.
   * @param niter number of outer iterations
   */
  public void setNiter(int niter) {
    _niter = niter;
  }

  /**
   * Set the number of linear iterations.
   * @param liter number of linear iterations
   */
  public void setLiter(int liter) {
    _liter = liter;
  }

  /**
   * Set the initial in-line dip value.
   * @param p0 initial in-line dip
   */
  public void setP0(float p0){
    _p0 = p0;
  }

  /**
   * Set the initial cross-line dip value.
   * @param q0 initial cross-line dip
   */
  public void setQ0(float q0){
    _q0 = q0;
  }

  /**
   * Set the minimum and maximum in-line dip values.
   * @param pmin minimum in-line dip
   * @param pmax maximum in-line dip
   */
  public void setPBounds(float pmin, float pmax) {
    _pmin = pmin;
    _pmax = pmax;
  }

  /**
   * Set the minimum and maximum cross-line dip values.
   * @param qmin minimum cross-line dip
   * @param qmax maximum cross-line dip
   */
  public void setQBounds(float qmin, float qmax) {
    _qmin = qmin;
    _qmax = qmax;
  }

  /**
   * Set the both flag (either "y" or "n"), if "y" computes both left
   * and right predictions.
   * @param both "y" or "n" flag to compute both left
   * and right predictions or not
   */
  public void setBoth(String both) {
    _both = both;
  }

  /**
   * Set the verbosity flag (either "y" or "n").
   * @param verb verbosity flag
   */
  public void setVerb(String verb) {
    _verb = verb;
  }

  /**
   * Uses Madagascar to run sfdip to find the slopes in seismic data.
   * @param x the array[n2][n1] of inputs
   * @param y the array[n2][n1] of outputs 
   */
  public void findSlopes(float[][] x, float[][] y) {
    int n2 = x.length;
    int n1 = x[0].length;

    String in_dat = "this_file_in.dat";
    String out_dat = "this_file_out.dat";
    writeBinary(x,in_dat);

    String new_rsf = "new.rsf";
    String out_rsf = "this_file_out.rsf";
    String gom_rsf = "gom.rsf";
    String file_cmd1 = "sfspike < "+gom_rsf+" > "+new_rsf;
    String file_cmd2 = "sfspike < "+gom_rsf+" > "+out_rsf;
    String file_cmd3 = "echo in="+in_dat+" >> "+new_rsf;
    String dip_cmd = "sfdip < "+new_rsf+" > "+out_rsf+" both="+_both+
      " n4="+_n4+" niter="+_niter+" liter="+_liter+" rect1="+_rect1+
      " rect2="+_rect2+" rect3="+_rect3+" p0="+_p0+" q0="+_q0+
      " order="+_order+" nj1="+_nj1+" nj2="+_nj2+" verb="+_verb+
      " pmin="+_pmin+" pmax="+_pmax+" qmin="+_qmin+" qmax="+_qmax;
    String file2array_cmd = "mv /var/tmp/"+out_rsf+"@ "+out_dat;

    String[] cmd1 = {"bash","-c",file_cmd1};
    String[] cmd2 = {"bash","-c",file_cmd2};
    String[] cmd3 = {"bash","-c",file_cmd3};
    String[] cmd4 = {"bash","-c",dip_cmd};
    String[] cmd5 = {"bash","-c",file2array_cmd};
    pb(cmd1);
    pb(cmd2);
    pb(cmd3);
    pb(cmd4);
    pb(cmd5);
    readImage(out_dat,y);
  }

  /*********************************Private*********************************/

  /**
   * A method that allows the user to call it multiple times when command line
   * commands must be run one after the other.
   * @param cmd the array of string commands to be run on the command line
   */
  private void pb(String[] cmd) {
    try {
      Process pb = new ProcessBuilder(cmd).start();
      InputStreamReader isri = new InputStreamReader(pb.getInputStream());
      BufferedReader stdInput = new BufferedReader(isri);

      InputStreamReader isre = new InputStreamReader(pb.getErrorStream());
      BufferedReader stdError = new BufferedReader(isre);

      System.out.println("Here is the standard output of the command:\n");
      String s1;
      String s2;
      while ((s1=stdInput.readLine()) !=null) {
        System.out.println(s1);
      }

      System.out.println("Here is the standard error of the command:\n");
      while ((s2=stdError.readLine()) !=null) {
        System.out.println(s2);

      }
      //pb.waitFor();
    } catch (IOException e) {
      throw new RuntimeException(e);
    }
  }

  /**
   * Reads a binary file.
   * @param fileName the name of the file to be read
   * @param x the array[n2][n1] of output data read from the file
   */
  private static void readImage(String fileName, float[][] x) {
    try {
      ArrayInputStream ais = new ArrayInputStream(fileName);
      ais.readFloats(x);
      ais.close();
    } catch (IOException e) {
      throw new RuntimeException(e);
    }
  }

  /**
   * Writes seismic data to binary file.
   * @param x array[n2][n1] of data to write to the binary file
   * @param fileName name of output binary file
   */
  private static void writeBinary(float[][] x, String fileName) {
    ByteOrder byteOrder = ByteOrder.BIG_ENDIAN;
    try {
      ArrayOutputStream aos = new ArrayOutputStream(fileName,byteOrder);
      aos.writeFloats(x);
      aos.close();
    } catch (IOException e) {
      throw new RuntimeException(e);
    }
  }

  /****************************Private Variables****************************/
  private String _both="n";
  private String _verb="n";
  private int _n4=1,_nj1=1,_nj2=1,_rect1=1,_rect2=1,_rect3=1,_order=1;
  private int _niter=5,_liter=20;
  private float _p0=0.0f, _q0=0.0f;
  private float _pmin, _pmax;
  private float _qmin=-100.0f, _qmax=-100.0f;
}
