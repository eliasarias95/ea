import java.io.*;

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
  public Sfdip() {
    _both = "n";
    _verb = "n";
    _n4 = 2;
    _nj1 = 1;
    _nj2 = 1;
    _rect1 = 1;
    _rect2 = 1;
    _rect3 = 1;
    _order = 1;
    _niter = 5;
    _liter = 20;
    _p0 = 0.0f;
    _q0 = 0.0f;
    _pmin = -100;
    _pmax = 100;
    _qmin = -100;
    _qmax = 100;
  }

  /**
   * Constructor when given only the amount of information necessary to
   * run sfdip on a 2D data set.
   * @param both string options "y" and "n", if "y" compute both left
   * and right predictions
   * @param verb string options "y" and "n", verbosity flag...?
   * @param n4 what to compute in 3D, 0: in-line, 1: cross-line, 2: both
   * @param nj1 in-line antialiasing
   * @param rect1 dip smoothness in the 1st dimension
   * @param rect2 dip smoothness in the 2nd dimension
   * @param order accuracy order
   * @param niter number of outer iterations
   * @param liter number of linear iterations
   * @param p0 initial in-line dip
   * @param pmin minimum in-line dip
   * @param pmax maximum in-line dip
   */
  public Sfdip(String both, String verb, int n4, int nj1, 
      int rect1, int rect2, int order, int niter, int liter, 
      float p0, float pmin, float pmax) {
    _both = both;
    _verb = verb;
    _n4 = n4;
    _nj1 = nj1;
    _nj2 = 1;
    _rect1 = rect1;
    _rect2 = rect2;
    _rect3 = 1;
    _order = order;
    _niter = niter;
    _liter = liter;
    _p0 = p0;
    _q0 = 0.0f;
    _pmin = pmin;
    _pmax = pmax;
    _qmin = -100;
    _qmax = 100;
  }

  /**
   * Constructor when given all information.
   * @param both string options "y" and "n", if "y" compute both left
   * and right predictions
   * @param verb string options "y" and "n", verbosity flag...?
   * @param n4 what to compute in 3D, 0: in-line, 1: cross-line, 2: both
   * @param nj1 in-line antialiasing value
   * @param nj2 cross-line antialiasing value
   * @param rect1 dip smoothness in the 1st dimension
   * @param rect2 dip smoothness in the 2nd dimension
   * @param rect3 dip smoothness in the 3rd dimension
   * @param order accuracy order
   * @param niter number of outer iterations
   * @param liter number of linear iterations
   * @param p0 initial in-line dip
   * @param q0 initial cross-line dip
   * @param pmin minimum in-line dip
   * @param pmax maximum in-line dip
   * @param qmin minimum cross-line dip
   * @param qmax maximum cross-line dip
   */
  public Sfdip(String both, String verb, int n4, int nj1, int nj2, 
      int rect1, int rect2, int rect3, int order, int niter, int liter, 
      float p0, float q0, float pmin, float pmax, float qmin, float qmax) {
    _both = both;
    _verb = verb;
    _n4 = n4;
    _nj1 = nj1;
    _nj2 = nj2;
    _rect1 = rect1;
    _rect2 = rect2;
    _rect3 = rect3;
    _order = order;
    _niter = niter;
    _liter = liter;
    _p0 = p0;
    _q0 = q0;
    _pmin = pmin;
    _pmax = pmax;
    _qmin = qmin;
    _qmax = qmax;
  }

  /**
   * Set the in-line antialiasing parameter value.
   * @param nj1 in-line antialiasing value
   */
  public static void setNj(int nj1) {
    _nj1 = nj1;
    _nj2 = 1;
  }

  /**
   * Set the in-line and cross-line antialiasing parameter values.
   * @param nj1 in-line antialiasing value
   * @param nj2 cross-line antialiasing value
   */
  public static void setNj(int nj1, int nj2) {
    _nj1 = nj1;
    _nj2 = nj2;
  }

  /**
   * Set the dip smoothness in the 1st dimension.
   * @param rect1 dip smoothness in the 1st dimension
   */
  public static void setRect(int rect1) {
    _rect1 = rect1;
    _rect2 = 1;
    _rect3 = 1;
  }

  /**
   * Set the dip smoothness in the 1st and 2nd dimensions.
   * @param rect1 dip smoothness in the 1st dimension
   * @param rect2 dip smoothness in the 2nd dimension
   */
  public static void setRect(int rect1, int rect2) {
    _rect1 = rect1;
    _rect2 = rect2;
    _rect3 = 1;
  }

  /**
   * Set the dip smoothness in the 1st, 2nd, and 3rd dimensions.
   * @param rect1 dip smoothness in the 1st dimension
   * @param rect2 dip smoothness in the 2nd dimension
   * @param rect3 dip smoothness in the 3rd dimension
   */
  public static void setRect(int rect1, int rect2, int rect3) {
    _rect1 = rect1;
    _rect2 = rect2;
    _rect3 = rect3;
  }

  /**
   * Set the order of accuracy parameter value.
   * @param order accuracy order
   */
  public static void setOrder(int order) {
    _order = order;
  }

  /**
   * Set the number of outer iterations.
   * @param niter number of outer iterations
   */
  public static void setNiter(int niter) {
    _niter = niter;
  }

  /**
   * Set the number of linear iterations.
   * @param liter number of linear iterations
   */
  public static void setLiter(int liter) {
    _liter = liter;
  }

  /**
   * Set the initial in-line dip value.
   * @param p0 initial in-line dip
   */
  public static void setP0(float p0){
    _p0 = p0;
  }

  /**
   * Set the initial cross-line dip value.
   * @param q0 initial cross-line dip
   */
  public static void setQ0(float q0){
    _q0 = q0;
  }

  /**
   * Set the minimum and maximum in-line dip values.
   * @param pmin minimum in-line dip
   * @param pmax maximum in-line dip
   */
  public static void setPBounds(float pmin, float pmax) {
    _pmin = pmin;
    _pmax = pmax;
  }

  /**
   * Set the minimum and maximum cross-line dip values.
   * @param qmin minimum cross-line dip
   * @param qmax maximum cross-line dip
   */
  public static void setQBounds(float qmin, float qmax) {
    _qmin = qmin;
    _qmax = qmax;
  }

  /**
   * Set the both flag (either "y" or "n"), if "y" computes both left
   * and right predictions.
   * @param both "y" or "n" flag to compute both left
   * and right predictions or not
   */
  public static void setBoth(String both) {
    _both = both;
  }

  /**
   * Set the verbosity flag (either "y" or "n").
   * @param verb verbosity flag
   */
  public static void setVerb(String verb) {
    _verb = verb;
  }

  /**
   * Uses Madagascar to run sfdip to find the slopes in seismic data.
   * @param in the name of the input file
   * @param out the name of the output file
   */
  public static void findSlopes(String in, String out) {
    try {
      String cmd = "sfdip < "+in+" > "+out+" both="+_both+" n4="+_n4+ 
        " niter="+_niter+" liter="+_liter+" rect1="+_rect1+
        " rect2="+_rect2+" rect3="+_rect3+" p0="+_p0+" q0="+_q0+
        " order="+_order+" nj1="+_nj1+" nj2="+_nj2+" verb="+_verb+
        " pmin="+_pmin+" pmax="+_pmax+" qmin="+_qmin+
        " qmax="+_qmax;

      /**String cmd = "sfdip < "+in+" > "+out+" both="+_both+" n4="+_n4+ 
        " niter="+_niter+" liter="+_liter+" rect1="+_rect1+
        " rect2="+_rect2+" p0="+_p0+
        " order="+_order+" nj1="+_nj1+" verb="+_verb+
        " pmin="+_pmin+" pmax="+_pmax;
        */             
      //String cmd = "sfdip < "+in+" > "+out;

      Process pb = new ProcessBuilder("bash","-c",cmd).start();
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
    } catch (IOException e) {
      throw new RuntimeException(e);
    }
  }

  public static void main(String[] args) {
    
  } 

  /**Private Variables*/
  private String _both, _verb;
  private int _n4, _nj1, _nj2, _rect1, _rect2, _rect3, _order, _niter, _liter;
  private float _p0, _q0, _pmin, _pmax, _qmin, _qmax;
}
