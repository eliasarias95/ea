package slopes;

import edu.mines.jtk.dsp.Sampling;
import edu.mines.jtk.util.Check;

/**
 * Java class to allow for more java-like use of the sfdip code in Madagascar.
 * @author Elias Arias, Colorado School of Mines, CWP
 * @version 25.3.2015
 */
public class Sfdip {

  /**
   * Default constructor that initializes all fields to the default values 
   * used in Madagascar.
   */
  public Sfdip(float pmin, float pmax) {
    _pmin = pmin;
    _pmax = pmax;
    _qmin = pmin;
    _qmax = pmax;
  }

  /**
   * What to compute in 3D, 0: in-line, 1: cross-line, 2: both.
   * @param n4 value either 0, 1, or 2.
   */
  public void setN4(int n4) {
    Check.argument(n4==0 || n4==1 || n4==2, "must be an integer 0, 1 or 2");
    _n4 = n4;
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
   * @param s1 sampling of 1st dimension
   * @param s2 sampling of 2nd dimension
   * @param f the array[n2][n1] of inputs
   * @param p the array[n2][n1] of output slopes 
   */
  public void findSlopes(Sampling s1, Sampling s2, float[][] f, float[][] p) {
    RsfFilter rf = new RsfFilter("sfdip","both="+_both,
      "n4="+_n4,"niter="+_niter,"liter="+_liter,"rect1="+_rect1,
      "rect2="+_rect2,"rect3="+_rect3,"p0="+_p0,"q0="+_q0,
      "order="+_order,"nj1="+_nj1,"nj2="+_nj2,"verb="+_verb,
      "pmin="+_pmin,"pmax="+_pmax,"qmin="+_qmin,"qmax="+_qmax);
    rf.apply(s1,s2,f,p);
  }

  /**
   * Uses Madagascar to run sfdip to find the slopes in seismic data.
   * @param s1 sampling of 1st dimension
   * @param s2 sampling of 2nd dimension
   * @param s3 sampling of 3rd dimension
   * @param f the array[n3][n2][n1] of inputs
   * @param p2 the array[n3][n2][n1] of output in-line slopes 
   * @param p3 the array[n3][n2][n1] of output cross-line slopes 
   */
  public void findSlopes(Sampling s1, Sampling s2, Sampling s3, 
      float[][][] f, float[][][] p2, float[][][] p3) {
    RsfFilter rf = new RsfFilter("sfdip","both="+_both,
      "n4="+_n4,"niter="+_niter,"liter="+_liter,"rect1="+_rect1,
      "rect2="+_rect2,"rect3="+_rect3,"p0="+_p0,"q0="+_q0,
      "order="+_order,"nj1="+_nj1,"nj2="+_nj2,"verb="+_verb,
      "pmin="+_pmin,"pmax="+_pmax,"qmin="+_qmin,"qmax="+_qmax);
    rf.apply(s1,s2,s3,f,p2,p3);
  }

  /**
   * Uses Madagascar to run sfdip to find the slopes in seismic data.
   * @param s1 sampling of 1st dimension
   * @param s2 sampling of 2nd dimension
   * @param f the array[n2][n1] of inputs
   * @param pi the initial slope values
   * @param p the array[n2][n1] of output slopes 
   */
  public void findSlopes(Sampling s1, Sampling s2, float[][] f,
     float[][] pi, float[][] p) {
    RsfFilter rf = new RsfFilter("sfdip","both="+_both,
      "n4="+_n4,"niter="+_niter,"liter="+_liter,"rect1="+_rect1,
      "rect2="+_rect2,"rect3="+_rect3,"p0="+_p0,"q0="+_q0,
      "order="+_order,"nj1="+_nj1,"nj2="+_nj2,"verb="+_verb,
      "pmin="+_pmin,"pmax="+_pmax,"qmin="+_qmin,"qmax="+_qmax);
    rf.apply(s1,s2,f,p,pi);
  }

  /****************************Private Variables****************************/
  private String _both="n";
  private String _verb="n";
  private int _n4=1,_nj1=1,_nj2=1,_rect1=1,_rect2=1,_rect3=1,_order=1;
  private int _niter=5,_liter=20;
  private float _p0=0.0f, _q0=0.0f;
  private float _pmin, _pmax, _qmin, _qmax;
}
