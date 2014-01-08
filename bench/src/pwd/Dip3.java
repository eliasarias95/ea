package pwd;

/**
 * 3D dip estimation.
 * This is a Java implementation of the following file: dip3.c
 * Located here:
 * http://sourceforge.net/p/rsf/code/HEAD/tree/trunk/user/pwd/dip3.c
 *
 * @author Elias Arias, Colorado School of Mines, CWP
 * @version 09.12.2013
 */
public class Dip3 {
  private float eps;
  private float[] u1,u2,dp,p0;
  private int n, n1, n2, n3; 
  private int[] nn = new int[3];

  /**
   * Dip constructor, dip3_init
   * @param m1
   * @param m2
   * @param m3
   * @param rect 
   */
  public Dip3(int m1, int m2, int m3, int[] rect, 
              int niter, float eps1, boolean verb) {
    n1=m1;
    n2=m2;
    n3=m3;
    n = n1*n2*n3;
    eps = eps1;

    u1 = new float[n];
    u2 = new float[n];   
    dp = new float[n];
    p0 = new float[n];

    nn[0]=n1;
    nn[1]=n2;
    nn[2]=n3;

    //sf_divn_init (3, n, nn, rect, niter, verb);
  }
}
