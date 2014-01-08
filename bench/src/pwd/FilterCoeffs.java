package pwd;

/**
 * All-pass plane-wave destruction filter coefficients. 
 * This is a Java implementation of the following file: apfilt.c
 * Located here:
 * http://sourceforge.net/p/rsf/code/HEAD/tree/trunk/user/pwd/apfilt.c#l4
 *
 * @author Elias Arias, Colorado School of Mines, CWP
 * @version 09.12.2013
 */
public class FilterCoeffs {

  /**
   * FilterCoeffs constructor that initializes n and b, apfilt_init.
   * @param nw the filter order
   */
  public FilterCoeffs(int nw) {
    double bk;
    _n = nw*2;
    _b = new double[_n+1];

    for (int k=0; k<_n+1; ++k) {
	    bk = 1.0;
	    for (int j=0; j<_n; ++j) {
	      if (j < _n-k) {
		      bk *= (k+j+1.0)/(2.0*(2.0*j+1.0)*(j+1.0));
	      } else {
		      bk *= 1.0/(2.0*(2.0*j+1.0));
	      }
	    }
	    _b[k] = bk;
    }        
  }

  /**
   * Finds the coefficients of the filter, passfilter.
   * @param p the slope
   * @param a the output filter [n+1]
   */
  public void passFilter(float p, float[] a) {
    double ak;
    
    for (int k=0; k<_n+1; ++k) {
	    ak = _b[k];
	    for (int j=0; j<_n; ++j) {
	      if (j < _n-k) {
		      ak *= (_n-j-p);
	      } else {
		      ak *= (p+j+1);
	      }
	    }
	    a[k] = (float)ak;
    }
  }

  /**
   * Find the coefficients of the filter's derivative, aderfilter.
   * @param p the slope
   * @param a the output filter [n+1]
   */
  public void derPassFilter(float p, float[] a) {
    double ak, ai;
    
    for (int k=0; k<_n+1; ++k) {
    	ak = 0.0;
	    for (int i=0; i<_n; ++i) {
	      ai = -1.0;
	        for (int j=0; j<_n; ++j) {
		        if (j != i) {			
		          if (j < _n-k) {
			          ai *= (_n-j-p);
		          } else {
			            ai *= (p+j+1.0);
		            }
		        } else if (j < _n-k) {
		            ai *= (-1.0);
		          }
	        }
	        ak += ai;
	    }
	    a[k] = (float)(ak*_b[k]);
    }
  }

/***********************Private Fields***********************/
  private int _n;
  private double[] _b;
}
