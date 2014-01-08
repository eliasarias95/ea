package pwd;

import edu.mines.jtk.util.Check;

/**
 * 3D plane-wave destruction filter. 
 * This is a Java implementation of the following file: allp3.c
 * Located here:
 * http://sourceforge.net/p/rsf/code/HEAD/tree/trunk/user/pwd/allp3.c#l242
 * @author Elias Arias, Colorado School of Mines, CWP
 * @version 09.12.2013
 */
public class Filter3 {
  
  /**
   * Filter3 constructor that initializes the following variables.
   * @param n1 data size in 1st-dimension
   * @param n2 data size in 2nd-dimension
   * @param n3 data size in 3rd-dimension
   * @param nw filter size
   * @param nj filter step
   * @param pp dip [n1*n2*n3]
   */
  public Filter3(int n1, int n2, int n3, int nw, int nj, float[] pp) {
    _n1 = n1;
    _n2 = n2;
    _n3 = n3;
    _nw = nw;
    _nj = nj;
    _pp = pp;
    _flt = new float[2*_nw+1];
    _fc = new FilterCoeffs(_nw);
  }  

  /**
   * Inline plane-wave destruction, allpass1
   * @param left left or right prediction
   * @param der derivative flag
   * @param x input
   * @param y output
   */
  public void inlineDestructor(boolean left, boolean der, 
                               float[] x, float[] y) {
    int i1, i2, i3, iw, is, i, n1, n2, n3, j1, j2, ip;

    n1 = _n1;
    n2 = _n2;
    n3 = _n3;

    if (left) {
	    j1=1;
      j2=n2;
      ip=-n1;
    } else {
	    j1=0;
      j2=n2-1;
      ip=n1;
    }

    for (i3=0; i3<n3; ++i3) {
	    for (i2=0; i2<n2; ++i2) {
	      for (i1=0; i1<n1; ++i1) {
		      i = i1 + n1*(i2+n2*i3);
		      y[i] = 0.0f;
	      }
	    }
    }
  
    for (i3=0; i3<n3; ++i3) {
	    for (i2=j1; i2<j2; ++i2) {
	      for (i1=_nw*_nj; i1<n1-_nw*_nj; ++i1) {
		      i = i1 + n1*(i2+n2*i3);
		      if (der) {
		        _fc.derPassFilter(_pp[i],_flt);
		      } else {
		        _fc.passFilter(_pp[i],_flt);
		        }
		        for (iw=0; iw<=2*_nw; ++iw) {
		          is = (iw-_nw)*_nj;
		          y[i] += (x[i+is+ip] - x[i-is])*_flt[iw];
		        }
	      }
	    }
    }
  }

  /**
   * Inline plane-wave destruction, allpass2
   * @param left left or right prediction
   * @param der derivative flag
   * @param x input
   * @param y output
   */
  public void crosslineDestructor(boolean left, boolean der, 
                                  float[] x, float[] y) {
    int i1, i2, i3, iw, j1, j2, ip, is, i;
    int n1, n2, n3;

    n1 = _n1;
    n2 = _n2;
    n3 = _n3;

    if (left) {
	    j1=1; 
      j2=n3;   
      ip=-n1*n2;
    } else {
	    j1=0; 
      j2=n3-1; 
      ip=n1*n2;
    }
    
    for (i3=0; i3<n3; ++i3) {
	    for (i2=0; i2<n2; ++i2) {
	      for (i1=0; i1<n1; ++i1) {
		      i = i1 + n1*(i2+n2*i3);
		      y[i] = 0.0f;
	      }
	    }
    }
    
    for (i3=j1; i3<j2; ++i3) {
	    for (i2=0; i2<n2; ++i2) {
	      for (i1=_nw*_nj; i1<n1-_nw*_nj; ++i1) {
		      i = i1 + n1*(i2+n2*i3);
		      if (der) {
		        _fc.derPassFilter(_pp[i], _flt);
		      } else {
		        _fc.passFilter(_pp[i], _flt);
		      }
		        for (iw=0; iw<=2*_nw; ++iw) {
		          is = (iw-_nw)*_nj;
		          y[i] += (x[i+is+ip] - x[i-is])*_flt[iw];
		        }
	        }
	      }
      }
    }

  /**
   * Plane-wave destructor as linear operator, allpass3_lop
   * @param adj
   * @param add
   * @param m1
   * @param m2
   * @param x
   * @param y
   * @param f
   */
  public void linopDestructor(boolean adj, boolean add, int m1, int m2, 
                              float[] x, float[] y, Filter3 f){
    int i, i1, i2, i3, iw, is, n1, n2, n3, nw, nj;
 
    Check.argument(m2 != 2*m1,
                   "size mismatch: m2 != 2*m1");
    n1 = _n1;
    n2 = _n2;
    n3 = _n3;
    nw = _nw;
    nj = _nj;
 
    Check.argument(n1*n2*n3 != m1,
                   "size mismatch: n1*n2*n3 != m1");

    for (i3=0; i3<n3; ++i3) {
      for (i2=0; i2<n2-1; ++i2) {
        for (i1=nw*nj; i1<n1-nw*nj; ++i1) {
          i = i1 + n1*(i2 + n2*i3);
          _fc.passFilter(_pp[i], _flt);
          for (iw=0; iw<=2*nw; ++iw) {
            is = (iw-nw)*nj;
            if (adj) {
              x[i+n1+is] += y[i] * _flt[iw];
              x[i-is]    -= y[i] * _flt[iw];
            } else {
              y[i] += (x[i+n1+is] - x[i-is]) * _flt[iw];
            }
          }
        }
      }
    }

    n1 = f._n1;
    n2 = f._n2;
    n3 = f._n3;
    nw = f._nw;
    nj = f._nj;

    Check.argument(n1*n2*n3 != m1,
                   "size mismatch: n1*n2*n3 != m1");

    for (i3=0; i3<n3-1; ++i3) {
      for (i2=0; i2<n2; ++i2) {
        for (i1=nw*nj; i1<n1-nw*nj; ++i1) {
          i = i1 + n1*(i2 + n2*i3);
          f._fc.passFilter(f._pp[i], f._flt);
          for (iw = 0; iw <= 2*nw; ++iw) {
            is = (iw-nw)*nj;
            if (adj) {
              x[i+n1*n2+is] += y[i+m1] * f._flt[iw];
              x[i-is]       -= y[i+m1] * f._flt[iw];
            } else {
              y[i+m1] += (x[i+n1*n2+is] - x[i-is]) * f._flt[iw];
            }
          }
        }
      }
    }
  }

  /***********************Private Fields***********************/
  private int _n1,_n2,_n3,_nw,_nj;
  private float[] _flt,_pp;
  private FilterCoeffs _fc;
}
