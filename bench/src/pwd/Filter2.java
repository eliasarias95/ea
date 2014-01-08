package pwd;

/**
 * 2D Plane-wave destruction filter.
 * This is a Java implementation of the following file: allp2.c
 * Located here: 
 * http://sourceforge.net/p/rsf/code/HEAD/tree/trunk/user/pwd/allp2.c#l1
 * @author Elias Arias, Colorado School of Mines, CWP
 * @version 09.12.2013
 */

public class Filter2 {
  /**
   * A constructor for a 2D plane-wave destruction filter.
   * @param n1 data size in the 1st-dimension
   * @param n2 data size in the 2nd-dimension
   * @param nw filter order
   * @param nj filter step
   * @param pp dip [ny][nx]
   */
  public Filter2(int n1, int n2, int nw, int nj, float[][] pp) {
    _n1 = n1;
    _n2 = n2;
    _nw = nw;
    _nj = nj;
    _pp = pp;
    _flt = new float[2*_nw+1];
    _fc = new FilterCoeffs(_nw);
  }

  /**
   * Plane-wave destructor as linear operator, allpass21_lop.
   * @param adj
   * @param add
   * @param n1
   * @param n2
   * @param x
   * @param y
   * @param f
   */
  public void linopDestructor(boolean adj, boolean add, int n1, int n2,
      float[] x, float[] y, Filter2 f) {

    int j;

    if(add) return;

    if(adj) {
      for (j=0; j<n1; ++j) {
        x[j] = 0.0f;
      }
    } else {
      for (j=0; j<n2; ++j) {
        y[j] = 0.0f;
      }
    }

    int i, i1, i2, iw, is;

    for (i2=0; i2<_n2-1; ++i2) {
      for (i1=_nw*_nj; i1<_n1-_nw*_nj; ++i1) {
        _fc.passFilter(_pp[i2][i1], _flt);
        i = i1 + i2*_n1;
        for (iw=0; iw<=2*_nw; ++iw) {
          is = (iw-_nw)*_nj;
          if (adj) {
            x[i+is+_n1] += y[i]*_flt[iw];
            x[i-is]     -= y[i]*_flt[iw];
          } else {
            y[i] += (x[i+is+_n1] - x[i-is]) * _flt[iw];
          }
        }
      }
    }
  }

  /**
   * Plane-wave destruction, allpass21
   * @param der derivative flag
   * @param x input 
   * @param y output 
   */
  public void destructor(boolean der, float[][] x, float[][] y) {
    int i1, i2, iw, is;

    for (i2=0; i2<_n2-1; ++i2) {
      for (i1=_nw*_nj; i1<_n1-_nw*_nj; ++i1) {
        if (der) {
          _fc.derPassFilter(_pp[i2][i1],_flt);
        } else {
          _fc.passFilter(_pp[i2][i1],_flt);
        }
        for (iw=0; iw<=2*_nw; ++iw) {
          is = (iw-_nw)*_nj;
          y[i2][i1] += (x[i2+1][i1+is] - 
                        x[i2  ][i1-is])*_flt[iw];
        }
      }
    }
  }

  /***********************Private Fields***********************/
  private int _n1,_n2,_nw,_nj;
  private float[] _flt;
  private float[][] _pp;
  private FilterCoeffs _fc;
}
