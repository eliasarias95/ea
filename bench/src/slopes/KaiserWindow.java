/****************************************************************************
Copyright (c) 2004, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package slopes;

import static java.lang.Math.*;

/**
 * @author Dave Hale, Colorado School of Mines
 * @version 2005.08.04
 */
public class KaiserWindow {

  /**
   * Returns a Kaiser window with specified transition width and window length.
   * The product width*length cannot be less than one.
   * @param width the transition width
   * @param length the two-sided window length.
   * @return the window.
   */
  public static KaiserWindow fromWidthAndLength(double width, double length) {
    double d = width*length;
    double a = 14.36*d+7.95;
    double error = pow(10.0,-a/20.0);
    return new KaiserWindow(error,width,length);
  }

  /**
   * Returns the value of this Kaiser window function w(x) for specified x.
   * @param x the argument for which to evaluate w(x).
   * @return the value w(x).
   */
  public double evaluate(double x) {
    double xx = x*x;
    return (xx<=_xxmax)?_scale*ino(_alpha*sqrt(1.0-xx/_xxmax)):0.0;
  }

  /**
   * Gets the maximum absolute error.
   * @return the maximum absolute error.
   */
  public double getError() {
    return _error;
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private double _error;
  private double _width;
  private double _length;
  private double _alpha;
  private double _scale;
  private double _xxmax;
  private static final double DBL_EPSILON = 2.2204460492503131e-16d;

  private KaiserWindow(double error, double width, double length) {
    _error = error;
    _width = width;
    _length = length;
    double a = -20.0*log10(_error);
    if (a<=21.0) {
      _alpha = 0.0;
    } else if (a<=50.0) {
      _alpha = 0.5842*pow(a-21.0,0.4)+0.07886*(a-21.0);
    } else {
      _alpha = 0.1102*(a-8.7);
    }
    _scale = 1.0/ino(_alpha);
    _xxmax = 0.25*_length*_length;
  }

  private double ino(double x) {
    double s = 1.0;
    double ds = 1.0;
    double d = 0.0;
    do {
      d += 2.0;
      ds *= (x*x)/(d*d);
      s += ds;
    } while (ds>s*DBL_EPSILON);
    return s;
  }
}
