#include <kaiser_window.h>

/**************************************************************************
  Copyright (c) 2004, Colorado School of Mines and others. All rights reserved.
  This program and accompanying materials are made available under the terms of
  the Common Public License - v1.0, which accompanies this distribution, and is
  available at http://www.eclipse.org/legal/cpl-v10.html
 **************************************************************************/
/**
 *@author Dave Hale, Colorado School of Mines
 *C++ translation by Elias Arias, Colorado School of Mines
 *@version 2005.08.04
 *@version 2015.07.23
 */
kaiser_window* kaiser_window::fromWidthAndLength(
    double width, double length) {
  double d = width*length;
  double a = 14.36*d+7.95;
  double error = pow(10.0,-a/20.0);
  return new kaiser_window(error,width,length);
}

double kaiser_window::evaluate(double x) {
  double xx = x*x;
  return (xx <= _xxmax)?_scale*ino(_alpha*sqrt(1.0-xx/_xxmax)):0.0;
}

double kaiser_window::getError() {
  return _error;
}

kaiser_window::kaiser_window(double error, double width, double length) {
  _error = error;
  _width = width;
  _length = length;
  double a = -20.0*log10(_error);
  if (a<=21.0) {
    _alpha = 0.0;
  }
  else if (a<=50.0) {
    _alpha = 0.5842*pow(a-21.0,0.4) + 0.07886*(a-21.0);
  }
  else {
    _alpha = 0.1102*(a-8.7);
  }
  _scale = 1.0/ino(_alpha);
  _xxmax = 0.25*_length*_length;
}

double kaiser_window::ino(double x) {
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
