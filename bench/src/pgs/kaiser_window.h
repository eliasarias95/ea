#ifndef KAISER_WINDOW_H
#define KAISER_WINDOW_H
#include <cmath>

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
	class kaiser_window {

  public:
	  static kaiser_window *fromWidthAndLength(double width, double length);
	  virtual double evaluate(double x);
	  virtual double getError();

  private:
	  double _error, _width, _length, _alpha, _scale, _xxmax;
	  static const double DBL_EPSILON = 2.2204460492503131e-16;

	  kaiser_window(double error, double width, double length);
	  double ino(double x);
	};
#endif
