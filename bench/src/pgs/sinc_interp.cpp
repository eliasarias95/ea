#include <sinc_interp.h>

/**************************************************************************
  Copyright (c) 2012, Colorado School of Mines and others. All rights reserved.
  This program and accompanying materials are made available under the terms of
  the Common Public License - v1.0, which accompanies this distribution, and is
  available at http://www.eclipse.org/legal/cpl-v10.html
 ***************************************************************************/

/**
 *@author Dave Hale, Colorado School of Mines
 *@author Bill Harlan, Landmark Graphics
 *C++ translation by Elias Arias, Colorado School of Mines 
 *@version 2012.12.21
 *@version 2015.07.23
 */
public:
enum class Extrapolation {
  ZERO,
  CONSTANT,
};

static sinc_interp *sinc_interp::fromErrorAndLength(double emax, int lmax) {
  return new sinc_interp(emax,0.0,lmax);
}

static sinc_interp *sinc_interp::fromErrorAndFrequency(
    double emax, double fmax) {
  return new sinc_interp(emax,fmax,0);
}

static sinc_interp *sinc_interp::fromFrequencyAndLength(
    double fmax, int lmax) {
  return new sinc_interp(0.0,fmax,lmax);
}

sinc_interp() : sinc_interp(0.0,0.3,8) {}

virtual double sinc_interp::getMaximumError() {
  return _table->design->emax;
}

virtual double sinc_interp::getMaximumFrequency() {
  return _table->design->fmax;
}

virtual int sinc_interp::getMaximumLength() {
  return _table->design->lmax;
}

virtual long long sinc_interp::getTableBytes() {
  long long nbytes = 4LL;
  nbytes *= _table->lsinc;
  nbytes *= _table->nsinc;
  return nbytes;
}

virtual Extrapolation sinc_interp::getExtrapolation() {
  return _extrap;
}

virtual void sinc_interp::setExtrapolation(Extrapolation extrap) {
  _extrap = extrap;
}

virtual float sinc_interp::interpolate(
    int nxu, double dxu, double fxu, float yu[], double xi) {
  double xscale = 1.0/dxu;
  double xshift = _lsinc-fxu*xscale;
  int nxum = nxu-_lsinc;
  return interpolate(xscale,xshift,nxum,nxu,yu,xi);
}

virtual void sinc_interp::interpolate(
    int nxu, double dxu, double fxu, float yu[],
    int nxi, float xi[], float yi[]) {
  double xscale = 1.0/dxu;
  double xshift = _lsinc-fxu*xscale;
  int nxum = nxu-_lsinc;
  for (int ixi=0; ixi<nxi; ++ixi) {
    yi[ixi] = interpolate(xscale,xshift,nxum,nxu,yu,xi[ixi]);
  }
}

virtual void sinc_interp::interpolate(
    int nxu, double dxu, double fxu, float yu[], 
    int nxi, double dxi, double fxi, float yi[]) {
  if (dxu == dxi) {
    shift(nxu,dxu,fxu,yu,nxi,fxi,yi);
  }
  else {
    double xscale = 1.0/dxu;
    double xshift = _lsinc-fxu*xscale;
    int nxum = nxu-_lsinc;
    for (int ixi=0; ixi<nxi; ++ixi) {
      double xi = fxi+ixi*dxi;
      yi[ixi] = interpolate(xscale,xshift,nxum,nxu,yu,xi);
    }
  }
}

virtual float sinc_interp::interpolate(axis *sxu, float yu[], double xi) {
  return interpolate(sxu->n,sxu->d,sxu->o,yu,xi);
}

virtual void sinc_interp::interpolate(
    axis *sxu, float yu[], axis *sxi, float yi[]) {
  int nxu = sxu->n;
  int nxi = sxi->n;
  double xscale = 1.0/sxu->d;
  double xshift = _lsinc-sxu->o*xscale;
  int nxum = nxu-_lsinc;
  for (int ixi=0; ixi<nxi; ++ixi) {
    double xi = sxi->get_val(ixi);
    yi[ixi] = interpolate(xscale,xshift,nxum,nxu,yu,xi);
  }
}

virtual void sinc_interp::accumulate(
    double xa, float ya, int nxu, double dxu, double fxu, float yu[]) {
  double xscale = 1.0 / dxu;
  double xshift = _lsinc - fxu*xscale;
  int nxum = nxu - _lsinc;
  accumulate(xscale,xshift,nxum,xa,ya,nxu,yu);
}

virtual void sinc_interp::accumulate(
    int nxa, float xa[], float ya[], 
    int nxu, double dxu, double fxu, float yu[]) {
  double xscale = 1.0 / dxu;
  double xshift = _lsinc - fxu*xscale;
  int nxum = nxu - _lsinc;
  for (int ixa = 0; ixa < nxa; ++ixa) {
    accumulate(xscale,xshift,nxum,xa[ixa],ya[ixa],nxu,yu);
  }
}

virtual float **sinc_interp::getTable() {
  return copy(_table->asinc);
}

virtual int sinc_interp::getNumberInTable() {
  return _table->asinc->length;
}

virtual int sinc_interp::getLengthInTable() {
  return _table->asinc[0].length;
}

sinc_interp::sinc_interp(double emax, double fmax, int lmax) {
  _table = getTable(emax,fmax,lmax);
  _lsinc = _table->lsinc;
  _nsinc = _table->nsinc;
  _nsincm1 = _table->nsincm1;
  _ishift = _table->ishift;
  _dsinc = _table->dsinc;
  _asinc = _table->asinc;
}

private:
class Design {
  public:
    double emax = 0;
    double fmax = 0;
    int lmax = 0;
    Design(double emax, double fmax, int lmax) {
      this->emax = emax;
      this->fmax = fmax;
      this->lmax = lmax;
    }
};

private:
class Table {
  public:
    Design *design; // here, all three design parameters are non-zero
    int lsinc = 0, nsinc = 0, nsincm1 = 0, ishift = 0;
    double dsinc = 0;
    //need to delete somewhere below (where appropriate).
    float **asinc;
};

private:
static Table *sinc_interp::makeTable(Design *design) {
  double emax = design->emax;
  double fmax = design->fmax;
  int lmax = design->lmax;

  // The Kaiser window transition width is twice the difference 
  // between the Nyquist frequency 0.5 and the maximum frequency.
  double wwin = 2.0*(0.5 - fmax);

  // The Kaiser window accounts for a hard-wired fraction of the maximum 
  // interpolation error. The other error will be due to table lookup.
  double ewin = emax*EWIN_FRAC;
  kaiser_window *kwin = kaiser_window::fromWidthAndLength(wwin,lmax);
  ewin = 3.0*kwin->getError();
  emax = ewin / EWIN_FRAC;
  double etabMin = 1.1*PI*fmax / (NTAB_MAX - 1);
  double emaxMin = etabMin / (1.0 - EWIN_FRAC);
  if (emax < emaxMin) {
    emax = emaxMin;
    ewin = emax*EWIN_FRAC;
  }

  double etab = emax - ewin;
  double dsinc = (fmax > 0.0)?etab / (PI*fmax):1.0;
  int nsincMin = 1 + static_cast<int>(ceil(1.0 / dsinc));
  int nsinc = 2;
  while (nsinc < nsincMin) {
    nsinc *= 2;
  }
  ++nsinc;
  int lsinc = lmax;
  Table *table = makeTable(nsinc,lsinc,kwin);
  table->design = new Design(emax,fmax,lmax);
  _tables[design] = table; // key is design with one zero parameter
  return table;
}

static Table *sinc_interp::makeTable(
    int nsinc, int lsinc, kaiser_window *kwin) {
  float asinc[nsinc][lsinc];
  int nsincm1 = nsinc - 1;
  int ishift = -lsinc - lsinc / 2 + 1;
  double dsinc = 1.0 / (nsinc - 1);

  // The first and last interpolators are shifted unit impulses.
  // Handle these two cases exactly, with no rounding errors.
  for (int j = 0; j < lsinc; ++j) {
    asinc[0      ][j] = 0.0f;
    asinc[nsinc-1][j] = 0.0f;
  }
  asinc[0      ][lsinc/2-1] = 1.0f;
  asinc[nsinc-1][lsinc/2  ] = 1.0f;

  // Other interpolators are sampled Kaiser-windowed sinc functions.
  for (int isinc=1; isinc<nsinc-1; ++isinc) {
    double x = -lsinc / 2 + 1 - dsinc*isinc;
    for (int i = 0; i < lsinc; ++i,x += 1.0) {
      asinc[isinc][i] = static_cast<float>(sinc(x)*kwin->evaluate(x));
    }
  }
  Table *table = new Table();
  table->lsinc = lsinc;
  table->nsinc = nsinc;
  table->nsincm1 = nsincm1;
  table->ishift = ishift;
  table->dsinc = dsinc;
  table->asinc = asinc;
  return table;
}
static double sinc_interp::sinc(double x) {
  return (x != 0.0)?sin(PI*x) / (PI*x):1.0;
}

static Table *sinc_interp::getTable(double emax, double fmax, int lmax) {
  Design *design = new Design(emax,fmax,lmax);
  Table *table = _tables[design];
  if (table == nullptr) {
    table = makeTable(design);
  }
  return table;
}

float sinc_interp::interpolate(
    double xscale, double xshift, int nxum, int nxu, float yu[], double x) {
  // Which uniform samples?
  double xn = xshift + x*xscale;
  int ixn = static_cast<int>(xn);
  int kyu = _ishift + ixn;

  // Which sinc approximation?
  double frac = xn - ixn;
  if (frac < 0.0) {
    frac += 1.0;
  }
  int ksinc = static_cast<int>(frac*_nsincm1 + 0.5);
  //need to delete somewhere below (where appropriate).
  float *asinc = _asinc[ksinc];

  // If no extrapolation is necessary, use a fast loop.
  // Otherwise, extrapolate uniform samples, as necessary.
  float yr = 0.0f;
  if (kyu >= 0 && kyu <= nxum) {
    for (int isinc = 0; isinc < _lsinc; ++isinc,++kyu) {
      yr += yu[kyu]*asinc[isinc];
    }
  }
  else if (_extrap == Extrapolation::ZERO) {
    for (int isinc = 0; isinc < _lsinc; ++isinc,++kyu) {
      if (0 <= kyu && kyu < nxu) {
        yr += yu[kyu]*asinc[isinc];
      }
    }
  }
  else if (_extrap == Extrapolation::CONSTANT) {
    for (int isinc = 0; isinc < _lsinc; ++isinc,++kyu) {
      int jyu = (kyu < 0)?0:(nxu <= kyu)?nxu - 1:kyu;
      yr += yu[jyu]*asinc[isinc];
    }
  }
  return yr;
}

void sinc_interp::shift(
    int nxu, double dxu, double fxu, float yu[],
    int nxi,             double fxi, float yi[]) {
  double lxu = fxu + (nxu - 1)*dxu;
  double xscale = 1.0 / dxu;
  double xshift = _lsinc - fxu*xscale;
  int nxum = nxu - _lsinc;

  // Which output samples are near beginning and end of uniform sequence?
  double dx = dxu;
  double x1 = fxu + dxu*_lsinc / 2;
  double x2 = lxu - dxu*_lsinc / 2;
  double x1n = (x1 - fxi) / dx;
  double x2n = (x2 - fxi) / dx;
  int ix1 = max(0,min(nxi,static_cast<int>(x1n) + 1));
  int ix2 = max(0,min(nxi,static_cast<int>(x2n) - 1));

  // Interpolate output samples near beginning of uniform sequence.
  for (int ixi = 0; ixi < ix1; ++ixi) {
    double xi = fxi + ixi*dx;
    yi[ixi] = interpolate(xscale,xshift,nxum,nxu,yu,xi);
  }

  // Interpolate output samples near end of uniform sequence.
  for (int ixi = ix2; ixi < nxi; ++ixi) {
    double xi = fxi + ixi*dx;
    yi[ixi] = interpolate(xscale,xshift,nxum,nxu,yu,xi);
  }

  // Now we ignore the ends, and use a single sinc approximation.

  // Which uniform samples?
  double xn = xshift + (fxi + ix1*dx)*xscale;
  int ixn = static_cast<int>(xn);
  int kyu = _ishift + ixn;

  // Which sinc approximation?
  double frac = xn - ixn;
  if (frac < 0.0) {
    frac += 1.0;
  }
  int ksinc = static_cast<int>(frac*_nsincm1 + 0.5);
  //need to delete somewhere below (where appropriate).
  float *asinc = _asinc[ksinc];

  // Interpolate for output indices ix1 <= ix <= ix2.
  for (int ix = ix1; ix < ix2; ++ix,++kyu) {
    float yr = 0.0f;
    for (int isinc = 0,jyu = kyu; isinc < _lsinc; ++isinc,++jyu) {
      yr += yu[jyu]*asinc[isinc];
    }
    yi[ix] = yr;
  }
}

void sinc_interp::accumulate(
    double xscale, double xshift, int nxum, 
    double x, float y, int nxu, float yu[]) {
  // Which uniform samples?
  double xn = xshift + x*xscale;
  int ixn = static_cast<int>(xn);
  int kyu = _ishift + ixn;

  // Which sinc approximation?
  double frac = xn - ixn;
  if (frac < 0.0) {
    frac += 1.0;
  }
  int ksinc = static_cast<int>(frac*_nsincm1 + 0.5);
  //need to delete somewhere below (where appropriate).
  float *asinc = _asinc[ksinc];

  // If no extrapolation is necessary, use a fast loop.
  // Otherwise, extrapolate uniform samples, as necessary.
  if (kyu >= 0 && kyu <= nxum) {
    for (int isinc = 0; isinc < _lsinc; ++isinc,++kyu) {
      yu[kyu] += y*asinc[isinc];
    }
  }
  else if (_extrap == Extrapolation::ZERO) {
    for (int isinc = 0; isinc < _lsinc; ++isinc,++kyu) {
      if (0 <= kyu && kyu < nxu) {
        yu[kyu] += y*asinc[isinc];
      }
    }
  }
  else if (_extrap == Extrapolation::CONSTANT) {
    for (int isinc = 0; isinc < _lsinc; ++isinc,++kyu) {
      int jyu = (kyu < 0)?0:(nxu <= kyu)?nxu - 1:kyu;
      yu[jyu] += y*asinc[isinc];
    }
  }
}
