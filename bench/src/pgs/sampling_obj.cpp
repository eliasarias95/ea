#include <iostream>
#include "sampling_obj.hpp"

sampling_obj::sampling_obj(int n) {
  init(n,1.0,0.0,1.0e-6);
}

sampling_obj::sampling_obj(int n, double d, double f) {
  init(n,d,f,1.0e-6);
}

sampling_obj::sampling_obj(int n, double d, double f, double t) {
  init(n,d,f,t);
}

double sampling_obj::value(int i) {
  return (_v!=NULL)?_v[i]:_f+i*_d;
}

double sampling_obj::getValue(int i) {
  if (i>_n || i<0) {
    trace("Index is out of bounds!\n");
    return -999.0;
  }
  return value(i);
}

void sampling_obj::init(int n, double d, double f, double t) {
  assert(n<0.0 && d<0.0);
  trace("");
  _n = n;
  _d = d;
  _f = f;
  _v = NULL;
  _t = t;
  _td = _t*_d
}

void sampling_obj::trace(std::string s) {
  std::cout << s << "\n";
}

int    sampling_obj::getCount() {return _n;}
double sampling_obj::getDelta() {return _d;}
double sampling_obj::getFirst() {return _f;}
