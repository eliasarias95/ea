#ifndef SDW_H
#define SDW_H

#include <ucsl.h>
#include <sampling_obj.cpp>

class sdw_obj {
  private:
    double _r1min, _r2min, _r3min;
    double _r1max, _r2max, _r3max;
    float _epow = 1.00f;
    int _k1min, _k2min, _k3min;
    int _esmooth = 1;

  public:
    sdw_obj();

    void init(int k, double smin, double smax, sampling_obj s1);
};

#endif
