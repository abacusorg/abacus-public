#pragma once

#include "avxdirect.h"

class AVXDirectFloatNR {
public:
    AVXDirectFloatNR(int maxsrc);
    ~AVXDirectFloatNR();
    void compute(int nsrc, ThreeVector<float> *psrc, int nsink, ThreeVector<float> *psink, float3 &delta, float eps,
                 ThreeVector<float> *pacc);

private:
  void storejpdata(int nsrc, ThreeVector<float> *psrc);

  jpstruct<float> *jpdata;
  
#ifdef DIRECTCUBICSPLINE
  ipstruct<float,8> *ipdata;
  ipstruct<float,8> *deltas;
  apstruct<float,8> *accdata;

  void KernelAccPot(ipstruct<float,8> *ipdata, jpstruct<float> *jpdata, int nsrc,
                      ipstruct<float,8> *deltas, float eps2, apstruct<float,8> *accdata);
#else
  ipstruct<float,4> *ipdata;
  ipstruct<float,4> *deltas;
  apstruct<float,4> *accdata;
  
  void KernelAccPot(ipstruct<float,4> *ipdata, jpstruct<float> *jpdata, int nsrc,
                      ipstruct<float,4> *deltas, float eps2, apstruct<float,4> *accdata);
#endif

  int maxsrc;
};


