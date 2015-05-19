class AVXDirectFloatNR {
public:
    AVXDirectFloatNR(int maxsrc);
    ~AVXDirectFloatNR();
    void compute(int nsrc, ThreeVector<float> *psrc, int nsink, ThreeVector<float> *psink, float3 &delta, float eps,
                 ThreeVector<float> *pacc);

private:
  void storejpdata(int nsrc, ThreeVector<float> *psrc);

  void KernelAccPot(ipstruct<float> *ipdata, jpstruct<float> *jpdata, int nsrc,
                      ipstruct<float> *deltas, float eps2, apstruct<float> *accdata);

  int maxsrc;
  jpstruct<float> *jpdata;
  ipstruct<float> *ipdata;
  ipstruct<float> *deltas;
  apstruct<float> *accdata;
};


