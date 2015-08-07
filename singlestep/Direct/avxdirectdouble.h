class AVXDirectDouble {
public:
    AVXDirectDouble(int maxsrc);
    ~AVXDirectDouble();
    void compute(int nsrc, ThreeVector<double> *psrc, int nsink, ThreeVector<double> *psink, double3 &delta, double eps,
                 ThreeVector<double> *pacc);

private:
    void storejpdata(int nsrc, ThreeVector<double> *psrc);

  void KernelAccPot(ipstruct<double,4> *ipdata, jpstruct<double> *jpdata, int nsrc,
                      ipstruct<double,4> *deltas, double eps2, apstruct<double,4> *accdata);

    int maxsrc;
    jpstruct<double> *jpdata;
    ipstruct<double,4> *ipdata;
    ipstruct<double,4> *deltas;
    apstruct<double,4> *accdata;
};

