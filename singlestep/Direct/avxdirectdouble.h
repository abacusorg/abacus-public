class AVXDirectDouble {
public:
    AVXDirectDouble(int maxsrc);
    ~AVXDirectDouble();
    void compute(int nsrc, ThreeVector<double> *psrc, int nsink, ThreeVector<double> *psink, double3 &delta, double eps,
                 ThreeVector<double> *pacc);

private:
    void storejpdata(int nsrc, ThreeVector<double> *psrc);

  void KernelAccPot(ipstruct<double> *ipdata, jpstruct<double> *jpdata, int nsrc,
                      ipstruct<double> *deltas, double eps2, apstruct<double> *accdata);

    int maxsrc;
    jpstruct<double> *jpdata;
    ipstruct<double> *ipdata;
    ipstruct<double> *deltas;
    apstruct<double> *accdata;
};

