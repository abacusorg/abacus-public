#ifndef TAYLOR
#define TAYLOR

#ifdef AVXMULTIPOLES
#include "externaltaylor.h"
#endif

#include "basemultipoles.h"


class Taylor : public basemultipoles {
public:
    Taylor(int order);
    ~Taylor(void);

    void AnalyticEvaluateTaylor(double *CT, FLOAT3 expansioncenter, int np, FLOAT3 *ps, FLOAT3 *a);
    void      ASMEvaluateTaylor(double *CT, FLOAT3 expansioncenter, int np, FLOAT3 *ps, FLOAT3 *a);


#ifdef AVXMULTIPOLES
    d4 *cx[128];
    d4 *cy[128];
    d4 *cz[128];

    d4 *ax[128];
    d4 *ay[128];
    d4 *az[128];

    d4 *px[128];
    d4 *py[128];
    d4 *pz[128];

    d4 *Qx[128];
    d4 *Qy[128];
    d4 *Qz[128];
#endif

};

#endif
