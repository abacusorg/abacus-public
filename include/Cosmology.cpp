// Copyright 2012-2025 The Abacus Developers
// SPDX-License-Identifier: GPL-3.0-or-later

#include "Cosmology.h"

#include <fmt/base.h>

Cosmology::Cosmology(double a_initial, MyCosmology& _C) {
    /// This sets epoch current and next to a_initial, 
    /// as well as initializing the 'early' epoch.
    C = _C;
    n = (sqrt(25.0-24.0*C.Omega_smooth/C.Omega_m)-1.0)/4.0;

    InitializeEarly(0.000001);
    BuildEpoch(early, today, 1.0);
    BuildEpoch(today, today, 1.0);

    BuildEpoch(early, current, a_initial);
    BuildEpoch(current, current, a_initial);
    CopyEpoch(next, current);

    return;
}

void Cosmology::CopyEpoch(Epoch& dest, const Epoch& source) {
    /// standard "copy(dest, source)"
    Epoch tmp = source; // need to make tmp copy to allow source<-source copies
    dest = tmp;
}

void Cosmology::printepoch(Epoch& E) {
    /// This does a screen dump of an epoch, for debugging purposes.
    FILE *fp = stdout;
    fmt::print(fp, "\n");
    fmt::print(fp, "         Scale factor = {:15.13f}\n", E.a);
    fmt::print(fp, "             Redshift = {:15.11f}\n", E.z);
    fmt::print(fp, "                    H = {:15.9le}\n", E.H);
    fmt::print(fp, "                 Time = {:15.9le}\n", E.t);
    fmt::print(fp, "               H*Time = {:15.13f}\n", E.H*E.t);
    fmt::print(fp, "                 EtaK = {:15.9le}\n", E.etaK);
    fmt::print(fp, "             a*H*EtaK = {:15.13f}\n", E.H*E.etaK*E.a);
    fmt::print(fp, "                 EtaD = {:15.9le}\n", E.etaD);
    fmt::print(fp, "           a^2*H*EtaD = {:15.13f}\n", E.H*E.etaD*E.a*E.a);
    fmt::print(fp, "                    w = {:15.12f}\n", E.w);
    fmt::print(fp, "           OmegaHat_m = {:15.12f}\n", E.OmegaHat_m);
    fmt::print(fp, "           OmegaHat_X = {:15.12f}\n", E.OmegaHat_X);
    fmt::print(fp, "           OmegaHat_K = {:15.12f}\n", E.OmegaHat_K);
    double total = E.OmegaHat_m + E.OmegaHat_X + E.OmegaHat_K;
    fmt::print(fp, "           Omega_m(a) = {:15.12f}\n", E.OmegaHat_m/total);
    fmt::print(fp, "           Omega_X(a) = {:15.12f}\n", E.OmegaHat_X/total);
    fmt::print(fp, "           Omega_K(a) = {:15.12f}\n", E.OmegaHat_K/total);
    fmt::print(fp, "               growth = {:15.13f}\n", E.growth);
    fmt::print(fp, "           growth/a^n = {:15.13f}\n", E.growth/pow(E.a,n));
    fmt::print(fp, "                    n = {:15.13f}\n", n);
    fmt::print(fp, "             f_growth = {:15.13f}\n", E.f_growth);
    fmt::print(fp, "f_growth/Omega_m^0.55 = {:15.13f}\n", E.f_growth/pow(E.OmegaHat_m/total,0.55));
}

void Cosmology::InitializeEarly(double a_early) {
    /// This starts with the EdS solution at very early times
    /// TODO: This version is ok for things near LCDM, but will be inaccurate
    /// for models in which dark energy isn't negligible at redshift million.
    double vars[5];

    vars[0] = 2.0/3.0;       // = Ht
    vars[1] = 2.0;           // = a*H*etaK
    vars[2] = -2.0;          // = a*a*H*etaD
    vars[3] = 1.0;           // = D/a^n, where n = (sqrt(25-24*Omega_smooth/Omega_m)-1)/4
    vars[4] = 0.0;           // = d(D/a^n)/dlna

    Unpack(early,log(a_early),vars);
    return;
}

void Cosmology::BuildEpoch(const Epoch& start, Epoch& end, double a) {
    /// This will construct epoch 'end' at scale factor 'a', starting
    /// the integration from 'start'.  The epoch 'start' is not altered
    /// by this routine.  It is legal to use end==start, in which case
    /// the epoch is updated to the new one.
    double  vars0[5],  vars1[5],  vars2[5],  vars3[5], varsfinal[5];
    double dvars0[5], dvars1[5], dvars2[5], dvars3[5];
    double step, delta_lna, lna0, lna1, lna2, lna3;
    int nstep, j, k;

    /// We are doing all integrations with small steps in ln(a) and
    /// with scaled variables that hold EdS constant.
    const double DLNA = 0.0025;
    const int MIN_STEPS = 1;

    delta_lna = log(a/start.a);
    step = delta_lna/MIN_STEPS;
    nstep = MIN_STEPS;
    if (fabs(step)>DLNA) {
        nstep = floor(fabs(delta_lna)/DLNA)+1;
        step = delta_lna/nstep;
    }

    CopyEpoch(end, start);

    for (j=0;j<nstep;j++) {
        Pack(end, &lna0, vars0);

	// Fourth-order Runga-Kutta
        Deriv(lna0, vars0, dvars0);
        for (k=0; k<5; k++) vars1[k] = vars0[k] + dvars0[k]*0.5*step; 
        lna1 = lna0 + 0.5*step;
        Deriv(lna1, vars1, dvars1);
        for (k=0; k<5; k++) vars2[k] = vars0[k] + dvars1[k]*0.5*step; 
        lna2 = lna0 + 0.5*step;
        Deriv(lna2, vars2, dvars2);
        for (k=0; k<5; k++) vars3[k] = vars0[k] + dvars2[k]*step; 
        lna3 = lna0 + step;
        Deriv(lna3, vars3, dvars3);
        for (k=0; k<5; k++) 
            varsfinal[k] = vars0[k]+
                step*(dvars0[k]/6 + dvars1[k]/3 + dvars2[k]/3 + dvars3[k]/6);

        Unpack(end, lna3, varsfinal);
    }
}


void Cosmology::Pack(Epoch& epoch, double *lna, double *vars) {
    /// This converts normal cosmological quantities into the EdS-scaled quantities 
    /// that we want to use in our integrals, for stability reasons.

    *lna = log(epoch.a);
    vars[0] = (epoch.H)*(epoch.t);
    vars[1] = (epoch.a)*(epoch.H)*(epoch.etaK);
    vars[2] = (epoch.a)*(epoch.a)*(epoch.H)*(epoch.etaD);
    vars[3] = (epoch.growth)/pow(epoch.a,n);
    vars[4] = vars[3]*(epoch.f_growth-n);

    return;
}

void Cosmology::Unpack(Epoch& epoch, double lna, double *vars) {
    /// This converts EdS-scaled quantities back into the normal cosmological variables.

    epoch.a = exp(lna);
    double a2 = epoch.a*epoch.a;
    double a3 = epoch.a*a2;

    epoch.z = 1/epoch.a - 1;
    epoch.OmegaHat_X = C.OmegaX(epoch.a);
    epoch.OmegaHat_m = C.Omega_m/a3;
    epoch.OmegaHat_K = C.Omega_K/a2;
    epoch.w = C.w(epoch.a);
    epoch.H = C.H0*sqrt(epoch.OmegaHat_m + epoch.OmegaHat_K + epoch.OmegaHat_X);

    epoch.t        = vars[0]/epoch.H;
    epoch.etaK     = vars[1]/epoch.H/epoch.a;
    epoch.etaD     = vars[2]/epoch.H/a2;
    epoch.growth   = vars[3]*pow(epoch.a,n);
    epoch.f_growth = vars[4]/vars[3] + n;

    return;
}

void Cosmology::Deriv(double lna, double *vars, double *dvars) {
    /// Compute the derivative dvars/dx (where x = ln(a)) at the given time
    /// TODO: this is recomputing things that are done in Unpack()
    /// and taking Omega values from Cosmology rather than Epoch.
    /// For now, this is fine, but it is a bug risk in the future.
    /// Omega_m is used for the H(z) expansion history, but we decrease it 
    /// by Omega_smooth when sourcing the growth function.
    double a, a2, a3, OmX, E2, dlnH2dx;
    a = exp(lna);
    a2 = a*a;
    a3 = a*a2;
    OmX = C.OmegaX(a);
    E2 = C.Omega_m/a3 + C.Omega_K/a2 + OmX;
    dlnH2dx = ( -3*C.Omega_m/a3 - 2*C.Omega_K/a2 + C.dOmegaXdlna(a) ) / E2;

    dvars[0] = 1 + 0.5*vars[0]*dlnH2dx;
    dvars[1] = 1 + 0.5*vars[1]*dlnH2dx + vars[1];
    dvars[2] = 1 + 0.5*vars[2]*dlnH2dx + vars[2]*2;
    dvars[3] = vars[4];
    dvars[4] = -vars[4]*(2.0+2.0*n + 0.5*dlnH2dx) - vars[3]*(n*n + 2.0*n + n*0.5*dlnH2dx - 1.5*(C.Omega_m-C.Omega_smooth)/a3/E2);

    return;
}

double Cosmology::KickFactor(double a, double da) {
    /// This reports the KickFactor required in the time stepper
    //    CopyEpoch(next, early);  // this increases execution time w/o increasing accuracy
    Epoch tmp;
    BuildEpoch(current,tmp,a);
    double kfl = tmp.etaK;
    BuildEpoch(current,tmp,a+da);
    return tmp.etaK - kfl;
}

double Cosmology::DriftFactor(double a, double da) {
    /// This reports the DriftFactor required in the time stepper
    //    CopyEpoch(next, early);  // this increases execution time w/o increasing accuracy
    Epoch tmp;
    BuildEpoch(current,tmp,a);
    double dfl = tmp.etaD;
    BuildEpoch(current,tmp,a+da);
    return tmp.etaD - dfl;
}

double Cosmology::a2t(double a) {
    /// Return the time for a given scale factor
    //    CopyEpoch(next, early);  // this increases execution time w/o increasing accuracy
    Epoch tmp;
    BuildEpoch(current,tmp,a);
    return tmp.t;
}

double Cosmology::H(double a) {
    /// Return the Hubble parameter for a given scale factor
    Epoch tmp;
    BuildEpoch(current,tmp,a);
    return tmp.H;
}

#define sign(x) ((x)<0?-1:1)
double Cosmology::t2a(double t) {
    /// Return the scale factor for a given time, which requires a search
    /// We do the search by Newton-Raphson.
    CopyEpoch(search,next);
    double a;
    int n = 0;

    a = 0.5;

    double f, fp, da;
    double dela = 1e-6;
    n = 0;
    do {
        BuildEpoch(search,search,a);
        f = search.t - t;
	// STDLOG(1,"t2a(): a={:12.8f}   t={:12.8f}\n", a, search.t);
        BuildEpoch(search,search,a-dela);
        fp = search.t;
        BuildEpoch(search,search,a+dela);
        fp = (search.t-fp)/(2*dela);
        da = - f/fp;
        a += da;
        n++;
	// STDLOG(1,"t2a(): a={:12.8f}    da={:12.4e}\n", a, da);
        //    } while(fabs(f)/t>1.0e-15 || fabs(da)>1e-13);
    } while(fabs(da)>1e-14);
    BuildEpoch(next,search,a);
    return a;

}
#undef sign
