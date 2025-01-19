#ifndef __COSMOLOGY_H
#define __COSMOLOGY_H

#include <cmath>
#include <cstdio>

class MyCosmology {
// The cosmology:
// Must include Omega_m, Omega_smooth, Omega_K, Omega_DE, H0
// and the functions OmegaX(a) and dOmegaXdlna(a)
// Omega_smooth is the portion of Omega_m that is being forced to be homogeneous.
// That is, Omega_m determines H(z), but the linear perturbations are carried by 
// a component of density Omega_m-Omega_smooth 
public:
    double Omega_m, Omega_smooth, Omega_K, Omega_DE, H0;

    double w0, wa;
    double w(double a) {
        return w0+wa*(1-a);
    }
    double OmegaX(double a) {
        return Omega_DE*pow(a,-3.0*(1.0+w0+wa))*exp(3.0*wa*(a-1.0));
    }
    double dOmegaXdlna(double a) {
        return -3*(1+w(a))*OmegaX(a);
    }
};

/// This structure contains the cosmological variables at any single epoch.
typedef struct {
        double a, z, t, etaK, etaD, growth, f_growth;
        double w, H, OmegaHat_X, OmegaHat_m, OmegaHat_K;
} Epoch;

class Cosmology {
/// This class carries out the cosmological integrations.
/// One instantiates with a given cosmology and initial redshift.
/// One then can integrate to any chosen epoch.
/// There are functions for the usual time-step applications.
public:

    Cosmology(double a_initial, MyCosmology& cos);
    void InitializeEarly(double a_early);
    //void SetupTimeStep(double next_a);
    void BuildEpoch(const Epoch& start, Epoch& end, double a);  
    double DriftFactor(double a, double da);
    double KickFactor(double a, double da);
    double a2t(double a);
    double t2a(double t);
    double H(double a);
    void printepoch(Epoch& E);

    Epoch early, today, current, next, search;
    /// In normal usage: 
    ///    today will be a=1
    ///    early will be some very early epoch, e.g., z=1e6
    ///    current will be the previous time step
    ///    next will be the next time step
    ///    search holds the result of bisection searches for a(t)
    MyCosmology C;

    double n;   // This is the EdS growth function exponent in the presence of smooth components.
private:
    void Pack(Epoch& epoch, double *lna, double *vars);
    void Unpack(Epoch& epoch, double lna, double *vars);
    void Deriv(double lna, double *vars, double *dvars);
    void CopyEpoch(Epoch& dest, const Epoch& source);
};

#endif
