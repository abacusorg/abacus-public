#ifdef DOUBLEPRECISION
#define RSQRT 1.0/sqrt
#define MIN fmin
#else
#define RSQRT 1.0f/sqrt
#define MIN fminf
#endif

#ifdef DIRECTCUBICSPLINE
template<int R>
inline void directkernel(
        FLOAT3   &sink,
        FLOAT3   &source,
        FLOAT3   &a,
        FLOAT &eps_inv){

    FLOAT drx, dry, drz, rinv,u,spf;
    drx = source.x - sink.x;
    dry = source.y - sink.y;
    drz = source.z - sink.z;

    rinv = RSQRT( drx*drx + dry*dry + drz*drz);
    if(R <= 1){
        u = eps_inv/rinv;
        //u = isfinite(u)? u : 0;  // this can never happen

        // Cubic spline from Hernquist & Katz, ApjS Vol 70 p 419, 1989
        if(u >= (FLOAT)(2.0))
            spf = rinv*rinv*rinv;
        else if (u >= (FLOAT)(1.0)){
            spf = rinv*rinv*rinv * ( (FLOAT)(-1.0/15.0)              +
                    (FLOAT)(8.0/3.0)*   u*u*u       -
                    (FLOAT)(3.0) *      u*u*u*u     +
                    (FLOAT)(6.0/5.0)*   u*u*u*u*u   -
                    (FLOAT)(1.0/6.0)*   u*u*u*u*u*u);
        }else{
            spf = eps_inv*eps_inv*eps_inv*(
                    (FLOAT)(4.0/3.0)            -
                    (FLOAT)(6.0/5.0)*   u*u     +
                    (FLOAT)(1.0/2.0)*   u*u*u);
        }
    }
    else spf = rinv*rinv*rinv;

    a.x -= spf*drx;
    a.y -= spf*dry;
    a.z -= spf*drz;
}

#elif defined DIRECTSINGLESPLINE
// Warning: the following does not return 1/r^2; it still needs to be multiplied by 1/eps^3
template<int R>
inline void directkernel(
        FLOAT3   &sink,
        FLOAT3   &source,
        FLOAT3   &a,
        FLOAT &inv_eps2){

    FLOAT dx, dy, dz, f, dr2;
    dx = source.x - sink.x;
    dy = source.y - sink.y;
    dz = source.z - sink.z;
    
    dr2 = (dx*dx + dy*dy + dz*dz)*inv_eps2 + (FLOAT)1e-32;
    f = RSQRT(dr2);
    
    if(R <= 1){
        f = MIN(f*f*f, (-15*f + 6)*dr2 + 10);
    }
    else{
       f = f*f*f;
    }

    a.x -= f*dx;
    a.y -= f*dy;
    a.z -= f*dz;
}

#elif defined DIRECTCUBICPLUMMER
#define TAU2 ((FLOAT)(1e-16))

template <int R>
inline void directkernel(
        FLOAT3   &sink,
        FLOAT3   &source,
        FLOAT3   &a,
        FLOAT &eps3){

    FLOAT drx, dry, drz, r;
    drx = source.x - sink.x;
    dry = source.y - sink.y;
    drz = source.z - sink.z;

    r = drx*drx + dry*dry + drz*drz + TAU2;
    r *= r*RSQRT(r);  //r^3

    r+=eps3;
    r=(FLOAT)(1.0)/r;

    a.x -= r * drx;
    a.y -= r * dry;
    a.z -= r * drz;
}

#else
template <int R>
inline void directkernel(
        FLOAT3   &sink,
        FLOAT3   &source,
        FLOAT3   &a,
        FLOAT &eps2){

    FLOAT drx, dry, drz, r;
    drx = source.x - sink.x;
    dry = source.y - sink.y;
    drz = source.z - sink.z;

    r = RSQRT( drx*drx + dry*dry + drz*drz  + eps2);
    r *=r*r;

    a.x -= r * drx;
    a.y -= r * dry;
    a.z -= r * drz;

}
#endif
