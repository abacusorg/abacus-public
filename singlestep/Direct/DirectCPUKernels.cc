#ifdef DOUBLEPRECISION
#define RSQRT 1.0/sqrt
#else
#define RSQRT 1.0f/sqrt
#endif

#ifdef DIRECTSPLINE
template<int R>
inline void directkernel(
        FLOAT3   &sink,
        FLOAT3   &source,
        FLOAT3   &a,
        FLOAT &eps2){

    FLOAT drx, dry, drz, r,u,spf;
    drx = source.x - sink.x;
    dry = source.y - sink.y;
    drz = source.z - sink.z;

    r = RSQRT( drx*drx + dry*dry + drz*drz);
    FLOAT eps_inv = RSQRT(eps2);
    u = eps_inv/r;
    u = isfinite(u)? u : 0;

    r *=r*r;

    // Cubic spline from Hernquist & Katz, ApjS Vol 70 p 419, 1989?
    if(u >= (FLOAT)(2.0))
        spf = r;
    else if (u >= (FLOAT)(1.0)){
        spf = r * ( (FLOAT)(-1.0/15.0)              +
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
    a.x -= spf*drx;
    a.y -= spf*dry;
    a.z -= spf*drz;
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
