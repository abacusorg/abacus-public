#ifdef DIRECTSPLINE
template<int R>
__device__ __inline__ void direct(
        FLOAT   &sinkx,     FLOAT &sinky,    FLOAT &sinkz,
        FLOAT   &sourcex,   FLOAT &sourcey,  FLOAT &sourcez,
        FLOAT &ax, FLOAT &ay, FLOAT &az,
        FLOAT &eps2){

    FLOAT drx, dry, drz, r,u,spf;
    drx = sourcex - sinkx;
    dry = sourcey - sinky;
    drz = sourcez - sinkz;

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

    ax -= spf*drx;
    ay -= spf*dry;
    az -= spf*drz;
}

#else
template <int R>
__device__ __inline__ void direct(
        FLOAT   &sinkx,     FLOAT &sinky,    FLOAT &sinkz,
        FLOAT   &sourcex,   FLOAT &sourcey,  FLOAT &sourcez,
        FLOAT &ax, FLOAT &ay, FLOAT &az,
        FLOAT &eps2){

    FLOAT drx, dry, drz, r;
    drx = sourcex - sinkx;
    dry = sourcey - sinky;
    drz = sourcez - sinkz;

    r = RSQRT( drx*drx + dry*dry + drz*drz  + eps2);
    r *=r*r;//*source.w;// * r * r * r;

    ax -= r * drx;
    ay -= r * dry;
    az -= r * drz;

}
#endif

