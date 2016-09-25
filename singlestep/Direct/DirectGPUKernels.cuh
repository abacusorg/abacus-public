#ifdef DIRECTCUBICSPLINE
template<int R>
__device__ __inline__ void direct(
        FLOAT   &sinkx,     FLOAT &sinky,    FLOAT &sinkz,
        FLOAT   &sourcex,   FLOAT &sourcey,  FLOAT &sourcez,
        FLOAT &ax, FLOAT &ay, FLOAT &az,
        FLOAT &eps_inv){

    FLOAT drx, dry, drz, rinv,u,spf;
    drx = sourcex - sinkx;
    dry = sourcey - sinky;
    drz = sourcez - sinkz;

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

    ax -= spf*drx;
    ay -= spf*dry;
    az -= spf*drz;
}

#elif defined DIRECTSINGLESPLINE
template<int R>
__device__ __inline__ void direct(
        FLOAT   &sinkx,     FLOAT &sinky,    FLOAT &sinkz,
        FLOAT   &sourcex,   FLOAT &sourcey,  FLOAT &sourcez,
        FLOAT &ax, FLOAT &ay, FLOAT &az,
        FLOAT &inv_eps2){

    FLOAT dx, dy, dz, f, dr2;
    dx = sourcex - sinkx;
    dy = sourcey - sinky;
    dz = sourcez - sinkz;
    
    dr2 = (dx*dx + dy*dy + dz*dz)*inv_eps2 + (FLOAT)1e-32;
    f = RSQRT(dr2);
    
    if(R <= 1){
        f = MIN(f*f*f, (-15*f + 6)*dr2 + 10);
    }
    else{
       f = f*f*f;
    }

    ax -= f*dx;
    ay -= f*dy;
    az -= f*dz;
}

#elif defined DIRECTCUBICPLUMMER

#define TAU2 ((FLOAT)(1e-16))

template <int R>
__device__ __inline__ void direct(
        FLOAT   &sinkx,     FLOAT &sinky,    FLOAT &sinkz,
        FLOAT   &sourcex,   FLOAT &sourcey,  FLOAT &sourcez,
        FLOAT &ax, FLOAT &ay, FLOAT &az,
        FLOAT &eps3){

    FLOAT drx, dry, drz, r;
    drx = sourcex - sinkx;
    dry = sourcey - sinky;
    drz = sourcez - sinkz;

    r = drx*drx + dry*dry + drz*drz + TAU2;
    r *= r*RSQRT(r);  //r^3

    r+=eps3;
    r=(FLOAT)(1.0)/r;

    ax -= r * drx;
    ay -= r * dry;
    az -= r * drz;

}
#else
// Default (plummer kernel)
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
