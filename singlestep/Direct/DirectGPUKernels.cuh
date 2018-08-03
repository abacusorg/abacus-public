//* \file This file contains the pairwise direct kernel.
 * 
 * We have various choices for softening, selected by pre-processor #ifdef
 * 
 *  This is set up to be input twice, once with COMPUTE_FOF_DENSITY_SET
 *  and once without.  With it, the function name is 
 *  direct_density() and it takes the extra &aw argument.
 *  Without, it's called direct().
 */

#ifdef COMPUTE_FOF_DENSITY_SET
#define DirectKernelName direct_density
#else
#define DirectKernelName direct
#endif

#ifdef DIRECTCUBICSPLINE
template<int R>
__device__ __inline__ void DirectKernelName(
        FLOAT   &sinkx,     FLOAT &sinky,    FLOAT &sinkz,
        FLOAT   &sourcex,   FLOAT &sourcey,  FLOAT &sourcez,
        FLOAT &ax, FLOAT &ay, FLOAT &az, 
	#ifdef COMPUTE_FOF_DENSITY_SET
	FLOAT &aw,
	#endif
        FLOAT &eps_inv, FLOAT &b2){

    FLOAT drx, dry, drz, rinv,u,spf;
    drx = sourcex - sinkx;
    dry = sourcey - sinky;
    drz = sourcez - sinkz;

    rinv = drx*drx + dry*dry + drz*drz;
    #ifdef COMPUTE_FOF_DENSITY_SET
    if (rinv<b2) { aw = aw+b2-rinv; }
    #endif
    rinv = RSQRT( rinv );
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
__device__ __inline__ void DirectKernelName(
        FLOAT   &sinkx,     FLOAT &sinky,    FLOAT &sinkz,
        FLOAT   &sourcex,   FLOAT &sourcey,  FLOAT &sourcez,
        FLOAT &ax, FLOAT &ay, FLOAT &az, 
	#ifdef COMPUTE_FOF_DENSITY_SET
	FLOAT &aw,
	#endif
        FLOAT &inv_eps2, FLOAT &b2){

    FLOAT dx, dy, dz, f, dr2;
    dx = sourcex - sinkx;
    dy = sourcey - sinky;
    dz = sourcez - sinkz;
    
    dr2 = (dx*dx + dy*dy + dz*dz);
    #ifdef COMPUTE_FOF_DENSITY_SET
    if (dr2<b2) { aw = aw+ b2-dr2; }
    #endif
    dr2 = dr2*inv_eps2+(FLOAT)1e-32;
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
__device__ __inline__ void DirectKernelName(
        FLOAT   &sinkx,     FLOAT &sinky,    FLOAT &sinkz,
        FLOAT   &sourcex,   FLOAT &sourcey,  FLOAT &sourcez,
        FLOAT &ax, FLOAT &ay, FLOAT &az, 
	#ifdef COMPUTE_FOF_DENSITY_SET
	FLOAT &aw,
	#endif
        FLOAT &eps3, FLOAT &b2){

    FLOAT drx, dry, drz, r;
    drx = sourcex - sinkx;
    dry = sourcey - sinky;
    drz = sourcez - sinkz;

    r = drx*drx + dry*dry + drz*drz;
    #ifdef COMPUTE_FOF_DENSITY_SET
    if (r<b2) { aw = aw+ b2-r; }
    #endif
    r += TAU2;
    r *= r*RSQRT(r);  //r^3

    r+=eps3;
    r=(FLOAT)(1.0)/r;

    ax -= r * drx;
    ay -= r * dry;
    az -= r * drz;

}
#elif defined DIRECTPLUMMER
// Default (plummer kernel)
template <int R>
__device__ __inline__ void DirectKernelName(
        FLOAT   &sinkx,     FLOAT &sinky,    FLOAT &sinkz,
        FLOAT   &sourcex,   FLOAT &sourcey,  FLOAT &sourcez,
        FLOAT &ax, FLOAT &ay, FLOAT &az, 
	#ifdef COMPUTE_FOF_DENSITY_SET
	FLOAT &aw,
	#endif
        FLOAT &eps2, FLOAT &b2){

    FLOAT drx, dry, drz, r;
    drx = sourcex - sinkx;
    dry = sourcey - sinky;
    drz = sourcez - sinkz;

    r = drx*drx + dry*dry + drz*drz;
    #ifdef COMPUTE_FOF_DENSITY_SET
    if (r<b2) { aw = aw+ b2-r; }
    #endif
    r += eps2;

    r = RSQRT(r);
    r *=r*r;//*source.w;// * r * r * r;

    ax -= r * drx;
    ay -= r * dry;
    az -= r * drz;

}

#else
    #error "Must have defined a DIRECT* technique!"
#endif

#undef DirectKernelName
