#!/usr/bin/env python3
'''
This program generates assembly code for the Taylors and Multipoles on
IBM POWER architectures using AltiVec/VSX intrinsics.
You shouldn't have to run this program manually; it will
be invoked from the Makefile.  But to see the manual usage, run:

$ ./generateCartesianVSX.py --help

The compiler can't figure out how to unroll some of the nested
loops, so the purpose of this program is to "manually" unroll them.
'''

from multipoles_metacode_utils \
    import emit_dispatch_function, Writer, default_metacode_argparser, cmapper


def emit_VSX_Multipoles(orders, fn='CM_VSX.cpp'):
    w = Writer(fn)

    w('''
        #include "threevector.hh"
        #include "header.cpp"

        #ifdef VSXMULTIPOLES
        #include <altivec.h>

        #include "assert.h"

        #define VSX_NVEC_DOUBLE 2
        typedef vector double VSX_DOUBLES;

        template <int Order>
        void MultipoleVSXKernel(FLOAT3 *particles, int n, double3 center, double *CM);
        ''')

    for order in orders:
        cml = (order+1)*(order+2)*(order+3)//6

        w(f'''
            template <>
            void MultipoleVSXKernel<{order}>(FLOAT3 *particles, int n, double3 center, double *CM){{''')
        w.indent()

        w(f'''
            VSX_DOUBLES cx = vec_splats(center.x);
            VSX_DOUBLES cy = vec_splats(center.y);
            VSX_DOUBLES cz = vec_splats(center.z);

            VSX_DOUBLES CMvec[{cml}];
            for(int k = 0; k < {cml}; k++)
                CMvec[k] = vec_splats(0.);

            int n_left = n % VSX_NVEC_DOUBLE;
            int n_aligned = n - n_left;
            ''')

        w('for(int i = 0; i < n_aligned; i += VSX_NVEC_DOUBLE){')
        w.indent()

        w('''
            VSX_DOUBLES px;
            VSX_DOUBLES py;
            VSX_DOUBLES pz;

            for(int j = 0; j < VSX_NVEC_DOUBLE; j++){
                px[j] = particles[i+j].x;
                py[j] = particles[i+j].y;
                pz[j] = particles[i+j].z;
            }

            VSX_DOUBLES fi, fij, fijk;
            fi = vec_splats(1.);
            ''')

        w('''
            VSX_DOUBLES deltax, deltay, deltaz;
            deltax = px - cx;
            deltay = py - cy;
            deltaz = pz - cz;
            ''')

        i = 0
        nflop = 3  # deltas
        for a in range(order+1):
            w('fij = fi;')
            for b in range(order-a+1):
                w('fijk = fij;')
                for c in range(order-a-b+1):
                    w(f'CMvec[{i}] += fijk;')
                    i += 1
                    nflop += 1
                    if c < order-a-b:
                        w('fijk *= deltaz;')
                        nflop += 1
                if b < order-a:
                    w('fij *= deltay;')
                    nflop += 1
            if a < order:
                w('\nfi *= deltax;')
                nflop += 1

        w.dedent()
        w('}')  # particle loop

        print(f'nflop: {nflop}; cml: {cml}')

        w('if (n_left != 0){')
        w.indent()
        w('''
            // plain C, only one element!
            FLOAT3 p = particles[n_aligned];
            
            double fi, fij, fijk;
            fi = 1.0;

            double deltax, deltay, deltaz;
            deltax = p.x - center.x;
            deltay = p.y - center.y;
            deltaz = p.z - center.z;
            ''')

        cml = (order+1)*(order+2)*(order+3)//6

        i = 0
        for a in range(order+1):
            w('fij = fi;')
            for b in range(order-a+1):
                w('fijk = fij;')
                for c in range(order-a-b+1):
                    w(f'CM[{i}] += fijk;')
                    i += 1
                    nflop += 1
                    if c < order-a-b:
                        w('fijk *= deltaz;')
                        nflop += 1
                if b < order-a:
                    w('fij *= deltay;')
                    nflop += 1
            if a < order:
                w('\nfi *= deltax;')
                nflop += 1
        w.dedent()
        w('}')  # remainder block

        w(f'''
            for(int k = 0; k < {cml}; k++)
                for(int j = 0; j < VSX_NVEC_DOUBLE; j++)
                    CM[k] += CMvec[k][j];
            ''')

        w.dedent()
        w('}\n')  # Kernel

    emit_dispatch_function(w, 'MultipoleVSXKernel(FLOAT3 *particles, int n, double3 center, double *CM)', orders)

    w('''
        #undef VSX_NVEC_DOUBLE
        #endif''')  # VSXMULTIPOLES


def emit_VSX_Multipoles_interleaved(orders, fn='CM_VSX.cpp'):
    w = Writer(fn)

    w('''
        #include "threevector.hh"
        #include "header.cpp"

        #ifdef VSXMULTIPOLES
        #include <altivec.h>

        #include "assert.h"

        #define VSX_NVEC_DOUBLE 2
        typedef vector double VSX_DOUBLES;

        template <int Order>
        void MultipoleVSXKernel(FLOAT3 *particles, int n, double3 center, double *CM);
        ''')

    for order in orders:
        cml = (order+1)*(order+2)*(order+3)//6

        w(f'''
            template <>
            void MultipoleVSXKernel<{order}>(FLOAT3 *particles, int n, double3 center, double *CM){{''')
        w.indent()

        w(f'''
            VSX_DOUBLES cx = vec_splats(center.x);
            VSX_DOUBLES cy = vec_splats(center.y);
            VSX_DOUBLES cz = vec_splats(center.z);

            VSX_DOUBLES CMvec[{cml}], CMvec2[{cml}];
            for(int k = 0; k < {cml}; k++){{
                CMvec[k] = vec_splats(0.);
                CMvec2[k] = vec_splats(0.);
            }}

            int n_left = n % (2*VSX_NVEC_DOUBLE);
            int n_aligned = n - n_left;
            ''')

        w('for(int i = 0; i < n_aligned; i += 2*VSX_NVEC_DOUBLE){')
        w.indent()

        w('''
            VSX_DOUBLES px, py, pz;
            VSX_DOUBLES px2, py2, pz2;

            for(int j = 0; j < VSX_NVEC_DOUBLE; j++){
                px[j] = particles[i+j].x;
                px2[j] = particles[i+j+VSX_NVEC_DOUBLE].x;
                py[j] = particles[i+j].y;
                py2[j] = particles[i+j+VSX_NVEC_DOUBLE].y;
                pz[j] = particles[i+j].z;
                pz2[j] = particles[i+j+VSX_NVEC_DOUBLE].z;
            }

            VSX_DOUBLES fi, fij, fijk;
            VSX_DOUBLES fi2, fij2, fijk2;
            fi = vec_splats(1.);
            fi2 = vec_splats(1.);
            ''')

        w('''
            VSX_DOUBLES deltax, deltay, deltaz;
            VSX_DOUBLES deltax2, deltay2, deltaz2;
            deltax = px - cx;
            deltax2 = px2 - cx;
            deltay = py - cy;
            deltay2 = py2 - cy;
            deltaz = pz - cz;
            deltaz2 = pz2 - cz;
            ''')

        i = 0
        nflop = 3  # deltas
        for a in range(order+1):
            w('fij = fi;')
            w('fij2 = fi2;')
            for b in range(order-a+1):
                w('fijk = fij;')
                w('fijk2 = fij2;')
                for c in range(order-a-b+1):
                    w(f'CMvec[{i}] += fijk;')
                    w(f'CMvec2[{i}] += fijk2;')
                    i += 1
                    nflop += 1
                    if c < order-a-b:
                        w('fijk *= deltaz;')
                        w('fijk2 *= deltaz2;')
                        nflop += 1
                if b < order-a:
                    w('fij *= deltay;')
                    w('fij2 *= deltay2;')
                    nflop += 1
            if a < order:
                w('\nfi *= deltax;')
                w('\nfi2 *= deltax2;')
                nflop += 1

        w.dedent()
        w('}')  # particle loop

        print(f'nflop: {nflop}; cml: {cml}')

        w('if (n_left != 0){')
        w.indent()
        w('''
            assert(0);  // TODO
            // plain C, only one element!
            FLOAT3 p = particles[n_aligned];
            
            double fi, fij, fijk;
            fi = 1.0;

            double deltax, deltay, deltaz;
            deltax = p.x - center.x;
            deltay = p.y - center.y;
            deltaz = p.z - center.z;
            ''')

        cml = (order+1)*(order+2)*(order+3)//6

        i = 0
        for a in range(order+1):
            w('fij = fi;')
            for b in range(order-a+1):
                w('fijk = fij;')
                for c in range(order-a-b+1):
                    w(f'CM[{i}] += fijk;')
                    i += 1
                    nflop += 1
                    if c < order-a-b:
                        w('fijk *= deltaz;')
                        nflop += 1
                if b < order-a:
                    w('fij *= deltay;')
                    nflop += 1
            if a < order:
                w('\nfi *= deltax;')
                nflop += 1
        w.dedent()
        w('}')  # remainder block

        w(f'''
            for(int k = 0; k < {cml}; k++)
                for(int j = 0; j < VSX_NVEC_DOUBLE; j++){{
                    CM[k] += CMvec[k][j] + CMvec2[k][j];
                }}
            ''')

        w.dedent()
        w('}\n')  # Kernel

    emit_dispatch_function(w, 'MultipoleVSXKernel(FLOAT3 *particles, int n, double3 center, double *CM)', orders)

    w('''
        #undef VSX_NVEC_DOUBLE
        #endif''')  # VSXMULTIPOLES


def emit_VSX_Multipoles_FMA(orders, fn='CM_VSX.cpp', max_zk=None):
    assert max_zk is None, "max_zk implementation not finished"
    if max_zk is None:
        max_zk = max(orders) + 1

    w = Writer(fn)

    w('''
        #include "threevector.hh"
        #include "header.cpp"

        #ifdef VSXMULTIPOLES
        #include <altivec.h>

        #include "assert.h"

        #define VSX_NVEC_DOUBLE 2
        typedef vector double VSX_DOUBLES;

        template <int Order>
        void MultipoleVSXKernel(FLOAT3 *particles, int n, double3 center, double *CM);
        ''')

    for order in orders:
        cml = (order+1)*(order+2)*(order+3)//6
        this_max_zk = min(order+1, max_zk)

        w(f'''
            template <>
            void MultipoleVSXKernel<{order}>(FLOAT3 *particles, int n, double3 center, double *CM){{''')
        w.indent()

        w(f'''
            VSX_DOUBLES cx = vec_splats(center.x);
            VSX_DOUBLES cy = vec_splats(center.y);
            VSX_DOUBLES cz = vec_splats(center.z);

            VSX_DOUBLES CMvec[{cml}];
            for(int k = 0; k < {cml}; k++)
                CMvec[k] = vec_splats(0.);

            int n_left = n % VSX_NVEC_DOUBLE;
            int n_aligned = n - n_left;
            ''')

        w('for(int i = 0; i < n_aligned; i += VSX_NVEC_DOUBLE){')
        w.indent()

        w(f'''
            VSX_DOUBLES px;
            VSX_DOUBLES py;
            VSX_DOUBLES pz;

            for(int j = 0; j < VSX_NVEC_DOUBLE; j++){{
                px[j] = particles[i+j].x;
                py[j] = particles[i+j].y;
                pz[j] = particles[i+j].z;
            }}

            VSX_DOUBLES fi, fij;
            fi = vec_splats(1.);
            
            VSX_DOUBLES deltax, deltay, deltaz;
            deltax = px - cx;
            deltay = py - cy;
            deltaz = pz - cz;

            // precompute powers of deltaz^k
            VSX_DOUBLES zk[{this_max_zk}];
            zk[0] = vec_splats(1.);
            ''')
        for c in range(1,this_max_zk):
            w(f'zk[{c}] = deltaz*zk[{c-1}];')

        i = 0
        nflop = 3 + this_max_zk-1
        for a in range(order+1):
            w('fij = fi;')
            for b in range(order-a+1):
                for c in range(order-a-b+1):
                    w(f'CMvec[{i}] += fij*zk[{c}];')
                    i += 1
                    nflop += 2
                if b < order-a:
                    w('fij *= deltay;')
                    nflop += 1
            if a < order:
                w('\nfi *= deltax;')
                nflop += 1

        w.dedent()
        w('}')  # particle loop

        print(f'nflop: {nflop}; cml: {cml}')

        w('if (n_left != 0){')
        w.indent()
        w(f'''
            // plain C, only one element!
            FLOAT3 p = particles[n_aligned];
            
            double fi, fij;
            fi = 1.0;

            double deltax, deltay, deltaz;
            deltax = p.x - center.x;
            deltay = p.y - center.y;
            deltaz = p.z - center.z;

            double zk[{this_max_zk}];
            zk[0] = 1.;
            ''')
        for c in range(1,this_max_zk):
            w(f'zk[{c}] = deltaz*zk[{c-1}];')

        cml = (order+1)*(order+2)*(order+3)//6

        i = 0
        for a in range(order+1):
            w('fij = fi;')
            for b in range(order-a+1):
                for c in range(order-a-b+1):
                    w(f'CM[{i}] += fij*zk[{c}];')
                    i += 1
                if b < order-a:
                    w('fij *= deltay;')
            if a < order:
                w('\nfi *= deltax;')
        w.dedent()
        w('}')  # remainder block

        w(f'''
            for(int k = 0; k < {cml}; k++)
                for(int j = 0; j < VSX_NVEC_DOUBLE; j++)
                    CM[k] += CMvec[k][j];
            ''')

        w.dedent()
        w('}\n')  # Kernel

    emit_dispatch_function(w, 'MultipoleVSXKernel(FLOAT3 *particles, int n, double3 center, double *CM)', orders)

    w('''
        #undef VSX_NVEC_DOUBLE
        #endif''')  # VSXMULTIPOLES


def emit_VSX_Multipoles_FMA_interleaved(orders, fn='CM_VSX.cpp', max_zk=None):
    assert max_zk is None, "max_zk implementation not finished"
    if max_zk is None:
        max_zk = max(orders) + 1

    w = Writer(fn)

    w('''
        #include "threevector.hh"
        #include "header.cpp"

        #ifdef VSXMULTIPOLES
        #include <altivec.h>

        #include "assert.h"

        #define VSX_NVEC_DOUBLE 2
        typedef vector double VSX_DOUBLES;

        template <int Order>
        void MultipoleVSXKernel(FLOAT3 *particles, int n, double3 center, double *CM);
        ''')

    for order in orders:
        cml = (order+1)*(order+2)*(order+3)//6
        this_max_zk = min(order+1, max_zk)

        w(f'''
            template <>
            void MultipoleVSXKernel<{order}>(FLOAT3 *particles, int n, double3 center, double *CM){{''')
        w.indent()

        w(f'''
            VSX_DOUBLES cx = vec_splats(center.x);
            VSX_DOUBLES cy = vec_splats(center.y);
            VSX_DOUBLES cz = vec_splats(center.z);

            VSX_DOUBLES CMvec[{cml}], CMvec2[{cml}];
            for(int k = 0; k < {cml}; k++){{
                CMvec[k] = vec_splats(0.);
                CMvec2[k] = vec_splats(0.);
            }}

            int n_left = n % (2*VSX_NVEC_DOUBLE);
            int n_aligned = n - n_left;
            ''')

        w('for(int i = 0; i < n_aligned; i += 2*VSX_NVEC_DOUBLE){')
        w.indent()

        w(f'''
            VSX_DOUBLES px, py, pz;
            VSX_DOUBLES px2, py2, pz2;

            for(int j = 0; j < VSX_NVEC_DOUBLE; j++){{
                px[j] = particles[i+j].x;
                px2[j] = particles[i+j+VSX_NVEC_DOUBLE].x;
                py[j] = particles[i+j].y;
                py2[j] = particles[i+j+VSX_NVEC_DOUBLE].y;
                pz[j] = particles[i+j].z;
                pz2[j] = particles[i+j+VSX_NVEC_DOUBLE].z;
            }}

            VSX_DOUBLES fi, fij;
            VSX_DOUBLES fi2, fij2;
            fi = vec_splats(1.);
            fi2 = vec_splats(1.);
            
            VSX_DOUBLES deltax, deltay, deltaz;
            VSX_DOUBLES deltax2, deltay2, deltaz2;
            deltax = px - cx;
            deltax2 = px2 - cx;
            deltay = py - cy;
            deltay2 = py2 - cy;
            deltaz = pz - cz;
            deltaz2 = pz2 - cz;

            // precompute powers of deltaz^k
            VSX_DOUBLES zk[{this_max_zk}], zk2[{this_max_zk}];
            zk[0] = vec_splats(1.);
            zk2[0] = vec_splats(1.);
            ''')
        for c in range(1,this_max_zk):
            w(f'zk[{c}] = deltaz*zk[{c-1}];')
            w(f'zk2[{c}] = deltaz2*zk2[{c-1}];')

        i = 0
        nflop = 3 + this_max_zk-1
        for a in range(order+1):
            w('fij = fi;')
            w('fij2 = fi2;')
            for b in range(order-a+1):
                for c in range(order-a-b+1):
                    w(f'CMvec[{i}] += fij*zk[{c}];')
                    w(f'CMvec2[{i}] += fij2*zk2[{c}];')
                    i += 1
                    nflop += 2
                if b < order-a:
                    w('fij *= deltay;')
                    w('fij2 *= deltay2;')
                    nflop += 1
            if a < order:
                w('\nfi *= deltax;')
                w('\nfi2 *= deltax2;')
                nflop += 1

        w.dedent()
        w('}')  # particle loop

        print(f'nflop: {nflop}; cml: {cml}')

        w('if (n_left != 0){')
        w.indent()
        w(f'''
            assert(0);  // TODO
            // plain C, only one element!
            FLOAT3 p = particles[n_aligned];
            
            double fi, fij;
            fi = 1.0;

            double deltax, deltay, deltaz;
            deltax = p.x - center.x;
            deltay = p.y - center.y;
            deltaz = p.z - center.z;

            double zk[{this_max_zk}];
            zk[0] = 1.;
            ''')
        for c in range(1,this_max_zk):
            w(f'zk[{c}] = deltaz*zk[{c-1}];')

        cml = (order+1)*(order+2)*(order+3)//6

        i = 0
        for a in range(order+1):
            w('fij = fi;')
            for b in range(order-a+1):
                for c in range(order-a-b+1):
                    w(f'CM[{i}] += fij*zk[{c}];')
                    i += 1
                if b < order-a:
                    w('fij *= deltay;')
            if a < order:
                w('\nfi *= deltax;')
        w.dedent()
        w('}')  # remainder block

        w(f'''
            for(int k = 0; k < {cml}; k++)
                for(int j = 0; j < VSX_NVEC_DOUBLE; j++){{
                    CM[k] += CMvec[k][j] + CMvec2[k][j];
                }}
            ''')

        w.dedent()
        w('}\n')  # Kernel

    emit_dispatch_function(w, 'MultipoleVSXKernel(FLOAT3 *particles, int n, double3 center, double *CM)', orders)

    w('''
        #undef VSX_NVEC_DOUBLE
        #endif''')  # VSXMULTIPOLES


def emit_VSX_Taylors(orders, fn='ET_VSX.cpp'):
    w = Writer(fn)

    w('''
        #include "threevector.hh"
        #include "header.cpp"
        
        #ifdef VSXMULTIPOLES
        #include <altivec.h>

        #include "assert.h"

        #define VSX_NVEC_DOUBLE 2
        typedef vector double VSX_DOUBLES;

        template <int Order>
        void TaylorVSXKernel(FLOAT3 *particles, int n, double3 center, double *CT, float3 *acc);
        ''')

    for order in orders:
        cml_orderm1 = (order)*(order+1)*(order+2)//6
        nflop = 0
        nflop_fixed = 0

        cmap = cmapper(order)

        w(f'''
            template <>
            void TaylorVSXKernel<{order}>(FLOAT3 *particles, int n, double3 center, double *CT, float3 *acc){{

                VSX_DOUBLES Qx[{cml_orderm1}], Qy[{cml_orderm1}], Qz[{cml_orderm1}];
            ''')
        w.indent()

        # Precompute Qxyz
        i = 0
        for a in range(order):
            for b in range(order-a):
                for c in range(order-a-b):
                    w(f'''
                        Qx[{i}] = vec_splats({a+1}*CT[{cmap(a+1, b  , c  )}]);
                        Qy[{i}] = vec_splats({b+1}*CT[{cmap(a  , b+1, c  )}]);
                        Qz[{i}] = vec_splats({c+1}*CT[{cmap(a  , b  , c+1)}]);
                    ''')
                    i += 1
                    nflop_fixed += 3

        # Now compute the accelerations
        w(f'''
                VSX_DOUBLES cx = vec_splats(center.x);
                VSX_DOUBLES cy = vec_splats(center.y);
                VSX_DOUBLES cz = vec_splats(center.z);

                int n_left = n % VSX_NVEC_DOUBLE;
                int n_aligned = n - n_left;

                for(int i = 0; i < n_aligned; i += VSX_NVEC_DOUBLE){{
            ''')
        w.indent()

        w(f'''
            VSX_DOUBLES px;
            VSX_DOUBLES py;
            VSX_DOUBLES pz;

            for(int j = 0; j < VSX_NVEC_DOUBLE; j++){{
                px[j] = particles[i+j].x;
                py[j] = particles[i+j].y;
                pz[j] = particles[i+j].z;
            }}

            VSX_DOUBLES fi, fij, fijk;
            fi = vec_splats(1.);
            
            VSX_DOUBLES deltax, deltay, deltaz;
            deltax = px - cx;
            deltay = py - cy;
            deltaz = pz - cz;

            VSX_DOUBLES ax, ay, az;
            ax = vec_splats(0.); 
            ay = vec_splats(0.);
            az = vec_splats(0.);
            ''')  # TODO: could make first loop set instead of accumulate
        nflop += 3 

        i = 0
        for a in range(order):
            w('fij = fi;')
            for b in range(order-a):
                w('fijk = fij;')
                for c in range(order-a-b):
                    w(f'''
                        ax -= Qx[{i}] * fijk;
                        ay -= Qy[{i}] * fijk;
                        az -= Qz[{i}] * fijk;
                        ''')
                    nflop += 6
                    i += 1
                    if c < order-a-b-1:
                        w('fijk *= deltaz;')
                        nflop += 1
                if b < order-a-1:
                    w('fij *= deltay;')
                    nflop += 1
            if a < order-1:
                w('\nfi *= deltax;')
                nflop += 1

        w(f'''
            for(int j = 0; j < VSX_NVEC_DOUBLE; j++){{
                acc[i+j].x = ax[j];
                acc[i+j].y = ay[j];
                acc[i+j].z = az[j];
            }}
            ''')
        w.dedent()
        w('}\n')  # particle loop

        w.dedent()
        w('}\n')  # Kernel

        print(f'Order {order}: {nflop} flop per particle ({nflop_fixed} fixed flop)')

    emit_dispatch_function(w, 'TaylorVSXKernel(FLOAT3 *particles, int n, double3 center, double *CT, float3 *acc)', orders)

    w('#endif')  # VSXMULTIPOLES


def emit_VSX_Taylors_noQ(orders, fn='ET_VSX.cpp'):
    w = Writer(fn)

    w('''
        #include "threevector.hh"
        #include "header.cpp"
        
        #ifdef VSXMULTIPOLES
        #include <altivec.h>

        #include "assert.h"

        #define VSX_NVEC_DOUBLE 2
        typedef vector double VSX_DOUBLES;

        template <int Order>
        void TaylorVSXKernel(FLOAT3 *particles, int n, double3 center, double *CT, float3 *acc);
        ''')

    for order in orders:
        cml_orderm1 = (order)*(order+1)*(order+2)//6
        nflop = 0
        nflop_fixed = 0

        cmap = cmapper(order)

        w(f'''
            template <>
            void TaylorVSXKernel<{order}>(FLOAT3 *particles, int n, double3 center, double *CT, float3 *acc){{
            ''')
        w.indent()

        # Now compute the accelerations
        w(f'''
                VSX_DOUBLES cx = vec_splats(center.x);
                VSX_DOUBLES cy = vec_splats(center.y);
                VSX_DOUBLES cz = vec_splats(center.z);

                int n_left = n % VSX_NVEC_DOUBLE;
                int n_aligned = n - n_left;

                for(int i = 0; i < n_aligned; i += VSX_NVEC_DOUBLE){{
            ''')
        w.indent()

        w(f'''
            VSX_DOUBLES px;
            VSX_DOUBLES py;
            VSX_DOUBLES pz;

            for(int j = 0; j < VSX_NVEC_DOUBLE; j++){{
                px[j] = particles[i+j].x;
                py[j] = particles[i+j].y;
                pz[j] = particles[i+j].z;
            }}

            VSX_DOUBLES fi, fij, fijk;
            fi = vec_splats(1.);
            
            VSX_DOUBLES deltax, deltay, deltaz;
            deltax = px - cx;
            deltay = py - cy;
            deltaz = pz - cz;

            VSX_DOUBLES ax, ay, az;
            ''')
        nflop += 3 

        i = 0
        op = '= -'  # set on the first iteration
        nflop -= 3
        for a in range(order):
            w('fij = fi;')
            for b in range(order-a):
                w('fijk = fij;')
                for c in range(order-a-b):
                    w(f'''
                        ax {op} {a+1}*CT[{cmap(a+1, b  , c  )}] * fijk;
                        ay {op} {b+1}*CT[{cmap(a  , b+1, c  )}] * fijk;
                        az {op} {c+1}*CT[{cmap(a  , b  , c+1)}] * fijk;
                        ''')
                    op = '-='
                    nflop += 9
                    i += 1
                    if c < order-a-b-1:
                        w('fijk *= deltaz;')
                        nflop += 1
                if b < order-a-1:
                    w('fij *= deltay;')
                    nflop += 1
            if a < order-1:
                w('\nfi *= deltax;')
                nflop += 1

        w(f'''
            for(int j = 0; j < VSX_NVEC_DOUBLE; j++){{
                acc[i+j].x = ax[j];
                acc[i+j].y = ay[j];
                acc[i+j].z = az[j];
            }}
            ''')
        w.dedent()
        w('}\n')  # particle loop

        w(f'''
            if(n_left){{
                double deltax, deltay, deltaz;
                deltax = particles[n-1].x - center.x;
                deltay = particles[n-1].y - center.y;
                deltaz = particles[n-1].z - center.z;

                double fi, fij, fijk;
                fi = 1.;

                double ax, ay, az;
            ''')
        w.indent()

        i = 0
        op = '= -'  # set on the first iteration
        for a in range(order):
            w('fij = fi;')
            for b in range(order-a):
                w('fijk = fij;')
                for c in range(order-a-b):
                    w(f'''
                        ax {op} {a+1}*CT[{cmap(a+1, b  , c  )}] * fijk;
                        ay {op} {b+1}*CT[{cmap(a  , b+1, c  )}] * fijk;
                        az {op} {c+1}*CT[{cmap(a  , b  , c+1)}] * fijk;
                        ''')
                    op = '-='
                    i += 1
                    if c < order-a-b-1:
                        w('fijk *= deltaz;')
                if b < order-a-1:
                    w('fij *= deltay;')
            if a < order-1:
                w('\nfi *= deltax;')

        w('''
            acc[n-1].x = ax;
            acc[n-1].y = ay;
            acc[n-1].z = az;
            ''')

        w.dedent()
        w('}')  # remainder loop
        
        w.dedent()
        w('}\n')  # Kernel

        print(f'Order {order}: {nflop} flop per particle ({nflop_fixed} fixed flop)')

    emit_dispatch_function(w, 'TaylorVSXKernel(FLOAT3 *particles, int n, double3 center, double *CT, float3 *acc)', orders)

    w('#endif')  # VSXMULTIPOLES


def emit_VSX_Taylors_interleaved(orders, fn='ET_VSX.cpp'):
    w = Writer(fn)

    w('''
        #include "threevector.hh"
        #include "header.cpp"
        
        #ifdef VSXMULTIPOLES
        #include <altivec.h>

        #include "assert.h"

        #define VSX_NVEC_DOUBLE 2
        typedef vector double VSX_DOUBLES;

        template <int Order>
        void TaylorVSXKernel(FLOAT3 *particles, int n, double3 center, double *CT, float3 *acc);
        ''')

    for order in orders:
        cml_orderm1 = (order)*(order+1)*(order+2)//6
        nflop = 0
        nflop_fixed = 0

        cmap = cmapper(order)

        w(f'''
            template <>
            void TaylorVSXKernel<{order}>(FLOAT3 *particles, int n, double3 center, double *CT, float3 *acc){{

                VSX_DOUBLES Qx[{cml_orderm1}], Qy[{cml_orderm1}], Qz[{cml_orderm1}];
                VSX_DOUBLES Qx2[{cml_orderm1}], Qy2[{cml_orderm1}], Qz2[{cml_orderm1}];
            ''')
        w.indent()

        # Precompute Qxyz
        i = 0
        for a in range(order):
            for b in range(order-a):
                for c in range(order-a-b):
                    w(f'''
                        Qx[{i}] = vec_splats({a+1}*CT[{cmap(a+1, b  , c  )}]);
                        Qx2[{i}] = Qx[{i}];
                        Qy[{i}] = vec_splats({b+1}*CT[{cmap(a  , b+1, c  )}]);
                        Qy2[{i}] = Qy[{i}];
                        Qz[{i}] = vec_splats({c+1}*CT[{cmap(a  , b  , c+1)}]);
                        Qz2[{i}] = Qz[{i}];
                    ''')
                    i += 1
                    nflop_fixed += 3

        # Now compute the accelerations
        w(f'''
                VSX_DOUBLES cx = vec_splats(center.x);
                VSX_DOUBLES cy = vec_splats(center.y);
                VSX_DOUBLES cz = vec_splats(center.z);

                int n_left = n % (2*VSX_NVEC_DOUBLE);
                int n_aligned = n - n_left;

                for(int i = 0; i < n_aligned; i += 2*VSX_NVEC_DOUBLE){{
            ''')
        w.indent()

        w(f'''
            VSX_DOUBLES px,py,pz;
            VSX_DOUBLES px2,py2,pz2;

            for(int j = 0; j < VSX_NVEC_DOUBLE; j++){{
                px[j] = particles[i+j].x;
                px2[j] = particles[i+j+VSX_NVEC_DOUBLE].x;
                py[j] = particles[i+j].y;
                py2[j] = particles[i+j+VSX_NVEC_DOUBLE].y;
                pz[j] = particles[i+j].z;
                pz2[j] = particles[i+j+VSX_NVEC_DOUBLE].z;
            }}

            VSX_DOUBLES fi, fij, fijk;
            VSX_DOUBLES fi2, fij2, fijk2;
            fi = vec_splats(1.);
            fi2 = vec_splats(1.);
            
            VSX_DOUBLES deltax, deltay, deltaz;
            VSX_DOUBLES deltax2, deltay2, deltaz2;
            deltax = px - cx;
            deltax2 = px2 - cx;
            deltay = py - cy;
            deltay2 = py2 - cy;
            deltaz = pz - cz;
            deltaz2 = pz2 - cz;

            VSX_DOUBLES ax, ay, az;
            VSX_DOUBLES ax2, ay2, az2;
            ax = vec_splats(0.);
            ax2 = vec_splats(0.);  
            ay = vec_splats(0.);
            ay2 = vec_splats(0.);
            az = vec_splats(0.);
            az2 = vec_splats(0.);
            ''')  # TODO: could make first loop set instead of accumulate
        nflop += 3 

        i = 0
        for a in range(order):
            w('fij = fi;')
            w('fij2 = fi2;')
            for b in range(order-a):
                w('fijk = fij;')
                w('fijk2 = fij2;')
                for c in range(order-a-b):
                    w(f'''
                        ax -= Qx[{i}] * fijk;
                        ax2 -= Qx2[{i}] * fijk2;
                        ay -= Qy[{i}] * fijk;
                        ay2 -= Qy2[{i}] * fijk2;
                        az -= Qz[{i}] * fijk;
                        az2 -= Qz2[{i}] * fijk2;
                        ''')
                    nflop += 6
                    i += 1
                    if c < order-a-b-1:
                        w('fijk *= deltaz;')
                        w('fijk2 *= deltaz2;')
                        nflop += 1
                if b < order-a-1:
                    w('fij *= deltay;')
                    w('fij2 *= deltay2;')
                    nflop += 1
            if a < order-1:
                w('\nfi *= deltax;')
                w('\nfi2 *= deltax2;')
                nflop += 1

        w(f'''
            for(int j = 0; j < VSX_NVEC_DOUBLE; j++){{
                acc[i+j].x = ax[j];
                acc[i+j+VSX_NVEC_DOUBLE].x = ax2[j];
                acc[i+j].y = ay[j];
                acc[i+j+VSX_NVEC_DOUBLE].y = ay2[j];
                acc[i+j].z = az[j];
                acc[i+j+VSX_NVEC_DOUBLE].z = az2[j];
            }}
            ''')
        w.dedent()
        w('}\n')  # particle loop

        w.dedent()
        w('}\n')  # Kernel

        print(f'Order {order}: {nflop} flop per particle ({nflop_fixed} fixed flop)')

    emit_dispatch_function(w, 'TaylorVSXKernel(FLOAT3 *particles, int n, double3 center, double *CT, float3 *acc)', orders)

    w('#endif')  # VSXMULTIPOLES

if __name__ == '__main__':
    parser = default_metacode_argparser(doc=__doc__)
    args = parser.parse_args()

    if args.onlyorder:
        orders = [args.onlyorder]
    else:
        orders = list(range(1,args.maxorder+1))

    #emit_VSX_Multipoles(orders)
    #emit_VSX_Multipoles_interleaved(orders)
    emit_VSX_Multipoles_FMA(orders)
    #emit_VSX_Multipoles_FMA_interleaved(orders)
    #emit_VSX_Taylors(orders)
    emit_VSX_Taylors_noQ(orders)  # currently fastest on Summit for all sizes
    #emit_VSX_Taylors_interleaved(orders)
