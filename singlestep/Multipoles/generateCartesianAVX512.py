#!/usr/bin/env python3
'''
This program generates AVX-512 assembly for Taylors and Multipoles
evaluation. You shouldn't have to run this program manually; it will
be invoked from the Makefile.  But to see the manual usage, run:

$ ./generateCartesianAVX512.py --help

The compiler can't figure out how to unroll some of the nested
loops, so the purpose of this program is to "manually" unroll them.

Curiously, the "un-unrolled" version of the Taylors seem faster.
That's implemented in EvaluateTaylors.cpp.
'''

# TODO: recently reformatted, test on hal

from multipoles_metacode_utils \
    import emit_dispatch_function, Writer, default_metacode_argparser, cmapper

def emit_AVX512_Multipoles(orders, fn='CMAVX512.cpp'):
    w = Writer(fn)

    w('''
        #include "config.h"
        #ifdef AVX512MULTIPOLES

        #include "header.cpp"
        #include "threevector.hh"

        #include "avx512_calls.h"
        #include "assert.h"
        
        template <int Order>
        void Multipole512Kernel(FLOAT3 *xyz, int n, FLOAT3 center, double *CM);''')

    for order in orders:
        cml = (order+1)*(order+2)*(order+3)//6
        
        w(f'''
            template <>
            void Multipole512Kernel<{order}>(FLOAT3 *xyz, int n, FLOAT3 center, double *_CM){{

                AVX512_DOUBLES cx = AVX512_SET_DOUBLE(center.x);
                AVX512_DOUBLES cy = AVX512_SET_DOUBLE(center.y);
                AVX512_DOUBLES cz = AVX512_SET_DOUBLE(center.z);

                AVX512_DOUBLES CM[{cml}];
                for(int i = 0; i < {cml}; i++)
                    CM[i] = AVX512_SETZERO_DOUBLE();

                int n_aligned = n - (n % AVX512_NVEC_DOUBLE);
                int nleft = n - n_aligned;

                for(int k=0; k <= n_aligned-AVX512_NVEC_DOUBLE; k += AVX512_NVEC_DOUBLE) {{
                    // Load 8 DOUBLE3s as List3s
                    AVX512_DOUBLES px, py, pz;
                    for(int j = 0; j < AVX512_NVEC_DOUBLE; j++){{
                        px[j] = xyz[k+j].x;
                        py[j] = xyz[k+j].y;
                        pz[j] = xyz[k+j].z;
                    }}

                    AVX512_DOUBLES fi = AVX512_SET_DOUBLE(1.);
                    AVX512_DOUBLES fij, fijk;
                    AVX512_DOUBLES deltax = AVX512_SUBTRACT_DOUBLES(px, cx);
                    AVX512_DOUBLES deltay = AVX512_SUBTRACT_DOUBLES(py, cy);
                    AVX512_DOUBLES deltaz = AVX512_SUBTRACT_DOUBLES(pz, cz);
        ''')
        w.indent(2)

        i = 0
        for a in range(order+1):
            w('fij = fi;')
            for b in range(order-a+1):
                w('fijk = fij;')
                for c in range(order-a-b+1):
                    w(f'CM[{i}] = AVX512_ADD_DOUBLES(CM[{i}], fijk);')
                    i += 1
                    if(c < order-a-b):
                        w('fijk = AVX512_MULTIPLY_DOUBLES(fijk, deltaz);')
                if(b < order-a):
                    w('\nfij = AVX512_MULTIPLY_DOUBLES(fij, deltay);')
            if(a < order):
                w('\nfi = AVX512_MULTIPLY_DOUBLES(fi, deltax);')

        w.dedent()
        w('}')  # particle loop

        # remainder loop
        w(f'''
            if(n_aligned < n){{
                // Load nleft DOUBLE3s as List3s
                AVX512_DOUBLES px = AVX512_SETZERO_DOUBLE();
                AVX512_DOUBLES py = AVX512_SETZERO_DOUBLE();
                AVX512_DOUBLES pz = AVX512_SETZERO_DOUBLE();
                for(int j = 0; j < nleft; j++){{
                    px[j] = xyz[n_aligned+j].x;
                    py[j] = xyz[n_aligned+j].y;
                    pz[j] = xyz[n_aligned+j].z;
                }}

                AVX512_DOUBLES deltax = AVX512_SUBTRACT_DOUBLES(px, cx);
                AVX512_DOUBLES deltay = AVX512_SUBTRACT_DOUBLES(py, cy);
                AVX512_DOUBLES deltaz = AVX512_SUBTRACT_DOUBLES(pz, cz);

                // Evaluate these 8 particles
                AVX512_DOUBLES fi,fij,fijk;
                fi = AVX512_SET_DOUBLE(1.0);
                AVX512_MASK_DOUBLE left_mask = masks_per_misalignment_value_double[nleft];
                ''')
        w.indent()

        i = 0
        ops = 0
        for a in range(order+1):
            w('fij = fi;')
            for b in range(order-a+1):
                w('fijk = fij;')
                for c in range(order-a-b+1):
                    w(f'CM[{i}] = AVX512_MASK_ADD_DOUBLES(CM[{i}], left_mask, CM[{i}], fijk);')
                    ops += 1
                    i += 1
                    if(c < order-a-b):
                        w('fijk = AVX512_MULTIPLY_DOUBLES(fijk, deltaz);')
                        ops += 1
                        
                if(b < order-a):
                    w('\nfij = AVX512_MULTIPLY_DOUBLES(fij, deltay);')
                    ops += 1
            if(a < order):
                w('\nfi = AVX512_MULTIPLY_DOUBLES(fi, deltax);')
                ops += 1
        w.dedent()
        w('}')

        w(f'''
            for(int i = 0; i < {cml}; i++){{
                _CM[i] = AVX512_HORIZONTAL_SUM_DOUBLES(CM[i]);
            }}
            ''')

        w.dedent()
        w('}')  # Kernel
        
        #print(f'Order {order}: {ops} FLOP')

    emit_dispatch_function(w, '''Multipole512Kernel(FLOAT3 *xyz, int n, FLOAT3 center, double *CM)''', orders)

    w('#endif')  # AVX512MULTIPOLES

# this is only marginally faster than the non-FMA version
# TODO: finish bringing particle loop inside
def emit_AVX512_Multipoles_FMA(orders, fn='CMAVX512.cpp', max_zk=None, unroll=4, naccum=1):
    assert max_zk is None, "max_zk probably not faster; implementation not finished"

    if max_zk is None:
        max_zk = max(orders) + 1

    w = Writer(fn)

    w('''
        #include "config.h"
        #ifdef AVX512MULTIPOLES

        #include "header.cpp"
        #include "threevector.hh"

        #include "avx512_calls.h"
        #include "assert.h"
        
        template <int Order>
        void Multipole512Kernel(FLOAT3 *xyz, int n, FLOAT3 center, double *CM);
        ''')

    for order in orders:
        this_max_zk = min(order+1, max_zk)

        cml = (order+1)*(order+2)*(order+3)//6
        
        w(f'''
            template <>
            void Multipole512Kernel<{order}>(FLOAT3 *xyz, int n, FLOAT3 center, double *_CM){{

                AVX512_DOUBLES cx = AVX512_SET_DOUBLE(center.x);
                AVX512_DOUBLES cy = AVX512_SET_DOUBLE(center.y);
                AVX512_DOUBLES cz = AVX512_SET_DOUBLE(center.z);
                
                AVX512_DOUBLES CM[{cml}][{naccum}];
                for(int i = 0; i < {cml}; i++)
                    for(int ii = 0; ii < {naccum}; ii++)
                        CM[i][ii] = AVX512_SETZERO_DOUBLE();
                
                int n_partial = n % ({unroll}*AVX512_NVEC_DOUBLE);
                int n_full = n - n_partial;
                int n_left = n_partial % AVX512_NVEC_DOUBLE;
                n_partial -= n_left;
                
                for(int k=0; k < n_full; k += {unroll}*AVX512_NVEC_DOUBLE) {{
            ''')
        w.indent(2)

        def _emit_unrolled(start, unroll, mask=False):
            w(f'''
                // Load 8 DOUBLE3s as List3s
                AVX512_DOUBLES px[{unroll}], py[{unroll}], pz[{unroll}];

                AVX512_DOUBLES fi[{unroll}];
                AVX512_DOUBLES fij[{unroll}];
                AVX512_DOUBLES deltax[{unroll}];
                AVX512_DOUBLES deltay[{unroll}];
                AVX512_DOUBLES deltaz[{unroll}];
                ''')
            for c in range(this_max_zk):
                w(f'AVX512_DOUBLES zk{c}[{unroll}];')
            
            ops = 0
            nload = mask or 'AVX512_NVEC_DOUBLE'
            for u in range(unroll):
                w(f'''
                    fi[{u}] = AVX512_SET_DOUBLE(1.);
                    for(int j = 0; j < {nload}; j++){{
                        px[{u}][j] = xyz[{start}+j+{u}*AVX512_NVEC_DOUBLE].x;
                        py[{u}][j] = xyz[{start}+j+{u}*AVX512_NVEC_DOUBLE].y;
                        pz[{u}][j] = xyz[{start}+j+{u}*AVX512_NVEC_DOUBLE].z;
                    }}
                    zk0[{u}] = AVX512_SET_DOUBLE(1.0);
                    for(int j = {nload}; j < AVX512_NVEC_DOUBLE; j++){{
                        px[{u}][j] = center.x;
                        py[{u}][j] = center.y;
                        pz[{u}][j] = center.z;
                        zk0[{u}][j] = 0.;
                    }}

                    deltax[{u}] = AVX512_SUBTRACT_DOUBLES(px[{u}], cx);
                    deltay[{u}] = AVX512_SUBTRACT_DOUBLES(py[{u}], cy);
                    deltaz[{u}] = AVX512_SUBTRACT_DOUBLES(pz[{u}], cz);
                ''')
            
            ops += 3
            
            for c in range(1,this_max_zk):
                for u in range(unroll):
                    w(f'zk{c}[{u}] = AVX512_MULTIPLY_DOUBLES(zk{c-1}[{u}], deltaz[{u}]);')
                ops += 1

            i = 0
            for a in range(order+1):
                for u in range(unroll):
                    w(f'fij[{u}] = fi[{u}];')
                for b in range(order-a+1):
                    for c in range(order-a-b+1):
                        if c % this_max_zk == 0 and c > 0:
                            for u in range(unroll):
                                w(f'fij[{u}] = AVX512_MULTIPLY_DOUBLES(fij[{u}], zk{this_max_zk-1}[{u}]);')
                            ops += 1
                        for u in range(unroll):
                            w(f'CM[{i}][{u%naccum}] = AVX512_FMA_ADD_DOUBLES(fij[{u}], zk{c % this_max_zk}[{u}], CM[{i}][{u%naccum}]);')
                        ops += 2
                        i += 1
                    if b < order-a:
                        for u in range(unroll):
                            w(f'\nfij[{u}] = AVX512_MULTIPLY_DOUBLES(fij[{u}], deltay[{u}]);')
                        ops += 1
                if a < order:
                    for u in range(unroll):
                        w(f'\nfi[{u}] = AVX512_MULTIPLY_DOUBLES(fi[{u}], deltax[{u}]);')
                    ops += 1
            
            #print(f'Order {order}: {ops} FLOP')
        _emit_unrolled('k', unroll)
        w.dedent()
        w('}')  # Particle loop
        
        w('''
        switch(n_partial){
        ''')
        w.indent()
        for u in range(1,unroll):
            w(f'case {u}*AVX512_NVEC_DOUBLE: {{')
            w.indent()
            _emit_unrolled('n_full', u)
            w('break;\n}')
            w.dedent()
        w.dedent()
        w('}\n')  # switch
        
        w('if(n_left){')
        w.indent()
        _emit_unrolled('n_full+n_partial', 1, mask='n_left')
        w.dedent()
        w('}')  # mask iteration
        
        # All done! Sum the vectors.
        w(f'''
            for(int i = 0; i < {cml}; i++){{
        ''')
        w.indent(1)
        for u in range(1,naccum):
            w(f'CM[i][0] = AVX512_ADD_DOUBLES(CM[i][0], CM[i][{u}]);')
            
        w('_CM[i] = AVX512_HORIZONTAL_SUM_DOUBLES(CM[i][0]);')
        w.dedent()
        w('}')  # cml
        
        w.dedent()
        w('}')  # Kernel

    emit_dispatch_function(w, '''Multipole512Kernel(FLOAT3 *xyz, int n, FLOAT3 center, double *CM)''', orders)

    w('#endif')  # AVX512MULTIPOLES


def emit_AVX512_Taylors(orders, fn='ETAVX512.cpp'):
    w = Writer(fn)

    w('''
        #include "config.h"
        
        #ifdef AVX512MULTIPOLES
        #include "header.cpp"
        #include "threevector.hh"
        
        #include "avx512_calls.h"
        #include "assert.h"

        template <int Order>
        void Taylor512Kernel(double *CT, FLOAT3 center, int n, FLOAT3 *xyz, FLOAT3 *acc);
        ''')

    for order in orders:
        # taylors use a=0..order, multipoles use a=0..order+1
        cml_orderm1 = (order)*(order+1)*(order+2)//6

        cmap = cmapper(order)

        w(f'''
            template <>
            void Taylor512Kernel<{order}>(double *CT, FLOAT3 center, int n, FLOAT3 *xyz, FLOAT3 *acc) {{

                AVX512_DOUBLES Qx[{cml_orderm1}], Qy[{cml_orderm1}], Qz[{cml_orderm1}];
            ''')
        w.indent()

        i = 0
        for a in range(order):
            for b in range(order-a):
                for c in range(order-a-b):
                    w(f'''
                        Qx[{i}] = AVX512_SET_DOUBLE(-{a+1}*CT[{cmap(a+1,b  ,c  )}]);
                        Qy[{i}] = AVX512_SET_DOUBLE(-{b+1}*CT[{cmap(a  ,b+1,c  )}]);
                        Qz[{i}] = AVX512_SET_DOUBLE(-{c+1}*CT[{cmap(a  ,b  ,c+1)}]);
                        ''')
                    i += 1

        # Now compute the accelerations
        w(f'''
            AVX512_DOUBLES cx = AVX512_SET_DOUBLE(center.x);
            AVX512_DOUBLES cy = AVX512_SET_DOUBLE(center.y);
            AVX512_DOUBLES cz = AVX512_SET_DOUBLE(center.z);

            int n_aligned = n - (n % AVX512_NVEC_DOUBLE);
            int nleft = n - n_aligned;

            for(int k=0; k < n_aligned; k += AVX512_NVEC_DOUBLE) {{
            ''')
        w.indent()

        w(f'''
            // Load 8 DOUBLE3s as List3s
            AVX512_DOUBLES px, py, pz;
            for(int j = 0; j < AVX512_NVEC_DOUBLE; j++){{
                px[j] = xyz[k+j].x;
                py[j] = xyz[k+j].y;
                pz[j] = xyz[k+j].z;
            }}

            AVX512_DOUBLES fi = AVX512_SET_DOUBLE(1.);
            AVX512_DOUBLES fij, fijk;
            AVX512_DOUBLES deltax = AVX512_SUBTRACT_DOUBLES(px, cx);
            AVX512_DOUBLES deltay = AVX512_SUBTRACT_DOUBLES(py, cy);
            AVX512_DOUBLES deltaz = AVX512_SUBTRACT_DOUBLES(pz, cz);

            AVX512_DOUBLES ax, ay, az;
            ax = AVX512_SETZERO_DOUBLE();
            ay = AVX512_SETZERO_DOUBLE();
            az = AVX512_SETZERO_DOUBLE();
            ''')

        i = 0
        for a in range(order):
            w('fij = fi;')
            for b in range(order-a):
                w('fijk = fij;')
                for c in range(order-a-b):
                    w(f'''
                        ax = AVX512_FMA_ADD_DOUBLES(Qx[{i}], fijk, ax);
                        ay = AVX512_FMA_ADD_DOUBLES(Qy[{i}], fijk, ay);
                        az = AVX512_FMA_ADD_DOUBLES(Qz[{i}], fijk, az);
                        ''')
                    i += 1
                    if(c < order-a-b-1):
                        w('fijk = AVX512_MULTIPLY_DOUBLES(fijk, deltaz);')
                if(b < order-a-1):
                    w('\nfij = AVX512_MULTIPLY_DOUBLES(fij, deltay);')
            if(a < order-1):
                w('\nfi = AVX512_MULTIPLY_DOUBLES(fi, deltax);')

        w(f'''
            for(int j = 0; j < AVX512_NVEC_DOUBLE; j++){{
                acc[k+j].x = ax[j];
                acc[k+j].y = ay[j];
                acc[k+j].z = az[j];
            }}
            ''')

        w.dedent()
        w('}')  # particle loop

        # Remainder loop
        w(f'''
            if (nleft) {{
                // Load nleft DOUBLE3s as List3s
                // Init these to zero to avoid FP signals
                // There are other ways to do this, such as masked operations
                AVX512_DOUBLES px = AVX512_SETZERO_DOUBLE();
                AVX512_DOUBLES py = AVX512_SETZERO_DOUBLE();
                AVX512_DOUBLES pz = AVX512_SETZERO_DOUBLE();

                for(int j = 0; j < nleft; j++){{
                    px[j] = xyz[n_aligned+j].x;
                    py[j] = xyz[n_aligned+j].y;
                    pz[j] = xyz[n_aligned+j].z;
                }}

                AVX512_DOUBLES fi = AVX512_SET_DOUBLE(1.);
                AVX512_DOUBLES fij, fijk;
                AVX512_DOUBLES deltax = AVX512_SUBTRACT_DOUBLES(px, cx);
                AVX512_DOUBLES deltay = AVX512_SUBTRACT_DOUBLES(py, cy);
                AVX512_DOUBLES deltaz = AVX512_SUBTRACT_DOUBLES(pz, cz);

                AVX512_DOUBLES ax, ay, az;
                ax = AVX512_SETZERO_DOUBLE();
                ay = AVX512_SETZERO_DOUBLE();
                az = AVX512_SETZERO_DOUBLE();
            ''')
        w.indent()

        i = 0
        for a in range(order):
            w('fij = fi;')
            for b in range(order-a):
                w('fijk = fij;')
                for c in range(order-a-b):
                    w(f'''
                        ax = AVX512_FMA_ADD_DOUBLES(Qx[{i}], fijk, ax);
                        ay = AVX512_FMA_ADD_DOUBLES(Qy[{i}], fijk, ay);
                        az = AVX512_FMA_ADD_DOUBLES(Qz[{i}], fijk, az);
                        ''')
                    i += 1
                    if(c < order-a-b-1):
                        w('fijk = AVX512_MULTIPLY_DOUBLES(fijk, deltaz);')
                if(b < order-a-1):
                    w('\nfij = AVX512_MULTIPLY_DOUBLES(fij, deltay);')
            if(a < order-1):
                w('\nfi = AVX512_MULTIPLY_DOUBLES(fi, deltax);')

        w(f'''
            for(int j = 0; j < nleft; j++){{
                acc[n_aligned+j].x = ax[j];
                acc[n_aligned+j].y = ay[j];
                acc[n_aligned+j].z = az[j];
            }}
            ''')
        w.dedent()
        w('}')  # remainder loop

        w.dedent()
        w('}')  # kernel

    emit_dispatch_function(w, 'Taylor512Kernel(double *CT, FLOAT3 center, int n, FLOAT3 *xyz, FLOAT3 *acc)', orders)

    w('#endif')  # AVX512MULTIPOLES


if __name__ == '__main__':
    parser = default_metacode_argparser(doc=__doc__)
    args = parser.parse_args()

    if args.onlyorder:
        orders = [args.onlyorder]
    else:
        orders = list(range(1,args.maxorder+1))

    emit_AVX512_Taylors(orders)
    #emit_AVX512_Multipoles(orders)
    emit_AVX512_Multipoles_FMA(orders)
