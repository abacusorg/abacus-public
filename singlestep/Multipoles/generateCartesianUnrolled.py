#!/usr/bin/env python3
'''
This program generates plain C code for the Taylors and Multipoles
evaluation. You shouldn't have to run this program manually; it will
be invoked from the Makefile.  But to see the manual usage, run:

$ ./generateCartesianUnrolled.py --help

The compiler can't figure out how to unroll some of the nested
loops, so the purpose of this program is to "manually" unroll them.
'''

from multipoles_metacode_utils \
    import emit_dispatch_function, Writer, default_metacode_argparser, cmapper


def emit_unrolled_Multipoles(orders, fn='CM_unrolled.cpp', unroll=1):
    w = Writer(fn)

    w('''
        #include "threevector.hh"
        #include "header.cpp"

        #ifdef UNROLLEDMULTIPOLES

        #include "assert.h"

        template <int Order>
        void MultipoleUnrolledKernel(FLOAT3 *particles, int n, double3 center, double *CM);
        ''')

    for order in orders:
        cml = (order+1)*(order+2)*(order+3)//6

        w(f'''
            template <>
            void MultipoleUnrolledKernel<{order}>(FLOAT3 *particles, int n, double3 center, double *CM){{
                double CMlocal[{cml}] = {{}};''')
        w.indent()

        w(f'for(int i = 0; i < n; i+={unroll}){{')
        w.indent()

        w(f'''
            FLOAT3 p[{unroll}];
            double fi[{unroll}], fij[{unroll}], fijk[{unroll}];
            ''')
        for m in range(unroll):
            w(f'p[{m}] = particles[i+{m}];')
            w(f'fi[{m}] = 1.;')

        w(f'double deltax[{unroll}], deltay[{unroll}], deltaz[{unroll}];')
        for m in range(unroll):
            w(f'''
                deltax[{m}] = p[{m}].x - center.x;
                deltay[{m}] = p[{m}].y - center.y;
                deltaz[{m}] = p[{m}].z - center.z;
                ''')

        i = 0
        nflop = 3  # deltas
        for a in range(order+1):
            for m in range(unroll):
                w(f'fij[{m}] = fi[{m}];')
            for b in range(order-a+1):
                for m in range(unroll):
                    w(f'fijk[{m}] = fij[{m}];')
                for c in range(order-a-b+1):
                    for m in range(unroll):
                        w(f'CMlocal[{i}] += fijk[{m}];')
                    i += 1
                    nflop += 1
                    if c < order-a-b:
                        for m in range(unroll):
                            w(f'fijk[{m}] *= deltaz[{m}];')
                        nflop += 1
                if b < order-a:
                    for m in range(unroll):
                        w(f'fij[{m}] *= deltay[{m}];')
                    nflop += 1
            if a < order:
                for m in range(unroll):
                    w(f'\nfi[{m}] *= deltax[{m}];')
                nflop += 1

        #print(f'nflop: {nflop}; cml: {cml}')

        w.dedent()
        w('}')  # particle loop

        w(f'''
            for(int i = 0; i < {cml}; i++)
                CM[i] = CMlocal[i];
            ''')
        w.dedent()
        w('}\n')  # Kernel

    emit_dispatch_function(w, 'MultipoleUnrolledKernel(FLOAT3 *particles, int n, double3 center, double *CM)', orders)

    w('#endif')  # UNROLLEDMULTIPOLES


def emit_unrolled_Multipoles_FMA(orders, fn='CM_unrolled.cpp', max_zk=None):
    assert max_zk is None, "max_zk implementation not finished"
    if max_zk is None:
        max_zk = max(orders) + 1

    w = Writer(fn)

    w('''
        #include "threevector.hh"
        #include "header.cpp"

        #ifdef UNROLLEDMULTIPOLES

        #include "assert.h"

        template <int Order>
        void MultipoleUnrolledKernel(FLOAT3 *particles, int n, double3 center, double *CM);
        ''')

    for order in orders:
        this_max_zk = min(order+1, max_zk)
        cml = (order+1)*(order+2)*(order+3)//6

        w(f'''
            template <>
            void MultipoleUnrolledKernel<{order}>(FLOAT3 *particles, int n, double3 center, double *CM){{''')
        w.indent()

        w('for(int i = 0; i < n; i++){')
        w.indent()

        w('''
            FLOAT3 p = particles[i];
            double fi, fij;
            fi = 1.0;
            ''')

        w(f'''
            double deltax, deltay, deltaz;
            deltax = p.x - center.x;
            deltay = p.y - center.y;
            deltaz = p.z - center.z;

            double zk[{this_max_zk}];
            zk[0] = 1.;
            ''')
        for c in range(1,this_max_zk):
            w(f'zk[{c}] = deltaz*zk[{c-1}];')

        i = 0
        nflop = 3  # deltas
        for a in range(order+1):
            w('fij = fi;')
            for b in range(order-a+1):
                for c in range(order-a-b+1):
                    w(f'CM[{i}] = CM[{i}] + fij*zk[{c}];')
                    i += 1
                    nflop += 2
                if b < order-a:
                    w('fij *= deltay;')
                    nflop += 1
            if a < order:
                w('\nfi *= deltax;')
                nflop += 1

        #print(f'nflop: {nflop}; cml: {cml}')

        w.dedent()
        w('}')  # particle loop
        w.dedent()
        w('}\n')  # Kernel

    emit_dispatch_function(w, 'MultipoleUnrolledKernel(FLOAT3 *particles, int n, double3 center, double *CM)', orders)

    w('#endif')  # UNROLLEDMULTIPOLES


def emit_unrolled_Taylors(orders, fn='ET_unrolled.cpp'):
    w = Writer(fn)

    w('''
        #include "threevector.hh"
        #include "header.cpp"
        
        #ifdef UNROLLEDMULTIPOLES
        
        #include "assert.h"

        template <int Order>
        void TaylorUnrolledKernel(FLOAT3 *particles, int n, double3 center, double *CT, float3 *acc);

        ''')

    for order in orders:
        cml_orderm1 = (order)*(order+1)*(order+2)//6
        cmap = cmapper(order)

        w(f'''
            template <>
            void TaylorUnrolledKernel<{order}>(FLOAT3 *particles, int n, double3 center, double *CT, float3 *acc){{

                double3 Q[{cml_orderm1}];
            ''')
        w.indent()

        # Precompute Qxyz
        i = 0
        for a in range(order):
            for b in range(order-a):
                for c in range(order-a-b):
                    w(f'''
                        Q[{i}].x = {a+1}*CT[{cmap(a+1, b  , c  )}];
                        Q[{i}].y = {b+1}*CT[{cmap(a  , b+1, c  )}];
                        Q[{i}].z = {c+1}*CT[{cmap(a  , b  , c+1)}];
                    ''')
                    i += 1

        w('''
            #pragma GCC ivdep
            for(int j = 0; j < n; j++){
            ''')
        w.indent()

        w('''
            FLOAT3 p = particles[j];
            double3 thisacc(0.);
            ''')

        w('''
            double fi, fij, fijk;
            fi = 1.0;
            ''')

        w('''
            double deltax, deltay, deltaz;
            deltax = p.x - center.x;
            deltay = p.y - center.y;
            deltaz = p.z - center.z;
            ''')

        i = 0
        for a in range(order):
            w('fij = fi;')
            for b in range(order-a):
                w('fijk = fij;')
                for c in range(order-a-b):
                    w(f'''
                        thisacc -= Q[{i}] * fijk;
                        ''')
                    i += 1
                    if c < order-a-b-1:
                        w('fijk *= deltaz;')
                if b < order-a-1:
                    w('fij *= deltay;')
            if a < order-1:
                w('\nfi *= deltax;')

        w('acc[j] = thisacc;')
        w.dedent()
        w('}\n')  # particle loop

        w.dedent()
        w('}\n')  # Kernel

    emit_dispatch_function(w, 'TaylorUnrolledKernel(FLOAT3 *particles, int n, double3 center, double *CT, float3 *acc)', orders)

    w('#endif')  # UNROLLEDMULTIPOLES


def emit_unrolled_Taylors_noQ(orders, fn='ET_unrolled.cpp'):
    w = Writer(fn)

    w('''
        #include "threevector.hh"
        #include "header.cpp"
        
        #ifdef UNROLLEDMULTIPOLES
        
        #include "assert.h"

        template <int Order>
        void TaylorUnrolledKernel(FLOAT3 *particles, int n, double3 center, double *CT, float3 *acc);

        ''')

    for order in orders:
        cml_orderm1 = (order)*(order+1)*(order+2)//6
        cmap = cmapper(order)

        w(f'''
            template <>
            void TaylorUnrolledKernel<{order}>(FLOAT3 *particles, int n, double3 center, double *CT, float3 *acc){{
            ''')
        w.indent()

        w('''
            #pragma GCC ivdep
            for(int j = 0; j < n; j++){
            ''')
        w.indent()

        w('''
            FLOAT3 p = particles[j];
            double3 thisacc(0.);
            ''')

        w('''
            double fi, fij, fijk;
            fi = 1.0;
            ''')

        w('''
            double deltax, deltay, deltaz;
            deltax = p.x - center.x;
            deltay = p.y - center.y;
            deltaz = p.z - center.z;
            ''')

        i = 0
        for a in range(order):
            w('fij = fi;')
            for b in range(order-a):
                w('fijk = fij;')
                for c in range(order-a-b):
                    w(f'''
                        thisacc.x -= {a+1}*CT[{cmap(a+1, b  , c  )}] * fijk;
                        thisacc.y -= {b+1}*CT[{cmap(a  , b+1, c  )}] * fijk;
                        thisacc.z -= {c+1}*CT[{cmap(a  , b  , c+1)}] * fijk;
                        ''')
                    i += 1
                    if c < order-a-b-1:
                        w('fijk *= deltaz;')
                if b < order-a-1:
                    w('fij *= deltay;')
            if a < order-1:
                w('\nfi *= deltax;')

        w('acc[j] = thisacc;')
        w.dedent()
        w('}\n')  # particle loop

        w.dedent()
        w('}\n')  # Kernel

    emit_dispatch_function(w, 'TaylorUnrolledKernel(FLOAT3 *particles, int n, double3 center, double *CT, float3 *acc)', orders)

    w('#endif')  # UNROLLEDMULTIPOLES


if __name__ == '__main__':
    parser = default_metacode_argparser(doc=__doc__)
    args = parser.parse_args()

    if args.onlyorder:
        orders = [args.onlyorder]
    else:
        orders = list(range(1,args.maxorder+1))

    emit_unrolled_Multipoles(orders)
    #emit_unrolled_Multipoles_FMA(orders)
    #emit_unrolled_Taylors(orders)
    emit_unrolled_Taylors_noQ(orders)  # a little faster
