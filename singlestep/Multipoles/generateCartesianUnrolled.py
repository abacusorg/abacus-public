#!/usr/bin/env python
'''
This program generates plain C code for the Taylors and Multipoles
evaluation. You shouldn't have to run this program manually; it will
be invoked from the Makefile.  But to see the manual usage, run:

$ ./generateCartesianUnrolled.py --help

The compiler can't figure out how to unroll some of the nested
loops, so the purpose of this program is to "manually" unroll them.
'''

from multipoles_metacode_utils \
    import emit_dispatch_function, Writer, default_metacode_argparser


def emit_unrolled_Multipoles(orders, fn='CM_unrolled.cpp'):
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
        w(f'''
            template <>
            void MultipoleUnrolledKernel<{order}>(FLOAT3 *particles, int n, double3 center, double *CM){{''')
        w.indent()

        w('for(int i = 0; i < n; i++){')
        w.indent()

        w('''
            FLOAT3 p = particles[i];
            double fi, fij, fijk;
            fi = 1.0;
            ''')

        w('''
            double deltax, deltay, deltaz;
            deltax = p.x - center.x;
            deltay = p.y - center.y;
            deltaz = p.z - center.z;
            ''')

        cml = (order+1)*(order+2)*(order+3)//6

        i = 0
        nflop = 3  # deltas
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

        #print(f'nflop: {nflop}; cml: {cml}')

        w.dedent()
        w('}')  # particle loop
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
        void TaylorUnrolledKernel(FLOAT3 *particles, int n, double3 center, double3 *Q, float3 *acc);

        ''')

    for order in orders:
        w(f'''
            template <>
            void TaylorUnrolledKernel<{order}>(FLOAT3 *particles, int n, double3 center, double3 *Q, float3 *acc){{''')
        w.indent()

        cml_orderm1 = (order)*(order+1)*(order+2)//6

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

    emit_dispatch_function(w, 'TaylorUnrolledKernel(FLOAT3 *particles, int n, double3 center, double3 *Q, float3 *acc)', orders)

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
    emit_unrolled_Taylors(orders)
