#!/usr/bin/env python
'''
This program generates assembly code for the Taylors and Multipoles on
IBM POWER architectures using AltiVec/VSX intrinsics.
You shouldn't have to run this program manually; it will
be invoked from the Makefile.

The compiler can't figure out how to unroll some of the nested
loops, so the purpose of this program is to "manually" unroll them.
'''

import textwrap
import argparse

class Writer:
    def __init__(self, fn):
        self.indent_level = 0
        self.fp = open(fn, 'w')

    def __call__(self, txt, *args,**kwargs):
        txt = textwrap.dedent(txt)
        txt = textwrap.indent(txt, ' '*4*self.indent_level)
        self.fp.write(txt, *args,**kwargs)
        self.fp.write('\n')

    def __del__(self):
        self.fp.close()

    def indent(self):
        self.indent_level += 1

    def dedent(self):
        self.indent_level -= 1


def emit_VSX_Multipoles(maxorder=16, onlyorder=None, fn='CM_VSX.cpp'):
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

    for order in range(maxorder+1):
        cml = (order+1)*(order+2)*(order+3)//6

        w(f'''
            template <>
            void MultipoleVSXKernel<{order}>(FLOAT3 *particles, int n, double3 center, double *CM){{''')
        w.indent()

        if order == 0 or (onlyorder is not None and order != onlyorder):
            w('assert(0); // should never be called')
            w.dedent()
            w('}\n')
            continue

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
    w('''
        #undef VSX_NVEC_DOUBLE
        #endif''')  # VSXMULTIPOLES


def emit_VSX_Multipoles_FMA(maxorder=16, onlyorder=None, fn='CM_VSX.cpp', max_zk=None):
    assert max_zk is None, "max_zk implementation not finished"
    if max_zk is None:
        max_zk = maxorder + 1

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

    for order in range(maxorder+1):
        cml = (order+1)*(order+2)*(order+3)//6
        this_max_zk = min(order+1, max_zk)

        w(f'''
            template <>
            void MultipoleVSXKernel<{order}>(FLOAT3 *particles, int n, double3 center, double *CM){{''')
        w.indent()

        if order == 0 or (onlyorder is not None and order != onlyorder):
            w('assert(0); // should never be called')
            w.dedent()
            w('}\n')
            continue

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
    w('''
        #undef VSX_NVEC_DOUBLE
        #endif''')  # VSXMULTIPOLES


def emit_VSX_Taylors(maxorder=16, fn='ET_VSX.cpp', onlyorder=None):
    w = Writer(fn)

    w('''
        #include "threevector.hh"
        #include "header.cpp"
        #ifdef VSXMULTIPOLES
        #include "assert.h"
        extern "C" {
        ''')

    for order in range(maxorder+1):
        w(f'void TaylorVSXKernel{order}(FLOAT3 *particles, int n, double3 center, double3 *Q, float3 *acc){{')
        w.indent()

        if order == 0 or (onlyorder != None and order != onlyorder):
            w('assert(0); // should never be called')
            w.dedent()
            w('}\n')
            continue

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
    w('}')  # extern "C"
    w('#endif')  # VSXMULTIPOLES

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--maxorder', help='Maximum Multipole/Taylor series order to output', type=int, default=16)
    parser.add_argument('--onlyorder', help='Stub all but the given order (useful for fast compilation/debugging)', type=int)
    args = parser.parse_args()
    args = vars(args)

    #emit_VSX_Multipoles(**args)
    emit_VSX_Multipoles_FMA(**args)
    emit_VSX_Taylors(**args)
