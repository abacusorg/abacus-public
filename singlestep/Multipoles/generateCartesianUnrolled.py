#!/usr/bin/env python
'''
This program generates plain C code for the Taylors and Multipoles
evaluation. You shouldn't have to run this program manually; it will
be invoked from the Makefile.

The compiler can't figure out how to unroll some of the nested
loops, so the purpose of this program is to "manually" unroll them.
'''

import textwrap

class Writer:
    def __init__(self, fn):
        self.indent_level = 0
        self.fp = open(fn, 'w')

    def __call__(self, txt, *args,**kwargs):
        txt = textwrap.dedent(txt)
        txt = textwrap.indent(txt, '    '*self.indent_level)
        #txt = txt.format(**d)
        self.fp.write(txt, *args,**kwargs)
        self.fp.write('\n')

    def __del__(self):
        self.fp.close()

    def indent(self):
        self.indent_level += 1

    def dedent(self):
        self.indent_level -= 1


def emit_unrolled_Multipoles(maxorder=16, fn='CM_unrolled.cpp'):
    w = Writer(fn)

    w('''
        #include "threevector.hh"
        #include "header.cpp"
        #ifdef UNROLLEDMULTIPOLES
        #include "assert.h"
        extern "C" {
        ''')

    for order in range(maxorder+1):
        w('void MultipoleUnrolledKernel{order}(FLOAT3 p, double3 center, double *CM){{'.format(order=order))
        w.indent()

        if order == 0:
            w('assert(0); // should never be called')
            w.dedent()
            w('}\n')
            continue

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

        cml = (order+1)*(order+2)*(order+3)//6

        i = 0
        for a in range(order+1):
            w('fij = fi;')
            for b in range(order-a+1):
                w('fijk = fij;')
                for c in range(order-a-b+1):
                    w('CM[{i}] += fijk;'.format(i=i))
                    i += 1
                    if(c < order-a-b):
                        w('fijk *= deltaz;')
                if(b < order-a):
                    w('fij *= deltay;')
            if(a < order):
                w('\nfi *= deltax;')

        w.dedent()
        w('}\n')  # Kernel
    w('}')  # extern "C"
    w('#endif')  # UNROLLEDMULTIPOLES


# try bringing particle loop inside?
def emit_unrolled_Taylors(maxorder=16, fn='ET_unrolled.cpp'):
    w = Writer(fn)

    w('''
        #include "threevector.hh"
        #include "header.cpp"
        #ifdef UNROLLEDMULTIPOLES
        #include "assert.h"
        extern "C" {
        ''')

    for order in range(maxorder+1):
        w('void TaylorUnrolledKernel{order}(FLOAT3 *particles, int n, double3 center, double3 *Q, float3 *acc){{'.format(order=order))
        w.indent()

        if order == 0:
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
                    w('''
                        thisacc -= Q[{i}] * fijk;
                        '''.format(i=i))
                    i += 1
                    if(c < order-a-b-1):
                        w('fijk *= deltaz;')
                if(b < order-a-1):
                    w('fij *= deltay;')
            if(a < order-1):
                w('\nfi *= deltax;')

        w('acc[j] = thisacc;')
        w.dedent()
        w('}\n')  # particle loop

        w.dedent()
        w('}\n')  # Kernel
    w('}')  # extern "C"
    w('#endif')  # UNROLLEDMULTIPOLES

if __name__ == '__main__':
    emit_unrolled_Multipoles()
    emit_unrolled_Taylors()
