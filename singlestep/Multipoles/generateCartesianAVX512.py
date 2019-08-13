#!/usr/bin/env python
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
    import emit_dispatch_function, Writer, default_metacode_argparser

def emit_AVX512_Multipoles(orders, fn='CMAVX512.cpp'):
    w = Writer(fn)

    w('''
        #include "config.h"
        #ifdef AVX512MULTIPOLES
        #include "avx512_calls.h"
        #include "assert.h"
        
        template <int Order>
        void Multipole512Kernel(AVX512_DOUBLES &px, AVX512_DOUBLES &py, AVX512_DOUBLES &pz,
            AVX512_DOUBLES &cx, AVX512_DOUBLES &cy, AVX512_DOUBLES &cz,
            AVX512_DOUBLES *CM );''')

    for order in orders:
        w(f'''
            template <>
            void Multipole512Kernel<{order}>(AVX512_DOUBLES &px, AVX512_DOUBLES &py, AVX512_DOUBLES &pz,
                                     AVX512_DOUBLES &cx, AVX512_DOUBLES &cy, AVX512_DOUBLES &cz,
                                     AVX512_DOUBLES *CM ){{

                AVX512_DOUBLES fi = AVX512_SET_DOUBLE(1.);
                AVX512_DOUBLES fij, fijk;
                AVX512_DOUBLES deltax = AVX512_SUBTRACT_DOUBLES(px, cx);
                AVX512_DOUBLES deltay = AVX512_SUBTRACT_DOUBLES(py, cy);
                AVX512_DOUBLES deltaz = AVX512_SUBTRACT_DOUBLES(pz, cz);
        ''')
        w.indent()

        cml = (order+1)*(order+2)*(order+3)/6

        i = 0
        for a in range(order+1):
            w('fij = fi;')
            for b in range(order-a+1):
                w('fijk = fij;')
                for c in range(order-a-b+1):
                    w('*CM = AVX512_ADD_DOUBLES(*CM, fijk);')
                    if i < cml-1:
                        w('CM++;')
                    i += 1
                    if(c < order-a-b):
                        w('fijk = AVX512_MULTIPLY_DOUBLES(fijk, deltaz);')
                if(b < order-a):
                    w('\nfij = AVX512_MULTIPLY_DOUBLES(fij, deltay);')
            if(a < order):
                w('\nfi = AVX512_MULTIPLY_DOUBLES(fi, deltax);')

        w.dedent()
        w('}')  # Kernel

    emit_dispatch_function(w, '''Multipole512Kernel(
        AVX512_DOUBLES &px, AVX512_DOUBLES &py, AVX512_DOUBLES &pz,
        AVX512_DOUBLES &cx, AVX512_DOUBLES &cy, AVX512_DOUBLES &cz,
        AVX512_DOUBLES *CM)''', orders)

    w('#endif')  # AVX512MULTIPOLES

# this is only marginally faster than the non-FMA version
def emit_AVX512_Multipoles_FMA(orders, fn='CMAVX512.cpp', max_zk=None):
    assert max_zk is None, "max_zk probably not faster; implementation not finished"

    if max_zk is None:
        max_zk = max(orders) + 1

    w = Writer(fn)

    w('''
        #include "config.h"
        #ifdef AVX512MULTIPOLES
        #include "avx512_calls.h"
        #include "assert.h"
        
        template <int Order>
        void Multipole512Kernel(AVX512_DOUBLES &px, AVX512_DOUBLES &py, AVX512_DOUBLES &pz,
            AVX512_DOUBLES &cx, AVX512_DOUBLES &cy, AVX512_DOUBLES &cz,
            AVX512_DOUBLES *CM );
        ''')

    for order in orders:
        this_max_zk = min(order+1, max_zk)

        w(f'''
            template <>
            void Multipole512Kernel<{order}>(AVX512_DOUBLES &px, AVX512_DOUBLES &py, AVX512_DOUBLES &pz,
                                     AVX512_DOUBLES &cx, AVX512_DOUBLES &cy, AVX512_DOUBLES &cz,
                                     AVX512_DOUBLES *CM ){{

                AVX512_DOUBLES fi = AVX512_SET_DOUBLE(1.);
                AVX512_DOUBLES fij;
                AVX512_DOUBLES deltax = AVX512_SUBTRACT_DOUBLES(px, cx);
                AVX512_DOUBLES deltay = AVX512_SUBTRACT_DOUBLES(py, cy);
                AVX512_DOUBLES deltaz = AVX512_SUBTRACT_DOUBLES(pz, cz);
            ''')

        for c in range(this_max_zk):
            w(f'''
                AVX512_DOUBLES zk{c};
                zk0 = AVX512_SET_DOUBLE(1.0);
                ''')
        for c in range(1,this_max_zk):
            w(f'zk{c} = AVX512_MULTIPLY_DOUBLES(zk{c-1}, deltaz);')

        cml = (order+1)*(order+2)*(order+3)//6

        i = 0
        for a in range(order+1):
            w('fij = fi;')
            for b in range(order-a+1):
                for c in range(order-a-b+1):
                    if c % this_max_zk == 0 and c > 0:
                        w(f'''
                            fij = AVX512_MULTIPLY_DOUBLES(fij, zk{this_max_zk-1});
                            CM[{i}] = AVX512_FMA_ADD_DOUBLES(fij, zk{c % this_max_zk}, CM[{i}]);
                            ''')
                    i += 1
                if(b < order-a):
                    w('\nfij = AVX512_MULTIPLY_DOUBLES(fij, deltay);')
            if(a < order):
                w('\nfi = AVX512_MULTIPLY_DOUBLES(fi, deltax);')

        w.dedent()
        w('}')  # Kernel

    emit_dispatch_function(w, '''Multipole512Kernel(
        AVX512_DOUBLES &px, AVX512_DOUBLES &py, AVX512_DOUBLES &pz,
        AVX512_DOUBLES &cx, AVX512_DOUBLES &cy, AVX512_DOUBLES &cz,
        AVX512_DOUBLES *CM)''', orders)

    w('#endif')  # AVX512MULTIPOLES


def emit_AVX512_Taylors(orders, fn='ETAVX512.cpp'):
    w = Writer(fn)

    w('''
        #include "config.h"
        
        #ifdef AVX512MULTIPOLES
        
        #include "avx512_calls.h"
        #include "assert.h"

        template <int Order>
        void Taylor512Kernel(AVX512_DOUBLES &px, AVX512_DOUBLES &py, AVX512_DOUBLES &pz,
                             AVX512_DOUBLES &cx, AVX512_DOUBLES &cy, AVX512_DOUBLES &cz,
                             AVX512_DOUBLES *tx, AVX512_DOUBLES *ty, AVX512_DOUBLES *tz,
                             AVX512_DOUBLES &ax, AVX512_DOUBLES &ay, AVX512_DOUBLES &az);
        ''')

    for order in orders:
        w(f'''
            template <>
            void Taylor512Kernel<{order}>(AVX512_DOUBLES &px, AVX512_DOUBLES &py, AVX512_DOUBLES &pz,
                                     AVX512_DOUBLES &cx, AVX512_DOUBLES &cy, AVX512_DOUBLES &cz,
                                     AVX512_DOUBLES *tx, AVX512_DOUBLES *ty, AVX512_DOUBLES *tz,
                                     AVX512_DOUBLES &ax, AVX512_DOUBLES &ay, AVX512_DOUBLES &az) {{

                AVX512_DOUBLES fi = AVX512_SET_DOUBLE(1.);
                AVX512_DOUBLES fij, fijk;
                AVX512_DOUBLES deltax = AVX512_SUBTRACT_DOUBLES(px, cx);
                AVX512_DOUBLES deltay = AVX512_SUBTRACT_DOUBLES(py, cy);
                AVX512_DOUBLES deltaz = AVX512_SUBTRACT_DOUBLES(pz, cz);
            ''')
        w.indent()

        # taylors use a=0..order, multipoles use a=0..order+1 for some reason
        cml_orderm1 = (order)*(order+1)*(order+2)/6
        i = 0
        w('int i = 0;')
        for a in range(order):
            w('fij = fi;')
            for b in range(order-a):
                w('fijk = fij;')
                for c in range(order-a-b):
                    w('''
                        ax = AVX512_FMA_ADD_DOUBLES(tx[i], fijk, ax);
                        ay = AVX512_FMA_ADD_DOUBLES(ty[i], fijk, ay);
                        az = AVX512_FMA_ADD_DOUBLES(tz[i], fijk, az);
                        ''')
                    if i < cml_orderm1-1:
                        w('i++;')
                    i += 1
                    if(c < order-a-b-1):
                        w('fijk = AVX512_MULTIPLY_DOUBLES(fijk, deltaz);')
                if(b < order-a-1):
                    w('\nfij = AVX512_MULTIPLY_DOUBLES(fij, deltay);')
            if(a < order-1):
                w('\nfi = AVX512_MULTIPLY_DOUBLES(fi, deltax);')

        w.dedent()
        w('}')  # Taylor512Kernel

    emit_dispatch_function(w, '''Taylor512Kernel(
        AVX512_DOUBLES &px, AVX512_DOUBLES &py, AVX512_DOUBLES &pz,
        AVX512_DOUBLES &cx, AVX512_DOUBLES &cy, AVX512_DOUBLES &cz,
        AVX512_DOUBLES *tx, AVX512_DOUBLES *ty, AVX512_DOUBLES *tz,
        AVX512_DOUBLES &ax, AVX512_DOUBLES &ay, AVX512_DOUBLES &az)''', orders)

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
