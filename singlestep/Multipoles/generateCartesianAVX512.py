#!/usr/bin/env python
'''
This program generates AVX-512 assembly for Taylors and Multipoles
evaluation. You shouldn't have to run this program manually; it will
be invoked from the Makefile.

The compiler can't figure out how to unroll some of the nested
loops, so the purpose of this program is to "manually" unroll them.

Curiously, the "un-unrolled" version of the Taylors seem faster.
That's implemented in EvaluateTaylors.cpp.
'''

def emit_AVX512_Multipoles(maxorder=16, fn='CMAVX512.cpp'):
    with open(fn, 'w') as fp:
        def w(txt, *args,**kwargs):
            d = {}
            try:
                d['order'] = order
            except NameError:
                pass
            txt = txt.format(**d)
            fp.write(txt, *args,**kwargs)
            fp.write('\n')

        w('#include "config.h"')
        w('#ifdef AVX512MULTIPOLES')
        w('#include "avx512_calls.h"')
        w('#include "assert.h"')
        w('''
            template <int Order>
            void Multipole512Kernel(AVX512_DOUBLES &px, AVX512_DOUBLES &py, AVX512_DOUBLES &pz,
                AVX512_DOUBLES &cx, AVX512_DOUBLES &cy, AVX512_DOUBLES &cz,
                AVX512_DOUBLES *CM );''')

        for order in range(maxorder+1):
            w('void Multipole512Kernel<{order}>(AVX512_DOUBLES &px, AVX512_DOUBLES &py, AVX512_DOUBLES &pz,')
            w('                         AVX512_DOUBLES &cx, AVX512_DOUBLES &cy, AVX512_DOUBLES &cz,')
            w('                         AVX512_DOUBLES *CM ){{')

            if order == 0:
                w('    assert(0); // should never be called\n}} \n')
                continue

            w('    AVX512_DOUBLES fi = AVX512_SET_DOUBLE(1.);')
            w('    AVX512_DOUBLES fij, fijk;')
            w('    AVX512_DOUBLES deltax = AVX512_SUBTRACT_DOUBLES(px, cx);')
            w('    AVX512_DOUBLES deltay = AVX512_SUBTRACT_DOUBLES(py, cy);')
            w('    AVX512_DOUBLES deltaz = AVX512_SUBTRACT_DOUBLES(pz, cz);')

            cml = (order+1)*(order+2)*(order+3)/6

            i = 0
            for a in range(order+1):
                w('        fij = fi;')
                for b in range(order-a+1):
                    w('        fijk = fij;')
                    for c in range(order-a-b+1):
                        w('        *CM = AVX512_ADD_DOUBLES(*CM, fijk);')
                        if i < cml-1:
                            w('        CM++;')
                        i += 1
                        if(c < order-a-b):
                            w('        fijk = AVX512_MULTIPLY_DOUBLES(fijk, deltaz);')
                    if(b < order-a):
                        w('\n        fij = AVX512_MULTIPLY_DOUBLES(fij, deltay);')
                if(a < order):
                    w('\n        fi = AVX512_MULTIPLY_DOUBLES(fi, deltax);')

            w('}}')  # Kernel
        w('#endif')  # AVX512MULTIPOLES

# this is only marginally faster than the non-FMA version
def emit_AVX512_Multipoles_FMA(maxorder=16, fn='CMAVX512.cpp', max_zk=None):
    assert max_zk is None, "max_zk probably not faster; implementation not finished"

    if max_zk is None:
        max_zk = maxorder + 1

    with open(fn, 'w') as fp:
        def w(txt, *args,**kwargs):
            d = {}
            try:
                d['order'] = order
            except NameError:
                pass
            txt = txt.format(**d)
            fp.write(txt, *args,**kwargs)
            fp.write('\n')

        w('#include "config.h"')
        w('#ifdef AVX512MULTIPOLES')
        w('#include "avx512_calls.h"')
        w('#include "assert.h"')
        w('''
            template <int Order>
            void Multipole512Kernel(AVX512_DOUBLES &px, AVX512_DOUBLES &py, AVX512_DOUBLES &pz,
                AVX512_DOUBLES &cx, AVX512_DOUBLES &cy, AVX512_DOUBLES &cz,
                AVX512_DOUBLES *CM );''')

        for order in range(maxorder+1):
            this_max_zk = min(order+1, max_zk)

            w('void Multipole512Kernel<{order}>(AVX512_DOUBLES &px, AVX512_DOUBLES &py, AVX512_DOUBLES &pz,')
            w('                         AVX512_DOUBLES &cx, AVX512_DOUBLES &cy, AVX512_DOUBLES &cz,')
            w('                         AVX512_DOUBLES *CM ){{')

            if order == 0:
                w('    assert(0); // should never be called\n}} \n')
                continue

            w('    AVX512_DOUBLES fi = AVX512_SET_DOUBLE(1.);')
            w('    AVX512_DOUBLES fij;')
            w('    AVX512_DOUBLES deltax = AVX512_SUBTRACT_DOUBLES(px, cx);')
            w('    AVX512_DOUBLES deltay = AVX512_SUBTRACT_DOUBLES(py, cy);')
            w('    AVX512_DOUBLES deltaz = AVX512_SUBTRACT_DOUBLES(pz, cz);')

            for c in range(this_max_zk):
                w('    AVX512_DOUBLES zk{c};'.format(c=c))
            w('    zk0 = AVX512_SET_DOUBLE(1.0);')
            for c in range(1,this_max_zk):
                w('    zk{c} = AVX512_MULTIPLY_DOUBLES(zk{cm1}, deltaz);'.format(c=c,cm1=c-1))

            cml = (order+1)*(order+2)*(order+3)/6

            i = 0
            for a in range(order+1):
                w('        fij = fi;')
                for b in range(order-a+1):
                    for c in range(order-a-b+1):
                        if c % this_max_zk == 0 and c > 0:
                            w('        fij = AVX512_MULTIPLY_DOUBLES(fij, zk{c});'.format(c=this_max_zk-1))
                        w('        CM[{i}] = AVX512_FMA_ADD_DOUBLES(fij, zk{c}, CM[{i}]);'.format(c=c % this_max_zk, i=i))
                        i += 1
                    if(b < order-a):
                        w('\n        fij = AVX512_MULTIPLY_DOUBLES(fij, deltay);')
                if(a < order):
                    w('\n        fi = AVX512_MULTIPLY_DOUBLES(fi, deltax);')

            w('}}')  # Kernel
        w('#endif')  # AVX512MULTIPOLES

def emit_AVX512_Taylors(maxorder=16, fn='ETAVX512.cpp'):
    with open(fn, 'w') as fp:
        def w(txt, *args,**kwargs):
            d = {}
            try:
                d['order'] = order
            except NameError:
                pass
            txt = txt.format(**d)
            fp.write(txt, *args,**kwargs)
            fp.write('\n')

        w('#include "config.h"')
        w('#ifdef AVX512MULTIPOLES')
        w('#include "avx512_calls.h"')
        w('#include "assert.h"')
        w('extern "C" {{')

        for order in range(maxorder+1):
            w('void Taylor512Kernel{order}(AVX512_DOUBLES &px, AVX512_DOUBLES &py, AVX512_DOUBLES &pz,')
            w('                         AVX512_DOUBLES &cx, AVX512_DOUBLES &cy, AVX512_DOUBLES &cz,')
            w('                         AVX512_DOUBLES *tx, AVX512_DOUBLES *ty, AVX512_DOUBLES *tz,')
            w('                         AVX512_DOUBLES &ax, AVX512_DOUBLES &ay, AVX512_DOUBLES &az) {{')

            if order == 0:
                w('    assert(0); // should never be called\n}} \n')
                continue

            w('    AVX512_DOUBLES fi = AVX512_SET_DOUBLE(1.);')
            w('    AVX512_DOUBLES fij, fijk;')
            w('    AVX512_DOUBLES deltax = AVX512_SUBTRACT_DOUBLES(px, cx);')
            w('    AVX512_DOUBLES deltay = AVX512_SUBTRACT_DOUBLES(py, cy);')
            w('    AVX512_DOUBLES deltaz = AVX512_SUBTRACT_DOUBLES(pz, cz);')

            # taylors use a=0..order, multipoles use a=0..order+1 for some reason
            cml_orderm1 = (order)*(order+1)*(order+2)/6
            i = 0
            w('    int i = 0;')
            for a in range(order):
                w('        fij = fi;')
                for b in range(order-a):
                    w('        fijk = fij;')
                    for c in range(order-a-b):
                        w('        ax = AVX512_FMA_ADD_DOUBLES(tx[i], fijk, ax);')
                        w('        ay = AVX512_FMA_ADD_DOUBLES(ty[i], fijk, ay);')
                        w('        az = AVX512_FMA_ADD_DOUBLES(tz[i], fijk, az);')
                        if i < cml_orderm1-1:
                            w('        i++;')
                        i += 1
                        if(c < order-a-b-1):
                            w('        fijk = AVX512_MULTIPLY_DOUBLES(fijk, deltaz);')
                    if(b < order-a-1):
                        w('\n        fij = AVX512_MULTIPLY_DOUBLES(fij, deltay);')
                if(a < order-1):
                    w('\n        fi = AVX512_MULTIPLY_DOUBLES(fi, deltax);')

            w('}}')  # Taylor512Kernel
        w('}}')  # extern "C"
        w('#endif')  # AVX512MULTIPOLES


if __name__ == '__main__':
    emit_AVX512_Taylors()
    #emit_AVX512_Multipoles()
    emit_AVX512_Multipoles_FMA()
