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


def emit_VSX_Multipoles_FMA_interleaved(orders, fn='CM_VSX.cpp', use_zk=False, verbose=False, unroll=4, explicit_fma=True, nCMvec=1):
    if use_zk in (None, True):
        use_zk = max(orders) + 1

    nCMvec = eval(str(nCMvec))
    VSX_NVEC_DOUBLE = 2

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
        this_max_zk = min(order+1, use_zk)

        w(f'''
            template <>
            void MultipoleVSXKernel<{order}>(FLOAT3 *particles, int n, double3 center, double *CM){{''')
        w.indent()

        w(f'''
            VSX_DOUBLES cx = vec_splats(center.x);
            VSX_DOUBLES cy = vec_splats(center.y);
            VSX_DOUBLES cz = vec_splats(center.z);

            VSX_DOUBLES CMvec[{cml}][{nCMvec}];''')

        w(f'for(int k = 0; k < {cml}; k++){{')
        w.indent()
        for m in range(nCMvec):
            w(f'CMvec[k][{m}] = vec_splats(0.);')
        w.dedent()
        w('}')

        w(f'''
            int n_left = n % ({unroll}*VSX_NVEC_DOUBLE);
            int n_aligned = n - n_left;
            int odd = n_left % 2;
            n_left -= odd;
            ''')

        w(f'for(int i = 0; i < n_aligned; i += {unroll}*VSX_NVEC_DOUBLE){{')
        w.indent()

        def _emit_unrolled(start, unroll):
            w(f'''
                VSX_DOUBLES px[{unroll}], py[{unroll}], pz[{unroll}];

                for(int j = 0; j < VSX_NVEC_DOUBLE; j++){{''')
            w.indent()

            for m in range(unroll):
                w(f'px[{m}][j] = particles[{start}+j+{m}*VSX_NVEC_DOUBLE].x;')
            for m in range(unroll):
                w(f'py[{m}][j] = particles[{start}+j+{m}*VSX_NVEC_DOUBLE].y;')
            for m in range(unroll):
                w(f'pz[{m}][j] = particles[{start}+j+{m}*VSX_NVEC_DOUBLE].z;')

            w.dedent()
            w('}')

            w(f'''
                VSX_DOUBLES fi[{unroll}], fij[{unroll}];
                VSX_DOUBLES deltax[{unroll}], deltay[{unroll}], deltaz[{unroll}];
                ''')
            if use_zk:
                w(f'VSX_DOUBLES zk[{this_max_zk}][{unroll}];')
            else:
                w(f'VSX_DOUBLES fijk[{unroll}];')

            for m in range(unroll):
                w(f'''
                    fi[{m}] = vec_splats(1.);

                    deltax[{m}] = vec_sub(px[{m}], cx);
                    deltay[{m}] = vec_sub(py[{m}], cy);
                    deltaz[{m}] = vec_sub(pz[{m}], cz);
                    ''')

            if use_zk:
                # precompute powers of deltaz^k
                for m in range(unroll):
                    w(f'zk[0][{m}] = vec_splats(1.);')

                for c in range(1,this_max_zk):
                    for m in range(unroll):
                        w(f'zk[{c}][{m}] = vec_mul(deltaz[{m}], zk[{c-1}][{m}]);')
                w('\n')

            i = 0
            nflop = 3 + this_max_zk-1
            nflop_fma = 0
            for a in range(order+1):
                for m in range(unroll):
                    w(f'fij[{m}] = fi[{m}];')
                for b in range(order-a+1):
                    if not use_zk:
                        for m in range(unroll):
                            w(f'fijk[{m}] = fij[{m}];')
                    for c in range(order-a-b+1):
                        for m in range(unroll):

                            if use_zk:
                                zkterms = []
                                _c = c
                                while _c > 0:
                                    step = min(_c, this_max_zk-1)
                                    zkterms += [f'zk[{step}][{m}]']
                                    _c -= step

                                if len(zkterms) > 1:
                                    zkterm = ''.join(f'vec_mul({_zk}, ' for _zk in zkterms[:-1])
                                    zkterm += zkterms[-1] + ')'*(len(zkterms)-1)
                                elif len(zkterms) > 0:
                                    zkterm = zkterms[0]
                                #zkterm = '*'.join(zkterms)

                                if c == 0:
                                    # No sense in generating instructions to multiply by 1!
                                    if explicit_fma:
                                        w(f'CMvec[{i}][{m%nCMvec}] = vec_add(fij[{m}], CMvec[{i}][{m%nCMvec}]);')
                                    else:
                                        w(f'CMvec[{i}][{m%nCMvec}] += fij[{m}];')
                                else:
                                    if explicit_fma:
                                        w(f'CMvec[{i}][{m%nCMvec}] = vec_madd(fij[{m}], {zkterm}, CMvec[{i}][{m%nCMvec}]);')
                                    else:
                                        w(f'CMvec[{i}][{m%nCMvec}] += fij[{m}]*{zkterm};')
                                nflop += 2
                                nflop_fma += 2
                            else:
                                w(f'CMvec[{i}][{m%nCMvec}] += fijk[{m}];')

                        if not use_zk:
                            if c < order-a-b:
                                for m in range(unroll):
                                    w(f'fijk[{m}] = vec_mul(deltaz[{m}], fijk[{m}]);')
                        i += 1
                            
                    if b < order-a:
                        for m in range(unroll):
                            w(f'fij[{m}] = vec_mul(deltay[{m}], fij[{m}]);')
                        nflop += 1
                if a < order:
                    for m in range(unroll):
                        w(f'fi[{m}] = vec_mul(deltax[{m}], fi[{m}]);')
                    nflop += 1
                    w('\n')

        _emit_unrolled('i', unroll)

        w.dedent()
        w('}')  # particle loop

        if verbose:
            print(f'nflop: {nflop} ({nflop_fma/nflop*100:.3g}% FMA); cml: {cml}')

        # Generate unrolled versions of 2, 4, 6 particles
        w('switch(n_left){')
        w.indent()
        for m in range(1,unroll):
            w(f'case {VSX_NVEC_DOUBLE*m}: {{')
            w.indent()
            _emit_unrolled('n_aligned', m)
            w('break;\n}')
            w.dedent()
        w.dedent()
        w('}\n')  # switch

        # Now the last particle
        w(f'''
            if(odd){{
                double fi = 1., fij, fijk;
                double deltax = particles[n-1].x - center.x;
                double deltay = particles[n-1].y - center.y;
                double deltaz = particles[n-1].z - center.z;
            ''')
        w.indent()

        if use_zk:
            # precompute powers of deltaz^k
            w(f'''
                double zk[{this_max_zk}];
                zk[0] = 1.;
                ''')
            for c in range(1,this_max_zk):
                w(f'zk[{c}] = deltaz*zk[{c-1}];')
            w('\n')

        i = 0
        for a in range(order+1):
            w(f'fij = fi;')
            for b in range(order-a+1):
                if not use_zk:
                    w(f'fijk = fij;')
                for c in range(order-a-b+1):
                    if use_zk:
                        zkterms = []
                        _c = c
                        while _c > 0:
                            step = min(_c, this_max_zk-1)
                            zkterms += [f'zk[{step}]']
                            _c -= step
                        zkterm = '*'.join(zkterms)

                        if c == 0:
                            # No sense in generating instructions to multiply by 1!
                            w(f'CM[{i}] += fij;')
                        else:
                            w(f'CM[{i}] += fij*{zkterm};')
                    else:
                        w(f'CM[{i}] += fijk;')

                    if not use_zk:
                        if c < order-a-b:
                            w(f'fijk *= deltaz;')
                    i += 1
                        
                if b < order-a:
                    w(f'fij *= deltay;')
            if a < order:
                w(f'fi *= deltax;')
                w('\n')

        w('}')  # remainder loop

        w(f'''
            for(int k = 0; k < {cml}; k++)
                for(int j = 0; j < VSX_NVEC_DOUBLE; j++)''')
        w.indent(2)
        w('CM[k] += ' + '+'.join(f'CMvec[k][{m}][j]' for m in range(nCMvec)) + ';')
        w.dedent(2)

        w.dedent()
        w('}\n')  # Kernel

    emit_dispatch_function(w, 'MultipoleVSXKernel(FLOAT3 *particles, int n, double3 center, double *CM)', orders)

    w('''
        #undef VSX_NVEC_DOUBLE
        #endif''')  # VSXMULTIPOLES


def emit_VSX_Taylors_interleaved(orders, fn='ET_VSX.cpp', verbose=False, unroll=4, explicit_fma=True, noQ=False):
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
        nflop, nflop_fma, nflop_fixed = 0,0,0
        VSX_NVEC_DOUBLE = 2

        cmap = cmapper(order)

        w(f'''
            template <>
            void TaylorVSXKernel<{order}>(FLOAT3 *particles, int n, double3 center, double *CT, float3 *acc){{
            ''')
        w.indent()

        if not noQ:
            for m in range(unroll):
                w(f'VSX_DOUBLES Qx{m}[{cml_orderm1}], Qy{m}[{cml_orderm1}], Qz{m}[{cml_orderm1}];')

            # Precompute Qxyz
            i = 0
            for a in range(order):
                for b in range(order-a):
                    for c in range(order-a-b):
                        w(f'''
                            Qx0[{i}] = vec_splats({a+1}*CT[{cmap(a+1, b  , c  )}]);
                            Qy0[{i}] = vec_splats({b+1}*CT[{cmap(a  , b+1, c  )}]);
                            Qz0[{i}] = vec_splats({c+1}*CT[{cmap(a  , b  , c+1)}]);
                            ''')
                        for m in range(1,unroll):
                            w(f'''
                                Qx{m}[{i}] = Qx0[{i}];
                                Qy{m}[{i}] = Qy0[{i}];
                                Qz{m}[{i}] = Qz0[{i}];
                                ''')
                        i += 1
                        nflop_fixed += 3

        # Now compute the accelerations
        w(f'''
                VSX_DOUBLES cx = vec_splats(center.x);
                VSX_DOUBLES cy = vec_splats(center.y);
                VSX_DOUBLES cz = vec_splats(center.z);

                int n_left = n % ({unroll}*VSX_NVEC_DOUBLE);
                int n_aligned = n - n_left;
                int odd = n_left % 2;
                n_left -= odd;

                for(int i = 0; i < n_aligned; i += {unroll}*VSX_NVEC_DOUBLE){{
            ''')
        w.indent()

        def _emit_unrolled(start, unroll):
            nflop, nflop_fma = 0,0

            w(f'VSX_DOUBLES px[{unroll}],py[{unroll}],pz[{unroll}];')

            w(f'for(int j = 0; j < VSX_NVEC_DOUBLE; j++){{')
            w.indent()
            # could initialize better? does it matter?
            for m in range(unroll):
                w(f'px[{m}][j] = particles[{start}+j+{m}*VSX_NVEC_DOUBLE].x;')
            for m in range(unroll):
                w(f'py[{m}][j] = particles[{start}+j+{m}*VSX_NVEC_DOUBLE].y;')
            for m in range(unroll):
                w(f'pz[{m}][j] = particles[{start}+j+{m}*VSX_NVEC_DOUBLE].z;')
            w.dedent()
            w(f'}}')

            w(f'VSX_DOUBLES fi[{unroll}], fij[{unroll}], fijk[{unroll}];')
            for m in range(unroll):
                w(f'fi[{m}] = vec_splats(1.);')

            w(f'VSX_DOUBLES deltax[{unroll}], deltay[{unroll}], deltaz[{unroll}];')

            for m in range(unroll):
                w(f'''
                    deltax[{m}] = vec_sub(px[{m}], cx);
                    deltay[{m}] = vec_sub(py[{m}], cy);
                    deltaz[{m}] = vec_sub(pz[{m}], cz);
                    ''')
            nflop += 3

            w(f'VSX_DOUBLES ax[{unroll}], ay[{unroll}], az[{unroll}];')
            for m in range(unroll):
                w(f'''
                    ax[{m}] = vec_splats(0.);
                    ay[{m}] = vec_splats(0.);
                    az[{m}] = vec_splats(0.);
                    ''')  # Ccould make first loop set instead of accumulate, but doesn't seem to matter
            
            i = 0
            for a in range(order):
                for m in range(unroll):
                    w(f'fij[{m}] = fi[{m}];')
                for b in range(order-a):
                    for m in range(unroll):
                        w(f'fijk[{m}] = fij[{m}];')
                    for c in range(order-a-b):
                        if explicit_fma:
                            if not noQ:
                                for m in range(unroll):
                                    w(f'''
                                        ax[{m}] = vec_nmsub(Qx{m}[{i}], fijk[{m}], ax[{m}]);
                                        ay[{m}] = vec_nmsub(Qy{m}[{i}], fijk[{m}], ay[{m}]);
                                        az[{m}] = vec_nmsub(Qz{m}[{i}], fijk[{m}], az[{m}]);
                                        ''')
                                nflop_fma += 6
                                nflop += 6
                            else:
                                for m in range(unroll):
                                    w(f'''
                                        ax[{m}] = vec_nmsub(vec_splats({a+1}*CT[{cmap(a+1, b  , c  )}]), fijk[{m}], ax[{m}]);
                                        ay[{m}] = vec_nmsub(vec_splats({b+1}*CT[{cmap(a  , b+1, c  )}]), fijk[{m}], ay[{m}]);
                                        az[{m}] = vec_nmsub(vec_splats({c+1}*CT[{cmap(a  , b  , c+1)}]), fijk[{m}], az[{m}]);
                                        ''')
                                nflop_fma += 6
                                nflop += 7.5
                        else:
                            raise NotImplementedError
                            #w(f'''
                            #    ax -= Qx[{i}] * fijk;
                            #    ax2 -= Qx2[{i}] * fijk2;
                            #    ay -= Qy[{i}] * fijk;
                            #    ay2 -= Qy2[{i}] * fijk2;
                            #    az -= Qz[{i}] * fijk;
                            #    az2 -= Qz2[{i}] * fijk2;
                            #    ''')
                        i += 1
                        if c < order-a-b-1:
                            for m in range(unroll):
                                w(f'fijk[{m}] = vec_mul(deltaz[{m}], fijk[{m}]);')
                            nflop += 1
                    if b < order-a-1:
                        for m in range(unroll):
                            w(f'fij[{m}] = vec_mul(deltay[{m}], fij[{m}]);')
                        nflop += 1
                if a < order-1:
                    w('\n')
                    for m in range(unroll):
                        w(f'fi[{m}] = vec_mul(deltax[{m}], fi[{m}]);')
                    nflop += 1

            w(f'for(int j = 0; j < VSX_NVEC_DOUBLE; j++){{')
            w.indent()
            for m in range(unroll):
                w(f'acc[{start}+j+{m}*VSX_NVEC_DOUBLE].x = ax[{m}][j];')
            for m in range(unroll):
                w(f'acc[{start}+j+{m}*VSX_NVEC_DOUBLE].y = ay[{m}][j];')
            for m in range(unroll):
                w(f'acc[{start}+j+{m}*VSX_NVEC_DOUBLE].z = az[{m}][j];')
            
            w.dedent()
            w('}')

            return nflop, nflop_fma

        nflop, nflop_fma = _emit_unrolled('i', unroll)

        w.dedent()
        w('}\n')  # particle loop

        # Generate unrolled versions of 2, 4, 6 particles
        w('switch(n_left){')
        w.indent()
        for m in range(1,unroll):
            w(f'case {VSX_NVEC_DOUBLE*m}: {{')
            w.indent()
            _emit_unrolled('n_aligned', m)
            w('break;\n}')
            w.dedent()
        w.dedent()
        w('}\n')  # switch

        # Now the last particle
        w(f'''
            if(odd){{
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

        if verbose:
            print(f'Order {order}: {nflop} flop per particle ({nflop_fma/nflop*100:.3g}% FMA, {nflop_fixed} per cell)')

    emit_dispatch_function(w, 'TaylorVSXKernel(FLOAT3 *particles, int n, double3 center, double *CT, float3 *acc)', orders)

    w('#endif')  # VSXMULTIPOLES

if __name__ == '__main__':
    parser = default_metacode_argparser(doc=__doc__)
    args = parser.parse_args()

    if args.onlyorder:
        orders = [args.onlyorder]
    else:
        orders = list(range(1,args.maxorder+1))

    emit_VSX_Multipoles_FMA_interleaved(orders, verbose=args.verbose)
    emit_VSX_Taylors_interleaved(orders, verbose=args.verbose, noQ=True)
