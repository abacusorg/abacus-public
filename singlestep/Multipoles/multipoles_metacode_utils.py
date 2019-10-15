'''
This module contains a few functions related to the loop-unrolling
meta-code for the multipoles/Taylors (generateCartesianUnrolled.py
and the like).
'''

import textwrap
import argparse
import re

import numpy as np

from Abacus.Tools import ArgParseFormatter

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

    def indent(self, n=1):
        self.indent_level += n

    def dedent(self, n=1):
        self.indent_level -= n


def emit_dispatch_function(w, func_sig, orders):
    # Parse the function signature
    sig_re = re.compile(r'(?P<func>\w+)\s*\((?P<args>(.|\n)+)\)')
    match = sig_re.match(func_sig)
    name = match.group('func')
    args = match.group('args')
    argnames = [a.strip() for a in args.split(',')]
    name_re = re.compile(r'\w+\s*\W*\s*(?P<name>\w+)')
    matches = [name_re.match(a) for a in argnames]
    argnames = [m.group('name') for m in matches]
    
    argnames = ', '.join(argnames)

    w(f'''
        void Dispatch{name}(int order, {args}){{
            switch(order){{''')
    w.indent(2)
    for order in orders:
        w(f'''
            case {order}:
                {name}<{order}>({argnames});
                break;''')
    w(f'''
        default:
            assert(0 && "No kernel for order >{max(orders)}.  Try ./configure --with-max-order=MAXORDER?");
            break;''')
    w.dedent(2)
    w('    }\n}')


def default_metacode_argparser(doc):
    parser = argparse.ArgumentParser(description=doc, formatter_class=ArgParseFormatter)
    parser.add_argument('--maxorder', help='Maximum Multipole/Taylor series order to output', type=int, default=8)
    parser.add_argument('--onlyorder', help='Stub all but the given order (useful for fast compilation/debugging)', type=int)

    return parser


class cmapper:
    def __init__(self, order):
        self.order = order
        self.cml = (order+1)*(order+2)*(order+3)//6
        self._cmap = np.empty((order+1)**3, dtype=int)
        self._lmap = lambda a,b,c: (order+1)**2*a + (order+1)*b + c

        i = 0
        for a in range(order+1):
            for b in range(order+1-a):
                for c in range(order+1-a-b):
                    self._cmap[self._lmap(a,b,c)] = i
                    i += 1

    def __call__(self, a, b, c):
        return self._cmap[self._lmap(a,b,c)]
