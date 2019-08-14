'''
This module contains a few functions related to the loop-unrolling
meta-code for the multipoles/Taylors (generateCartesianUnrolled.py
and the like).
'''

import textwrap
import argparse
import re

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
