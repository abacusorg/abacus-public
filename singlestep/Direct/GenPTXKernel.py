#!/usr/bin/env python

#see ASMKernel.ptx for a code breakdown/explanation

import sys


if len(sys.argv) != 4:
    print("Usage: {} [n_sources] [ilp] [unroll]".format(sys.argv[0]))
    sys.exit(1)

n_sources = int(sys.argv[1])
ilp = int(sys.argv[2])
unroll = int(sys.argv[3])

preamble = '''
.version 4.2
.target sm_35
.address_size 64

.visible .func FullDirectTile (
                        .param .b64 pcache_x,
                        .param .b64 pcache_y,
                        .param .b64 pcache_z,
                        .param .b64 pa_x,
                        .param .b64 pa_y,
                        .param .b64 pa_z,
                        .param .b64 psink_x,
                        .param .b64 psink_y,
                        .param .b64 psink_z,
                        .param .b64 peps2)
{{
'''


single_declarations = '''
.reg .s32 count;

.reg .b64 cache_x;
.reg .b64 cache_y;
.reg .b64 cache_z;

.reg .f32 a_x;
.reg .f32 a_y;
.reg .f32 a_z;

.reg .f32 sink_x;
.reg .f32 sink_y;
.reg .f32 sink_z;

.reg .f32 eps2;

.reg .b64 addr;

.reg .pred loop_predicate;
'''
regcount = 16

prologue = '''
ld.param.b64 addr,[pa_x];
ld.f32 a_x,[addr];

ld.param.b64 addr,[pa_y];
ld.f32 a_x,[addr];

ld.param.b64 addr,[pa_z];
ld.f32 a_x,[addr];

ld.param.b64 cache_x,[pcache_x];
ld.param.b64 cache_y,[pcache_y];
ld.param.b64 cache_z,[pcache_z];

ld.param.b64 addr,[psink_x];
ld.f32 sink_x,[addr];

ld.param.b64 addr,[psink_y];
ld.f32 sink_y,[addr];

ld.param.b64 addr,[psink_z];
ld.f32 sink_z,[addr];

ld.param.b64 addr,[peps2];
ld.f32 eps2,[addr];

mov.u32 count,-{loop_count:d};

DirectKernel:
'''

src_declaration = '''
.reg .f32 src_x_{i:d};
.reg .f32 src_y_{i:d};
.reg .f32 src_z_{i:d};
.reg .f32 r_{i:d};
.reg .f32 r2_{i:d};
'''

regcount += ilp*5

src_load = '''
ld.f32 src_x_{i:d}, [cache_x+{src_offset:d}];
ld.f32 src_y_{i:d}, [cache_y+{src_offset:d}];
ld.f32 src_z_{i:d}, [cache_z+{src_offset:d}];
'''

src_delta = '''
sub.f32 src_x_{i:d}, src_x_{i:d}, sink_x;
sub.f32 src_y_{i:d}, src_y_{i:d}, sink_y;
sub.f32 src_z_{i:d}, src_z_{i:d}, sink_z;
'''

pot1 = "fma.rn.f32 r_{i:d},src_x_{i:d},src_x_{i:d},eps2;\n"
pot2 = "fma.rn.f32 r_{i:d},src_y_{i:d},src_y_{i:d},r_{i:d};\n"
pot3 = "fma.rn.f32 r_{i:d},src_z_{i:d},src_z_{i:d},r_{i:d};\n"
pot4 = "rsqrt.approx.f32 r_{i:d},r_{i:d};\n"
pot5 = "mul.f32 r2_{i:d},r_{i:d},r_{i:d};\n"
pot6 = "mul.f32 r_{i:d},r_{i:d},r2_{i:d};\n"

ax = '''
mul.f32 r2_{i:d},r_{i:d},src_x_{i:d};
sub.f32 a_x,a_x,r2_{i:d};
'''

ay = '''
mul.f32 r2_{i:d},r_{i:d},src_y_{i:d};
sub.f32 a_y,a_y,r2_{i:d};
'''

az = '''
mul.f32 r2_{i:d},r_{i:d},src_z_{i:d};
sub.f32 a_z,a_z,r2_{i:d};
'''

loop = '''
add.s64 cache_x,cache_x,{loop_offset:d};
add.s64 cache_y,cache_y,{loop_offset:d};
add.s64 cache_z,cache_z,{loop_offset:d};

add.s32 count,count,{loopsize:d};
setp.ne.s32 loop_predicate,count,0;
@loop_predicate bra.uni DirectKernel;
'''

epilogue = '''
ld.param.b64 addr,[pa_x];
st.f32 [addr],a_x;

ld.param.b64 addr,[pa_x];
st.f32 [addr],a_x;

ld.param.b64 addr,[pa_x];
st.f32 [addr],a_x;

ret;
}
'''

outfile = open("GeneratedASMKernel.ptx",'wb')

src_per_loop = ilp*unroll

outfile.write(preamble.format(loop_count = n_sources))

outfile.write(single_declarations.format(loop_count = n_sources))

for s in range(ilp):
    outfile.write(src_declaration.format(i = s))

outfile.write(prologue.format(loop_count = n_sources))

for l in range(unroll):
    for s in range(ilp):
        outfile.write(src_load.format(i=s,src_offset = 4*(s + ilp*l)))
    for s in range(ilp):
        outfile.write(src_delta.format(i=s))
    for s in range(ilp):
        outfile.write(pot1.format(i=s))
    for s in range(ilp):
        outfile.write(pot2.format(i=s))
    for s in range(ilp):
        outfile.write(pot3.format(i=s))
    for s in range(ilp):
        outfile.write(pot4.format(i=s))
    for s in range(ilp):
        outfile.write(pot5.format(i=s))
    for s in range(ilp):
        outfile.write(pot6.format(i=s))
    for s in range(ilp):
        outfile.write(ax.format(i=s))
    for s in range(ilp):
        outfile.write(ay.format(i=s))
    for s in range(ilp):
        outfile.write(az.format(i=s))
if src_per_loop != n_sources:
    outfile.write(loop.format(loop_offset = 4*src_per_loop,loopsize=src_per_loop))

outfile.write(epilogue)
outfile.close()
print("Used {:d} registers.".format(regcount))
