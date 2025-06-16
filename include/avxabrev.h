/*
 * Copyright 2012-2025 The Abacus Developers
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#define VZEROALL                   asm("vzeroall");
#define VLOADPS(mem, dst)          asm("vmovaps %0, %" dst::"m"(mem));
#define VLOADPD(mem, dst)          asm("vmovapd %0, %" dst::"m"(mem));
#define VSTORPS(src, mem)          asm("vmovaps %" src ", %0"::"m"(mem));
#define VSTORPD(src, mem)          asm("vmovapd %" src ", %0"::"m"(mem));
#define VMOVAPD(src,dst)           asm("vmovapd " src "," dst );

#define VXORPS( src1,src2, dst)    asm("vxorps " src1 "," src2 "," dst);

#define VADDPS(src1, src2, dst)    asm("vaddps " src1 "," src2 "," dst);
#define VADDPS_M(mem, src, dst)    asm("vaddps %0, %" src ", %" dst""::"m"(mem));
#define VADDPD(src1, src2, dst)    asm("vaddpd " src1 "," src2 "," dst);
#define VADDPD_M(mem, src, dst)    asm("vaddpd %0, %" src ", %" dst""::"m"(mem));

#define VSUBPS(src1, src2, dst)    asm("vsubps " src1 "," src2 "," dst);
#define VSUBPS_M(mem, src, dst)    asm("vsubps %0, %" src ", %" dst""::"m"(mem));
#define VSUBPD(src1, src2, dst)    asm("vsubpd " src1 "," src2 "," dst);
#define VSUBPD_M(mem, src, dst)    asm("vsubpd %0, %" src ", %" dst""::"m"(mem));

#define VMULPS(src1, src2, dst)    asm("vmulps " src1 "," src2 "," dst);
#define VMULPS_M(mem, src, dst)    asm("vmulps %0,  %" src",  %" dst""::"m"(mem));
#define VMULPD(src1, src2, dst)    asm("vmulpd " src1 "," src2 "," dst);
#define VMULPD_M(mem, src, dst)    asm("vmulpd %0,  %" src",  %" dst""::"m"(mem));

#define VHADDPD(src1, src2, dst)   asm("vhaddpd " src1 "," src2 "," dst);
#define VRSQRTPS(src, dst)         asm("vrsqrtps " src "," dst);
#define VCVTPD2PS(src, dst)        asm("vcvtpd2ps " src "," dst);
#define VCVTPS2PD(src, dst)        asm("vcvtps2pd " src "," dst);
#define VMERGE(src1, src2, dst)    asm("vperm2f128 %0,  %" src1 ", %" src2 ", %" dst""::"g"(2));
#define VUP2LOW(src, dst)          asm("vextractf128 %0,  %" src ", %" dst""::"g"(1));
#define VANDPS(src1, src2, dst)    asm("vandps " src1 "," src2 "," dst);
#define VCMPNEQPS(src1, src2, dst) asm("vcmpps %0, %" src1 ", %" src2 ", %" dst""::"g"(4));
#define PREFETCH(mem)              asm("prefetcht0 %0"::"m"(mem))

#define VBROADCASTSD(mem, reg)     asm("vbroadcastsd %0, %" reg::"m"(mem));
#define VBROADCASTSS(mem, reg)     asm("vbroadcastss %0, %" reg::"m"(mem));
#define VSQRTPD(src, dst)          asm("vsqrtpd ", src, "," dst);

#define SWAPHILOW(src1, dest)      asm("vperm2f128 $0x01," src1 "," src1 "," dest );

#define ALIGN16 __attribute__ ((aligned(16)))
#define ALIGN64 __attribute__ ((aligned(64)))
