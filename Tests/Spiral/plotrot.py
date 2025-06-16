#!/usr/bin/env python3
# Copyright 2012-2025 The Abacus Developers
# SPDX-License-Identifier: GPL-3.0-or-later

import time, sys
from numpy import *
from matplotlib.font_manager import FontProperties
from pylab import *
from math import acos, atan

from matplotlib import rcParams

rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['DroidSans']
rcParams['text.usetex'] = True
rcParams['font.size'] = 12.0
rcParams['font.weight'] = 'bold'
limits = array([[  0,   5, -225, 225],
                [2.3, 2.7, -210, 210],
                [2.47, 2.53,  -37,  37]])


#get figure width from latex...
fig_width = 5
fig_height = 3
fig = plt.figure(figsize=(fig_width,fig_height))
fig.subplots_adjust(bottom=0.2,wspace=0.4)

H0 = 50
fac = 1/sqrt(3)   # norm of k-vector
Lbox = 10/fac     # size of box for one wavelength in diagonal of 10 Mpc box
mid = Lbox*fac/2
print("mid = ", mid)
print(mid-0.4, mid+0.4)

limits = array([[  0,   Lbox*fac, -225, 225],
                [mid-0.05*Lbox*fac, mid+0.05*Lbox*fac, -225, 225],
                [mid-0.007*Lbox*fac, mid+0.007*Lbox*fac,  -37,  37]])

ax = []
ax.append(subplot(131))
filename = "spiral.rot"
ld = loadtxt(filename)
xd = ld[:,0]
vd = ld[:,1]

filename = "analyticspiral"
l = loadtxt(filename)
x = (l[:,0]+0.5)*Lbox/sqrt(3)
v = l[:,1]*H0*Lbox/sqrt(3)

ax[0].plot(x,v,linewidth=1)
ax[0].plot(xd,vd,'r.',linewidth=1,markersize=2)
ax[0].axis(limits[0]);
tickpos = [int(i) for i in Lbox*fac*frange(0,1,.2)]
ax[0].set_xticks(tickpos)
minorx0 = MultipleLocator(2)
minory0 = MultipleLocator(20)
ax[0].xaxis.set_minor_locator(minorx0)
ax[0].yaxis.set_minor_locator(minory0)
ax[0].set_yticks((-200,-100,0,100,200))
xlabel('x [Mpc]',size=12)
ylabel('v [km s$^{-1}$]',size=12)

ax.append(subplot(132))
ax[1].plot(x,v,linewidth=1)
ax[1].plot(xd,vd,'r.',linewidth=1,markersize=2)
ax[1].axis(limits[1]);
ax[1].set_xticks((4.6,5.0,5.4))
minorx1 = MultipleLocator(.05)
minory1 = MultipleLocator(20)
ax[1].xaxis.set_minor_locator(minorx1)
ax[1].yaxis.set_minor_locator(minory1)
ax[1].set_yticks((-200,-100,0,100,200))
xlabel('x [Mpc]',size=12)

ax.append(subplot(133))
ax[2].plot(x,v,linewidth=1)
ax[2].plot(xd,vd,'r.',linewidth=1,markersize=2)
ax[2].axis(limits[2]);
ax[2].set_xticks((4.94,5.0,5.06))
minorx2 = MultipleLocator(.01)
minory2 = MultipleLocator(5)
ax[2].xaxis.set_minor_locator(minorx2)
ax[2].yaxis.set_minor_locator(minory2)
ax[2].set_yticks((-30,0,30))
xlabel('x [Mpc]',size=12)

# set axis number font 
font = FontProperties(family = 'Helvetica',size=12);
for p in range(3):
    ax[p].set_yticklabels(ax[p].get_yticks(), fontproperties = font)
    ax[p].set_xticklabels(ax[p].get_xticks(), fontproperties = font)

#savefig('figa.eps',format='eps',bbox_inches='tight')
show()
