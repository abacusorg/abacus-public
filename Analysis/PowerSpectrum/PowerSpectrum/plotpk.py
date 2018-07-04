#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as p
import sys
import scipy.interpolate

p.xscale("log")
p.yscale("log")
def Fg(a,omega_m):
    x = (1-omega_m)/omega_m *a**3
    return a * (1+x)**.5 * (1 + 1.175*x + 0.3064*x**2 + 0.005355*x**3)/(1 + 1.857*x + 1.021*x**2 + 0.1530*x**3)
rescale = 0
for i in range(1,len(sys.argv)):
	fn = sys.argv[i]
        #rescale = 0
	if ".npz" in fn:
		l = np.load(fn)
		k = l["k"]
		print(len(k))
		P = l["P"]
		p_k = scipy.interpolate.interp1d(k,P)
		if rescale ==0:
			rescale = p_k(.1)
			print(("Rescale: {}".format(rescale)))
                else: P = P* rescale/p_k(.1)
		sigmaP = 0.  #P * 1.0/np.sqrt(l["nb"])
		p.errorbar(k,P,yerr=sigmaP,fmt=".",label="Observed")
		continue
	kP = np.loadtxt(fn)
	if kP.shape[1]==4:
		p.plot(kP[:,1],kP[:,2],label=fn)
		continue
        p_k = scipy.interpolate.interp1d(kP[:,0],kP[:,1])
        if rescale ==0:
		rescale = p_k(.1)
		print(("Rescale: {:f}".format(1.0/rescale)))
	if kP.shape[1] == 2:
		p.plot(kP[:,0],kP[:,1] *  rescale/p_k(.1),label=fn)
	else:
		p.plot(kP[:,0],kP[:,1] *  rescale/p_k(.1),label=fn+" bao")
		p.plot(kP[:,0],kP[:,2] *  rescale/p_k(.1),label=fn+"smooth")
for i in range(1,2):
	vx = np.array([i*2.0*np.pi/3176,i*2.0*np.pi/3176])
	vy = np.array([1e-10,1e10])
	p.plot(vx,vy,"--",label = "")
p.legend(loc="lower left")
p.ylim(1e4,1e6)
p.xlim(2e-2,3e-1)
p.show()
