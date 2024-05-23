#!/usr/bin/env python3
'''
Runs the spiral test and compares the result to the analytic answer

'''

import sys
import os
from os.path import join as pjoin
import shutil
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as p
import subprocess

from Abacus import abacus
from Abacus import GenParam
from Abacus import InputFile
abacuspath = abacus.abacuspath

def run():
    kvec = (1,0,0)
    phase = (np.pi,0,0)
    n1d = 256
    BoxSize = 17.3205080756888
    ainitial = 0.09
    across = 0.1
    astop =  0.11
    sf = 0.1*BoxSize/n1d

    #make the spiral initial conditions
    params = GenParam.makeInput('spiral.par', 'spiral.par2')
    try:
        shutil.rmtree(params["InitialConditionsDirectory"])
    except FileNotFoundError:
        pass
    os.makedirs(params["InitialConditionsDirectory"])
    subprocess.run([pjoin(os.path.curdir,'makespiralics') ,str(n1d), str(ainitial), str(across),
                     str(kvec[0]),str(kvec[1]),str(kvec[2]),
                     str(phase[0]),str(phase[1]),str(phase[2]),
                     str(params["Omega_Smooth"]),
                     params["InitialConditionsDirectory"] + "/ic_0"],
                     check=True)
    for i in range(1,params["CPD"]):
        f = open(params["InitialConditionsDirectory"]+"/ic_"+str(i),"wb")
        f.close()

    #run the problem
    abacus.run('spiral.par2', clean=True, erase_ic=False,
                param_kwargs=dict(NP = n1d**3, nTimeSlice=1, TimeSliceRedshifts=1/astop -1, SofteningLength=sf, InitialRedshift=1/ainitial - 1, BoxSize=BoxSize)
              )
    
    #plot the results and check the answer
    #writestate = InputFile.InputFile("write/state")
    ReadState = InputFile.InputFile(pjoin(params['WorkingDirectory'], "read", "state"))
    #laststep = writestate.FullStepNumber
    
    timeslice = "final.ts"
    os.system(abacuspath+"/util/phdata " + "%s/slice%5.3f/%s.z%5.3f.*"%(params["OutputDirectory"], ReadState.Redshift, params["SimName"],ReadState.Redshift) + " > " +timeslice)
    data = np.fromfile(timeslice,dtype = np.float64)


    subprocess.call([abacuspath+"/Tests/Spiral/makeanalytic",str(ainitial),str(across),str(astop),str(params["Omega_Smooth"])])
    analytic = np.fromfile("./analytic",sep = " ")
    analytic = np.reshape(analytic,(-1,2))
    # The makeanalytic program returns canonical velocities, whereas abacus returns ZSpace velocities.
    # Let's standardize on ZSpace, because that connects to the displacements
    analytic[:,1] /= ReadState.VelZSpace_to_Canonical

    xv = np.reshape(data, (-1,6))
    assert len(xv) == n1d**3
    print(xv.shape)
    print(analytic.shape)
    #p.plot(analytic[:,0], analytic[:,1])
    #p.scatter(xv[:,0],xv[:,3],marker=".",s = 10*(4 * np.abs(xv[0]) + 1), c= "r",linewidth=0 )
    p.xlabel("X")
    p.ylabel("Vx")
    p.savefig("spiral.png", dpi=500)
    p.figure()
    p.xlim(-.1,.1)
    p.ylim(-.05,.05)
    #p.plot(analytic[:,0], analytic[:,1])
    #p.scatter(xv[:,0],xv[:,3],marker=".",s = 10*(4 * np.abs(xv[0]) + 1), c= "r",linewidth=0 )
    p.xlabel("X")
    p.ylabel("Vx")
    p.savefig("spiral-zoomed.png",dpi=500)

    if astop<across:
        # There's a skewer of particles that always starts at -0.50, but it might have wrapped to +0.5
        sel = np.where(np.logical_and(abs(xv[:,1])>0.5-1e-3, abs(xv[:,2])>0.5-1e-3))
        print("Selecting one skewer of %d points."%(len(sel[0])))
        print("The following will only make sense before shell crossing:")
        print("The deviation is relative to Zel'dovich, so one expects breakdowns approaching shell crossing")
        print("Grid position     X_final    VX_final     (X-VX-grid) deviation")
        xsel = xv[sel[0],0]
        vsel = xv[sel[0],3]
        #print xsel, vsel
        index = np.argsort(xsel)
        #print index
        xsel = xsel[index]
        vsel = vsel[index]
        #print xsel, vsel
        trackback = (xsel-vsel)*n1d
        startguess = np.round(trackback)
        residual = trackback - startguess
        #print startguess, residual
        np.set_printoptions(precision=4)
        print(np.transpose([startguess,xsel,vsel,residual]))
        print()
        print("RMS of deviations from Zeldovich %e compared to rms Vx of %f." % (np.std(residual), np.std(xv[:,3])))
        print()

    print("Vx: rms %f, max %f"%( np.std(xv[:,3]), np.max(np.absolute(xv[:,3])) ))
    print("Vy: rms %e, max %e"%( np.std(xv[:,4]), np.max(np.absolute(xv[:,4])) ))
    print("Vz: rms %e, max %e"%( np.std(xv[:,5]), np.max(np.absolute(xv[:,5])) ))
    print("Ratio of max velocity (analytic/computed): %f"%(np.max(analytic[:,1])/np.max(xv[:,3])))

    check_pids(params, n1d)


from Abacus import ReadAbacus
def check_pids(params, n1d):
    particles = ReadAbacus.from_dir(pjoin(params.get('LocalWorkingDirectory', params['WorkingDirectory']), 'read'), pattern='position_*', return_pid=True, format='state')
    pids = particles['pid']
    pids.sort()
    assert (np.diff(pids) > 0).all()
    assert len(pids) == n1d**3
    print('All particles present with no duplicates.')


def animate_analytic(astart,astop,across,npoints):
    os.chdir("/tmp/")
    a_out = np.linspace(astart,astop,npoints)
    ReadState = InputFile.InputFile("/mnt/raid/doug/spiral/read/state")
    i = 0
    for a in a_out:
        subprocess.call([abacuspath+"/Tests/Spiral/makeanalytic",str(astart),str(across),str(a)])
        analytic = np.fromfile("./analytic",sep = " ")
        analytic = np.reshape(analytic,(-1,2))
        analytic[:,1] /= ReadState.VelZSpace_to_Canonical
        p.plot(analytic[:,0], analytic[:,1])
        p.xlabel("X")
        p.ylabel("Vx")
        p.title("a={:4.3f}".format(a))
        p.savefig("spiral{:05d}.png".format(i))
        p.cla()
        i+=1


if __name__ == '__main__':
    args = sys.argv
    if len(args) == 1:
        run()
    else:
        print("Usage: ./SpiralTest.py <directory to run test in>")
        sys.exit(1)
