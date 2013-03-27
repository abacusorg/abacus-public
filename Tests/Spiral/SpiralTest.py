#!/usr/bin/python
'''
'Runs the spiral test and compares the result to the analytic answer
'
'''

import sys
import os
import shutil
import numpy as np
import matplotlib.pyplot as p
from mpl_toolkits.mplot3d import Axes3D
import ctypes as ct
import subprocess

    

#we have to do a little bit of trickery to import things from the abacus/python directory
abacuspath = os.getenv("ABACUS","NONE")
if abacuspath == "NONE":
    print("Error: Please define $ABACUS to be the absolute path to your abacus distribution")
    sys.exit(1)
abacuspythondir = abacuspath+"/python/"
if abacuspythondir not in sys.path:
    sys.path.insert(0, abacuspythondir)


import abacus
import GenParam
import InputFile
import Tools

def run(basedir = "NONE"):
    tmpdir = ""
    if basedir == "NONE":
        tmpdir = os.getenv("ABACUS_TMP","NONE")
        if tmpdir == "NONE":
            tmpdir = "/tmp/abacus/"
        basedir = tmpdir +"/spiral/"
        
    if os.path.exists(basedir):
    
        erase = raw_input("Test Directory exists! Erase? (y/n)")
        if erase =="y":   
            shutil.rmtree(basedir)
            print "Erased previous test directory"
    if not os.path.exists(basedir):
        os.makedirs(basedir)
        
    kvec = (1,0,0)
    phase = (np.pi,0,0)
    n1d = 16
    ainitial = 0.1
    across = 4
    astop =  1.0
    sf = .1/n1d#7.5e-03
    
    
    #check if we are done
    if not os.path.exists(basedir+"write/state"):
    
        params = GenParam.makeInput(basedir+"spiral.par", abacuspath +"/Tests/Spiral/spiral.par2", NP = n1d**3,nTimeSlice = 1, TimeSlicez = 1/astop -1, SofteningLength = sf,InitialRedshift = 1/ainitial -1,CPD = 35,BoxSize = 17.3205080756888)
        os.makedirs(params["InitialConditionsDirectory"])
        #make the spiral initial conditions
        subprocess.call([abacuspath+"/Tests/Spiral/makespiralics",str(n1d), str(ainitial),str(across),
                         str(kvec[0]),str(kvec[1]),str(kvec[2]),
                         str(phase[0]),str(phase[1]),str(phase[2]),
                         params["InitialConditionsDirectory"] + "/ic_0"])
        for i in range(1,params["CPD"]):
            f = open(params["InitialConditionsDirectory"]+"/ic_"+str(i),"wb")
            f.close()
        
        #run the problem
        os.chdir(basedir)
        abacus.run("spiral.par",1000,1)
    else:
        os.chdir(basedir)
        params = GenParam.parseInput("spiral.par") 
    #plot the results and check the answer
    writestate = InputFile.InputFile("write/state")
    ReadState = InputFile.InputFile("read/state")
    laststep = writestate.FullStepNumber
    
    timeslice = "final.ts"
    os.system(abacuspath+"/util/phdata " + "%s/slice%5.3f/%s.z%5.3f.*"%(params["OutputDirectory"], ReadState.Redshift, params["SimName"],ReadState.Redshift) + " > " +timeslice)
    data = np.fromfile(timeslice,dtype = np.float64)
    
    
    subprocess.call([abacuspath+"/Tests/Spiral/makeanalytic",str(ainitial),str(across),str(astop)])
    analytic = np.fromfile("./analytic",sep = " ")
    analytic = np.reshape(analytic,(-1,2))
    
    xv = np.reshape(data, (-1,6))
    print xv.shape
    p.plot(xv[:,0],xv[:,3],".")
    p.plot(analytic[:,0], analytic[:,1])
    p.xlabel("X")
    p.ylabel("Vx")
    p.show()


if __name__ == '__main__':
    
    args = sys.argv
    if len(args) == 1:
        run()
    elif len(args) == 2:
        run(basedir = sys.argv[1])
    else:
        print("Usage: ./SpiralTest.py <directory to run test in>")
        sys.exit(1)
