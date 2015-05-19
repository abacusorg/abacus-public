#/usr/bin/python
'''
'Template for all tests. Put any functionality in the *run* method
'''
#put general imports here
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

#put imports from the abacus/python directory here
import abacus
import GenParam
import InputFile
import Tools

#you can define any functions needed for the test here


def run(basedir = "NONE"):
    tmpdir = ""
    if basedir == "NONE":
        tmpdir = os.getenv("ABACUS_TMP","NONE")
        if tmpdir == "NONE":
            tmpdir = "/tmp/abacus/"
        #change this to reflect the name of your test
        basedir = tmpdir +"/TEST/"
        
    if os.path.exists(basedir):
    
        erase = raw_input("Test Directory exists! Erase? (y/n)")
        if erase =="y":   
            shutil.rmtree(basedir)
            print "Erased previous test directory"
    if not os.path.exists(basedir):
        os.makedirs(basedir)
    #You should put the test functionality here.
    
    


if __name__ == '__main__':
    
    args = sys.argv
    if len(args) == 1:
        run()
    elif len(args) == 2:
        run(basedir = sys.argv[1])
    else:
        #change this line to reflect the name of your test
        print("Usage: ./TEST.py <directory to run test in>")
        sys.exit(1)