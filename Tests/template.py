# Copyright 2012-2025 The Abacus Developers
# SPDX-License-Identifier: GPL-3.0-or-later

#/usr/bin/env python
'''
Template for all tests. Put any functionality in the *run* method
'''
#put general imports here
import sys
import os
import shutil
import numpy as np
import matplotlib.pyplot as plt

#put imports from the abacus/python directory here
from Abacus import abacus
from Abacus import GenParam
from Abacus import InputFile
from Abacus import Tools

#you can define any functions needed for the test here


def run(basedir=None):
    tmpdir = os.getenv("ABACUS_TMP","/tmp/abacus/")
    if not basedir:
        #change this to reflect the name of your test
        basedir = tmpdir + "/TEST/"
        
    if os.path.exists(basedir):
        erase = input("Test directory exists! Erase? (y/[n]) ")
        if erase == "y":
            shutil.rmtree(basedir)
            print("Erased previous test directory")
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
        print("Usage: {} [directory to run test in].format(args[0]))")
        sys.exit(1)