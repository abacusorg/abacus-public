# Copyright 2012-2025 The Abacus Developers
# SPDX-License-Identifier: GPL-3.0-or-later

'''
Created on Feb 2, 2012

@author: dferrer
Updates by Lehman Garrison

Wrapper for input files. Files consist of two parts: a parameters section that essentially looks like python code, and a binary section. When this class is instantiated
for a file, it attempts to 'register' all of the parameters by simply executing this portion of of the file as slightly modified python code. This is a sketchy way to
do it, but it avoids the difficulty of having to rewrite ParseHeader in python. These parameters can then be accessed as class variables. Right now, full functionality
is only available for scalar data. Any parameter type that isn't understood is turned into a string, which can be parsed elsewhere.

e.g if infile.txt has the lines

a = 1
b = 2 2 2
c = 35v kls a

then if infile = InputFile('infile.txt')

infile.a returns 1
infile.b returns [2,2,2]
infile.c returns "35v kls a"

EXERCISE EXTREME CAUTION: THIS CLASS ALLOWS ARBITRARY EXECUTION OF TEXT FILES AS CODE
'''

from io import StringIO
import shlex

class InputFile:
    def __init__(self, fn=None, str_source=None):
        """
        Construct an InputFile from filename `fn`, or a string containing the file contents `str_source`.
        """
        
        # Read the input as either a file or string
        if fn and str_source:
            raise ValueError('Cannot specify both `fn` = "{}" and `str_source` = "{}"!'.format(fn, str_source))
        if not fn and not str_source:
            raise ValueError('Must specify one of `fn` and `str_source`!')
            
        if fn:
            param = open(fn, "r")
        elif str_source:
            param = StringIO(str_source)
        else:
            raise RuntimeError('Invalid state. Should never be able to get here!')
            
        code = param.readlines()
        param.close()
        
        for line in code:
            line = line.strip()
            comment = line.find('#')
            if comment > -1:  # Trim everything after a comment
                line = line[:comment]
            equals = line.find('=')
            if equals == -1:
                continue
            key = line[:equals].strip()
            valstr = line[equals+1:].strip()
            try:
            	# valid as a tuple, if we replace spaces with commas?
                items = shlex.split(valstr)
                vec  = ','.join(items)
                value = eval('({})'.format(vec))
            except (SyntaxError,NameError):
                try:
                    # valid python as-is?
                    value = eval(valstr)
                except:
                    try:
                    	# valid as a string if we wrap the whole RHS in quotes?
                    	value = eval('"{}"'.format(valstr))
                    except:
                        raise ValueError('Error: Could not parse line: "{:s}"\n\tfrom file "{}"'.format(line, fn))
            setattr(self, key, value)
        
        # A common pathology is that BoxSize is interpreted as an int by ParseHeader
        # We should endeavor to always write "50." instead of "50" in the .par files
        for field in ['BoxSize', 'InitialRedshift', 'ZD_PLT_target_z', 'wa', 'w0']:
            if field in self:
                setattr(self, field, float(self[field]))

        # HDF5 can't write lists of unicode strings, but a single string is fine (strange!)
        # So, we store the code as a single string that can later be split again if desired
        self.code = ''.join(code)

    # Allows access as "params[key]"
    def __getitem__(self, key):
        if key not in vars(self):
            raise KeyError("InputFile has no field {}".format(repr(key)))
        return getattr(self, key)
    def __setitem__(self, key, value):
        return setattr(self, key, value)

    def keys(self):
        selfvars = vars(self).copy()
        try:
            del selfvars['code']
        except KeyError:
            pass
        return selfvars.keys()
        
    def get(self, key, default=None):
        try:
            return self[key]
        except KeyError:
            return default
    
    # Allows testing "key in params"
    def __contains__(self, key):
        return hasattr(self, key)
    
    # Don't return 'code' when asking for a string representation
    def __repr__(self):
        selfvars = vars(self).copy()
        try:
            del selfvars['code']
        except KeyError:
            pass
        return str(selfvars)
    