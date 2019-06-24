'''
Created on Jun 11, 2012

@author: dwferrer

Collection of methods to create an input file using some combination of user specified parameters and default values. Reads in a level 2  dynamic parameter file 
and writes a level 1 parameter file readable by parseheader 
'''

import csv
import time
import sys
import os
import os.path
import copy
import shlex

import Abacus.Tools

def getDefaults(filename):
    infile = csv.reader(open(filename,'r'),delimiter = ',',quoting = csv.QUOTE_NONE)
    kw = {}
    for line in infile:
        try:
            x = 0
            x = eval(line[1])
            kw[line[0]] = x
        except SyntaxError:
            kw[line[0]] = line[1]
    return kw

#find the first item looks like $VAR$ in input and replace it with the corresponding shell variable value
def shellreplace(line):
    pos1 = line.find("$")
    if pos1 == -1:
        return line
    pos2 = line.find("$",pos1+1)
    if pos2 == -1:
        print(("No matching $ in :%s"%(line)))
        raise SyntaxError
    
    varname = line[pos1+1:pos2]
    
    value = os.getenv(varname, "NONE")
    if value == "NONE":
        print(("Variable %s is not defined!"%(varname)))
        
        raise ValueError
    return line[:pos1] + value + line[pos2 +1:]
    
# Replace tokens like @VAR@ with their previously defined value
# Replace tokens like @_VAR@ indicates an optional parameter
def varreplace(line, values):
    pos1 = line.find("@")
    if pos1 == -1:
        return line
    pos2 = line.find("@",pos1+1)
    if pos2 == -1:
        print(("No matching @ in :%s"%(line)))
        raise SyntaxError
    
    varname = line[pos1+1:pos2]
    optional = varname[0] == '_'
    if optional:
        varname = varname[1:]
    if not (varname in values.keys()):
        if optional:
            # variable wasn't defined but was optional; do nothing
            return line[:pos1] + '""' + line[pos2+1:]
        else:
            print(("Variable %s is not defined!"%(varname)))
            raise SyntaxError
    return line[:pos1] + 'varreplace_values["{:s}"]'.format(varname) + line[pos2 +1:]

def parseInput(filename, values=None, fixvalues=None, varreplace_values=None):
    if not values:
        values = {}
    if not fixvalues:
        fixvalues = {}
    if not varreplace_values:
        varreplace_values = values

    with open(filename, 'r') as param:
        for line in param:
            if '\n' in line:
                # The end of header token is '^B\n'.
                break
            if not line.strip():
                continue
            #process shell variables
            while "$" in line:
                line = shellreplace(line)
            
            #process previously defined variables
            while "@" in line:
                line = varreplace(line,varreplace_values)
                            
            if line.startswith("#include"):
                #We are including the values from another parameter file. This will overwrite any already specified values
                qt1loc = line.find('"')
                qt2loc = line[qt1loc+1:].find('"')
                if qt1loc == -1 or qt2loc ==-1:
                    print("Bad syntax on include: %s"%(line))
                    raise SyntaxError
                
                fn = line[qt1loc+1:qt2loc+qt1loc+1]
                with Abacus.Tools.chdir(os.path.dirname(os.path.abspath(filename))):  # includes are specified relative to the containing file
                    parseInput(fn, values,fixvalues)
                continue
                           
            equals = line.find('=')
            comment = line.find('#')
            if equals ==-1:
                if line[:comment].strip(): 
                    print("Could not parse line: %s"%(line))
                    raise SyntaxError
                continue
            
            
            x = "I HAVE NO VALUE"
            key = line[:equals].strip()
            if key in fixvalues.keys():
                continue
            value = 0
            if comment != -1:
                value = line[equals+1:comment]
            else:
                value =line[equals+1:]
            # try to parse value
            try:
                items = value.split()
                vec  = ','.join(items)
                x = eval(f"({vec})")
            except (SyntaxError, NameError):
                #perhaps value is a vector
                try:
                    x = eval(value)
                except SyntaxError:
                    #fall back and try to parse as string
                    try:
                        x = eval('"{}"'.format(value))
                    except SyntaxError:
                        print("Error: Could not parse line:"+line)
                        raise

            values[key] = x
            #print("'%s': %s"%(str(key), str(x)))
        return values

'''Create a parameters file with name filename. Optionally, DefFilename as a list of parameters and default values to output.
 Following this is a list of specific keywords and values to set. If the keyword strict is True, then only keywords that appear exactly in the default
  list will be accepted. If other keywords are present, a ValueError will be raised. This is to prevent misspelled or mis-capitalized keywords from being
  set. If you want to set keywords not in the Defaults file, then set the strict keyword to false to disable this behavior
Example: makeInput("test.par",CPD = 45, BoxSize = 4096) will output an abacus formatted input file using default values for every parameter except
CPD and BoxSize. IF BoxSize does not appear in the file defaults.param, e.g. the keyword present there is 'Boxsize', a value error will be raised.
'''

def tostr(item): #recursively decompose the argument into a string we can print in the parameters file
    if isinstance(item, str):
        if not ('"' in item):
            return '"' + item + '"'
        return item
    else:
        try:
            it = iter(item)
        except TypeError:
            return repr(item)
        else:
            b = " ".join([tostr(i) for i in item])
            return b
         
        

def makeInput(outfn, inputfile, strict=False, del_keywords=[], **keywords):
    defaults = parseInput(inputfile)
    for keyword in keywords.keys():
        if strict and (keyword not in defaults.keys()):
            raise ValueError(keyword, "Keyword not in default parameter set, and `strict=True` was set.")

    res = parseInput(inputfile, keywords, copy.deepcopy(keywords))

    # Special behavior: do not write out keys whose value is None
    for keyword in list(res.keys()):
        if res[keyword] is None:
            del res[keyword]

    # Delete any keys listed in `del_keywords`
    for dk in del_keywords:
        try:
            del res[dk]
        except KeyError:
            pass
    
    # Now write the results to file
    with open(outfn,'w') as outfile:
        for keyword in sorted(res.keys()):
            outfile.write(keyword)
            outfile.write(' = ')
            value = res[keyword]
            if isinstance(value, str) and not ('"' in value):
                value = '"' +value +'"'
            outfile.write(tostr(value))
            outfile.write('\n')
        outfile.write('#created by GenParam from L2 file: {:s}    :: {:s}\n'.format(inputfile, time.asctime()))
        outfile.write('\n') #finalize the header
    return res
        

if __name__ == '__main__':
    params =  parseInput("test/test.par")
    for key in sorted(params.keys()):
        print("%s = %s"%(str(key),str(params[key])))
    
    #makeInput('test.par')
