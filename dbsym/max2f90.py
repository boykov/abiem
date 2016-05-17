#!/usr/bin/env python
# -*- coding: utf-8 -*-

# from sage.all import *
# from sage.interfaces.maple import maple
# from sage.interfaces.maxima import maxima

import optparse
import sys
import os
import re

def strtr(text, dic): 
    """ Replace in 'text' all occurences of any key in the given
    dictionary by its corresponding value.  Returns the new tring.""" 
    # http://code.activestate.com/recipes/81330/
    # Create a regulaur expression  from the dictionary keys
    import re
    regex = re.compile("(%s)" % "|".join(map(re.escape, dic.keys())))
    # For each match, look-up corresponding value in dictionary
    return regex.sub(lambda mo: str(dic[mo.string[mo.start():mo.end()]]), text)

onestring = {
    "&\n&"       : "",
    }

f90replace =    {
    "sin"       : "dsin",
    "cos"       : "dcos",
    "sqrt"       : "asqrt",
    "exp"        : "cdexp",
    "erf"        : "cderf",
    "pxe"        : "dexp",
    "'"          : "",
    "δ"          : "dn",
    "kron_delta" : "dn",    
    "a("         : "axes(",
    "βv"         : "bt",    
    "β"          : "bt",
    "φ"          : "ph",
    "ρ"          : "rh",
    "σ"          : "sigm",
    "%pi"        : "PI"
    }

if __name__ == "__main__":
    sys.path.append(os.getcwd())
    usage = 'Usage: %prog [options] <maxima_form.mpl> <fortran_form.f90>'
    parser =  optparse.OptionParser(usage)
    parser.add_option("-o", "--output", action="store", dest="filename", 
                      default="False",metavar="FILE", help="write output to FILE")


    (options, args) = parser.parse_args()
    maxima_form = args[0]
    fortran_form = args[1]
    pattern = open(fortran_form).read()

    values = open(maxima_form).read()

    ar = strtr(strtr(values,onestring),f90replace).split("\n")

    r = re.compile(r'\%\(.*\)s')

    names = map(lambda x: x[2:-2],r.findall(pattern))
    names.reverse()

    dct = dict(map(lambda i: [names[i],ar[i]],
                   map(lambda x: names.index(x),
                       names)))

    # dct = dict([['jacobian' , ar[0]],
    #             ['x' , ar[1]],
    #             ['beta' , ar[2]]])    
    
    print pattern % dct

