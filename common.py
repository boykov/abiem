#!/usr/bin/python
# -*- coding: utf-8 -*-
"""

"""

__author__ = "Evgeny Boykov. <artscan@list.ru>"

import unittest
from numpy import *
import sys,os
import logging

class common():
    types = {
        "dp"   : ["double precision"  , ""          , "=" , ""],
        "i"    : ["integer"           , ""          , "=" , ""],
        "dc"   : ["double complex"    , ""          , "=" , ""],
        "dp1d" : ["double precision"  , "(:)"       , "=>", ", pointer"],
        "dp2d" : ["double precision"  , "(:,:)"     , "=>", ", pointer"],
        "i2d"  : ["integer"           , "(:,:)"     , "=>", ", pointer"],
        "dc1d" : ["double complex"    , "(:)"       , "=>", ", pointer"],
        "dc2d" : ["double complex"    , "(:,:)"     , "=>", ", pointer"],
        "dc3d" : ["double complex"    , "(:,:,:)"   , "=>", ", pointer"],
        "dp4d" : ["double precision"  , "(:,:,:,:)" , "=>", ", pointer"]}

    names = {
        "PI"                  : ["3.14159265358979324" , 0, "dp", "phi"],
        "max_neighbors"       : ["100"                 , 0, "i", "phi"],
        "dim_3d"              : ["3"                   , 0, "i", "phi"],
        "k_wave"              : ["0"                   , 0, "dc", "phi"],
        "axes"                : ["zeros((3))"          , 0, "dp1d", "phi"],

        "numnodes"            : ["N"   , 10, "i", "phi"],
        "hval"                : ["h"   , 10, "dp", "phi"],
        "hval2"               : ["h*h" , 10, "dp", "phi"],

        "node_coordinates"    : ["zeros((self.numnodes,3), order = 'Fortran')"    , 11, "dp2d", "phi"],
        "normal_coordinates"  : ["zeros((self.numnodes,3), order = 'Fortran')"    , 11, "dp2d", "phi"],
        "node_neighbors1"     : ["""zeros((self.numnodes,self.max_neighbors),
                                           dtype = int32, order = 'Fortran')"""   , 11, "i2d", "phi"],
        "node_neighbors2"     : ["""zeros((self.numnodes,self.max_neighbors),
                                           dtype = int32, order = 'Fortran')"""   , 11, "i2d", "phi"],
        "nstroke_coordinates" : ["""zeros((self.numnodes,3),
                                           order = 'Fortran')"""                  , 11, "dp2d", "phi"],
        "intphi_over"         : ["zeros((self.numnodes))"                         , 11, "dp1d", "integ"],
        "intphi_under"        : ["zeros((self.numnodes))"                         , 11, "dp1d", "integ"],
        "area"                : ["zeros((1))"                                     , 11, "dp1d", "integ"],
        "counter"             : ["zeros((1))"                                     , 11, "dp1d", "integ"],

        "sigma"               : ["zeros((self.numnodes))"                         , 11, "dp1d", "integ"],
        "gauss"               : ["""zeros((self.numnodes,10),
                                           dtype = complex, order = 'Fortran')""" , 11, "dc2d", "integ"],

        "q_density"           : ["zeros((self.numnodes), dtype = complex)"        , 11, "dc1d", "integ"],

        "dim_quad"            : ["q", 20, "i", "integ"],

        "quadphi_over"        : ["zeros((self.dim_quad,2),order = 'Fortran')" , 21, "dp2d", "integ"],
        "quadphi_under"       : ["zeros((self.dim_quad,2),order = 'Fortran')" , 21, "dp2d", "integ"],
        "quadsingular"        : ["zeros((self.dim_quad,2),order = 'Fortran')" , 21, "dp2d", "integ"],
        "centres"             : ["zeros((self.dim_quad))"                     , 21, "dp1d", "integ"],
        "weights"             : ["zeros((self.dim_quad))"                     , 21, "dp1d", "integ"],
        "jacobian"            : ["""zeros((4,self.dim_quad,4*self.dim_quad),
                                dtype = complex, order = 'Fortran')"""        , 21, "dc3d", "integ"],
        "nodes"               : ["""zeros((4,self.dim_quad,4*self.dim_quad,3),
                                     order = 'Fortran')"""                    , 21, "dp4d", "integ"]}

    def create_set_module(self):
        s = "! module set_params\n"
        for t in self.types.keys():
            s = s + "\n"
            s = s + "subroutine set_" + t + "_ptr(name,val)\n"
            s = s + "  character(128) :: name\n"
            s = s + "  " + self.types[t][0] + ", intent(in), target :: val" + self.types[t][1] + "\n"
            for n in self.names.keys():
                if self.names[n][2]==t:
                    s = s + "  if (name .eq. \"" + n + "\") " + n + " " + self.types[t][2] + " val\n"
            s = s + "end subroutine set_" + t + "_ptr\n"
        return s

    def create_module(self):
        s = "module params\n"
        for n in self.names.keys():
            s = s + "  " + self.types[self.names[n][2]][0]  + self.types[self.names[n][2]][3] + " :: " + n + self.types[self.names[n][2]][1] + "\n"
        s = s + "end module params"
        return s

    def __init__(self):
        for n in self.names.keys():
            if self.names[n][1]==0:
                setattr(self, n, eval(self.names[n][0]))

    def level1(self, N, h):
        for n in self.names.keys():
            if self.names[n][1]==10:
                setattr(self, n, eval(self.names[n][0]))
        for n in self.names.keys():
            if self.names[n][1]==11:
                setattr(self, n, eval(self.names[n][0]))

    def level2(self, q):
        for n in self.names.keys():
            if self.names[n][1]==20:
                setattr(self, n, eval(self.names[n][0]))
        for n in self.names.keys():
            if self.names[n][1]==21:
                setattr(self, n, eval(self.names[n][0]))

    def setObjPhi(self,obj):
        for n in self.names.keys():
            if self.names[n][3]=="phi":
                s = "obj.set_" + self.names[n][2] + "_ptr(\"" + n + "\", self." + n + ")"
                eval(s)

    def setObjInteg(self, obj):
        self.setObjPhi(obj)
        for n in self.names.keys():
            if self.names[n][3]=="integ":
                s = "obj.set_" + self.names[n][2] + "_ptr(\"" + n + "\", self." + n + ")"
                eval(s)
