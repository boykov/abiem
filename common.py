#!/usr/bin/python
# -*- coding: utf-8 -*-
"""

"""

__author__ = "Evgeny Boykov. <artscan@list.ru>"

import unittest
from numpy import *
import sys,os
import logging

common_names = {
    "PI"                  : [lambda (s): 3.14159265358979324 , 0, "dp", "phi"],
    "max_neighbors"       : [lambda (s): 100                 , 0, "i", "phi"],
    "dim_3d"              : [lambda (s): 3                   , 0, "i", "phi"],
    "k_wave"              : [lambda (s): 0                   , 0, "dc", "phi"],
    "axes"                : [lambda (s): zeros((3))          , 0, "dp1d", "phi"],
    "gauss6"              : [lambda (s): False               , 0, "l", "phi"],
    "qbx"                 : [lambda (s): False               , 0, "l", "phi"],
    "matrixa6_p"          : [lambda (s): False               , 0, "l", "phi"],
    "use_int_neighbors_p" : [lambda (s): False               , 0, "l", "phi"],

    "numnodes"            : [lambda (s): s.N    , 10, "i", "phi"],
    "hval"                : [lambda (s): s.h    , 10, "dp", "phi"],
    "hval2"               : [lambda (s): s.h*s.h, 10, "dp", "phi"],

    "node_coordinates"    : [lambda (s): zeros((s.numnodes,3), order = 'Fortran') , 11, "dp2d", "phi"],
    "normal_coordinates"  : [lambda (s): zeros((s.numnodes,3), order = 'Fortran') , 11, "dp2d", "phi"],
    "int_neighbors1"      : [lambda (s): zeros((s.numnodes,s.max_neighbors),
                                               dtype = complex, order = 'Fortran'), 11, "dc2d", "phi"],
    "int_neighbors2"      : [lambda (s): zeros((s.numnodes,s.max_neighbors),
                                               dtype = complex, order = 'Fortran'), 11, "dc2d", "phi"],
    "node_neighbors1"     : [lambda (s): zeros((s.numnodes,s.max_neighbors),
                                               dtype = int32, order = 'Fortran')  , 11, "i2d", "phi"],
    "node_neighbors2"     : [lambda (s): zeros((s.numnodes,s.max_neighbors),
                                               dtype = int32, order = 'Fortran')  , 11, "i2d", "phi"],
    "nstroke_coordinates" : [lambda (s): zeros((s.numnodes,3), order = 'Fortran') , 11, "dp2d", "phi"],
    "intphi_over"         : [lambda (s): zeros((s.numnodes))                      , 11, "dp1d", "integ"],
    "intphi_under"        : [lambda (s): zeros((s.numnodes))                      , 11, "dp1d", "integ"],
    "area"                : [lambda (s): zeros((1))                               , 11, "dp1d", "integ"],
    "counter"             : [lambda (s): zeros((1))                               , 11, "dp1d", "integ"],

    "sigma"               : [lambda (s): zeros((s.numnodes))                      , 11, "dp1d", "integ"],
    "gauss"               : [lambda (s): zeros((s.numnodes,10),
                                               dtype = complex, order = 'Fortran'), 11, "dc2d", "integ"],

    "q_density"           : [lambda (s): zeros((s.numnodes), dtype = complex)     , 11, "dc1d", "integ"],
    "intG_x"              : [lambda (s): zeros((s.numnodes,6,11),
                                               dtype = complex, order = 'Fortran'), 11, "dc3d", "integ"],
    "valG_y"              : [lambda (s): zeros((s.numnodes,6,11),
                                               dtype = complex, order = 'Fortran'), 11, "dc3d", "integ"],

    "dim_quad"            : [lambda (s): s.q, 20, "i", "integ"],


    "quadphi_over"        : [lambda (s): zeros((s.dim_quad,2),order = 'Fortran')  , 21, "dp2d", "integ"],
    "quadphi_under"       : [lambda (s): zeros((s.dim_quad,2),order = 'Fortran')  , 21, "dp2d", "integ"],
    "quadsingular"        : [lambda (s): zeros((s.dim_quad,2),order = 'Fortran')  , 21, "dp2d", "integ"],
    "centres"             : [lambda (s): zeros((s.dim_quad))                      , 21, "dp1d", "integ"],
    "weights"             : [lambda (s): zeros((s.dim_quad))                      , 21, "dp1d", "integ"],
    "jacobian"            : [lambda (s): zeros((4,s.dim_quad,4*s.dim_quad),
                                               dtype = complex, order = 'Fortran'), 21, "dc3d", "integ"],
    "farr"                : [lambda (s): zeros((4,s.dim_quad,4*s.dim_quad),
                                               dtype = complex, order = 'Fortran'), 21, "dc3d", "integ"],
    "nodes"               : [lambda (s) :zeros((4,s.dim_quad,4*s.dim_quad,3),
                                               order = 'Fortran')                 , 21, "dp4d", "integ"]}

class common():
    types = {
        "l"    : ["logical"           , ""          , "=" , ""],
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

    def create_set_module(self):
        s = "! module set_params\n"
        for t in self.types.keys():
            s = s + "\n"
            s = s + "subroutine set_" + t + "_ptr(name,val)\n"
            s = s + "  character(128) :: name\n"
            s = s + "  " + self.types[t][0] + ", intent(in), target :: val" + self.types[t][1] + "\n"
            for n in common_names.keys():
                if common_names[n][2]==t:
                    s = s + "  if (name .eq. \"" + n + "\") " + n + " " + self.types[t][2] + " val\n"
            s = s + "end subroutine set_" + t + "_ptr\n"
        return s

    def create_module(self):
        s = "module params\n"
        for n in common_names.keys():
            s = s + (
                "  " +
                self.types[common_names[n][2]][0] +
                self.types[common_names[n][2]][3] +
                " :: " +
                n +
                self.types[common_names[n][2]][1] +
                "\n")
        s = s + "end module params"
        return s

    def write_module(self):
        file = open('params.f90', 'w')
        file.write(self.create_module())

    def write_set_module(self):
        file = open('set_params.f90', 'w')
        file.write(self.create_set_module())

    def __init__(self):
        for n in common_names.keys():
            if common_names[n][1]==0:
                setattr(self, n, common_names[n][0](self))

    def level1(self, N, h):
        self.N = N
        self.h = h
        for n in common_names.keys():
            if common_names[n][1]==10:
                setattr(self, n, common_names[n][0](self))
        for n in common_names.keys():
            if common_names[n][1]==11:
                setattr(self, n, common_names[n][0](self))

    def level2(self, q):
        self.q = q
        for n in common_names.keys():
            if common_names[n][1]==20:
                setattr(self, n, common_names[n][0](self))
        for n in common_names.keys():
            if common_names[n][1]==21:
                setattr(self, n, common_names[n][0](self))

    def setObjPhi(self,obj):
        for n in common_names.keys():
            if common_names[n][3]=="phi":
                s = "obj.set_" + common_names[n][2] + "_ptr(\"" + n + "\", self." + n + ")"
                eval(s)

    def setObjInteg(self, obj):
        self.setObjPhi(obj)
        for n in common_names.keys():
            if common_names[n][3]=="integ":
                s = "obj.set_" + common_names[n][2] + "_ptr(\"" + n + "\", self." + n + ")"
                eval(s)
