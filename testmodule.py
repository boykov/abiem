#!/usr/bin/python
# -*- coding: utf-8 -*-
"""

"""

__author__ = "Evgeny Boykov. <artscan@list.ru>"

import unittest
from ellipsoid import ellipsoid
from dataelement import *
from numpy import *
import cmath
import math
from utilities import *
from phi import modphi as phi
from integ import modinteg as integ

import sys,os
import slaeahmed
from memo import memoize

import logging
from testsql import *

max_neighbors = 100
PI = 3.14159265358979324
dim_3d = 3
axes = zeros((3))
numpoints = 400
axes[:] = [float(0.75),float(1.),float(0.5)]

k_wave = 0
data = DataElement(numpoints)
data.magic = 0.410
name_matrixa = 'integ.matrixa'
name_vectorb = 'integ.vectorb'
name_approximateu = 'integ.approximateu'
integ_places = 4
slae_tol = 0.0
slae_places = 0
orderquad = 20
eps_matgen = 1e-4
eps_aggl = eps_matgen
eps_gmres = eps_matgen*0.01
eta = 0.8
bmin = 15
rankmax = 1000
flagMemo = True

def setEllipsoid(axes, numpoints):
    global numnodes, hval, hval2
    global normal_coordinates, node_coordinates
    e = ellipsoid(axes, numpoints)
    numnodes = e.points.shape[0]
    hval = e.get_h()
    hval2 = hval * hval
    node_coordinates = zeros((numnodes,3), order = 'Fortran')
    normal_coordinates = zeros((numnodes,3), order = 'Fortran')
    node_coordinates[:,:] = e.points[:,:]
    normal_coordinates[:,:] = e.normalvectors[:,:]

def setObjPhi(obj):
    global dim_3d, PI, max_neighbors, numnodes, hval, hval2
    global normal_coordinates, node_coordinates, nstroke_coordinates
    global axes, node_neighbors1, node_neighbors2
    obj.set_dim_3d(dim_3d)
    obj.set_pi(PI)
    obj.set_max_neighbors(max_neighbors)
    obj.set_hval(hval)
    obj.set_numnodes(numnodes)

    obj.set_dp1d_ptr("axes", axes)
    node_neighbors1 = zeros((numnodes,max_neighbors), dtype = int32, order = 'Fortran')
    node_neighbors2 = zeros((numnodes,max_neighbors), dtype = int32, order = 'Fortran')
    nstroke_coordinates = zeros((numnodes,3), order = 'Fortran')

    obj.set_dp2d_ptr("node_coordinates", node_coordinates)
    obj.set_dp2d_ptr("normal_coordinates", normal_coordinates)
    obj.set_dp2d_ptr("nstroke_coordinates", nstroke_coordinates)
    obj.set_node_neighbors1(node_neighbors1)
    obj.set_node_neighbors2(node_neighbors2)
    
    phi.filter_neighbors(2*hval, node_neighbors2, numnodes)
    phi.filter_neighbors(  hval, node_neighbors1, numnodes)
    phi.normal_vector_stroke(numnodes, node_neighbors1)
