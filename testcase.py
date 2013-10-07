#!/usr/bin/python
# -*- coding: utf-8 -*-
"""

"""

__author__ = "Evgeny Boykov. <artscan@list.ru>"

import unittest
from ellipsoid import ellipsoid
from dataelement import *
from numpy import *
from utilities import *
from phi import modphi as phi
from integ import modinteg as integ

import sys,os
import petsc4py
from petsc4py import PETSc

class params():
    def __init__(self, numpoints = 400, axes = [float(0.75),float(1.),float(0.5)]):
        self.max_neighbors = 100
        self.PI = 3.14159265358979324
        self.nd = 3
        self.axes = zeros((3), order = 'Fortran')
        self.numpoints = numpoints
        self.axes[:] = [float(0.75),float(1.),float(0.5)]
        pass

    def initQuad(self, orderquad):
        data = DataElement(1)
        data.orderquad = orderquad
        data.setupquad()
        self.Nz = data.orderquad
        self.quadphi_over = zeros((self.Nz,2),order = 'Fortran')
        self.quadphi_over[:,:] = data.quadphi_over[:,:]

    def initEllipsoid(self):
        self.e = ellipsoid(self.axes, self.numpoints)
        self.numnodes = self.e.points.shape[0]
        self.h = self.e.get_h()
        self.h2 = self.h * self.h
        self.node_coordinates = zeros((self.numnodes,3), order = 'Fortran')
        self.normal_coordinates = zeros((self.numnodes,3), order = 'Fortran')
        self.node_coordinates[:,:] = self.e.points[:,:]
        self.normal_coordinates[:,:] = self.e.normalvectors[:,:]

    def initPhi(self):
        self.node_neighbors1 = zeros((self.numnodes,100), dtype = int32, order = 'Fortran')
        self.node_neighbors2 = zeros((self.numnodes,100), dtype = int32, order = 'Fortran')
        self.nstroke_coordinates = zeros((self.numnodes,3), order = 'Fortran')

        self.setObjPhi(phi)

        phi.filter_neighbors(2*self.h, self.node_neighbors2, self.numnodes)
        phi.filter_neighbors(  self.h, self.node_neighbors1, self.numnodes)
        phi.normal_vector_stroke(self.numnodes, self.node_neighbors1)


    def initInteg(self):
        self.intphi_over = zeros((self.numnodes), order = 'Fortran')

        self.centres = zeros((self.Nz), order = 'Fortran')
        self.C = zeros((self.Nz), order = 'Fortran')
        self.jacobian = zeros((4,self.Nz,4*self.Nz), dtype = complex, order = 'Fortran')
        self.nodes = zeros((4,self.Nz,4*self.Nz,3), order = 'Fortran')

        self.centres[:] = self.quadphi_over[:,0]
        self.C[:] = self.quadphi_over[:,1]

        self.area = zeros((1), order = 'Fortran')

        self.setObjInteg(integ)


    def setObjPhi(self,obj):
        obj.set_nd(self.nd)
        obj.set_pi(self.PI)
        obj.set_max_neighbors(self.max_neighbors)
        obj.set_h(self.h)
        obj.set_numnodes(self.numnodes)

        obj.set_axes(self.axes)

        obj.set_node_coordinates(self.node_coordinates)
        obj.set_normal_coordinates(self.normal_coordinates)
        obj.set_nstroke_coordinates(self.nstroke_coordinates)
        obj.set_node_neighbors1(self.node_neighbors1)
        obj.set_node_neighbors2(self.node_neighbors2)

    def setObjInteg(self, obj):
        self.setObjPhi(obj)

        obj.set_intphi_over(self.intphi_over)
        obj.set_h2(self.h2)
        obj.set_nz(self.Nz)

        obj.set_centres(self.centres)
        obj.set_c(self.C)
        obj.set_jacobian(self.jacobian)
        obj.set_nodes(self.nodes)

        obj.set_area(self.area)

class testBIE(unittest.TestCase):
    def setUp(self):
        pass

    def testEllipsoid(self):
        self.P = params(400)
        self.P.initQuad(20)
        self.P.initEllipsoid()
        self.assertAlmostEqual(center_points(self.P.node_coordinates), 1.0e-11, places = 10)

    def testPhi(self):
        self.P = params(400)
        self.P.initQuad(20)
        self.P.initEllipsoid()
        self.P.initPhi()
        self.assertAlmostEqual(phi.test_phi(gen_points(self.P.numnodes,self.P.axes)), self.P.numnodes, places = 12)

    def testInteg(self):
        self.P = params(400)
        self.P.initQuad(20)
        self.P.initEllipsoid()
        self.P.initPhi()
        self.P.initInteg()
        integ.calcomp()
        self.P.area[0] = sum(self.P.intphi_over)
        self.assertAlmostEqual(self.P.area[0],6.971610618375645,places = 6)
