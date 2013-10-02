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

class testBIE(unittest.TestCase):
    def setUp(self):
        pass

    def scenario(self):
        self.e = ellipsoid(self.data.axes, self.data.numpoints)
        self.data.setupquad()
        self.numnodes = self.e.points.shape[0]

        self.node_neighbors1 = zeros((self.numnodes,100), dtype = int32, order = 'Fortran')
        self.node_neighbors2 = zeros((self.numnodes,100), dtype = int32, order = 'Fortran')
        self.nstroke_coordinates = zeros((self.numnodes,3), order = 'Fortran')
        self.intphi_over = zeros((self.numnodes), order = 'Fortran')

        self.centres = zeros((self.data.orderquad), order = 'Fortran')
        self.C = zeros((self.data.orderquad), order = 'Fortran')
        self.jacobian = zeros((4,self.data.orderquad,4*self.data.orderquad), dtype = complex, order = 'Fortran')
        self.nodes = zeros((4,self.data.orderquad,4*self.data.orderquad,3), order = 'Fortran')

        self.centres[:] = self.data.quadphi_over[:,0]
        self.C[:] = self.data.quadphi_over[:,1]

    def setObj(self,obj):
        obj.set_nd(3)
        obj.set_pi(3.14159265358979324)
        obj.set_max_neighbors(100)
        obj.set_h(self.e.get_h())
        obj.set_numnodes(self.numnodes)

        obj.set_node_coordinates(self.e.points)
        obj.set_normal_coordinates(self.e.normalvectors)
        obj.set_nstroke_coordinates(self.nstroke_coordinates)
        obj.set_node_neighbors1(self.node_neighbors1)
        obj.set_node_neighbors2(self.node_neighbors2)
        obj.set_axes(self.data.axes)
        
    def setPhi(self):
        self.setObj(phi)    

        phi.filter_neighbors(2*self.e.get_h(), self.node_neighbors2, self.numnodes)
        phi.filter_neighbors(  self.e.get_h(), self.node_neighbors1, self.numnodes)
        phi.normal_vector_stroke(self.numnodes, self.node_neighbors1)

    def setInteg(self):
        self.setObj(integ)    

        integ.set_intphi_over(self.intphi_over)
        integ.set_h2(self.e.get_h()*self.e.get_h())
        integ.set_nz(self.data.orderquad)

        integ.set_centres(self.centres)
        integ.set_c(self.C)
        integ.set_jacobian(self.jacobian)
        integ.set_nodes(self.nodes)

    def testEllipsoid(self):
        self.data = DataElement(400)
        self.data.orderquad = 20
        self.scenario()
        self.assertAlmostEqual(center_points(self.e.points), 1.0e-11, places = 10)

    def testPhi(self):
        self.data = DataElement(400)
        self.data.orderquad = 20
        self.scenario()
        self.setPhi()
        self.assertAlmostEqual(phi.test_phi(gen_points(self.numnodes,self.data.axes)), self.numnodes, places = 12)

    def testInteg(self):
        self.data = DataElement(400)
        self.data.orderquad = 20
        self.scenario()
        self.setPhi()
        self.setInteg()
        self.assertAlmostEqual(integ.calcarea(),6.971610618375645,places = 6)
