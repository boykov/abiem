#!/usr/bin/python
# -*- coding: utf-8 -*-
"""

"""

__author__ = "Evgeny Boykov. <artscan@list.ru>"

import unittest
from ellipsoid import ellipsoid
from numpy import *
from utilities import *
from phi import modphi as phi

class testBIE(unittest.TestCase):
    def setUp(self):
        self.number_points = 400
        self.list_axes = zeros((3), order = 'Fortran')
        self.list_axes[:] = [float(0.75),float(1.),float(0.5)]
        pass

    def testEllipsoid(self):
        e = ellipsoid(self.list_axes, self.number_points)
        self.assertAlmostEqual(center_points(e.points), 1.0e-11, places = 10)

    def testPhi(self):
        self.e = ellipsoid(self.list_axes, self.number_points)
        
        numnodes = self.e.points.shape[0]
        node_neighbors1 = zeros((numnodes,100), dtype = int32, order = 'Fortran')
        node_neighbors2 = zeros((numnodes,100), dtype = int32, order = 'Fortran')
        nstroke_coordinates = zeros((numnodes,3), order = 'Fortran')
        
        phi.max_neighbors = 100
        phi.h = self.e.get_h()
        phi.numnodes = numnodes
        
        phi.set_node_coordinates(self.e.points)
        phi.set_normal_coordinates(self.e.normalvectors)
        phi.set_nstroke_coordinates(nstroke_coordinates)
        phi.set_node_neighbors1(node_neighbors1)
        phi.set_node_neighbors2(node_neighbors2)
        phi.set_axes(self.list_axes)
        
        phi.filter_neighbors(2*phi.h, node_neighbors2, phi.numnodes)
        phi.filter_neighbors(  phi.h, node_neighbors1, phi.numnodes)
        phi.normal_vector_stroke(phi.numnodes, node_neighbors1)

        self.assertAlmostEqual(phi.test_phi(gen_points(numnodes,self.list_axes)), numnodes, places = 12)
