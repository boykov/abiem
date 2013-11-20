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
from dbsym import dbsym

import sys,os
import slaeahmed

class params():
    def __init__(self, numpoints = 400, axes = [float(0.75),float(1.),float(0.5)]):
        self.max_neighbors = 100
        self.PI = 3.14159265358979324
        self.nd = 3
        self.axes = zeros((3))
        self.numpoints = numpoints
        self.axes[:] = [float(0.75),float(1.),float(0.5)]
        self.k = 0
        self.data = DataElement(numpoints)
        self.data.magic = 0.410
        pass

    def initAHMED(self):
        """
        """
        self.eps = 1e-6
        self.eta = 0.8
        self.bmin = 15
        self.rankmax = 1000
        self.points = zeros((self.numnodes,3)) # C order !!!
        self.points[:,:] = self.node_coordinates[:,:]
        self.q_ahmed = zeros((self.numnodes), dtype = complex)

        slaeahmed.set_Hmatrix(self.eps,
                              self.eta,
                              self.bmin,
                              self.rankmax)

        slaeahmed.set_points(self.points)
        slaeahmed.set_vectorb(integ.vectorb)
        slaeahmed.set_kernel(integ.matrixa)
        slaeahmed.solve_slae()
        slaeahmed.get_q(self.q_ahmed)

    def initQuad(self, orderquad):
        self.data.orderquad = orderquad
        self.data.setupquad()
        self.Nz = self.data.orderquad
        self.quadphi_over = zeros((self.Nz,2),order = 'Fortran')
        self.quadphi_over[:,:] = self.data.quadphi_over[:,:]

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
        self.intphi_over = zeros((self.numnodes))

        self.centres = zeros((self.Nz))
        self.C = zeros((self.Nz))
        self.jacobian = zeros((4,self.Nz,4*self.Nz), dtype = complex, order = 'Fortran')
        self.nodes = zeros((4,self.Nz,4*self.Nz,3), order = 'Fortran')

        self.centres[:] = self.quadphi_over[:,0]
        self.C[:] = self.quadphi_over[:,1]

        self.area = zeros((1))

        self.setObjInteg(integ)
        integ.calcomp()

        self.sigma = zeros((self.numnodes))
        self.sigma[:] = map(self.data.fsigma,self.intphi_over)[:]
        integ.set_sigma(self.sigma)
        integ.set_k(self.k)

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

class testBIE(object):
    @classmethod
    def setUpClass(self):
        self.P = params(self.numpoints)
        self.P.initQuad(20)
        self.P.initEllipsoid()
        self.P.initPhi()
        self.P.initInteg()
        self.P.initAHMED()
        integ.set_q(self.P.q_ahmed)

    def testEllipsoid(self):
        self.assertAlmostEqual(
            center_points(self.P.node_coordinates),
            1.0e-11,
            places = 10)

    def testPhi(self):
        self.assertAlmostEqual(
            phi.test_phi(
                gen_points(self.P.numnodes,self.P.axes)),
            self.P.numnodes,
            places = 12)

    def testInteg(self):
        self.P.area[0] = sum(self.P.intphi_over)
        self.assertAlmostEqual(
            self.P.area[0],
            6.971610618375645,
            places = self.integ_places)

    def testSLAE(self):
        self.assertAlmostEqual(
            self.P.data.criteria(self.P.axes,integ.approximateu,self.P.data.exactu),
            self.slae_tol, places = self.slae_places)

class testBIEsmall(testBIE, unittest.TestCase):
    numpoints = 200
    integ_places = 4
    slae_tol = 0.0003
    slae_places = 4

class testBIEmedium(testBIE, unittest.TestCase):
    numpoints = 3200
    integ_places = 6
    slae_tol = 0.00019
    slae_places = 5

class testBIEbig(testBIE, unittest.TestCase):
    numpoints = 12800
    integ_places = 7
    slae_tol = 0.00009
    slae_places = 5

class testBIEhuge(testBIE, unittest.TestCase):
    numpoints = 25600
    integ_places = 6
    slae_tol = 0.00006
    slae_places = 5

class testBIEgig(testBIE, unittest.TestCase):
    numpoints = 51200
    integ_places = 6
    slae_tol = 0.00004
    slae_places = 5

class testFOO(unittest.TestCase):
    def testFoo(self):
        self.assertAlmostEqual(6.9716106130103963, 6.971610618375645, places = 7)
