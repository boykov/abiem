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

class params():
    def __init__(self, numpoints = 400, axes = [float(0.75),float(1.),float(0.5)]):
        self.max_neighbors = 100
        self.PI = 3.14159265358979324
        self.dim_3d = 3
        self.axes = zeros((3))
        self.numpoints = numpoints
        self.axes[:] = axes[:]
        self.k_wave = 0
        self.data = DataElement(numpoints)
        self.data.magic = 0.410
        self.name_matrixa = 'integ.matrixa'
        self.name_vectorb = 'integ.vectorb'
        self.name_approximateu = 'integ.approximateu'
        self.integ_places = 4
        self.slae_tol = 0.0
        self.slae_places = 0
        self.orderquad = 20
        self.eps_matgen = 1e-4
        self.eps_aggl = self.eps_matgen
        self.eps_gmres = self.eps_matgen*0.01
        self.eta = 0.8
        self.bmin = 15
        self.rankmax = 1000
        self.flagMemo = True

    def initAHMED(self):
        """
        """
        self.points = zeros((self.numnodes,3)) # C order !!!
        self.points[:,:] = self.node_coordinates[:,:]
        self.q_ahmed = zeros((self.numnodes), dtype = complex)

        slaeahmed.set_Hmatrix(self.eps_matgen,
                              self.eps_aggl,
                              self.eps_gmres,
                              self.eta,
                              self.bmin,
                              self.rankmax)

        slaeahmed.set_points(self.points)
        slaeahmed.set_vectorb(eval(self.name_vectorb))
        slaeahmed.set_kernel(eval(self.name_matrixa))
        slaeahmed.solve_slae()
        slaeahmed.get_q(self.q_ahmed)

    def initQuad(self, orderquad):
        self.data.orderquad = orderquad
        self.data.setupquad()
        self.dim_quad = self.data.orderquad
        self.quadphi_over = zeros((self.dim_quad,2),order = 'Fortran')
        self.quadphi_over[:,:] = self.data.quadphi_over[:,:]
        self.quadphi_under = zeros((self.dim_quad,2),order = 'Fortran')
        self.quadphi_under[:,:] = self.data.quadphi_under[:,:]
        self.quadsingular = zeros((self.dim_quad,2),order = 'Fortran')
        self.quadsingular[:,:] = self.data.quadsingular[:,:]

    def initEllipsoid(self):
        self.e = ellipsoid(self.axes, self.numpoints)
        self.numnodes = self.e.points.shape[0]
        self.hval = self.e.get_h()
        self.hval2 = self.hval * self.hval
        self.node_coordinates = zeros((self.numnodes,3), order = 'Fortran')
        self.normal_coordinates = zeros((self.numnodes,3), order = 'Fortran')
        self.node_coordinates[:,:] = self.e.points[:,:]
        self.normal_coordinates[:,:] = self.e.normalvectors[:,:]

    def initPhi(self):
        self.node_neighbors1 = zeros(
            (self.numnodes,self.max_neighbors), dtype = int32, order = 'Fortran')
        self.node_neighbors2 = zeros(
            (self.numnodes,self.max_neighbors), dtype = int32, order = 'Fortran')
        self.nstroke_coordinates = zeros(
            (self.numnodes,3), order = 'Fortran')

        self.setObjPhi(phi)

        phi.filter_neighbors(2*self.hval, self.node_neighbors2, self.numnodes)
        phi.filter_neighbors(  self.hval, self.node_neighbors1, self.numnodes)
        phi.normal_vector_stroke(self.numnodes, self.node_neighbors1)

    def withWrapMemo(self,largs, body, larrs, flagMemo):
        @memoize
        def withWrapMemoInternal(largs,body):
            logging.info(body + ' is saved to memo')
            eval(body)
            return larrs
        if flagMemo:
            logging.info(body + ' uses memo')
            return withWrapMemoInternal(largs,body)
        else:
            logging.info(body + ' does not use memo')
            eval(body)
            return larrs

    def initInteg(self):
        self.intphi_over = zeros((self.numnodes))
        self.intphi_under = zeros((self.numnodes))

        self.centres = zeros((self.dim_quad))
        self.weights = zeros((self.dim_quad))
        self.jacobian = zeros((4,self.dim_quad,4*self.dim_quad), dtype = complex, order = 'Fortran')
        self.nodes = zeros((4,self.dim_quad,4*self.dim_quad,3), order = 'Fortran')

        self.centres[:] = self.quadphi_over[:,0]
        self.weights[:] = self.quadphi_over[:,1]

        self.area = zeros((1))

        self.setObjInteg(integ)
        [self.intphi_over[:]] = self.withWrapMemo([self.axes,
                                                   self.node_coordinates,
                                                   self.dim_quad],
                                                  'integ.calcomp()',
                                                  [self.intphi_over[:]],
                                                  self.flagMemo)

        self.sigma = zeros((self.numnodes))
        self.sigma[:] = map(self.data.fsigma,self.intphi_over)[:]
        integ.set_sigma(self.sigma)
        integ.set_k_wave(self.k_wave)

        self.gauss = zeros((self.numnodes,10), dtype = complex, order = 'Fortran')
        integ.set_gauss(self.gauss)
        integ.setgauss()

        self.centres[:] = self.quadsingular[:,0]
        self.weights[:] = self.quadsingular[:,1]
        [self.gauss[:,3]] = self.withWrapMemo([self.axes,
                                               self.node_coordinates,
                                               self.dim_quad,
                                               self.k_wave],
                                              'integ.calcsing()',
                                              [self.gauss[:,3]],
                                              self.flagMemo)

    def setObjPhi(self,obj):
        obj.set_dim_3d(self.dim_3d)
        obj.set_pi(self.PI)
        obj.set_max_neighbors(self.max_neighbors)
        obj.set_hval(self.hval)
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
        obj.set_intphi_under(self.intphi_under)
        obj.set_hval2(self.hval2)
        obj.set_dim_quad(self.dim_quad)

        obj.set_centres(self.centres)
        obj.set_weights(self.weights)
        obj.set_jacobian(self.jacobian)
        obj.set_nodes(self.nodes)

        obj.set_area(self.area)

class testBIE(object):
    @classmethod
    def setUpClass(self):
        logging.basicConfig(level=logging.DEBUG)
        self.P = self.tmpP
        self.P.data.k = self.P.k_wave # TODO cleanup
        self.P.initQuad(self.P.orderquad)
        self.P.initEllipsoid()
        self.P.initPhi()
        self.P.initInteg()
        self.P.initAHMED()
        integ.set_q_density(self.P.q_ahmed)

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
            places = self.P.integ_places)

    def testSLAE(self):
        self.assertAlmostEqual(
            self.P.data.criteria(self.P.axes,
                                 eval(self.P.name_approximateu),
                                 self.P.data.exactu),
            self.P.slae_tol, places = self.P.slae_places)

class testBIEtest_sigm(testBIE, unittest.TestCase):
    tmpP = params(800)
    tmpP.k_wave = 6
    tmpP.name_approximateu = 'integ.approximateu_sigm'
    tmpP.name_matrixa = 'integ.matrixa_sigm'
    tmpP.slae_tol = 0.02
    tmpP.slae_places = 3

class testBIEtest(testBIE, unittest.TestCase):
    tmpP = params(800)
    tmpP.k_wave = 6
    tmpP.slae_tol = 0.02
    tmpP.slae_places = 3

class testBIEinteg(testBIE, unittest.TestCase):
    tmpP = params(200)
    tmpP.orderquad = 160
    tmpP.integ_places = 9

class testBIEsmallNG(testBIE, unittest.TestCase):
    tmpP = params(200)
    tmpP.name_matrixa = 'integ.matrixa2'
    tmpP.name_vectorb = 'integ.vectorb2'
    tmpP.k_wave = 0.1
    tmpP.slae_tol = 0.009
    tmpP.slae_places = 3

class testBIEsmall(testBIE, unittest.TestCase):
    tmpP = params(200)
    tmpP.integ_places = 5
    tmpP.slae_tol = 0.003
    tmpP.slae_places = 3

class testBIEsmall3(testBIE, unittest.TestCase):
    tmpP = params(200)
    tmpP.name_matrixa = 'integ.matrixa3'
    tmpP.orderquad = 20
    tmpP.integ_places = 5
    tmpP.slae_tol = 0.006
    tmpP.slae_places = 3

class testBIEmedium(testBIE, unittest.TestCase):
    tmpP = params(3200)
    tmpP.integ_places = 6
    tmpP.slae_tol = 0.0002
    tmpP.slae_places = 4

class testBIEbig(testBIE, unittest.TestCase):
    tmpP = params(12800)
    tmpP.integ_places = 7
    tmpP.slae_tol = 0.00008
    tmpP.slae_places = 5

class testBIEhuge3(testBIE, unittest.TestCase):
    tmpP = params(25600)
    tmpP.name_matrixa = 'integ.matrixa3'
    tmpP.integ_places = 6
    tmpP.slae_tol = 0.00007
    tmpP.slae_places = 5

class testBIEhuge(testBIE, unittest.TestCase):
    tmpP = params(25600)
    tmpP.integ_places = 6
    tmpP.slae_tol = 0.00006
    tmpP.slae_places = 5

class testBIEgig(testBIE, unittest.TestCase):
    tmpP = params(51200)
    tmpP.integ_places = 6
    tmpP.slae_tol = 0.00004
    tmpP.slae_places = 5

class testFOO(unittest.TestCase):
    def testFoo(self):
        self.assertAlmostEqual(6.9716106130103963, 6.971610618375645, places = 7)
