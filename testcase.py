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
from common import common

class params(common):
    def __init__(self, numpoints = 400, axes = [float(0.75),float(1.),float(0.5)]):
        common.__init__(self)
        self.numpoints = numpoints
        self.axes[:] = axes[:]
        self.data = DataElement(numpoints)
        self.data.magic = 0.410
        self.name_matrixa = 'integ.matrixa'
        self.name_vectorb = 'integ.vectorb'
        self.name_approximateu = 'integ.approximateu'
        self.integ_places = 4
        self.slae_tol = 0.0
        self.slae_places = 0
        self.under_places = 3
        self.qbx_places = 5
        self.orderquad = 20
        self.eps_matgen = 1e-4
        self.eps_aggl = self.eps_matgen
        self.eps_gmres = self.eps_matgen*0.01
        self.steps_gmres = 20*self.numpoints
        self.eta = 0.8
        self.bmin = 15
        self.rankmax = 1000
        self.flagMemo = False
        self.flagTestUnder = False
        self.flagTestQBX = False

    def initAHMED(self):
        """
        """
        self.points = zeros((self.numnodes,3)) # C order !!!
        self.points[:,:] = self.node_coordinates[:,:]

        slaeahmed.set_Hmatrix(self.eps_matgen,
                              self.eps_aggl,
                              self.eps_gmres,
                              self.steps_gmres,
                              self.eta,
                              self.bmin,
                              self.rankmax)

        slaeahmed.set_points(self.points)
        slaeahmed.set_vectorb(eval(self.name_vectorb))
        slaeahmed.set_kernel(eval(self.name_matrixa))
        slaeahmed.solve_slae()
        slaeahmed.get_q(self.q_density)

    def initQuad(self, orderquad):
        import scipy.special.orthogonal as op
        self.data.orderquad = orderquad
        self.data.setupquad()
        self.level2(self.data.orderquad)

        self.quadphi_over[:,0] = op.j_roots(self.orderquad,0,1)[0]
        self.quadphi_over[:,1] = op.j_roots(self.orderquad,0,1)[1]
        self.quadphi_under[:,0] = op.j_roots(self.orderquad,0,0)[0]
        self.quadphi_under[:,1] = op.j_roots(self.orderquad,0,0)[1]
        self.quadsingular[:,:] = self.data.quadsingular[:,:] # TODO extract


    def initEllipsoid(self):
        self.e = ellipsoid(self.axes, self.numpoints)
        self.level1(self.e.points.shape[0], self.e.get_h())
        self.node_coordinates[:,:] = self.e.points[:,:]
        self.normal_coordinates[:,:] = self.e.normalvectors[:,:]

    def initPhi(self):
        self.setObjPhi(phi)

        phi.filter_neighbors(2*self.hval, self.node_neighbors2, self.numnodes)
        phi.filter_neighbors(  self.hval, self.node_neighbors1, self.numnodes)
        phi.normal_vector_stroke(self.numnodes, self.node_neighbors1)

    def calcomp(self):
        integ.setup_calcomp()
        for i in range(0,self.numnodes,1):
            # broken by interface f as argument
            integ.integrate(self.nstroke_coordinates[i,:],
                            self.node_coordinates[i,:],
                            i+1,
                            self.centres,self.weights,1)
            self.intphi_over[i] = integ.folding(i+1,i+1,self.dim_quad,1)

    def initInteg(self):
        self.setObjInteg(integ)

        integ.set_l_ptr("gauss6", (self.name_matrixa == 'integ.matrixa6' or self.flagTestUnder))
        integ.set_l_ptr("qbx", (self.flagTestQBX))
        integ.set_l_ptr("matrixa6_p", (self.name_matrixa == 'integ.matrixa6'))

        integ.calcomp()
        self.sigma[:] = map(self.data.fsigma,self.intphi_over)[:]

        integ.setgauss()
        integ.calcsing()

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
        logging.debug("counter = " + str(self.P.counter))

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

    def testUnder(self):
        if self.P.flagTestUnder:
            self.assertAlmostEqual(
                sum(self.P.gauss[:,5])/(4*math.pi),
                self.P.gauss[self.P.numnodes - 1,3],
                places = self.P.under_places)

    def testQBX(self):
        if self.P.flagTestQBX:
            self.assertAlmostEqual(
                self.P.gauss[30,5]/(4*math.pi),
                integ.foldingg(5,self.P.node_coordinates[self.P.numnodes-1,:],31,self.P.k_wave),
                places = self.P.qbx_places)

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
    tmpP.under_places = 5
    tmpP.flagTestUnder = True
    tmpP.slae_tol = 0.003
    tmpP.slae_places = 3

class testBIEsmall6(testBIE, unittest.TestCase):
    tmpP = params(200)
    tmpP.integ_places = 5
    tmpP.under_places = 4
    tmpP.k_wave = 0.1
    tmpP.name_approximateu = 'integ.approximateu4'
    tmpP.name_matrixa = 'integ.matrixa6'
    tmpP.qbx_places = 8
    tmpP.flagTestQBX = True
    tmpP.flagTestUnder = True
    tmpP.slae_tol = 0.00006
    tmpP.slae_places = 5

class testBIEmicro6(testBIE, unittest.TestCase):
    tmpP = params(50)
    tmpP.integ_places = 2
    tmpP.under_places = 3
    tmpP.k_wave = 0.1
    tmpP.name_approximateu = 'integ.approximateu4'
    tmpP.name_matrixa = 'integ.matrixa6'
    tmpP.qbx_places = 8
    tmpP.flagTestQBX = True
    tmpP.flagTestUnder = True
    tmpP.slae_tol = 0.01
    tmpP.slae_places = 2

class testBIEsmallQBX(testBIE, unittest.TestCase):
    tmpP = params(200)
    tmpP.integ_places = 5
    tmpP.under_places = 4
    tmpP.k_wave = 1
    tmpP.flagTestUnder = True
    tmpP.qbx_places = 8
    tmpP.flagTestQBX = True
    tmpP.slae_tol = 0.007
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
    tmpP.under_places = 6
    tmpP.flagTestUnder = True
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
