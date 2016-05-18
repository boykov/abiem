#!/usr/bin/python
# -*- coding: utf-8 -*-
"""

"""

__author__ = "Evgeny Boykov. <artscan@list.ru>"

import unittest
from datetime import datetime
from ellipsoid import ellipsoid
from numpy import *
import cmath
import math
from utilities import *
from phi import modphi as phi
from integ import modinteg as integ

import sys,os
import slaeahmed

import logging
from common import common
import savearrays

class params(common):
    def __init__(self, numpoints = 400):
        common.__init__(self)
        self.numpoints = numpoints
        self.axes[:] = self.listaxes[:]

    def eps_aggl(self):
        return self.eps_matgen

    def eps_gmres(self):
        return self.eps_matgen*10

    def initAHMED(self):
        """
        """
        self.points = zeros((self.numnodes,3)) # C order !!!
        self.points[:,:] = self.node_coordinates[:,:]

        slaeahmed.set_Hmatrix(self.eps_matgen,
                              self.eps_aggl(),
                              self.eps_gmres(),
                              self.steps_gmres,
                              self.eta,
                              self.bmin,
                              self.rankmax)

        slaeahmed.set_points(self.points)
        slaeahmed.set_vectorb(eval(self.name_vectorb))
        slaeahmed.set_kernel(eval(self.name_matrixa))
        slaeahmed.solve_slae()
        slaeahmed.get_q(self.q_density)

    def initPETSC(self):
        """
        """
        te = savearrays.TaskElement()
        te.initcover(self)
        te.save(te.name_savearrays)
        os.system("make -s solve")
        te = te.load(te.name_savearrays)
        self.q_density[:] = te.q[:]
        os.system("rm -f " + te.name_savearrays)

    def initQuad(self, orderquad):
        import scipy.special.orthogonal as op
        self.level2(orderquad)

        self.quadphi_over[:,0] = op.j_roots(self.orderquad,0,1)[0]
        self.quadphi_over[:,1] = op.j_roots(self.orderquad,0,1)[1]
        self.quadphi_under[:,0] = op.j_roots(self.orderquad,0,0)[0]
        self.quadphi_under[:,1] = op.j_roots(self.orderquad,0,0)[1]
        self.quadsingular[:,0] = op.j_roots(self.orderquad,2,-0.5)[0]
        self.quadsingular[:,1] = map(lambda y: y,
                                     op.j_roots(self.orderquad,2,-0.5)[1]/
                                     map(lambda x: (1-x)**2*(1+x)**(-0.5),
                                         op.j_roots(self.orderquad,2,-0.5)[0]))

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
        integ.set_l_ptr("qbx", (self.flagNeedQBX))
        integ.set_l_ptr("qbx_gauss6", (self.flagTestQBX_gauss6))
        integ.set_l_ptr("matrixa6_p", (self.name_matrixa == 'integ.matrixa6'))
        integ.set_l_ptr("use_int_neighbors_p", self.use_int_neighbors_p)

        integ.calcomp()
        self.sigma[:] = map(self.fsigma,self.intphi_over)[:]

        integ.setgauss()
        integ.calcsing()

    def fsigma(self,x):
        return math.sqrt(self.magic*x/(math.pi**2))

    def criteria(self,axes,f,g):
        koef = 0.8
        self.edata = ellipsoid(koef*axes,self.numpoints)
        shp = self.edata.points.shape
        adata = zeros(shp[0],dtype = complex)
        adata[:] = map(lambda(i):
                       self.epsilon(f,g,self.edata.points[i,:]),
                       range(0,shp[0],1))
        tmp = max(abs(adata[:]))
        return tmp

    def epsilon(self,f,g,x):
        # print "eps, ",g(x),f(x)
        return (g(x)-f(x))/abs(g(x))
    
class testBIE(object):
    @classmethod
    def setUpClass(self):
        tick = datetime.now()
        logging.basicConfig(level=logging.DEBUG)
        self.P = self.tmpP
        self.P.eps_gmres_ = self.P.eps_gmres()
        self.P.initQuad(self.P.orderquad)
        self.P.initEllipsoid()
        print "numnodes: ", self.P.numnodes
        self.P.initPhi()
        self.P.initInteg()
        if self.P.flagAHMED: self.P.initAHMED()
        if not self.P.flagAHMED: self.P.initPETSC()
        logging.debug("counter = " + str(self.P.counter/(self.P.numnodes**2)))
        tock = datetime.now()
        self.diff = tock - tick
        print "seconds: ", self.diff.seconds

    def testSeconds(self):
        self.assertLess(self.diff.seconds, self.P.test_seconds)

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
            self.P.criteria(self.P.axes,
                                 eval(self.P.name_approximateu),
                                 integ.exactu),
            self.P.slae_tol, places = self.P.slae_places)

    def testUnder(self):
        if self.P.flagTestUnder:
            self.assertAlmostEqual(
                sum(self.P.gauss[:,5])/(4*math.pi),
                self.P.gauss[self.P.numnodes - 1,3],
                places = self.P.under_places)

    def testQBX(self):
        if (self.P.flagNeedQBX and (not self.P.flagTestQBX_gauss6)):
            for nj in range(self.P.numnodes-1,self.P.numnodes):
                for i in range(0,self.P.numnodes):
                    for j in range(0,self.P.max_neighbors):
                        kj = self.P.node_neighbors2[nj,j]
                        if (kj-1 == i):
                            break
                    if ((kj == 0)):
                        self.assertAlmostEqual(
                                self.P.gauss[i,5]/(4*math.pi),
                                integ.foldingg(self.P.dim_intG,self.P.node_coordinates[nj,:],i+1,self.P.k_wave),
                                places = self.P.qbx_places)

class testBIEtest_sigm(testBIE, unittest.TestCase):
    tmpP = params(800)
    tmpP.k_wave = 6
    tmpP.name_approximateu = 'integ.approximateu_sigm'
    tmpP.name_matrixa = 'integ.matrixa_sigm'
    tmpP.slae_tol = 0.0036
    tmpP.slae_places = 3

class testBIEtest(testBIE, unittest.TestCase):
    tmpP = params(800)
    tmpP.k_wave = 6
    tmpP.slae_tol = 0.0036
    tmpP.slae_places = 3

class testBIEinteg(testBIE, unittest.TestCase):
    tmpP = params(200)
    tmpP.orderquad = 160
    tmpP.k_wave = 0.1
    tmpP.integ_places = 9
    tmpP.under_places = 8
    tmpP.flagTestUnder = True

class testBIEsmallNG(testBIE, unittest.TestCase):
    tmpP = params(200)
    tmpP.name_matrixa = 'integ.matrixa2'
    tmpP.name_vectorb = 'integ.vectorb2'
    tmpP.k_wave = 0.1
    tmpP.slae_tol = 0.009
    tmpP.slae_places = 3

class testBIEmediumQBX(testBIE, unittest.TestCase):
    tmpP = params(1600)
    tmpP.integ_places = 7
    tmpP.under_places = 7
    tmpP.k_wave = 0.1
    tmpP.orderquad = 30
    tmpP.dim_intG = 7
    tmpP.flagTestUnder = True
    tmpP.qbx_places = 7
    tmpP.flagNeedQBX = True
    tmpP.flagTestQBX_gauss6 = False
    tmpP.slae_tol = 0.0003
    tmpP.slae_places = 4

class testBIEsmall(testBIE, unittest.TestCase):
    tmpP = params(200)
    tmpP.integ_places = 5
    tmpP.under_places = 5
    tmpP.flagTestUnder = True
    tmpP.test_seconds = 15
    tmpP.slae_tol = 0.003
    tmpP.slae_places = 3

class testBIEsmall6(testBIE, unittest.TestCase):
    tmpP = params(200)
    tmpP.eps_matgen = 1e-7
    tmpP.test_seconds = 35
    tmpP.integ_places = 5
    tmpP.under_places = 4
    tmpP.k_wave = 0.1
    tmpP.name_approximateu = 'integ.approximateu4'
    tmpP.name_matrixa = 'integ.matrixa6'
    tmpP.qbx_places = 6
    tmpP.flagNeedQBX = True
    tmpP.flagTestUnder = True
    tmpP.flagTestQBX_gauss6 = True
    tmpP.slae_tol = 0.00005
    tmpP.slae_places = 5

class testBIEsmall6k(testBIE, unittest.TestCase):
    tmpP = params(200)
    tmpP.eps_matgen = 1e-4
    tmpP.test_seconds = 35
    tmpP.integ_places = 5
    tmpP.under_places = 4
    tmpP.k_wave = 8
    tmpP.name_approximateu = 'integ.approximateu4'
    tmpP.name_matrixa = 'integ.matrixa6'
    tmpP.qbx_places = 6
    tmpP.flagNeedQBX = True
    tmpP.flagTestUnder = True
    tmpP.flagTestQBX_gauss6 = True
    tmpP.slae_tol = 0.02
    tmpP.slae_places = 2

class testBIEmiddle6(testBIE, unittest.TestCase):
    tmpP = params(1600)
    tmpP.eps_matgen = 1e-8
    tmpP.integ_places = 7
    tmpP.under_places = 7
    tmpP.k_wave = 0.1
    tmpP.orderquad = 30
    tmpP.dim_intG = 7
    tmpP.name_approximateu = 'integ.approximateu4'
    tmpP.name_matrixa = 'integ.matrixa6'
    tmpP.qbx_places = 7
    tmpP.flagNeedQBX = True
    tmpP.flagTestUnder = True
    tmpP.flagTestQBX_gauss6 = True
    tmpP.slae_tol = 2e-6
    tmpP.slae_places = 6

class testBIEmiddle6k(testBIE, unittest.TestCase):
    tmpP = params(1600)
    tmpP.eps_matgen = 1e-5
    tmpP.test_seconds = 1600
    tmpP.integ_places = 7
    tmpP.under_places = 6
    tmpP.k_wave = 20
    tmpP.orderquad = 30
    tmpP.dim_intG = 7
    tmpP.name_approximateu = 'integ.approximateu4'
    tmpP.name_matrixa = 'integ.matrixa6'
    tmpP.qbx_places = 7
    tmpP.flagNeedQBX = True
    tmpP.flagTestUnder = True
    tmpP.flagTestQBX_gauss6 = True
    tmpP.slae_tol = 2e-2
    tmpP.slae_places = 2

class testBIEmedium6k(testBIE, unittest.TestCase):
    tmpP = params(3200)
    tmpP.eps_matgen = 1e-5
    tmpP.test_seconds = 3200
    tmpP.integ_places = 7
    tmpP.under_places = 6
    tmpP.k_wave = 28
    tmpP.orderquad = 30
    tmpP.dim_intG = 7
    tmpP.name_approximateu = 'integ.approximateu4'
    tmpP.name_matrixa = 'integ.matrixa6'
    tmpP.qbx_places = 7
    tmpP.flagNeedQBX = True
    tmpP.flagTestUnder = True
    tmpP.flagTestQBX_gauss6 = True
    tmpP.slae_tol = 2e-2
    tmpP.slae_places = 2

class testBIEbig6k(testBIE, unittest.TestCase):
    tmpP = params(12800)
    tmpP.eps_matgen = 1e-5
    tmpP.test_seconds = 22000
    tmpP.integ_places = 7
    tmpP.under_places = 7
    tmpP.k_wave = 56
    tmpP.orderquad = 30
    tmpP.dim_intG = 7
    tmpP.name_approximateu = 'integ.approximateu4'
    tmpP.name_matrixa = 'integ.matrixa6'
    tmpP.qbx_places = 7
    tmpP.flagNeedQBX = True
    tmpP.flagTestUnder = True
    tmpP.flagTestQBX_gauss6 = True
    tmpP.slae_tol = 2e-2
    tmpP.slae_places = 6

class testBIEmedium6(testBIE, unittest.TestCase):
    tmpP = params(3200)
    tmpP.eps_matgen = 1e-8
    tmpP.test_seconds = 3000
    tmpP.integ_places = 7
    tmpP.under_places = 7
    tmpP.k_wave = 0.1
    tmpP.orderquad = 30
    tmpP.dim_intG = 7
    tmpP.name_approximateu = 'integ.approximateu4'
    tmpP.name_matrixa = 'integ.matrixa6'
    tmpP.qbx_places = 7
    tmpP.flagNeedQBX = True
    tmpP.flagTestUnder = True
    tmpP.flagTestQBX_gauss6 = True
    tmpP.slae_tol = 8e-7
    tmpP.slae_places = 6

class testBIEbig6(testBIE, unittest.TestCase):
    tmpP = params(12800)
    tmpP.eps_matgen = 1e-8
    tmpP.test_seconds = 22000
    tmpP.integ_places = 7
    tmpP.under_places = 7
    tmpP.k_wave = 0.1
    tmpP.orderquad = 30
    tmpP.dim_intG = 7
    tmpP.name_approximateu = 'integ.approximateu4'
    tmpP.name_matrixa = 'integ.matrixa6'
    tmpP.qbx_places = 7
    tmpP.flagNeedQBX = True
    tmpP.flagTestUnder = True
    tmpP.flagTestQBX_gauss6 = True
    tmpP.slae_tol = 6e-8
    tmpP.slae_places = 6

class testBIEmicro6(testBIE, unittest.TestCase):
    tmpP = params(50)
    tmpP.integ_places = 2
    tmpP.under_places = 3
    tmpP.k_wave = 0.1
    tmpP.name_approximateu = 'integ.approximateu4'
    tmpP.name_matrixa = 'integ.matrixa6'
    tmpP.qbx_places = 8
    tmpP.flagNeedQBX = True
    tmpP.flagTestUnder = True
    tmpP.slae_tol = 0.01
    tmpP.slae_places = 2

class testBIEsmallQBX(testBIE, unittest.TestCase):
    tmpP = params(200)
    tmpP.integ_places = 5
    tmpP.under_places = 5
    tmpP.k_wave = 1
    tmpP.flagTestUnder = True
    tmpP.qbx_places = 6
    tmpP.flagNeedQBX = True
    tmpP.flagTestQBX_gauss6 = False
    tmpP.slae_tol = 0.007
    tmpP.slae_places = 3

class testBIEsmall3(testBIE, unittest.TestCase):
    tmpP = params(200)
    tmpP.name_matrixa = 'integ.matrixa3'
    tmpP.orderquad = 20
    tmpP.flagAHMED = False
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

class testBIEmedium3(testBIE, unittest.TestCase):
    tmpP = params(3200)
    tmpP.name_matrixa = 'integ.matrixa3'
    tmpP.integ_places = 6
    tmpP.test_seconds = 30
    tmpP.slae_tol = 0.0004
    tmpP.slae_places = 4

class testBIEmedium3k(testBIE, unittest.TestCase):
    tmpP = params(3200)
    tmpP.name_matrixa = 'integ.matrixa3'
    tmpP.k_wave = 7
    tmpP.integ_places = 6
    tmpP.test_seconds = 30
    tmpP.slae_tol = 0.007
    tmpP.slae_places = 3

class testBIEbig(testBIE, unittest.TestCase):
    tmpP = params(12800)
    tmpP.integ_places = 7
    tmpP.eps_matgen = 1e-6
    tmpP.slae_tol = 0.00008
    tmpP.slae_places = 5

class testBIEbig3k(testBIE, unittest.TestCase):
    tmpP = params(12800)
    tmpP.name_matrixa = 'integ.matrixa3'
    tmpP.k_wave = 14
    tmpP.integ_places = 7
    tmpP.slae_tol = 0.01
    tmpP.slae_places = 2

class testBIEhuge3k(testBIE, unittest.TestCase):
    tmpP = params(25600)
    tmpP.name_matrixa = 'integ.matrixa3'
    tmpP.k_wave = 21
    tmpP.eps_matgen = 1e-5
    tmpP.integ_places = 6
    tmpP.test_seconds = 450
    tmpP.slae_tol = 0.02
    tmpP.slae_places = 2

class testBIEhuge3(testBIE, unittest.TestCase):
    tmpP = params(25600)
    tmpP.name_matrixa = 'integ.matrixa3'
    tmpP.flagAHMED = False
    tmpP.eps_matgen = 1e-5
    tmpP.integ_places = 6
    tmpP.test_seconds = 450
    tmpP.slae_tol = 0.00007
    tmpP.slae_places = 5

class testBIEhuge(testBIE, unittest.TestCase):
    tmpP = params(25600)
    tmpP.integ_places = 6
    tmpP.eps_matgen = 1e-6
    tmpP.slae_tol = 0.00006
    tmpP.slae_places = 5

class testBIEgig(testBIE, unittest.TestCase):
    tmpP = params(51200)
    tmpP.integ_places = 6
    tmpP.eps_matgen = 1e-6
    tmpP.slae_tol = 0.00004
    tmpP.slae_places = 5

class testBIEgig3(testBIE, unittest.TestCase):
    tmpP = params(51200)
    tmpP.name_matrixa = 'integ.matrixa3'
    tmpP.test_seconds = 7000
    tmpP.flagAHMED = False
    tmpP.integ_places = 6
    tmpP.eps_matgen = 1e-6
    tmpP.slae_tol = 0.00003
    tmpP.slae_places = 5

class testFOO(unittest.TestCase):
    def testFoo(self):
        self.assertAlmostEqual(6.9716106130103963, 6.971610618375645, places = 7)
