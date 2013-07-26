#!/usr/bin/python
# -*- coding: utf-8 -*-
"""

"""

__author__ = "Evgeny Boykov. <artscan@list.ru>"

import unittest
from ellipsoid import ellipsoid
from numpy import *
from utilities import *

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
        from phi import modphi as phi
        self.e = ellipsoid(self.list_axes, self.number_points)
        # TODO add nneigh1, nneigh2
        shp = self.e.points.shape[0]
        nneigh1 = zeros((shp,100), order = 'Fortran')
        nneigh2 = zeros((shp,100), order = 'Fortran')
        phi.init_phi(self.e.points, self.e.normalvectors, self.e.get_h(), nneigh1, nneigh2, self.list_axes)
        self.assertAlmostEqual(phi.testphijac(), shp, places = 12)
