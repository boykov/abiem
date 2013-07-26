#!/usr/bin/python
# -*- coding: utf-8 -*-
"""

"""

__author__ = "Evgeny Boykov. <artscan@list.ru>"

import unittest
from ellipsoid import ellipsoid
from numpy import *
from utilities import *

class testEllipsoid(unittest.TestCase):
    def setup(self):
        pass

    def testEllipsoid(self):
        self.number_points = 400
        self.list_axes = zeros((3), order = 'Fortran')
        self.list_axes[:] = [float(0.75),float(1.),float(0.5)]
        e = ellipsoid(self.list_axes, self.number_points)
        self.assertAlmostEqual(center_points(e.points), 1.0e-11, places = 10)
