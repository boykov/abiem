#!/usr/bin/python
# -*- coding: utf-8 -*-
"""

"""

__author__ = "Evgeny Boykov. <artscan@list.ru>"

from numpy import *
import random
import math

def center_points(points):
    c = zeros((3))
    c[0] = sum(points[:,0])
    c[1] = sum(points[:,1])
    c[2] = sum(points[:,2])
    return sum(abs(c))

def test_phi(axes):
    a = [random.random(), random.random(), random.random()]
    a[0] = 2*(a[0] - 0.5)
    a[1] = 2*(a[1] - 0.5)*math.sqrt(1 - a[0]**2)
    a[2] = math.copysign(1.0, a[2] - 0.5)*math.sqrt(1 - a[0]**2 - a[1]**2)

    a[0] = a[0] * axes[0]
    a[1] = a[1] * axes[1]
    a[2] = a[2] * axes[2]

    return a

def gen_points(numnodes,axes):
    a = zeros((numnodes,3), order = 'Fortran')
    for i in range(0,numnodes):
        a[i,:] = test_phi(axes)[:]
    return a
