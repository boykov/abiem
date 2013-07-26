#!/usr/bin/python
# -*- coding: utf-8 -*-
"""

"""

__author__ = "Evgeny Boykov. <artscan@list.ru>"

from numpy import *

def center_points(points):
    c = zeros((3))
    c[0] = sum(points[:,0])
    c[1] = sum(points[:,1])
    c[2] = sum(points[:,2])
    return sum(abs(c))
