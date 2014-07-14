#!/usr/bin/python
# -*- coding: utf-8 -*-
"""

"""

__author__ = "Evgeny Boykov. <artscan@list.ru>"

import math
import cmath
import pickle
import shelve
import sys,os
from numpy import *

from integ import modinteg as integ

class CoverElement(object):
    def initcover(self, c):
        self.numnodes = c.numnodes
        self.hval = c.hval
        self.node_coordinates = c.node_coordinates
        self.intphi_over = c.intphi_over
        self.name_matrixa = c.name_matrixa
        self.eps_gmres_ = c.eps_gmres_
        self.k_wave = c.k_wave
        self.gauss = c.gauss
        self.sigma = c.sigma
        self.dim_3d = c.dim_3d

class TaskElement(CoverElement):
    name_savearrays = 'savearrays.dat'

    def save(self,fname):
        f = open(fname,'w')
        pickle.dump(self,f)
        f.close()

    @classmethod
    def load(self,fname):
        f = open(fname,'r')
        ce = pickle.load(f)
        f.close()
        return ce

    def setupkernel(self):
        """
        Linking with fortran object
        """
        self.k = integ

        self.k.set_i_ptr("numnodes", self.numnodes)
        self.k.set_i_ptr("dim_3d", self.dim_3d)
        self.k.set_dp_ptr("hval", self.hval)
        self.k.set_dc_ptr("k_wave", self.k_wave)
        self.k.set_dc2d_ptr("gauss", self.gauss)
        self.k.set_dp2d_ptr("node_coordinates", self.node_coordinates)
        self.k.set_dp1d_ptr("intphi_over", self.intphi_over)
        self.k.set_dp1d_ptr("sigma", self.sigma)

    def killfortran(self):
        if (hasattr(self,'k')):
            del self.k
