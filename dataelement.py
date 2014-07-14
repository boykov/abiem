#!/usr/bin/python
# -*- coding: utf-8 -*-

from numpy import *
from ellipsoid import ellipsoid
import math
import cmath

class DataBase:
    def __init__(self,numpoints,axes):
        self.numpoints = numpoints
        self.orderquad = 4
        self.ordersing = 4
        self.slaetol = 1.e-6
        self.matshell = False
        self.memoCalcArea = False
        self.memoCalcSing = False
        self.axes = zeros((3))
        self.axes[:] = axes[:]
        self.scen = ""
        self.addinfo = ""
        self.uid = ""
        self.magic = None
        self.k = None
        self.isng = False
        self.p_use_intf = False
    def setupquad(self):
        import scipy.special.orthogonal as op
        self.quadphi_over = zeros((self.orderquad,2),order = 'Fortran')
        self.quadsingular = zeros((self.ordersing,2),order = 'Fortran')
        self.quadphi_under = zeros((self.orderquad,2),order = 'Fortran')
        self.centres03 = op.j_roots(self.orderquad,0,1)[0]
        self.C03 = op.j_roots(self.orderquad,0,1)[1]

        self.centres2q1 = op.j_roots(self.ordersing,2,-0.5)[0]
        self.C2q1 = map(lambda y: y,
                        op.j_roots(self.ordersing,2,-0.5)[1]/
                        map(lambda x: (1-x)**2*(1+x)**(-0.5),
                            op.j_roots(self.ordersing,2,-0.5)[0]))

        self.quadsingular[:,0] = self.centres2q1
        self.quadsingular[:,1] = self.C2q1
        self.quadphi_over[:,0] = self.centres03
        self.quadphi_over[:,1] = self.C03
        self.quadphi_under[:,0] = op.j_roots(self.orderquad,0,0)[0]
        self.quadphi_under[:,1] = op.j_roots(self.orderquad,0,0)[1]

class DataElement(DataBase):
    def __init__(self,numpoints,axes = [float(0.75),float(1.),float(0.5)]):
        DataBase.__init__(self,numpoints,axes)
        self.magic = 0.41444
        self.k = complex(1,0)

    def fsigma2(self,x):
        return math.sqrt(self.magic/(math.pi))*x

    def fsigma(self,x):
        return math.sqrt(self.magic*x/(math.pi**2))

    def centerpoints(self,points):
        c = zeros((3))
        c[0] = sum(points[:,0])
        c[1] = sum(points[:,1])
        c[2] = sum(points[:,2])
        return sum(abs(c))

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

    def fb(self,x,ret):
        ret[0] = cmath.exp(complex(0,1)*self.k*x[2])
        return

    def exactu(self,x):
        return cmath.exp(complex(0,1)*self.k*x[2])

    def exactu2(self,x):
        return cmath.exp(complex(0,1)*self.k*x[2])
