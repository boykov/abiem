#!/usr/bin/python
# -*- coding: utf-8 -*-
"""

"""

__author__ = "Evgeny Boykov. <artscan@list.ru>"

import optparse
from max2f90 import *

def spec():
    spec_out = open("specfun.out").read()
    spherical_harmonic = open("spherical_harmonic.f90", 'w')
    spherical_bessel_y = open("spherical_bessel_y.f90", 'w')
    spherical_bessel_j = open("spherical_bessel_j.f90", 'w')
    ar = strtr(strtr(spec_out,onestring),f90replace).split("\n")
    spherical_harmonic.write(ar[2])
    spherical_bessel_y.write(ar[1])
    spherical_bessel_j.write(ar[0])

def xx_rhopy():
    xx_rho_out = open("xx_rho.out").read()
    xx_rho3 = open("xx_rho3.f90", 'w')
    xx_rho2 = open("xx_rho2.f90", 'w')
    xx_rho1 = open("xx_rho1.f90", 'w')
    ar = strtr(strtr(xx_rho_out,onestring),f90replace).split("\n")
    xx_rho3.write(ar[2])
    xx_rho2.write(ar[1])
    xx_rho1.write(ar[0])

def singularpy():
    singular_out = open("singular.out").read()
    singular = open("singular.f90", 'w')
    singular2 = open("singular2.f90", 'w')
    singular3 = open("singular3.f90", 'w')
    ar = strtr(strtr(singular_out,onestring),f90replace).split("\n")
    singular3.write(ar[2])
    singular2.write(ar[1])
    singular.write(ar[0])

def jacobianpy():
    jacobian_out = open("jacobian.out").read()
    beta = open("beta.f90", 'w')
    x = open("x.f90", 'w')
    jacobian = open("jacobian.f90", 'w')
    ar = strtr(strtr(jacobian_out,onestring),f90replace).split("\n")
    beta.write(ar[2])
    x.write(ar[1])
    jacobian.write(ar[0])

def jacobian2py():
    jacobian_out = open("jacobian2.out").read()
    beta = open("beta2.f90", 'w')
    x = open("x2.f90", 'w')
    jacobian = open("jacobian2.f90", 'w')
    ar = strtr(strtr(jacobian_out,onestring),f90replace).split("\n")
    beta.write(ar[2])
    x.write(ar[1])
    jacobian.write(ar[0])

def dirichlet_helmholtzpy():
    dirichlet_helmholtz_out = open("dirichlet-helmholtz.out").read()
    Amn = open("Amn.f90", 'w')
    Bmn = open("Bmn.f90", 'w')
    limdA = open("limdA.f90", 'w')
    dA = open("dA.f90", 'w')
    limA = open("limA.f90", 'w')
    A = open("A.f90", 'w')
    A_s = open("A_s.f90", 'w')
    cderf = open("cderf.f90", 'w')
    ar = strtr(strtr(dirichlet_helmholtz_out,onestring),f90replace).split("\n")
    cderf.write(ar[7])
    A_s.write(ar[6])
    Amn.write(ar[5])
    Bmn.write(ar[4])
    limdA.write(ar[3])
    dA.write(ar[2])
    limA.write(ar[1])
    A.write(ar[0])

parser =  optparse.OptionParser()
(options, args) = parser.parse_args()
scen = args[0]

exec(scen + '()')
