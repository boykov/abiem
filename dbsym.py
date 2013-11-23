#!/usr/bin/python
# -*- coding: utf-8 -*-
"""

"""

__author__ = "Evgeny Boykov. <artscan@list.ru>"

from max2f90 import *

singular_out = open("singular.out").read()

singular = open("singular.f90", 'w')
singular2 = open("singular2.f90", 'w')

ar = strtr(strtr(singular_out,onestring),f90replace).split("\n")

singular2.write(ar[1])
singular.write(ar[0])

jacobian_out = open("jacobian.out").read()

beta = open("beta.f90", 'w')
x = open("x.f90", 'w')
jacobian = open("jacobian.f90", 'w')

ar = strtr(strtr(jacobian_out,onestring),f90replace).split("\n")

beta.write(ar[2])
x.write(ar[1])
jacobian.write(ar[0])

jacobian_out = open("jacobian2.out").read()

beta = open("beta2.f90", 'w')
x = open("x2.f90", 'w')
jacobian = open("jacobian2.f90", 'w')

ar = strtr(strtr(jacobian_out,onestring),f90replace).split("\n")

beta.write(ar[2])
x.write(ar[1])
jacobian.write(ar[0])

dirichlet_helmholtz_out = open("dirichlet-helmholtz.out").read()

Amn = open("Amn.f90", 'w')
Bmn = open("Bmn.f90", 'w')
limdA = open("limdA.f90", 'w')
dA = open("dA.f90", 'w')
limA = open("limA.f90", 'w')
A = open("A.f90", 'w')

ar = strtr(strtr(dirichlet_helmholtz_out,onestring),f90replace).split("\n")

Amn.write(ar[5])
Bmn.write(ar[4])
limdA.write(ar[3])
dA.write(ar[2])
limA.write(ar[1])
A.write(ar[0])
