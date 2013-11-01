#!/usr/bin/python
# -*- coding: utf-8 -*-
"""

"""

__author__ = "Evgeny Boykov. <artscan@list.ru>"

from max2f90 import *

values = open("jacobian.out").read()

beta = open("beta.f90", 'w')
x = open("x.f90", 'w')
jacobian = open("jacobian.f90", 'w')

ar = strtr(strtr(values,onestring),f90replace).split("\n")

beta.write(ar[2])
x.write(ar[1])
jacobian.write(ar[0])
