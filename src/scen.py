#!/usr/bin/python
# -*- coding: utf-8 -*-
"""

"""

__author__ = "Evgeny Boykov. <artscan@list.ru>"

from testcase import *

class s(testBIE, unittest.TestCase):
    pass

l = {'A':'integ.matrixa3',
     'data': [
         [200,0.0062,4],
         [400,0.0029,4],
         # [800,0.0015,4],
         # [1600,0.00073,5],
         # [3200,0.00037,5],
         # [6400,0.0002,5],
         ]}

for n,t,p in l['data']:
    s.tmpP = params(n)
    s.tmpP.name_matrixa = l['A']
    s.tmpP.slae_tol = t
    s.tmpP.slae_places = p
    suite = unittest.TestSuite()
    suite.addTest(unittest.TestLoader().loadTestsFromTestCase(s))
    unittest.TextTestRunner(verbosity=2).run(suite)

