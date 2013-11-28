#!/usr/bin/python
# -*- coding: utf-8 -*-

import functools
import cPickle
import shelve
import inspect
# from ellipsoid import *
# from cover import cover
# from coverelement import DataElement, CoverElement
# from numpy import *


def memoize(fctn):
    memory = shelve.open("memo.dat")
    @functools.wraps(fctn)
    def memo(*args,**kwargs):
        haxh = cPickle.dumps((args, sorted(kwargs.iteritems())))

        if haxh not in memory:
            memory[haxh] = fctn(*args,**kwargs)

        return memory[haxh]

    if memo.__doc__:
        memo.__doc__ = "\n".join([memo.__doc__,"This function is memoized."])
    return memo

@memoize
def wrapellipsoid(axes,points):
    return ellipsoid(axes,points)

@memoize
def wrapcover(axes,points):
    d = DataElement(points)
    e = ellipsoid(axes,points)
    c = cover
    c.init_cover(e.points,e.normalvectors,e.get_h())
    c.set_intphi_over(d.quadphi_over,d.axes)
    c.calcarea(d.fb)
    ce = CoverElement()
    ce.initcover(c)
    return ce

if __name__ == "__main__":
    e = wrapcover(array([0.75,1,0.5]),100)
