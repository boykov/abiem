#!/usr/bin/python
# -*- coding: utf-8 -*-

import functools
import cPickle
import shelve
import inspect

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
