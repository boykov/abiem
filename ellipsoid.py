#!/usr/bin/python
# -*- coding: utf-8 -*-
"""Sampling ellipsoid

Filename: ellipsoid.py

Copyright (C) 2012  artscan

Author: Evgeny Boykov <artscan@list.ru>
Keywords: ellipsoid, sampling
URL: 
Commentary:

python-версия алгоритма дискретизации эллипсоида, ранее реализованного
в программах решения интегральных уравнений на языке fortran.
Изменения, в т.ч. и в проверке int(currentEllipseParts) % 2,
позволяющей использовать четное количество точек для симметрии.
Позволяет обеспечить несмещенное расположение центра тяжести.
"""

__author__ = "Evgeny Boykov. <artscan@list.ru>"


import math
from numpy import *

class ellipsoid(object):
    
    def __init__(self, axes, parts):
        # self.axes = zeros(3,order = 'Fortran')
        self.axes = axes
        self.parts = parts
        self.numpoints = self.samplingEllipsoid()
        self.normalVectors()
        self.radius_hat = 2 # 0.67
        
        
    def normalVectors(self):
        """Compute normal vector for each point from ellipsoid sampling points"""
        axesReverse = zeros(3)
        axesReverse[:] = 1/self.axes[:]
        self.normalvectors[0:self.numpoints,0:3] = 0.
        X = zeros(3)
        for i in range(0,self.numpoints,1):
            X[0:3] = self.points[i,0:3] * (axesReverse[0:3]) * (axesReverse[0:3])
            normXReverse = 1./ math.sqrt(X[0]*X[0]+X[1]*X[1]+X[2]*X[2])
            self.normalvectors[i,0:3] = X[0:3] * normXReverse
        return


    def get_h(self):
        sqrt_np = self.ellipseParts()
        return self.radius_hat * 2.*self.axes[1]*self.ep2p(1.-(self.axes[0]/self.axes[1])**2)/sqrt_np


    def ellipseParts(self):
        return int(0.5 + math.sqrt(self.parts))

    def divideOne(self,eccentricity,parts,angle):
        n = int(parts)
        addon = eccentricity * (parts - n)

        angle[0] = 0.
        for i in range(0,n,1):
            angle[i+1] = angle[i] + 1.+ (eccentricity - 1.)*abs((2.*(i+1)-1.)/parts-1.)
        
        coef = 1 / (angle[n] + addon)
        angle[0:n+1] = angle[0:n+1]*coef

    def completeArraySym(self,a,parts):
        rn = 0.5 * parts
        n = int(rn)
        j = n
        for i in range(int(parts) - n,int(parts),1):
            a[i,0] = - a[j,0]
            a[i,1] =   a[j,1]
            j = j - 1

    def ep2p(self,m):
        """
        Вычисление полного эллиптического интеграла 2-го рода
        """
        a  = array([0.44325141463e0,0.06260601220e0,0.04757383546e0,0.01736506451e0])
        b  = array([0.24998368310e0,0.09200180037e0,0.04069697526e0,0.00526449639e0])
        m1 = 1. - m
        c  = m1
        s1 = 1.
        s2 = 0.
        for i in range(0,4,1):
            s1 = s1 + a[i] * c
            s2 = s2 + b[i] * c
            c  = c * m1
        return s1 - log(m1) * s2

    def samplingEllipse(self,axes,parts,points):
        n = int(parts) + 1
        angle = zeros(n)

        eccentricity = axes.max()/axes.min()

        self.divideOne(eccentricity,parts,angle)

        points[1:n,0] = axes.min() * sin(math.pi * angle[1:n])
        points[1:n,1] = axes.max() * cos(math.pi * angle[1:n])        

        points[0,1] = axes.max()
        points[0,0] = 0

    def samplingEllipsoid(self):
            
        mainparts = self.ellipseParts()

        p5 = 5 * mainparts
        p7 = 7 * mainparts

        currentparts = zeros(mainparts + 1,integer)
        currentparts[0]         = 1
        currentparts[mainparts] = 1

        mainpoints = zeros((p5,2))
        fullpoints = zeros((p7,p5,2))
        currentpoints = zeros((p7,2))        

        mainaxes = array([self.axes[0],self.axes[1]])
        
        self.samplingEllipse(mainaxes,mainparts,mainpoints)

        currentaxes = zeros(2)
        
        for i in range(1,mainparts,1):
            angle = math.sqrt(1. - (mainpoints[i,1] / self.axes[1])**2)

            currentaxes[1] = self.axes[2] * angle
            currentaxes[0] = self.axes[0] * angle

            currentEllipseParts = float(int(0.5 + 2.* mainparts * angle * self.axes[0] * self.ep2p(1.- (self.axes[2] / self.axes[0]) ** 2) / (self.axes[1] * self.ep2p(1.- (self.axes[0] / self.axes[1]) ** 2))))

            if int(currentEllipseParts) % 2 <> 0: currentEllipseParts = currentEllipseParts - 1

            
            self.samplingEllipse(currentaxes,currentEllipseParts / 2,currentpoints)

            self.completeArraySym(currentpoints,currentEllipseParts)

            currentparts[i] = currentEllipseParts

            fullpoints[0,i,1] = currentaxes[0]
            fullpoints[1:currentEllipseParts,i,1] = currentpoints[1:currentEllipseParts,1]
            fullpoints[1:currentEllipseParts,i,0] = currentpoints[1:currentEllipseParts,0]

        totalParts = 0
        for i in range(0,mainparts+1,1):
            currentNumberPoints = currentparts[i]
            if currentNumberPoints == 0:
                continue
            for j in range(0,currentNumberPoints,1):
                totalParts = totalParts + 1

        self.points = zeros((totalParts,3),order = 'Fortran')
        self.normalvectors = zeros((totalParts,3),order = 'Fortran')        

        totalParts = 0
        for i in range(0,mainparts+1,1):
            currentNumberPoints = currentparts[i]
            if currentNumberPoints == 0:
                continue
            for j in range(0,currentNumberPoints,1):
                self.points[totalParts,0] = fullpoints[j,i,1]
                self.points[totalParts,1] = mainpoints[i,1]
                self.points[totalParts,2] = fullpoints[j,i,0]
                totalParts = totalParts + 1

        return totalParts
        
