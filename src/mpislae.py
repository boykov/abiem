#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys,os

from savearrays import TaskElement
import petsc4py

petsc4py.init(sys.argv)
from petsc4py import PETSc

class Del2Mat:

    def __init__(self):
        pass

    def create(self, A):
        pass

    def mult(self, A, x, y):
        "y <- A * x"

        scatter, x0 = PETSc.Scatter.toAll(x)
        scatter.scatter(x, x0, False, PETSc.ScatterMode.FORWARD)

        for m in range(0,len(y[...])):
            ind = y.getOwnershipRange()[0]
            if te.name_matrixa == "integ.matrixa3":
                y[...][m] = te.k.dota3(x0,ind+m)
            if te.name_matrixa == "integ.matrixa":
                y[...][m] = te.k.dota(x0,ind+m)

    def multTranspose(self, A, x, y):
        "y <- A' * x"
        pass

    def getDiagonal(self, A, D):
        "D[i] <- A[i,i]"
        pass

te = TaskElement.load(TaskElement.name_savearrays)
te.setupkernel()

A = PETSc.Mat().create()
A.setSizes([te.numnodes,te.numnodes])
A.setType('python')
shell = Del2Mat() # shell context
A.setPythonContext(shell)

A.setUp()

A.assemblyBegin()
A.assemblyEnd()

x,b = A.createVecs()

u = PETSc.Vec()
u.create(PETSc.COMM_WORLD)
u.setSizes(te.numnodes)
u.setType(PETSc.Vec.Type.MPI)

map(lambda i:u.setValue(i,te.k.vectorb(i+1)),range(0,te.numnodes,1))


none   = 1.0
b.axpy(none,u)


ksp = PETSc.KSP().create(PETSc.COMM_WORLD)
pc = ksp.getPC()
ksp.setType('gmres')
pc.setType('none')
ksp.setTolerances(te.eps_gmres_)

ksp.setOperators(A)
ksp.setFromOptions()

ksp.solve(b, x)

b1 = b.duplicate()
A.mult(x,b1)
none = - 1.0
b1.axpy(none,b)
norm = b1.norm(PETSc.NormType.NORM_2)
its = ksp.getIterationNumber()


scatter, U0 = PETSc.Scatter.toZero(x)
scatter.scatter(x, U0, False, PETSc.ScatterMode.FORWARD)

myrank = x.comm.rank
if myrank == 0:
    print "Norm of error ",norm,",iterations ",its
    te.q = U0.array
    te.killfortran()
    te.save(te.name_savearrays)

PETSc.COMM_WORLD.barrier()
scatter.destroy()
