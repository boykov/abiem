Introduction
============

ABIEM (Another Boundary Integral Equations Method) is a prototype
program built for the two methods of the solution of 3-D Helmholtz BIE
on smooth surfaces. These methods are described in
[kashirin2010](http://www.icita.org/2010/papers/04-ru-Kashirin.pdf)
and [ying2006](http://mrl.nyu.edu/~dzorin/papers/ying2006h3b.pdf).

Installation Instructions
=========================

We offer the following installation options:

## Using a Docker image of ABIEM

 * Download the [ABIEM Docker image](https://github.com/boykov/abiem/releases/download/v1.1/eab-abiem.tar.gz)

## Building ABIEM from scratch

 * All necessary commands for building the environment are
   [here](https://github.com/boykov/vagrant-box/blob/abeim/late_command.sh)
 * Requirements: [Maxima](maxima.sourceforge.net/), M. Bebendorf’s AHMED library,  [PETSc](https://www.mcs.anl.gov/petsc/) (optional)

Usage
=====

    cd /home/eab/git/difwave/abiem
    git pull
    make all

Docker
=====

    docker build . -t eab-abiem
    docker run -it --rm --privileged --hostname eab-abiem eab-abiem sh
    make test