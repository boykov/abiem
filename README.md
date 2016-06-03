Introduction
============

ABIEM (Another Boundary Integral Equations Method) is a prototype
program built for the two methods of the solution of 3-D Helmholtz BIE
on smooth sufaces. These methods are described in
[kashirin2010](http://www.icita.org/2010/papers/04-ru-Kashirin.pdf)
and [ying2006](http://mrl.nyu.edu/~dzorin/papers/ying2006h3b.pdf).

Installation Instructions
=========================

We offer the following installation options:

## Using a VirtualBox image of ABIEM

 * Download a recent version of [Oracle VirtualBox](https://www.virtualbox.org/) and install it. VirtualBox is a free and capable virtualization environment.
 * Download the [ABIEM VirtualBox image](https://github.com/boykov/abiem/releases/download/v1.0/abiem.ova)
 * Start VirtualBox and import the image using `File -> Import Appliance`.
 * Start the virtual machine on the main screen using the `Show` button. By holding `Shift` while clicking it the image will be started without an extra window appearing.
 * Once the machine is fully booted use the port `2222` in your ssh client to access the command line in the virtual image. The username is `eab` and the password is `eab`.

## Building ABIEM from scratch (coming soon)

Usage
=====

    cd /home/eab/git/difwave/abiem
    git pull
    make all
