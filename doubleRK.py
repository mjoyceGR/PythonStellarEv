#!/usr/bin/env python
##################################################################
#
# Author: Meridith Joyce
# Stellar Astrophysics, Dartmouth College
#
##################################################################
import math
import timeit
import numpy as np
import matplotlib.pyplot as plt

do_integrate=True
n=1.5

def rKN(x, theta, z, fx, gx, h):	#no iters
    k0=fx(x, theta, z)*h
    L0=gx(x, theta, z)*h

    k1=fx(x + 0.5*h, theta + 0.5*k0, z + 0.5*L0)*h
    L1=gx(x + 0.5*h, theta + 0.5*k0, z + 0.5*L0)*h

    k2=fx(x + 0.5*h, theta + 0.5*k1, z + 0.5*L1)*h
    L2=gx(x + 0.5*h, theta + 0.5*k1, z + 0.5*L1)*h

    k3=fx(x + h, theta + k2, z + L2)*h
    L3=gx(x + h, theta + k2, z + L2)*h

    x     = x     + h
    theta = theta + (k0 + 2*k1 + 2*k2 + k3)/6.0
    z     = z     + (L0 + 2*L1 + 2*L2 + L3)/6.0
    return x,theta,z

def LaneEmden_1(xi, theta, dtdxi):
    z = dtdxi
    #dzdxi = -(1.0/xi**2.0) * (2*xi*z + xi**2.0*theta**n)
    return z

def LaneEmden_2(xi, theta, dtdxi): #don't take n for now
    n=1.5
    z = dtdxi
    if xi == 0 or xi < 0 or theta < 0:
        dzdxi = 0
    else:
        dzdxi = -(1.0/xi**2.0) * (2*xi*z + xi**2.0*theta**n)
    return dzdxi

def LaneEmdenSolver(step,iters): 
    x1    = 0      ## ???
    theta1= 1      #likewise
    z1    = 0
    print >> outf,   "## xi\t\ttheta\t\tdtdxi "

    #for i in range(iters):
    while theta1 > theta_bound:
        x1, theta1, z1 = rKN(x1, theta1, z1, LaneEmden_1, LaneEmden_2, step)#,1,2] ## vectorized
        print >> outf, x1 ,"\t\t", theta1, "\t", z1
    #When done, print xi value corresponding to first theta ~ 0
    print "The first zero is: ", x1
    return x1

iters = 1000
theta_bound = 0
if do_integrate ==True:
    step  = 10**(-4)
    outf = open("theta_xi.dat","w")
    x1 = LaneEmdenSolver(step,iters)
    outf.close()
else:
    pass

xi, theta, dtdxi=np.loadtxt("theta_xi.dat",usecols=(0,1,2),unpack=True)
plt.plot(xi, theta,  "g-",label=r'$\theta(\xi)$')
plt.plot(xi, dtdxi,  "m-",label=r"$\frac{d\theta}{d\xi}$")
plt.axvline(x = x1, color='red',linestyle='--')
plt.xlabel(r'$\xi$')
plt.ylabel('Lane-Emden Solutions')
plt.legend(loc=1)
plt.savefig("LaneEmdenSolutions.png")
#plt.show()
plt.close()