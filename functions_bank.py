########################################################################
#
# Module functions_bank.py
#
# This module contains all of the numerical routines used in this project.
# It is imported in main.py and all routines are executed there.
# The routines are sorted with headers by which problem they address
#
# Author: Meridith Joyce
# Adv Stellar Astrophysics, Dartmouth College
#
########################################################################
import math
import numpy as np
from constants_bank import *
import os
import subprocess

########################################################################
#
# Routines for solving dimensionless Lane-Emden
#
########################################################################
## 2-D Runge-Kutta
def rKN(x, theta, z, fx, gx, h):
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

## First Lane-Emden reduced ODE
def LaneEmden_1(xi, theta, dtdxi):
    z = dtdxi
    return z

## Second Lane-Emden reduced ODE
def LaneEmden_2(xi, theta, dtdxi): 
    z = dtdxi
    if xi == 0 or xi < 0 or theta < 0:
        dzdxi = 0
    else:
        dzdxi = -(1.0/xi**2.0) * (2*xi*z + xi**2.0*theta**n)
    return dzdxi

## Routine that applies 2-D Runge-Kutta to the two Lane-Emden ODEs
def LaneEmdenSolver(step,outf,theta_bound): 
    print "\n\n>> Solving Lane-Emden <<\n"
    x1    = 0      
    theta1= 1     
    z1    = 0
    while theta1 > theta_bound:
        x1, theta1, z1 = rKN(x1, theta1, z1, LaneEmden_1, LaneEmden_2, step)
        print >> outf, x1 ,"\t\t", theta1, "\t", z1
    print "The first zero (xi_1) is: ", x1
    print "dtheta/dxi at xi_1 is: ", z1
    return x1, theta1, z1

## this calls rk2D.f90, which performs the integration in fortran instead. It is not currently 
## used in main.py because the python values are better, but the routine works!
def LaneEmden_fortran(x1_1, theta_1, z_1, stepsize, polyindex ):
    inf = open("fortran_input.txt","w")
    print >> inf, x1_1, theta_1, z_1, stepsize, polyindex
    inf.close()
    subprocess.call('gfortran rk2d.f90 -o LE.out',shell=True)
    subprocess.call('./LE.out',shell=True)
    return

########################################################################
#
# Physical Routines for pre-main sequence models
#
########################################################################
## converts dimensionless quantities to density
def Density(theta,rho_c):
    rho = rho_c*theta**n
    return rho
## converts dimensionless quantities to pressure
def Pressure(K,rho):
    P = K*rho**((n+1)/n)
    return P
## converts dimensionless quantities to temperature
def Temperature(P,rho):
    T = mu*mp*P/(kb*rho)
    return T   

########################################################################
#
# Quantum Gas Functions
#
########################################################################
## computes pressure from rho in the  non-relativistic regime
def PN(rho):
    PN = KN*rho**(5.0/3.0)
    return PN
## computes pressure from rho in the relativistic regime
def PR(rho):
    PR = KR*rho**(4.0/3.0)
    return PR
## the combined approximation using PN and PR
def Pcomb(rho, PN, PR):
    P = PN*PR/np.sqrt(PN**2.0 + PR**2.0)
    return P
## the combined approximation using rho_0
def Pcomb_otherform(rho): 
    P = KN*rho**(5.0/3.0)/np.sqrt(1.0 + (rho/rho_0)**(2.0/3.0) )
    return P

########################################################################
#
# White Dwarf Routines
#
########################################################################
## quantum gas analog to ODE 1
def wd_1(xi, theta, dtdxi):
    V = dtdxi
    return V
## quantum gas analog to ODE 2
def wd_2(xi, theta, dtdxi): 
    V = dtdxi
    if xi == 0 or xi < 0 or theta < 0:
        dzdxi = 0
    else:
        gamma = 1.0 + theta**(2.0/3.0)
        Bterm1 = (5.0/3.0)*theta**(-1.0/3.0)*(gamma)**(-0.5)
        Bterm2 = (1.0/3.0)*theta**(1.0/3.0)*(gamma)**(-3.0/2.0)
        B = Bterm1 - Bterm2
       
        Cterm1 = (5.0/9.0)*theta**(-4.0/3.0) *(gamma)**(-0.5) 
        Cterm2 = (2.0/3.0)*theta**(-2.0/3.0) *(gamma)**(-3.0/2.0) 
        Cterm3 = (1.0/3.0)*(gamma)**(-5.0/2.0)
        C = -Cterm1 - Cterm2 + Cterm3
        dzdxi = (-2.0/xi)*V - (C/B)* V**2.0 - (1.0/B)*theta
    return dzdxi

## 2-D Runge-Kutta including Boole's Rule for numerical integration in order to keep track of m-integral
def RK_withM(x, theta, z, fx, gx, h,m):  
    k0=fx(x, theta, z)*h
    L0=gx(x, theta, z)*h

    k1=fx(x + 0.5*h, theta + 0.5*k0, z + 0.5*L0)*h
    L1=gx(x + 0.5*h, theta + 0.5*k0, z + 0.5*L0)*h

    k2=fx(x + 0.5*h, theta + 0.5*k1, z + 0.5*L1)*h
    L2=gx(x + 0.5*h, theta + 0.5*k1, z + 0.5*L1)*h

    k3=fx(x + h, theta + k2, z + L2)*h
    L3=gx(x + h, theta + k2, z + L2)*h

    ## Boole's method for Numerical Integration
    b = x + h
    a = x 
    f0 = a
    f1 = a + (b-a)/4.0
    f2 = a + (b-a)/2.0
    f3 = a + 3.0*(b-a)/4.0
    f4 = b
    mnew = theta*((b-a)/90.0)*(7.0*f0**2.0 + 32.0*f1**2.0 + 12.0*f2**2.0 + 32.0*f3**2.0 + 7.0*f4**2.0)
    m = m + mnew

    x     = x     + h
    theta = theta + (k0 + 2*k1 + 2*k2 + k3)/6.0
    z     = z     + (L0 + 2*L1 + 2*L2 + L3)/6.0
    return x,theta,z, m

## Routine that applies Runge-Kutta + Boole to the two quantum gas analog ODEs
def wd_physical(theta_c, step,filename):
    print "\n\n>> generating white dwarf models <<\n"
    theta_0 = theta_c
    outf=open(filename,"w")
    print >> outf, "#m\t\txi\t\ttheta\t\tdtheta dxi"

    #Initial Conditions
    x1     = 0.0     
    theta1 = theta_c      
    z1     = 0.0
    m      = 0.0

    #stop when theta/theta_c = 0.001 (0.1%) to avoid numerical difficulties
    while theta1/theta_0 >= 0.001: 
        x1, theta1, z1, m = RK_withM(x1, theta1, z1, wd_1, wd_2, step,m)
        print >> outf,  m, "\t", x1, "\t", theta1, "\t", z1
    outf.close()

    print "m: ", m, "\tMsolar: ", Msolar, "\ttheta_c/theta: ", theta1/theta_0
    #alternative: totalM = 0.090*m*Msolar; they are equivalent
    totalM = 4.0*pi*rho_0*(alpha_wd**3.0)*m
    totalR = alpha_wd*x1
    return x1, theta1, z1, totalM, totalR