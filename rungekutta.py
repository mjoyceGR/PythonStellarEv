########################################################################
#
# Numerical Routines for Lane-Emden
#
########################################################################
import math
import numpy as np
from constants_bank import *

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
    return z

def LaneEmden_2(xi, theta, dtdxi): #don't take n for now
    z = dtdxi
    if xi == 0 or xi < 0 or theta < 0:
        dzdxi = 0
    else:
        dzdxi = -(1.0/xi**2.0) * (2*xi*z + xi**2.0*theta**n)
    return dzdxi

def LaneEmdenSolver(step,iters,outf,theta_bound): 
    x1    = 0      ## ???
    theta1= 1      #likewise
    z1    = 0
    while theta1 > theta_bound:
        x1, theta1, z1 = rKN(x1, theta1, z1, LaneEmden_1, LaneEmden_2, step)
        print >> outf, x1 ,"\t\t", theta1, "\t", z1
    print "The first zero (xi_1) is: ", x1
    print "dtheta/dxi at xi_1 is: ", z1
    return x1, theta1, z1

########################################################################
#
# Physical Routines for pre-main sequence models
#
########################################################################
def Density(theta,rho_c):
    rho = rho_c*theta**n
    return rho

def Pressure(K,rho):
    #Pressure(rho,T)
    #P = (kb/(mu*mp))*rho*T   
    P = K*rho**((n+1)/n)
    return P

def Temperature(P,rho):
    #Temperature(K,theta, rho_c)
    #T = ( K*mu*mp/kb) *theta * rho_c**(1.0/n)
    T = mu*mp*P/(kb*rho)
    return T   

########################################################################
#
# White Dwarf Routines
#
########################################################################
def wd_1(xi, theta, dtdxi):
    V = dtdxi
    return V

def wd_2(xi, theta, dtdxi): #take B and C?
    V = dtdxi
    if xi == 0 or xi < 0 or theta < 0:
        dzdxi = 0
    else:
        Bterm1 = 5.0*theta**(-1.0/3.0)*(1.0 + theta**(2.0/3.0))**(-0.5)/3.0
        Bterm2 = theta**(1.0/3.0)*(1+theta**(2.0/3.0))**(-3.0/2.0)/3.0
        B = Bterm1 - Bterm2
       
        Cterm1 = (5.0/9.0)*theta**(-4.0/3.0) *(1 + theta**(2.0/3.0))**(-0.5) 
        Cterm2 = (2.0/3.0)*theta**(-2.0/3.0) *(1 + theta**(2.0/3.0))**(-3.0/2.0) 
        Cterm3 = (1.0/3.0)*(1 + theta**(2.0/3.0))**(-5.0/2.0)
        C = -Cterm1 - Cterm2 + Cterm3

        dzdxi = (-2.0/xi)*V - (C/B)* V**2.0 - (1.0/B)*theta
    return dzdxi

def wdSolver(theta_c, step,outf,theta_bound): 
    x1    = 0      ## ???
    theta1 = theta_c      #likewise
    z1    = 0
    while theta1 > theta_bound:
        x1, theta1, z1 = rKN(x1, theta1, z1, wd_1, wd_2, step)
        print >> outf, x1 ,"\t\t", theta1, "\t", z1
    print "The first zero (xi_1) is: ", x1
    print "dtheta/dxi at xi_1 is: ", z1, "\n"
    return x1, theta1, z1  

def RK_withM(x, theta, z, fx, gx, h):    #no iters
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

    # #Boole's rule for Numerical Integration
    b = x
    a = x - h
    f0 = a
    f1 = a + (b-a)/4.0
    f2 = a + (b-a)/2.0
    f3 = a + 3.0*(b-a)/4.0
    f4 = b
    m = theta*((b-a)/90.0)*(7.0*f0**2.0 + 32.0*f1**2.0 + 12.0*f2**2.0 + 32.0*f3**2.0 + 7.0*f4**2.0)

    #trapezoid rule for numerical integration
    # m = h*theta*( (x-h)**2.0 + x**2.0 )/2.0 
    if m <0:
        m = 0
    return x,theta,z, m

def wd_physical(theta_c, step,theta_bound): 
    x1    = 0      ## ???
    theta1 = theta_c      #likewise
    z1    = 0
    mlist = []
    while theta1 > theta_bound:
        x1, theta1, z1, m = RK_withM(x1, theta1, z1, wd_1, wd_2, step)
        mlist.append(m)
        #print "m: ", m
    m = sum(np.array(mlist))
    print "m: ", m, "\tMsolar: ", Msolar
    #totalM = 0.090*m*Msolar
    totalM = 4.0*pi*rho_0*(alpha_wd**3.0)*m
    totalR = alpha_wd*x1
    return x1, theta1, z1, totalM, totalR