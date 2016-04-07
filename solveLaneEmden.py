#!/usr/bin/env python
import math
import timeit

iterations=1000


# RK Techinique with arrays: ------------------------------------------------------------
#x= x, fx = f(x), n= dimension, step = step size

def rKN(x, fx, n, step):	
    k1 = []
    k2 = []
    k3 = []
    k4 = []
    xk = []
    for i in range(n):  # n is the dimension NOT AN ITERATOR!!!!!
	# print "n: ", n
	# print "range(0,n): ",range(0,n)
	# print "fx[i]: ", fx[i]
	# print "fx[i](x): ", fx[i](x)
        k1.append(fx[i](x)*step)
    for i in range(n):
        xk.append(x[i] + k1[i]*0.5)
    for i in range(n):
        k2.append(fx[i](xk)*step)
    for i in range(n):
        xk[i] = x[i] + k2[i]*0.5
    for i in range(n):
        k3.append(fx[i](xk)*step)
    for i in range(n):
        xk[i] = x[i] + k3[i]
    for i in range(n):
        k4.append(fx[i](xk)*step)
    for i in range(n):
        x[i] = x[i] + (k1[i] + 2*(k2[i] + k3[i]) + k4[i])/6
    print "\n\nnRK ran\n\n"
    return x


def LaneEmden_1(xi, theta, dtdxi,n): ##Lane-Emden separated into two ODEs
    n = 1.5     #for this, always 
    z = dtdxi
    return z#, dzdxi

def LaneEmden_2(xi, theta, dtdxi,n):
    n=1.5
    z = dtdxi
    dzdxi = -(1.0/xi**2.0) * (2*xi*z + xi**2.0*theta**n)
    return dzdxi

def LaneEmdenSolver():
	fx=[LaneEmden] 			#RK takes function to be integrated as an argument
	w1=[0.25]				#likewise
	step = 0.1
	n=1 #1 dimension 
	for i in range(iterations):
		w1 = rKN(w1,fx,n,step) 	 	#0.25 #some value of w
		z1 = rKN(z1,fx,n,step)		#0.5 #some value of z
		print "w1: ", w1, "z1: ", z1

LaneEmdenSolver()