#!/usr/bin/env python
##################################################################
#
# executable main_withfortan.py
#
# Same as main.py but using fortran for numerical integration. This does not give
# as good of values as the python routines and it is completely tertiary 
#
# Author: Meridith Joyce
# Adv Stellar Astrophysics, Dartmouth College
#
##################################################################
# python intrinsic modules
import math
import numpy as np
import matplotlib.pyplot as plt

##################################################################
# modules I wrote
from constants_bank import *
from functions_bank import *
from plotmodule import *

### Operation Flags
do_integrate  = True       # Perform numerical integration of Lane-Emden
do_physical   = True       # Perform physical conversion to temperature, density, pressure
do_qm_gas     = True       # Plot the three quantum gas equations of state
do_whitedwarf = True       # Generate a group of white dwarf models 
do_wd_structure = True      # Save and plot the internal structure of white dwarf models with 1, 1.3 Msolar

### names of data files this program generates; edit at will
#LaneEmden_data_filename = "theta_xi.dat"
physical_soln_filename  = "physical.dat"
whitedwarf_data_filename= "wdstats.dat"
wd_structure_header     = "wd_structure"

plot_everything = True
##################################################################
# Lane-Emden Solver: Obtain xi_1 and dtheta/dxi at xi_1 using FORTRAN instead of python
if do_integrate ==True:
    x1 = 0.0
    theta_1 = 1.0
    z1 = 0
    stepsize = 1.0e-5
    LaneEmden_fortran(x1, theta_1, z1,stepsize,n)
    LaneEmdenplotter(filename="results.txt",x1= x1,outname="LaneEmdenSolutions") 

# Convert to Physical Space
if do_physical ==True:
    xis, thetas, zis =np.loadtxt("results.txt",usecols=(0,1,2), unpack=True)
    xi_1 = xis[-1]
    z_xi_1 = zis[-1] 
    theta_1 = 0.0

    alpha= R/xi_1
    rho_c = -M*xi_1/(z_xi_1 * 4*math.pi*R**3.0) 
    K = (4*math.pi)**(1/n) * (G/(n+1))*xi_1**(-(n+1)/n) * (-z_xi_1)**((1-n)/n) *M**((n-1)/n) *R**((3-n)/n) #My K

    xi, theta, dtdxi=np.loadtxt("results.txt", usecols=(0,1,2), unpack=True)
    
    physout = open(physical_soln_filename,"w") 
    print >>physout, "## mass\t\t\tPressure\t\tTemp\t\tdensity\t\tradius "

    i=0
    rho = rho_c
    while i <= len(xi)-1:
        rho = Density(theta[i],rho_c)
        P = Pressure(K, rho)
        T = Temperature(P, rho)
        print >> physout, P, "\t", T, "\t", rho, "\t", xi[i]*alpha
        i = i+1 
    physout.close()
    rmax = 0.8e12
    plotter(filename="physical.dat",rmax=rmax,outname="physical")


# quantum gas equations
if do_qm_gas == True:
    rho_exp = np.array(np.arange(3.0,12.0,0.001)) 
    rho = 10.0**(rho_exp) 
    P_nonrel = PN(rho)  
    P_rel    = PR(rho)
    P_comb   = Pcomb(rho, P_nonrel, P_rel)
    plot_xmax = 12.5
    qm_plot(rho, P_nonrel, P_rel, P_comb,plot_xmax,outname="QMgas_EOS")

# computes white dwarf models over the range 10^theta_c_exp, specified below
if do_whitedwarf==True:
    theta_c_exp = np.array(np.arange(3.0,12.0,0.3)) 
    theta_c = 10.0**(theta_c_exp) 
    rho_c = theta_c/rho_0 
    step = 10.0**(-5.0)     #step size h  

    wd_data = open(whitedwarf_data_filename,"w")
    print >> wd_data, "#stepsize \t rho_c \t\t theta1 \t\t x1 \t\t total M/Msolar\t100R/Rsolar\t" 
    for i in range(len(rho_c)):
        x1, theta1, z1, totalM, totalR =wd_physical(rho_c[i], step,"temp.dat") 
        print "xi_1: ", x1, "\ttotal M: ", totalM/Msolar, "\ttotal R: ", 100*totalR/Rsolar, "\n" 
        print "rho_c =",np.log10(rho_c[i]), "\ttheta_c = 10^",np.log10(theta_c[i])
        print >> wd_data, step, "\t\t", rho_c[i], "\t\t",theta1 ,"\t\t", x1,"\t\t", totalM/Msolar, "\t\t", 100*totalR/Rsolar 
    wd_data.close()

# Stores the structure for white dwarf models with 1.0 and 1.3 Msolar
if do_wd_structure==True:
    names = ["1Msolar","1.3Msolar"]
    theta_c_exp = np.array([7.4,8.4]) #theta exponents that give desired Msolar values (~1%)
    theta_c = 10.0**(theta_c_exp) 
    rho_c = theta_c/rho_0 
    step = 10.0**(-5.0)

    for i in range(len(rho_c)):
        x1, theta1, z1, totalM, totalR =wd_physical(rho_c[i], step,wd_structure_header+"_"+str(names[i]+".dat")) 
        print "xi_1: ", x1, "\ttotal M: ", totalM/Msolar, "\ttotal R: ", 100*totalR/Rsolar, "\n" 

# calls plotmodule.py
if plot_everything == True:
    plot_all()