#!/usr/bin/env python
##################################################################
#
# executable main.py
#
# Numerical Integration with Physical Quantities
#
# This is the main program. It requires the gfortran compiler and python version 2.7 or higher.
# This script utilizes the modules constants_bank.py, functions_bank.py, and plotmodule.py 
# which must be located in the same directory as this script. 
# To execute, type in the terminal
#   >$ chmod +x main.py
#   >$ ./main.py
#
# This script contains five subroutines, each marked with an operation flag (all set to "True" by default)
#   do_integrate solves the dimensionless Lane-Emden equation
#   do_physical performs the physical conversion to temperature, density, pressure
#       !! do_integrate must be performed FIRST before this can be run 
#       !! before each run, adjust the physical parameters to the ones you want in constants_bank.py
#   do_qm_gas computes and plots the three quantum gas equations of state   
#   do_whitedwarf generates a set of white dwarf models and plots them in physical units 
#   do_wd_structure saves the (large) data files for two white dwarf models that have 1 Msolar and 1.3 Msolar, respectively     
#
#   A plotting flag controls "plotmodule.py" whih can be used AFTER all of the data files have been generated
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
LaneEmden_data_filename = "theta_xi.dat"
physical_soln_filename  = "physical.dat"
whitedwarf_data_filename= "wdstats.dat"
wd_structure_header     = "wd_structure"

plot_everything = True
##################################################################
# Lane-Emden Solver: Obtain xi_1 and dtheta/dxi at xi_1
if do_integrate ==True:
    theta_bound = 0.0
    step  = 10**(-4.0)
    outf = open(LaneEmden_data_filename,"w")
    print >> outf,   "## xi\t\ttheta\t\tdtdxi "
    x1,t1, z1 = LaneEmdenSolver(step,outf,theta_bound)
    outf.close()

    print "\n\nSolutions: \nx1: ", x1, "\ntheta 1: ", t1 ,"\nz1: ", z1
    LaneEmdenplotter(filename=LaneEmden_data_filename,x1= x1,outname="LaneEmdenSolutions")

# Convert to Physical Space
if do_physical ==True:
    xi_1 = x1
    z_xi_1 = z1 
    theta_1 = 0.0
    alpha= R/xi_1
    rho_c = -M*xi_1/(z_xi_1 * 4*math.pi*R**3.0) 
    K = (4*math.pi)**(1/n) * (G/(n+1))*xi_1**(-(n+1)/n) * (-z_xi_1)**((1-n)/n) *M**((n-1)/n) *R**((3-n)/n) #My K

    xi, theta, dtdxi=np.loadtxt(LaneEmden_data_filename, usecols=(0,1,2), unpack=True)
    
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