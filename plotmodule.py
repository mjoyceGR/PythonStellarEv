#!/usr/bin/env python
##################################################################
#
# Module plotmodule.py
#
# Module for plotting results. These functions are called in main.py
# within the relevant subroutines. File names of the output figures 
# are hard-coded but may be manually edited. The names of data files from 
# which content is pulled (python function np.loadtxt()) are specified in main.py. 
# All of the data files and plots are generated in the directory from which
# main.py is executed.
#
# You can execute this as a stand-alone script, but you must change the input
# file names at the very bottom to the names of your files (if different from mine)
# To execute this as a stand-alone script, set flag "Runhere" to "True",
# make this executable with >$ chmod +x plotmodule.py
# and execute with >$ ./plotmodule.py
#
# Author: Meridith Joyce
# Adv Stellar Astrophysics, Dartmouth College
#
##################################################################
    
import math
import numpy as np
import matplotlib.pyplot as plt

from constants_bank import *
from functions_bank import *

def LaneEmdenplotter(filename, x1,outname):
	xi, theta, dtdxi=np.loadtxt(filename,usecols=(0,1,2),unpack=True) #theta_xi.dat
	plt.plot(xi, theta,  "g-",label=r'$\theta(\xi)$')
	plt.plot(xi, dtdxi,  "m-",label=r"$\frac{d\theta}{d\xi}$")
	plt.axvline(x = x1, color='red',linestyle='--',label=r'$\xi_1$')
	plt.xlabel(r'$\xi$')
	plt.ylabel('Lane-Emden Solutions')
	plt.legend(loc=6)
	plt.savefig(str(filename)+"_neq"+str(n)+".png")
	plt.close()
	return

def plotter(filename, rmax, outname):
	pres, temp,dens, radius=np.loadtxt(filename, usecols=(0,1,2,3), unpack=True)

	plt.plot(radius, pres,  "m-",label=r'P(r)')
	plt.xlim(0,rmax)
	plt.ylim(ymin=0)
	plt.xlabel(r'$r = \xi \alpha$')
	plt.ylabel(r'Pressure (Barye = g / cm s$^2$) ')
	plt.legend(loc=1)
	plt.savefig("pressure_neq"+str(n)+"_"+str(outname)+".png")
	plt.close()

	plt.plot(radius, temp,  "r-",label=r'T(r)')
	plt.xlim(0,rmax)
	plt.ylim(ymin=0)
	plt.xlabel(r'$r = \xi \alpha$')
	plt.ylabel('Temperature (K)')
	plt.legend(loc=1)
	plt.savefig("temperature_neq"+str(n)+"_"+str(outname)+".png")
	plt.close()

	plt.plot(radius, dens,  "b-",label=r'$\rho$(r)')
	plt.xlim(0,rmax)
	plt.ylim(ymin=0)
	plt.xlabel(r'$r = \xi \alpha$')
	plt.ylabel(r'Density (g/cm$^3$) ')
	plt.legend(loc=1)
	plt.savefig("density_neq"+str(n)+"_"+str(outname)+".png")
	plt.close()
	return 

def composition_plotter(file1, file2,rmax,outname):
	pres1, temp1, dens1, radius1=np.loadtxt(file1, usecols=(0,1,2,3), unpack=True)
	pres2, temp2, dens2, radius2=np.loadtxt(file2, usecols=(0,1,2,3), unpack=True)

	plt.plot(radius1, pres1,  "m--",linewidth=3,label=r'P(r) at solar composition')
	plt.plot(radius2, pres2,  "k-",linewidth=2, color="orange",label=r'P(r) at BBN composition')
	plt.xlim(0,rmax)
	plt.ylim(ymin=0)
	plt.xlabel(r'$r = \xi \alpha$')
	plt.ylabel(r'Pressure (Barye = g / cm s$^2$) ')
	plt.legend(loc=1)
	plt.savefig("pressure_"+str(outname)+".png")
	plt.close()

	plt.plot(radius1, temp1,  "r--",linewidth=3,label=r'T(r) at solar composition')
	plt.plot(radius2, temp2,  "m-",linewidth=2,label=r'T(r) at BBN composition')
	plt.xlim(0,rmax)
	plt.ylim(ymin=0)
	plt.xlabel(r'$r = \xi \alpha$')
	plt.ylabel('Temperature (K)')
	plt.legend(loc=1)
	plt.savefig("temperature_"+str(outname)+".png")
	plt.close()

	plt.plot(radius1, dens1,  "b--",linewidth=3,label=r'$\rho$(r) at solar composition')
	plt.plot(radius2, dens2,  "c-",linewidth=2,label=r'$\rho$(r) at BBN composition')
	plt.xlim(0,rmax)
	plt.ylim(ymin=0)
	plt.xlabel(r'$r = \xi \alpha$')
	plt.ylabel(r'Density (g/cm$^3$) ')
	plt.legend(loc=1)
	plt.savefig("density_"+str(outname)+".png")
	plt.close()
	return 


def qm_plot(rho, P1, P2, P3,xmax,outname):
	plt.plot(np.log10(rho), np.log10(P2), "m-",color="green",linewidth=3,label=r'P($\rho$) Relativistic')
	plt.plot(np.log10(rho), np.log10(P1), "m-",color="red",linewidth=3,label=r'P($\rho$) Non-relativistic')
	plt.plot(np.log10(rho), np.log10(P3), "m--",color="black",linewidth=3,label=r'P($\rho$) Combined approximation')
	plt.xlabel(r'Log$\rho$ (g/cm$^3$)')
	plt.ylabel(r'Log$P$ (dyne/cm$^2$)')
	plt.legend(loc=4)
	plt.savefig(str(outname)+".png")
	plt.close()

def wd_plot(filename, xmax,m,r,mr):
	density, M, R=np.loadtxt(filename, usecols=(1,4,5), unpack=True)

	plt.plot(np.log10(density*rho_0), M, "mo",color="orange",label=r'Mass vs Log($\rho$)')
	plt.axhline(y=Chandra_mass,color='red',linestyle='--', label=r'Chandrasekhar Limit')
	plt.xlabel(r'Log$\rho$ (g/cm$^3$)')
	plt.ylabel(r'Mass M$_{\odot}$')
	plt.legend(loc=4)
	plt.savefig(str(m)+".png")
	plt.close()

	plt.plot(np.log10(density*rho_0), R,  "mo" ,color="green",label=r'100 R vs Log($\rho$)')
	plt.xlabel(r'Log$\rho$ (g/cm$^3$)')
	plt.ylabel(r'100x Radius / R$_{\odot}$')
	plt.legend(loc=1)
	plt.savefig(str(r)+".png")
	plt.close()

	plt.plot( M,R,  "mo" ,color="purple",label=r'100 R vs Mass $M_{\odot}$')
	plt.axvline(x=Chandra_mass,color='red',linestyle='--', label=r'Chandrasekhar Limit')
	plt.xlabel(r'Mass M$_{\odot}$')
	plt.ylabel(r'100x Radius / R$_{\odot}$')
	plt.legend(loc=1)
	plt.savefig(str(mr)+".png")
	plt.close()
	return


def wd_rhovsR(filename,outfile,colorstr):
	m, xi, theta=np.loadtxt(filename, usecols=(0,1,2), unpack=True)
	density = theta*rho_0
	r = alpha_wd*xi
	plt.plot(r, np.log10(density),  "m-" ,color=str(colorstr),label=r'Log($\rho$) vs $r$')
		#+str(outfile.split('Msolar')[0]) +"$M_{\odot}$")
	plt.xlabel(r'Radial distance from center $r=\alpha\xi$ (cm) ')
	plt.ylabel(r'Log $\rho = \rho_0 \theta$ (g/cm$^3$)')
	plt.legend(loc=1)
	plt.savefig(str(outfile)+".png")
	plt.close()
	return


# Execute this as a stand-alone script AFTER all the information has been generated
# Filenames coded in this section are the ones I chose for my data/the defaults; update as necessary
def plot_all():
	LaneEmdenplotter(filename="theta_xi.dat",x1= x1_for_1p5,outname="LaneEmdenSolutions")

    ## Plot the various compositions
    #       !! sample data files have been included in the .tar for the plotting portion, 
    #       !! but to generate and plot your own data, you can change the name of
    #       !! <physout> and any other filenames right in this program
    #       !! physical_solar.dat, physical_BBN.dat, and physical_M1.3.dat are the names
    #       !! I chose for the output of the compositions assigned, e.g.
    # 		plotter(filename="physical.dat",rmax=rmax,outname="physical")
    # 		composition_plotter(file1="physical_solar.dat",file2="physical_BBN.dat",rmax=rmax,outname="BBN")
    # 		composition_plotter(file1="physical_solar.dat", file2="physical_M1.3.dat", rmax=rmax, outname="newMass")

	rmax = 0.8e12
	plotter(filename="physical.dat",rmax=rmax,outname="physical")
	#composition_plotter(file1="physical_solar.dat",file2="physical_BBN.dat",rmax=rmax,outname="BBN")
	#composition_plotter(file1="physical_solar.dat", file2="physical_M1.3.dat", rmax=rmax, outname="newMass")

	rho_exp = np.array(np.arange(3.0,11.0,0.001)) 
	rho = 10.0**(rho_exp) 
	P_nonrel = PN(rho)  # possibly use theta_c here
	P_rel    = PR(rho)
	P_comb   = Pcomb(rho, P_nonrel, P_rel)
	plot_xmax = 12.5
	qm_plot(rho, P_nonrel, P_rel, P_comb,plot_xmax,outname="QMgas_EOS")

	log_rho_max=4
	mass_vs_rho_outname= "wd_mass"
	radius_vs_rho_outname="wd_radius"
	radius_vs_mass_outname="wd_RvM"
	wd_plot("wdstats.dat",log_rho_max, mass_vs_rho_outname,  radius_vs_rho_outname,  radius_vs_mass_outname)

	wd_rhovsR("wd_structure_1Msolar.dat","1Msolar","magenta")
	wd_rhovsR("wd_structure_1.3Msolar.dat","1.3Msolar","darkgreen")
	return 