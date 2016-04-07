#!/usr/bin/env python
##################################################################
#
# Module constants_bank.py
#
# Module containing physical constants and trivial expressions
# derived from physical constants
#
# Edit values in this module to control structure and composition: 
# n, X, Y, Z, M, etc.
#
# Author: Meridith Joyce
# Adv Stellar Astrophysics, Dartmouth College
#
##################################################################
import math
import numpy as np

pi=math.pi

Rsolar = 6.995*10.0**10.0							# cm
Msolar = 1.989*10.0**33.0       					# number in cgs 
Lsolar = 3.839*10.0**33.0   						# number in cgs (erg/s)
sigma =  5.67*10.0**-5.0    						# Stefan-Boltzmann constant cgs

Z = 0.02                							# "metals"
Y = 0.27                							# Helium
X = 1.0 - (Y + Z)         							# Hydrogen
mu = 1.0/(2*X + (3.0/4.0)*Y + (1.0/2.0)*Z)          # ROUGHLY
mu_e = 2.0 											# adjusted chemical composition for high mass white dwarfs 

mp = 1.67*10**(-24) 				       			# mass of proton in grams
kb = 1.38*10**(-16)       							# boltzman constant in cgs

gamma = 5.0/3.0 									# unitless and possibly pointless
G = 6.674*10**-8.0        							# cgs

M = 1.3*Msolar         		 						# Mass   
Teff = 4000.0             							# degrees Kelvin
L = 10**1.5*Lsolar      							# Luminosity
R = math.sqrt(L/(4.0*math.pi*sigma*Teff**4.0))   	# Obtined using Stefan-Boltzmann 


n = 1.5 											# index
x1_for_1p5 = 3.6538 								# numerical solution for x1 at n=1; stored for plotting purposes only

KN = 3.166*10**(12.0) 								# dyne/cm^2
KR = 4.936*10**(14.0) 								# dyne/cm^2 

rho_0 = 3.789*10.0**6.0 							# g/cm^3; for white dwarf model
alpha_wd = 1.557*10.0**8.0 							# cm

Chandra_mass = 	1.45								# Chandrasekhar mass in solar masses