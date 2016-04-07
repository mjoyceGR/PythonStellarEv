#!/usr/bin/env python
import os 
import subprocess
import numpy as np 

compile_flag = True

if compile_flag ==True:
	subprocess.call("gfortran RungeKutta.f -o RK.out",shell=True)
	print "\bWarnings Ignored\ncompilation successful"
	subprocess.call("./RK.out",shell=True)
	print "execution complete"
else:
	subprocess.call("./RK.out",shell=True)
	print "execution complete"