#!/usr/bin/env python
"""This module contains the various subroutines
Usage import """
#######################################################################################

#####################  IMPORT STANDARD MODULES   ######################################   

import sys,os
import numpy as np #for numerical analysis
from scipy import integrate  #integration routine
from subprocess import call #for calling unix commands
from datetime import date  #to give a timestamp to output
import pdb	#for the debugger pdb.set_trace()
from matplotlib.pyplot import * # for the figures
from math import pi,exp,log,sqrt
from ConfigParser import SafeConfigParser
####################       IMPORT OWN MODULES     ######################################
sys.path.append('/users/Raj/CODES/REM-PYTHON')
import Rheo_forward
#####################  MODULES   ######################################################   


def convaspect(ifile,ofile,model,freq,N,gsopt=0):
	MODEL = readaspectfile(ifile)
	MODEL['P'][ MODEL['P'] < 0 ] = 0
	
# 	for ii in range(0,3):
# 		print	 Rheo_forward.PTd_VsQ(freq, MODEL['P'][ii], MODEL['T'][ii], MODEL['gs'][ii], MODEL['rho'][ii], model, 1, MODEL['Vs'][ii])
	
	Nreal = len(MODEL['P'])
	if N == -1:
		N = Nreal
		
		
	print "doing",N,"out of",Nreal,"points"
	Q = np.zeros((N))
	Vs = np.zeros((N))
	Vp = np.zeros((N))
	Vs_1hz = np.zeros((N))
	Vp_1hz = np.zeros((N))

	print "Looping through nodes, calculating Vs, Q"
	for ii in range(0,N):
# 		print ii
		if ii%(N/20)==1: 
			print 100*ii/N,"%"
# 	THIS IS NOW DONE IN THE DOWNSAMPLE MATLAB SCRIPT
# 		if MODEL['P'][ii] > 24.3e9:
# 			MODEL['gs'][ii]= MODEL['gs'][ii]/2.5
# 		
# 		if MODEL['gs'][ii] < 0:
# 			MODEL['gs'][ii]=1e-6
			
			
		# if MODEL['P'][ii] < 1e9: # only calculate for mantle pressures
# 			Vsi = MODEL['Vs'][ii]
# 			Qi  = 999.
# 		elif MODEL['T'][ii] < 900:
# 			Vsi = MODEL['Vs'][ii]
# 			Qi  = 999.
# 		else:
		if gsopt==0:
			Vsi, Qi = Rheo_forward.PTd_VsQ(freq,MODEL['P'][ii],MODEL['T'][ii],MODEL['gs'][ii],MODEL['rho'][ii],model,1,MODEL['Vs'][ii])
		else:
			Vsi, Qi = Rheo_forward.PTd_VsQ(freq,MODEL['P'][ii],MODEL['T'][ii],gsopt,MODEL['rho'][ii],model,1,MODEL['Vs'][ii])
			
		
		Q[ii] = Qi
		Vs[ii] = Vsi
		Vp[ii] = Vp0_to_VpF(MODEL['Vp'][ii],MODEL['Vs'][ii],Vs[ii])
		
		Vs_1hz[ii] = invKramers(Vs[ii],Q[ii],2*pi*freq)
		Vp_1hz[ii] = Vp0_to_VpF(MODEL['Vp'][ii],MODEL['Vs'][ii],Vs_1hz[ii])
		
	if gsopt!=0:
		print "THIS RUN WAS DONE WITH A UNIFORM GRAINSIZE OF ",gsopt
		ofile=ofile+"_gs"+str(gsopt)

	f= open(ofile,'w')
	for ii in range(0,N):
		f.write("%9.1f %9.1f %.3e %6.1f %f %e %7.2f %f %f %f %f %f \n" % (MODEL['X'][ii],MODEL['Y'][ii],MODEL['P'][ii],MODEL['T'][ii],MODEL['gs'][ii],MODEL['viscosity'][ii],MODEL['rho'][ii],Vp[ii],Vs[ii],Q[ii],Vp_1hz[ii],Vs_1hz[ii])) 
	
	return MODEL


def readaspectfile(ifile):
	"""Reads the splitting kernel file Modenjcoml.s. Returns kerarray"""
	print "Reading ASPECT output file",ifile

	mm=np.dtype([('P',np.float),('T',np.float),('gs',np.float),('viscosity',np.float),('rho',np.float),('Vp',np.float),('Vs',np.float),('X',np.float),('Y',np.float)])	
	f = open(ifile, 'r')
	count=0

	for line in f:
		ar=np.array(line.split(","))
		if count<1: # Don't read header
			count=count+1
		elif count==1:
			MODEL=np.array([(float(ar[3]),float(ar[4]),float(ar[5]),float(ar[6]),float(ar[7]),float(ar[12]),float(ar[13]),float(ar[16]),float(ar[17]))],dtype=mm)
			count=count+1	
		elif count>1:
				MODEL.resize(len(MODEL)+1,refcheck=0)
				count=count+1
				MODEL[len(MODEL)-1]=np.array([(float(ar[3]),float(ar[4]),float(ar[5]),float(ar[6]),float(ar[7]),float(ar[12]),float(ar[13]),float(ar[16]),float(ar[17]))],dtype=mm)
	
	MODEL['Vp'] = 1000*MODEL['Vp']
	MODEL['Vs'] = 1000*MODEL['Vs']
	return MODEL

def Vp0_to_VpF(Vp0,Vs0,VsF):
	Vp = sqrt( pow(Vp0,2) - (4/3)*( pow(Vs0,2) - pow(VsF,2) ) )
	return Vp
	
def invKramers(Vs_obs,Q_obs,omega_obs):
	Vs_1hz = Vs_obs/(1 + (1/(pi*Q_obs))*log(omega_obs/(2*pi)))
	return Vs_1hz
	


#########################  MAIN   ######################################################   	
def usage():
	"Displays the correct usage of translatePTd2VQ.py"
	print "Usage:  %s ifile ofile model freq N gsopt" % os.path.basename(sys.argv[0])

def main():
	if len(sys.argv) == 7:
		# Give name of aspect output file
		ifile=sys.argv[1]
		# Give name of translated output file
		ofile=sys.argv[2]
		#model="M11"
		model=sys.argv[3]
		#freq=1
		freq=float(sys.argv[4])
		#freq=1
		N=int(sys.argv[5])
		gsopt=float(sys.argv[6])
		convaspect(ifile,ofile,model,freq,N,gsopt)
	else:	
		usage()
		sys.exit(2)
		
	return

if __name__ == "__main__":
    main()