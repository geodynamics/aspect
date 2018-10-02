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
#####################  MODULES   ######################################################   


def PTd_VsQ(freq,P,T,gs,rho,model,highQapprox,Vs0=None):
	"""This subroutine calculates the values of Vs and Q at a given 
	PTd condition at a frequency f and density rho. If highQapprox is 1
	approximate values of Vs and Q can be calculated from 
	e.g. McCarthy et al. (2011) eqn. 24, else use their more explicit form eqn. B6"""
	
	######### NB this uses V_activation = sloping from 11 - 3 (e-6) in the LM ##########
		
	#### Value of the real and complex part of compliance			
	if Vs0==None:
		J1byJu=getJ1byJu(freq,P,T,gs,model)
		J2byJu=getJ2byJu(freq,P,T,gs,model)
		Ju=getJu(T,P,model)
	else:
		J1byJu=getJ1byJu(freq,P,T,gs,model,Vs0,rho)
		J2byJu=getJ2byJu(freq,P,T,gs,model,Vs0,rho)
		Ju=getJu(T,P,model,Vs0,rho)	

	J1=J1byJu*Ju
	J2=J2byJu*Ju
	
	# EVERYONE AGREES
	qs = J2/J1
	# VS DISAGREEMENT	
	if model=='M11':
		Vs = 1/(sqrt(rho*J1))
	elif model == 'P_M13':
		Vs = 1/(sqrt(rho*J1))
	elif model == 'JF10_eBurg':
		G = pow(J1**2 + J2**2,-0.5)
		Vs = sqrt(G/rho)
	elif model == 'Tak14':
# 		pdb.set_trace()
		G = pow(J1**2 + J2**2,-0.5)
		Vs = sqrt(G/rho)	

	#### When Qs^-1 << 1 (J2/J1 << 1), highQapprox can be 1. If Qs^-1 is not small, Qs^-1 does not equal Qmu^-1 ###### 
	if highQapprox == 0:
		factor=(1+sqrt(1+pow((J2/J1),2)))/2
		qs=qs/factor
		Vs=Vs/sqrt(factor)	

	Q=1/qs
	
# 	if Q > 1e9:
# 		pdb.set_trace()	
	
	return Vs,Q

def getJ1byJu(freq,P,T,gs,model,Vs0=None,rho=None):
	""" This subroutine gets the value of J1/Ju based on expressions in a paper """
	#	If model is McCarthy 2011
	if model == 'M11':
		if Vs0==None:
			fn=getfn(freq,P,T,gs,model)
		else:
			fn=getfn(freq,P,T,gs,model,Vs0,rho)
						
		#Getting J1 values from McCarthy et al., 2011
		a0=getconstant(model,'a0') 
		a1=getconstant(model,'a1') 
		a2=getconstant(model,'a2') 
		a3=getconstant(model,'a3') 
		a4=getconstant(model,'a4') 
		a5=getconstant(model,'a5') 
		a6=getconstant(model,'a6') 

		if fn <= 1.0e+13: 
			val=0
			val=val+a0*(pow(log(fn),0))
			val=val+a1*(pow(log(fn),1))
			val=val+a2*(pow(log(fn),2))
			val=val+a3*(pow(log(fn),3))
			val=val+a4*(pow(log(fn),4))
			val=val+a5*(pow(log(fn),5))
			val=val+a6*(pow(log(fn),6))
			J1byJu=1/val
		else:
			J1byJu=1

#	If model is Priestley and McKenzie 2013		
	elif model == 'P_M13':
		fn=getfn(freq,P,T,gs,model,Vs0,rho)
		#Getting J1 values from McCarthy et al., 2011
		a0=getconstant(model,'a0') 
		a1=getconstant(model,'a1') 
		a2=getconstant(model,'a2') 
		a3=getconstant(model,'a3') 
		a4=getconstant(model,'a4') 
		a5=getconstant(model,'a5') 
		a6=getconstant(model,'a6') 

		if fn <= 1.0e+13: 
			val=0
			val=val+a0*(pow(log(fn),0))
			val=val+a1*(pow(log(fn),1))
			val=val+a2*(pow(log(fn),2))
			val=val+a3*(pow(log(fn),3))
			val=val+a4*(pow(log(fn),4))
			val=val+a5*(pow(log(fn),5))
			val=val+a6*(pow(log(fn),6))
			J1byJu=1/val
		else:
			J1byJu=1
			
#	If model is extended Burgers from Jackson & Faul 2010
	elif model == 'JF10_eBurg':
		alpha=getconstant(model,'alpha')
		sigma=getconstant(model,'sigma')
		Delta=getconstant(model,'Delta')	
		tauHR=getconstant(model,'tauHR')	
		tauLR=getconstant(model,'tauLR')	
		tauPR=getconstant(model,'tauPR')
		DeltaP=getconstant(model,'DeltaP')
		ma=getconstant(model,'ma')
		
		omega=2*pi*freq	
		tauH=gettau(tauHR,ma,P,T,gs,model)
		tauL=gettau(tauLR,ma,P,T,gs,model)
		tauP=gettau(tauPR,ma,P,T,gs,model)
		
		j1b = (alpha*Delta)/(pow(tauH,alpha) - pow(tauL,alpha))
		j1p = DeltaP/(sigma*sqrt(2*pi))
		i1b = intgrate_j1b(tauL,tauH,alpha,omega,N=1000)
		i1p = integrate.quad(intgrnd_j1p, 0, np.inf, args=(omega,tauP,sigma))[0]

		J1byJu=1.0+(j1b*i1b)+(j1p*i1p)
				
#	If model is Takei 2014
	elif model == 'Tak14':
		m=getconstant(model,'m') # grainsize exponent
		eta0=getconstant(model,'eta0') # ref viscosity
		Ap=getconstant(model,'Ap') # pre-exponent for peak (normalisation factor)
		sigmap=getconstant(model,'sigmap') # controls peak width - see Takei 2014 eqn 17
		
		#### Value of the real and complex part of compliance			
		if Vs0==None:
			Ju=getJu(T,P,model)
		else:
			Ju=getJu(T,P,model,Vs0,rho)	
		
		omega=2*pi*freq	
		tauM=gettau(Ju*eta0,m,P,T,gs,model)
		
		i1 = integrate.quad(intgrnd_j1tak, 0, np.inf, args=(omega,tauM,Ap,sigmap))[0]

		J1byJu = 1.0 + i1
										
	return J1byJu

def getJ2byJu(freq,P,T,gs,model,Vs0=None,rho=None):
	""" This subroutine gets the value of J2 based on expressions in a paper """	
	
#	If model is McCarthy 2011 
	if model == 'M11':
		fn=getfn(freq,P,T,gs,model,Vs0,rho)
		taun=1/(2*pi*fn)
		#Getting J2 values from McCarthy et al., 2011
		Xn=mccarthyXn(taun)
		J2byJu=(pi/2*Xn)+taun/(2*pi)

#	If model is Priestley and McKenzie 2013		
	elif model == 'P_M13':
		fn=getfn(freq,P,T,gs,model,Vs0,rho)
		taun=1/(2*pi*fn)
		#Getting J2 values from McCarthy et al., 2011
		Xn=mccarthyXn(taun)
		J2byJu=(pi/2*Xn)+taun/(2*pi)

#	If model is extended Burgers from Jackson & Faul 2010
	elif model == 'JF10_eBurg':
		alpha=getconstant(model,'alpha')
		sigma=getconstant(model,'sigma')
		Delta=getconstant(model,'Delta')	
		tauHR=getconstant(model,'tauHR')	
		tauLR=getconstant(model,'tauLR')	
		tauPR=getconstant(model,'tauPR')
		tauMR=getconstant(model,'tauMR')
		DeltaP=getconstant(model,'DeltaP')
		ma=getconstant(model,'ma')
		mv=getconstant(model,'mv')
		
		omega=2*pi*freq	
		tauH=gettau(tauHR,ma,P,T,gs,model)
		tauL=gettau(tauLR,ma,P,T,gs,model)
		tauP=gettau(tauPR,ma,P,T,gs,model)
		tauM=gettau(tauMR,mv,P,T,gs,model)
		
		j2b = omega*(alpha*Delta)/(pow(tauH,alpha) - pow(tauL,alpha))
		j2p = omega*DeltaP/(sigma*sqrt(2*pi))
		
		tau_arr=np.linspace(tauL,tauH,100)
		i2b = intgrate_j2b(tauL,tauH,alpha,omega,N=1000)
		i2p = integrate.quad(intgrnd_j2p, 0, np.inf, args=(omega,tauP,sigma))[0]
		
		J2byJu = (j2b*i2b)+(j2p*i2p)+(1.0/(omega*tauM))
		
#	If model is Takei 2014
	elif model == 'Tak14':
		m=getconstant(model,'m') # grainsize exponent
		eta0=getconstant(model,'eta0') # ref viscosity
		Ap=getconstant(model,'Ap') # pre-exponent for peak (normalisation factor)
		sigmap=getconstant(model,'sigmap') # controls peak width - see Takei 2014 eqn 17

		#### Value of the real and complex part of compliance			
		if Vs0==None:
			Ju=getJu(T,P,model)
		else:
			Ju=1/(rho*pow(Vs0,2))
		
		omega=2*pi*freq	
		tauM=gettau(Ju*eta0,m,P,T,gs,model)
		
		i2 = integrate.quad(intgrnd_j2tak, 0, np.inf, args=(omega,tauM,Ap,sigmap))[0]
		J2byJu = i2 + (1.0/(omega*tauM))
				
	return J2byJu

def mccarthyXn(taun):
	"""This subroutine gets value of relaxation spectrum at a value of 
	normalized time scale (taun) from McCarthy et al. (2011) eqn. 25 """
	if taun >= 1.0e-11:
		Xn=0.32*(pow(taun,(0.39-0.28/(1+2.6*pow(taun,0.1)))))
	else:
		Xn=1853.0*sqrt(taun)
	return Xn
			
def getfn(freq,P,T,gs,model,Vs0 = None,rho=None):
	"""Get the normalized frequency from McCarthy et al. (2011) eqn. 19 and 22"""
	eta0=getconstant(model,'eta0')
	m=getconstant(model,'m')
	
	if Vs0==None:
		Ju=getJu(T,P,model)
	else:
		Ju=getJu(T,P,model,Vs0,rho)	
		
	tau0=Ju*eta0
	taur=gettau(tau0,m,P,T,gs,model)
	if model == "P_M13":
		taur = gettauM_P_M13(T,P,Ju)
	
	fn=freq*taur	
	return fn

def geteta(P,T,gs,model):
	"""Get the viscosity from a combination of power law and Arrhenius eqn. like 
	in McCarthy eqn. 8"""

	TR=getconstant(model,'TR') 
	gsR=getconstant(model,'gsR')
	eta0=getconstant(model,'eta0')
	if model ==  'M11' or model ==  'P_M13':
		m=getconstant(model,'m')
	elif model ==  'eBurgers':
		m=getconstant(model,'ma')
	E=getconstant(model,'E')
	R=getconstant('Common','R')

	eta=eta0*pow((gs/gsR),m)*exp((E/R)*(TR-T)/(TR*T))
	
	return eta
	
def getJu(T,P,model,Vs0 = None,rho = None ):
	"""This gets the unrelaxed modulus at T,P conditions based on constants from study"""

	if Vs0 is None:
		GUR=getconstant(model,'GUR')
		dGdT=getconstant(model,'dGdT') 
		dGdP=getconstant(model,'dGdP')
		TR=getconstant(model,'TR')
		PR=getconstant(model,'PR')
		Ju=1/(GUR + dGdT*(T-TR) + dGdP*(P-PR))
	else:
		Ju=1/(rho*pow(Vs0,2))	
	
	return Ju

def gettau(tau0,m,P,T,gs,model):
	"""calculate the relaxation time at P,T,gs, given model and reference(0) relaxation 
	time and grainsize exponent"""	
	PR=getconstant(model,'PR')
	TR=getconstant(model,'TR') 
	gsR=getconstant(model,'gsR')
	E=getconstant(model,'E')
	V=getconstant(model,'V')
	R=getconstant('Common','R')	
	
	if P > 24.3e9:
		V = 5e-6							# constant
# 		V = -7.2727e-17*P + 1.1770e-05		# 10 to 2
# 		V = -5.4545e-17*P + 1.13250e-05		# 10 to 4
# 		V = -7.2727e-17*P + 1.3767e-05		# 12 to 4
# 		Vact_lowM(P,Vtop,Vbot,Ptop,Pfold) where Pfold of 11e9 gives 1/2 LM die-out
# 		V = Vact_lowM(P,15e-6,4e-6,24.3e9,11e9)

	tau=tau0*pow((gs/gsR),m)*exp( (E/R)*(TR-T)/(TR*T) + (V/R)*(P*TR-PR*T)/(T*TR) )
	
	return tau

def gettauM_P_M13(T,P,Ju):
	"""This gets the Maxwell relaxation time at T,P conditions based on constants from study"""
	model='P_M13'
	TR=getconstant(model,'TR')
	PR=getconstant(model,'PR')
	E=getconstant(model,'E')
	V=getconstant(model,'V')
	R=getconstant('Common','R')	
	eta0=getconstant(model,'eta0')	
	
	if P > 24.3e9:
		V = 5e-6							# constant
# 		V = -7.2727e-17*P + 1.1770e-05		# 10 to 2
# 		V = -5.4545e-17*P + 1.13250e-05		# 10 to 4
# 		V = -7.2727e-17*P + 1.3767e-05		# 12 to 4
# 		Vact_lowM(P,Vtop,Vbot,Ptop,Pfold) where Pfold of 11e9 gives 1/2 LM die-out
# 		V = Vact_lowM(P,15e-6,4e-6,24.3e9,11e9)


	
	astar = exp((E - PR*V)/(R*TR))/exp((E-P*V)/(R*T))
	eta = eta0/astar
	tauM = Ju*eta
	
	return tauM
	
def getconstant(model,par):
	"""This subroutine gets the value of a constant relevant for a modeling 
	approach"""
	parser = SafeConfigParser()
	parser.read('constants_complete.ini')
	val=float(parser.get(model,par))
	
	return val
	
#####  Integrations for JF10  #####
def intgrate_j1b(limL,limH,alpha,omega,N=1000):
	xx = np.logspace(np.log10(limL),np.log10(limH),N)
	yy = np.zeros(len(xx))
	for ii in range(0,len(xx)):
		yy[ii] = pow(xx[ii],alpha-1)/(1 + pow(omega*xx[ii],2))
		
	I=integrate.cumtrapz(yy,xx)
	I = I[len(I)-1]
	return I
	
def intgrate_j2b(limL,limH,alpha,omega,N=1000):
	xx = np.logspace(np.log10(limL),np.log10(limH),N)
	yy = np.zeros(len(xx))
	for ii in range(0,len(xx)):
		yy[ii] = pow(xx[ii],alpha)/(1 + pow(omega*xx[ii],2))
		
	I=integrate.cumtrapz(yy,xx)
	I = I[len(I)-1]
	return I
	
def intgrate_j1p(limL,limH,omega,tauP,sigma,N=1000):
	xx = np.logspace(np.log10(limL),np.log10(limH),N)
	yy = np.zeros(len(xx))
	for ii in range(0,len(xx)):
		yy[ii] = (1/xx[ii])*(1/(1 + pow(omega*xx[ii],2)))*exp(-pow(log(xx[ii]/tauP),2)/(2*pow(sigma,2)))
		
	I=integrate.cumtrapz(yy,xx)
	I = I[len(I)-1]
	return I
	
# def intgrate_j2p(limL,limH,omega,tauP,sigma,N=1000):
# 	xx = np.logspace(np.log10(limL),np.log10(limH),N)
# 	yy = np.zeros(len(xx))
# 	for ii in range(0,len(xx)):
# 		yy[ii] = (1/(1 + pow(omega*xx[ii],2)))*exp(-pow(log(xx[ii]/tauP),2)/(2*pow(sigma,2)))
# 		
# 	I=integrate.cumtrapz(yy,xx)
# 	I = I[len(I)-1]
# 	return I
		
def intgrnd_j1b(tau,alpha,omega):
	return pow(tau,alpha-1)/(1. + pow(omega*tau,2.))

def intgrnd_j1p(tau,omega,tauP,sigma):
	return (1/tau)*(1/(1 + pow(omega*tau,2)))*exp(-pow(log(tau/tauP),2)/(2*pow(sigma,2)))

def intgrnd_j2b(tau,alpha,omega):
	return pow(tau,alpha)/(1. + pow(omega*tau,2.))
	
def intgrnd_j2p(tau,omega,tauP,sigma):
	return (1/(1 + pow(omega*tau,2)))*exp(-pow(log(tau/tauP),2)/(2*pow(sigma,2)))


#####  Integrations for Takei 2014	#####
def intgrnd_j1tak(tau,omega,tauM,Ap,sigmap):
	"""Integration part of J1byJu term (eqn 8 in Takei et al 2014)"""
	return (1/tau)*(1/(1 + pow(omega*tau,2)))*(0.444*pow(tau/tauM,0.38) + Ap*exp(-0.5*pow(((log(tau/tauM) + 8.1)/sigmap),2)))

def intgrnd_j2tak(tau,omega,tauM,Ap,sigmap):
	"""Integration part of J2byJu term (eqn 8 in Takei et al 2014)"""
	return (omega/(1 + pow(omega*tau,2)))*(0.444*pow(tau/tauM,0.38) + Ap*exp(-0.5*pow(((log(tau/tauM) + 8.1)/sigmap),2)))

def Vact_lowM(P,Vtop,Vbot,Ptop,Pfold):
	return Vbot + (Vtop-Vbot)*exp( (Ptop-P)/Pfold )


	
#########################  	PLOTTING   ############################################   	
