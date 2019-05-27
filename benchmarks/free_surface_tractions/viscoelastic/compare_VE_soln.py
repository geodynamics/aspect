# This script outputs analytical solutions for the surface
# deformation due to a cylindrical load applied over a
# viscoelastic half-space (from Nakiboglu and Lambeck,
# 1982). It outputs the solution for instantaneous loading
# (as 'soln_XX.txt'), linear unloading (as 'soln2_XX.txt')
# and instantaneous loading in the viscous limit (as 
# '../viscous/soln_XX.txt'), which approximates the
# solution of Haskell (1935a).

import math
import numpy as np
import matplotlib.pyplot as plt
import sys
from scipy.special import j0, j1

# Material properties
mu = 1e10 # shear modulus
A = 100e3 # load radius [m]
g = 9.8 # grav. acceleration [m/s^2]
eta = 3e20 # viscosity [Pa s]
rhos = 3300 # solid ("mantle") density [kg/m^3]
rhol = 900 # load ("ice") density [kg/m^3]
nu = 0.5 # Poisson's ratio

s2yr = 365.25*3600*24

# Heaviside function - chose H(0) = 0 (arg = -1)
# such that numerical results don't lag solution).
def heaviside(x,arg):
    if (arg==-1):
    	if (x>0):
	    return 1
    	else:
            return 0
    elif (arg==1):
        if (x>=0):
            return 1
        else:
            return 0
    elif (arg==0):
        if (x>0):
            return 1
	elif (x==0):
	    return 0.5
        else:
            return 0

# Alpha parameter
def alpha_fun(mu,A):
    alpha = rhos*g*A*(1-nu)/mu
    return alpha

# Beta parameter
def beta_fun(u,mu,eta,A):
   beta = alpha_fun(mu,A)/((eta/mu)*(u + alpha_fun(mu,A)))
   return beta

# Instantaneous loading (H1)
def W1h_fun(u,tt,mu,A,h0,eta,t0,t1): 
   W1 = h0*((1 - (u/(u+alpha_fun(mu,A)))*math.exp(-beta_fun(u,mu,eta,A)*(tt-t0)))*heaviside(tt-t0,-1) - (1 - (u/(u+alpha_fun(mu,A)))*math.exp(-beta_fun(u,mu,eta,A)*(tt-t1)))*heaviside(tt-t1,-1))
   return W1

# Linear loading (H2)
def W2h_fun(u,tt,mu,A,dhhdt,eta,t0,t1):
   tau = eta/mu
   W2 = dhhdt*tau*((tt/tau) - (u/alpha_fun(mu,A)) + (u/alpha_fun(mu,A))*(1 - t0*beta_fun(u,mu,eta,A))*math.exp(-beta_fun(u,mu,eta,A)*(tt - t0))*heaviside(tt-t0,-1)- ((tt/tau) - (u/alpha_fun(mu,A)) + (u/alpha_fun(mu,A))*(1 - t1*beta_fun(u,mu,eta,A))*math.exp(-beta_fun(u,mu,eta,A)*(tt - t1)))*heaviside(tt-t1,-1));
   return W2

# Total loading (H3)
def W3h_fun(u,tt,mu,A,h0,dhhdt,eta,t01,t11,t02,t12):
   W1t = W1h_fun(u,tt-t01,mu,A,h0,eta,t01,t11)*heaviside(tt-t01,1)
   W2t = W2h_fun(u,tt-t02,mu,A,dhhdt,eta,t02,t12)*heaviside(tt-t02,-1)
   W3 = W1t + W2t
   return W3

# Integrate H3.
def W3hu_fun(tt,r,mu,A,h,dhdt,eta,t01,t11,t02,t12):
   du = .1
   Ul = 1000/du
   W3 = np.zeros(int(Ul))
   for ui in range(0,int(Ul)):
       W3[ui] = (rhol*math.pow(rhos,-1))*j1(ui*du)*j0(ui*du*r/A)*W3h_fun(ui*du,tt,mu,A,h,dhdt,eta,t01,t11,t02,t12)	
   W3s = sum(W3)*du
   return W3s


# Time at which instantaneous loading/unloading ends.
t11 = s2yr*1e3
# Time at which instantaneous loading/unloading starts.
t01 = 0

# Time at which linear loading/unloading ends.
t12 = s2yr*1e3
# Time at which linear loading/unloading ends.
t02 = 0

# Load height (for H1).
h0 = 1e3
# Rate of loading (-ve is unloading).
dhdt = -h0/(t12-t02)

# radial profile values.
rs = np.arange(0,500e3,20e3)

W3i = np.zeros(len(rs))
W3i_v = np.zeros(len(rs))
W3i2 = np.zeros(len(rs))

# time values.
ts = np.arange(0,1500,100)

f0 = open('soln_z0.txt','w')
# can approximate Haskell solution by setting high shear modulus (mu).
mu_v = 1e12
f0_v = open('../viscous/soln_z0.txt','w')
f02 = open('soln2_z0.txt','w')


for ti in range(0,len(ts)):
    f = open(str('soln_%d.txt' % ti),'w')
    f_v = open(str('../viscous/soln_%d.txt' % ti),'w')
    f2 = open(str('soln2_%d.txt' % ti),'w')
    for ri in range(0,len(rs)):
	W3i[ri] = W3hu_fun(ts[ti]*s2yr,rs[ri],mu,A,-h0,0,eta,t01,t11,t02,t12)
	W3i_v[ri] = W3hu_fun(ts[ti]*s2yr,rs[ri],mu_v,A,-h0,0,eta,t01,t11,t02,t12)
	W3i2[ri] = W3hu_fun(ts[ti]*s2yr,rs[ri],mu,A,-h0,-dhdt,eta,t01,t11,t02,t12)

        f.write('%f\t%f\n' % (rs[ri],W3i[ri]))
	f_v.write('%f\t%f\n' % (rs[ri],W3i_v[ri]))
	f2.write('%f\t%f\n' % (rs[ri],W3i2[ri]))

    f0.write('%f\t%f\n' % (ts[ti],W3i[0]))
    f0_v.write('%f\t%f\n' % (ts[ti],W3i_v[0]))
    f02.write('%f\t%f\n' % (ts[ti],W3i2[0]))

    f.close()
    f_v.close()
    f2.close()

f0.close()
f0_v.close()
f02.close()
