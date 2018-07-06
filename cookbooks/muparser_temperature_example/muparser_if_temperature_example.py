#!/Users/billen/miniconda3/bin/python3
# Simple Python script to illustrate (?:) syntax for muparser

# Import modules from Python libraries
import numpy as np
from scipy.special import erfc
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter, MaxNLocator
from matplotlib.colors import BoundaryNorm

plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["font.size"] = "10"
########################################################################################
sec2yr = 60*60*24*365
kappa = 1e-6    # m^2/s
km2m = 1000

# Depth Profile
yminkm = 0
dykm = 10
ymaxkm = 1000
ykm = np.arange(yminkm,ymaxkm+dykm,dykm)
ny = ykm.shape[0]
# in meters
y = ykm*km2m
ymax = ymaxkm*km2m

xminkm = 0
dxkm = 50
xmaxkm = 5000
xkm = np.arange(xminkm,xmaxkm+dxkm,dxkm)
nx = xkm.shape[0]
# in meters
x = xkm*km2m

# Temperature 
Ts = 273
Tm = 1673
vsubcmyr = 2.5;
vsub = vsubcmyr/100/sec2yr
ageop = 30e6*sec2yr
xtrm = 2200*km2m
T = np.zeros((ny,nx))
# Note using (1-erfc() instead of erf() 
# because erfc() is defined in dealii muparser, but not erf()
for i in range(ny):  # Loop over each x position
	for j in range(nx):  # Loop over each y position
		if ((x[j]>0.0) and (x[j]<=xtrm)):	# subducting plate T (xm/vsub --> plate age)	
			T[i,j] = Ts + (Tm-Ts)*(1-erfc( (ymax-y[i])/(2*np.sqrt(kappa*(x[j]/vsub) ))))
		else:
			if (x[j]>xtrm): # overriding plate T (fixed age)
				T[i,j] = Ts + (Tm-Ts)*(1-erfc((ymax-y[i])/(2*np.sqrt(kappa*ageop))))
			else:
				T[i,j] = Tm # for x=0, otherwise would get divide by zero

print('ymax=',ymax,'xtrm=',xtrm,'vsub=',vsub)
print('ageop=',ageop,'Tm=',Tm,'Ts=',Ts,'kappa=',kappa)

# For muparser: 
# - no loops (done automatically by defining variable names)
# - remove T[i,j] =  
# - remove [i], [j]
# - replace xm and ym by x and y
# - replace 'and' by &&
# - remove np.
# - replace each if-else by (? : ) (delete if, replace : -> ?, else -> :)
# - use \ to continue lines
# - add some () to make expressions clearer
#
# subsection Initial temperature model
# 	set Model name = function
# 	subsection Function
#		set Variable names = x,y
#		set Function constants = ymax=1.0e6, xtrm= 2.200e6, vsub=7.927e-10, \
#			ageop=9.46e14, Tm=1673, Ts=273, kappa=1e-6 
#   	set Function expression = (((x>0.0) && (x<=xtrm)) ? \
#			   (Ts + (Tm-Ts)*(1-erfc((ymax-y)/(2*sqrt(kappa*(x/vsub)))))) : \
#			(x>xtrm) ? (Ts + (Tm-Ts)*(1-erfc((ymax-y)/(2*sqrt(kappa*ageop))))) :\
#              (Tm) ) 
#	end
# end 

# make a figure to show the temperature
levels = MaxNLocator(nbins=15).tick_values(T.min(), T.max())
cmap = plt.get_cmap('PiYG')
norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)

fig = plt.figure()
ax1 = fig.add_subplot(3,1,1)
im = ax1.pcolormesh(xkm,ykm,T)

fig.colorbar(im, ax=ax1)
ax1.set_xlabel('Distance (km)')
ax1.set_ylabel('Depth (km)')
ax1.set_title('Temperature')

pdffile = 'temperature_muparser_if.pdf'	
plt.tight_layout()
fig.savefig(pdffile,bbox_inches='tight')		
				