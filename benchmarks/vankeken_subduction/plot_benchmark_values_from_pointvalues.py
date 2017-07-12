
# coding: utf-8

# # Summary
# 
# This jupyter notebook compares results for the subduction benchmark in van Keken et al. 2008, case UM_1c with results calculated using ASPECT. The benchmark values for case 1c were provided by Peter Van Keken (pers. comm to Max Rudolph) below. 
# 
# On Mon, Sep 8, 2014 at 7:54 AM, Peter van Keken <keken@umich.edu> wrote:
# > Hi Max:
# >
# > sorry for the delay. The 6x6 files (same format as that requested in the
# > benchmark) for T, P, vx and vy are attached for UM_1c.
# >

# In[1]:
import sys
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

import glob
import tables
import scipy.interpolate as interpolate
from joblib import Parallel, delayed
import multiprocessing

num_cores = multiprocessing.cpu_count()

#get_ipython().magic('matplotlib notebook')
# Load Peter van Keken's benchmark values
vk_vx_file = 'vankeken_2008_values/vx.dat'
vk_vy_file = 'vankeken_2008_values/vy.dat'
vk_T_file = 'vankeken_2008_values/T.dat'

vk_vx = np.flipud(np.loadtxt(vk_vx_file))
vk_vy = np.flipud(np.loadtxt(vk_vy_file))
vk_T = np.flipud(np.loadtxt(vk_T_file))

output_dir = sys.argv[1]
print('Reading output from ',output_dir);

# In[2]:

# Make a list of points where output is desired.
nx = 111
ny = 101
x = np.linspace(0,6.6e5,nx)
y = np.linspace(0,6.0e5,ny)
n = nx*ny
xx,yy = np.meshgrid(x,y)
xxx = np.reshape(xx,(n,1))
yyy = np.reshape(yy,(n,1))

# parse the point_values file
pv_file = output_dir + '/point_values.txt'
pointvalues = np.loadtxt(pv_file)
# identify last time present
last_time = pointvalues[-1,0]
# isolate last timestep values
pointvalues = pointvalues[ pointvalues[:,0] == last_time , : ]
pvx = pointvalues[:,3]
pvy = pointvalues[:,4]
pvT = pointvalues[:,6]

a_T = np.reshape(pvT,(ny,nx))
a_vx = np.reshape(pvx,(ny,nx))
a_vy = np.reshape(pvy,(ny,nx))
    
dT = (a_T-273.0) - vk_T
#T_misfit.append(np.linalg.norm(dT,'fro'))

a_Tu = np.flipud(a_T)
T1111=( a_Tu[10,10] )
    
# Van Keken et al. 2008 define Twedge and Tslab as benchmark values
tmp = 0.0
for i in range(0,36):
    tmp += a_Tu[i,i]**2
    Tslab=np.sqrt(tmp/36.)
    
tmp=0.0
for i in range(9,21):
    for j in range(9,i+1):
        tmp+=a_Tu[j,i]**2
Twedge=np.sqrt(tmp/78.)

dT = (a_T-273.0) - vk_T
print(np.linalg.norm(dT,'fro'),'fro')

print('Making the graphical plots')
f, (ax1, ax2, ax3) = plt.subplots(1,3,figsize=(12,4))

cs1=ax1.pcolormesh(xx/1e3,yy/1e3,a_vx*100.-vk_vx,vmin=-5,vmax=5,cmap='bwr')
ax1.set_xlabel('(km)')
ax1.set_ylabel('(km)') 
ax1.set_title('x-velocity difference')
f.colorbar(cs1,ax=ax1)

cs2=ax2.pcolormesh(xx/1e3,yy/1e3,a_vy*100.-vk_vy,vmin=-5,vmax=5,cmap='bwr')
ax2.set_xlabel('(km)')
ax2.set_title('y-velocity difference')
f.colorbar(cs2,ax=ax2)

cs3=ax3.pcolormesh(xx/1e3,yy/1e3,a_T-273. - vk_T,vmin=-8.,vmax=8.,cmap='bwr')
ax3.set_xlabel('(km)')
ax3.set_title('Temperature difference')
f.tight_layout()
f.colorbar(cs3,ax=ax3)
f.savefig(output_dir + '/aspect_vankeken_difference_pv.eps')
#plt.show()


# In[12]:



print('Final values: \n','T_wedge=','{:f}'.format(Twedge-273.))
print('T_slab=','{:f}'.format(Tslab-273.))
print('T(11,11)=','{:f}'.format(T1111-273.))
