
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
plt.close('all')
num_cores = multiprocessing.cpu_count()

#get_ipython().magic('matplotlib notebook')
# Load Peter van Keken's benchmark values
vk_vx_file = 'vankeken_2008_values/vx.dat'
vk_vy_file = 'vankeken_2008_values/vy.dat'
vk_T_file = 'vankeken_2008_values/T.dat'

vk_vx = np.flipud(np.loadtxt(vk_vx_file))
vk_vy = np.flipud(np.loadtxt(vk_vy_file))
vk_T = np.flipud(np.loadtxt(vk_T_file))

#output_dir = sys.argv[1]
output_dir = 'output_case1'
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


# In[3]:

f, (ax1,ax2,ax3) = plt.subplots(1,3,figsize=(9,3))
cs=ax1.pcolormesh(xx/1e3,yy/1e3,vk_vx)
ax1.axis('tight')
plt.colorbar(cs,ax=ax1)
#plt.show()


# In[4]:

# Based on Ian Rose's code in ljhwang/ASPECT-jupyter

viz_files = sorted(glob.glob(output_dir + '/solution/solution*.h5'))
mesh_files = sorted(glob.glob(output_dir + '/solution/mesh*.h5'))

# Load the statistics file
stats_file = output_dir + '/statistics'
stats = np.loadtxt(stats_file,usecols=(0,1,2,4,12,15))
time = stats[:,1]
tavg = stats[:,4]
ncell = stats[:,3]
vrms = stats[:,5]
# load the output filenames from the statistics file
output_files=np.loadtxt(stats_file,usecols=(10),dtype=str)
output_times=[]
for i in range(len(output_files)):
    if output_files[i] != '""':
        output_times.append(time[i])
output_times = np.array(output_times)


T_misfit = []
mesh_file = mesh_files[0]
T1111 = []
Tslab = []
Twedge = []

#for viz_file in viz_files:
def process_file( viz_file ):
    viznum = int(viz_file.split('-')[1].split('.')[0])
    # find the mesh file with this or lower number
    for mf in mesh_files:
        meshnum = int(mf.split('-')[1].split('.')[0])
        if meshnum > viznum:
            break
        else:
            mesh_file = mf
    
    print('Loading ',viz_file)
    print('Loading ',mesh_file)
    
    mesh = tables.open_file(mesh_file,mode='r')
    fields = tables.open_file(viz_file,mode='r')
    nodes = mesh.root.nodes
    ele = mesh.root.cells
    
    a_T = getattr(fields.root,'T')
    a_v = getattr(fields.root,'velocity')
    a_P = getattr(fields.root,'p')
    
    coords = np.array([x for x in nodes])
    values = np.array([t for t in a_T])[:,0]
    p_values = np.array([p for p in a_P])[:,0] # extract the pressure field
    v = np.array([t for t in a_v])
    vx = v[:,0]
    vy = v[:,1]

    method = interpolate.LinearNDInterpolator
    #method = interpolate.NearestNDInterpolator
    
    fn = method(coords,values)
    a_T = fn(xx,yy)
    
    dT = (a_T-273.0) - vk_T
    T_misfit.append(np.linalg.norm(dT,'fro'))
    
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
    if viz_file == viz_files[-1]:
        fn = method(coords,vx)
        a_vx = fn(xx,yy)
        fn = method(coords,vy)
        a_vy = fn(xx,yy)
        fn = method(coords,p_values)
        a_P = fn(xx,yy)
    else:
        fields.close()
        mesh.close()
        a_vx=[]
        a_vy=[]
    return [viznum,a_T,a_vx,a_vy,a_P,T1111,Tslab,Twedge]

results = Parallel(n_jobs=num_cores)(delayed(process_file)(i) for i in viz_files)

Twedge = np.array([x[7] for x in results])
Tslab = np.array([x[6] for x in results])
T1111 = np.array([x[5] for x in results])

a_vx = results[-1][2]
a_vy = results[-1][3]
a_T = results[-1][1]
a_P = results[-1][4]

# In[5]:

dT = (a_T-273.0) - vk_T
print(np.linalg.norm(dT,'fro'),'fro')

#print('Plotting the mesh. This might take a while...')
## plot the final mesh
#f, ax = plt.subplots(1,1,figsize=(4,4))
#for cell in ele:
#    my_coords = coords[cell[[0,1,2,3,0]],:] # nodal coordinates
#    ax.plot(coords[:,0],coords[:,1],'k')
#f.tight_layout()
#f.savefig(output_dir + '/mesh.eps')


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

vlim = np.abs(a_T-273. - vk_T).max()
cs3=ax3.pcolormesh(xx/1e3,yy/1e3,a_T-273. - vk_T,vmin=-vlim,vmax=vlim,cmap='bwr')
ax3.set_xlabel('(km)')
ax3.set_title('Temperature difference')
f.tight_layout()
f.colorbar(cs3,ax=ax3)

f.savefig('figures/' + output_dir + '_aspect_vankeken_difference.png')
#f.savefig('figures/' + output_dir + '_aspect_vankeken_difference.eps')
# In[6]:

print('Making the graphical plots')
f, (ax1, ax2, ax3, ax4) = plt.subplots(1,4,figsize=(16,4))

cs1=ax1.pcolormesh(xx/1e3,yy/1e3,a_vx*100.,vmin=-5,vmax=5,cmap='bwr')
ax1.set_xlabel('(km)')
ax1.set_ylabel('(km)') 
ax1.set_title('x-velocity')
f.colorbar(cs1,ax=ax1)

cs2=ax2.pcolormesh(xx/1e3,yy/1e3,a_vy*100.,vmin=-5,vmax=5,cmap='bwr')
ax2.set_xlabel('(km)')
ax2.set_title('y-velocity')
f.colorbar(cs2,ax=ax2)

vlim = np.abs(a_T-273.).max()
cs3=ax3.pcolormesh(xx/1e3,yy/1e3,a_T-273.,vmin=0.,vmax=vlim,cmap='bwr')
ax3.set_xlabel('(km)')
ax3.set_title('Temperature difference')
f.tight_layout()
f.colorbar(cs3,ax=ax3)

vlim = np.abs(a_P).max()
vlim = 1e8
cs4=ax4.pcolormesh(xx/1e3,yy/1e3,a_P,vmin=-vlim,vmax=vlim,cmap='bwr')
ax4.set_xlabel('(km)')
ax4.set_title('Pressure and Streamlines')
ax4.streamplot(xx/1e3,yy/1e3,a_vx,a_vy)

f.colorbar(cs4,ax=ax4)

f.savefig('figures/' + output_dir + '_aspect_solution.png')
#f.savefig('figures/' + output_dir + '_aspect_solution.eps')
#plt.show()
plt.close()


# In[13]:
print('Making the benchmark value plots')
#get_ipython().magic('matplotlib inline')
f, (ax1,ax2, ax3, ax4) = plt.subplots(4,1,figsize=(4,12))
ax1.plot(time,ncell)
ax2.plot(time,vrms)
ax3.plot(time,tavg)
ax4.plot(T_misfit)
f.tight_layout()

f.savefig('figures' + '/' + output_dir + '_time-dependence.png')
f.savefig('figures' + '/' + output_dir + '_time-dependence.eps')

#plt.show()

f, (ax1,ax2,ax3) = plt.subplots(1,3,figsize=(12,4))
ax1.plot(output_times,Twedge-273.,'k')
ax1.plot([0, output_times[-1]],[850.50, 850.50],'b--')
ax1.plot([0, output_times[-1]],[854.99, 854.99],'r--')
ax1.set_ylabel('$T_{wedge}$ ($^\circ$C)')
ax1.set_xlabel('Time (yr)')

ax2.plot(output_times,Tslab-273.,'k')
ax2.plot([0, output_times[-1]],[503.04, 503.04],'b--')
ax2.plot([0, output_times[-1]],[511.09, 511.09],'r--')
ax2.set_ylabel('$T_{slab}$ ($^\circ$C)')
ax2.set_xlabel('Time (yr)')

ax3.plot(output_times,T1111-273.,'k')
ax3.plot([0, output_times[-1]],[387.78, 387.78],'b--')
ax3.plot([0, output_times[-1]],[397.55, 397.55],'r--')
ax3.set_ylabel('$T_{1111}$ ($^\circ$C)')
ax3.set_xlabel('Time (yr)')

ax4 = ax1.twinx()
ax4.plot(time,ncell,'g')
ax4.set_ylabel('# of cells',color='g')
ax4.tick_params('y',colors='g')
ax4.set_yscale('log')
f.tight_layout()

f.savefig('figures' + '/' + output_dir + '_benchmark_quantities.png')
f.savefig('figures' + '/' + output_dir + '_benchmark_quantities.eps')

# In[ ]:

print('Final values: \n','T_wedge=','{:f}'.format(Twedge[-1]-273.))
print('T_slab=','{:f}'.format(Tslab[-1]-273.))
print('T(11,11)=','{:f}'.format(T1111[-1]-273.))
