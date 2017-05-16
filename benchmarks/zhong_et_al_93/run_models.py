import numpy as np
import matplotlib.pyplot as plt
import os

wavelength = 1.
k = 2.*np.pi/wavelength
depth = 62./64.

# Analytic solution, given by equation 14 of Zhong et al, 1993
def analytic_solution(x):
  return np.cos(k*x)/np.sinh(k)/np.sinh(k) * ( k*(1.-depth) * np.sinh(k) * np.cosh(k*depth) - k * np.sinh(k*(1.-depth)) + np.sinh(k) * np.sinh(k*depth))

# Generate a new prm file for each refinement level.
# Also attempts to tighten up the delta function approximation
# as the mesh refinement increases.
def generate_prm( refine ):
  prmfile = open("zhong_case1.prm", "r")
  outfile = open("tmp.prm", "w")
  for l in prmfile.readlines():
    if 'Initial global refinement' in l:
      outfile.write('    set Initial global refinement = %i\n'%(refine) )
    elif 'Function constants' in l:
      outfile.write('    set Function constants  = lambda=1, pi=3.1415926536, L=1, depth=62, sigma = %f, N=64\n' % (2.*0.5**(refine-5)))
    else:
      outfile.write(l)
  prmfile.close()
  outfile.close()

# Used to inspect a topography solution.
def plot_topography():
  data = np.genfromtxt("output/dynamic_topography_surface.00000")
  x = data[:,0]
  y = data[:,1]
  topo = data[:,2]

  plt.plot(x, topo)
  plt.show()

# Get the point at the leftmost point, which is what is
# used in Zhong et al.
def get_midpoint():
  data = np.genfromtxt("output/dynamic_topography_surface.00000")
  x = data[:,0]
  y = data[:,1]
  topo = data[:,2]
  idx = abs(x-0.5).argmin()
  return topo[idx], x[idx]



# Run the model at different refinement levels
refinements = range(5,10)
h = np.array([np.power(0.5, r) for r in refinements])
convergence = []
for r in refinements:
  generate_prm(r)
  os.system('mpirun -n 4 ./aspect tmp.prm')
  topo, x = get_midpoint()
  convergence.append(np.abs(topo-analytic_solution(x)))

# Generate outputs.
np.savetxt('convergence.txt', np.vstack([h, np.array(convergence)]).T)
plt.loglog(h, convergence, '-*', basex=2)
plt.ylabel('log(Error)')
plt.xlabel('log(h)')
plt.title('Convergence')
plt.savefig('Zhong93_case1.png')
