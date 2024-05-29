#!/usr/bin/python

# This script can be used to plot errors for the
# 'transient rigid shear' benchmark described in Gassmoeller
# et al. (2023), "Benchmarking the accuracy of higher order
# particle methods in geodynamic models of transient flow"


import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

refinements = ["2","3","4","5","6","7"]
models = ["analytical_density", "compositional_field","continuous_compositional_field", \
          "higher_order_true_interpolation_bilinear_least_squares", \
          "higher_order_false_interpolation_bilinear_least_squares"]

labels = ["Density: Benchmark", "Density: FE field ($DGQ_2$)","Density: FE field ($Q_2$)", "Density: Particles (RK2)","Density: Particles (RK2 FOT)"]
errors = ["u_L2","p_L2","rho_L2"]
ylabels = [r"$\|\boldsymbol{u} - \boldsymbol{u}_h\|_{L_2}$",r"$\|p - p_h\|_{L_2}$",r"$\|\rho - \rho_h\|_{L_2}$"]

markers=['v','<','>','^','o','D','<','>','^','+','x']
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
plt.rcParams['lines.markersize'] = 8.5

h = []
for refinement in refinements:
    h.append(1/(2**int(refinement)))

def read_statistics(fname):
    """ Read the statistics file output by ASPECT
    
    return a pandas table, where names are taken from the statistics file.
    """
    # header:
    header = []
    header_read = True

    with open(fname) as f:
        while header_read :
            line = f.readline()
            if line[0] == '#':
                idx_start = line.find(":")
                header.append(line[idx_start+2:-1])
            else:
                header_read = False
                
    # data
    values = pd.read_csv(fname, skiprows=len(header), header=None, delim_whitespace=True, names=header)
    return values

def plot_error_over_time(statistics, output_file):
    figsize=(4,10)
    plt.rc('text', usetex=True)
    plt.rcParams['text.latex.preamble']='\\usepackage{relsize} \\usepackage{amsmath}'
    plt.rc('font', family='sanserif', size="16")
    plt.figure(figsize=figsize)

    for i_error in range(3):
        ax = plt.subplot(3,1,i_error+1)
        ax.set_ylabel(ylabels[i_error])
        plt.xlim(2e-3,1.1)

        for model,label,marker in zip(models,labels,markers):
            ax.loglog(statistics[model]["Time (seconds)"],statistics[model][errors[i_error]], label=label)

        ax.set_ybound(5e-8,1.0)

    plt.xlabel("Time $t$")
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
    plt.savefig(output_file, bbox_inches='tight',dpi=200)
    return None

def plot_error_over_resolution(statistics, timestep, output_file):
    figsize=(4,10)
    plt.rc('text', usetex=True)
    plt.rcParams['text.latex.preamble']='\\usepackage{relsize} \\usepackage{amsmath}'
    plt.rc('font', family='sanserif', size="16")
    figure=plt.figure(dpi=100,figsize=figsize)

    error_values = {}

    for model in models:
        error_values[model] = {}
        for error in errors:
            error_values[model][error] = []
            for refinement in refinements:
                error_values[model][error].append(statistics[refinement][model].iloc[timestep][error])

    scale_factors = [1.5,12.0,11.0]
    x = np.linspace(7e-3,0.25,100)
    y = 0.1 * x
    y2 = 0.1 * x*x
    y3 = 0.1 * x*x*x

    for i_error in range(3):
        ax = plt.subplot(3,1,i_error+1)
        ax.plot(x,scale_factors[i_error] * y, label='$h$', linestyle='--', color='grey')
        ax.plot(x,scale_factors[i_error] * y2, label='$h^2$', linestyle='-.', color='grey')
        ax.plot(x,scale_factors[i_error] * y3, label='$h^3$', linestyle=':', color='grey')
        ax.set_ylabel(ylabels[i_error])

        for i in range(len(models)):
            # dont plot density error for the analytical model (it is 0)
            if i_error == 2 and i == 0:
                continue
            ax.loglog(h,error_values[models[i]][errors[i_error]], marker=markers[i], label=labels[i], color=colors[i])

        ax.set_ybound(5e-8,1.0)

        if i_error == 0:
            ax.legend(ncol=2, bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)

    plt.xlabel("Cell size $h$")
    plt.savefig(output_file, bbox_inches='tight',dpi=200)
    return None

statistics = {}

for refinement in refinements:
    statistics[refinement] = {}
    for model in models:
        statistics[refinement][model] = read_statistics("refinement_" + refinement + "_" + model + "/statistics")

for refinement in refinements:
    plot_error_over_time(statistics[refinement], 'refinement_' + refinement + '_error_over_time.png')

plot_error_over_resolution(statistics, -1, "error_over_resolution_end.png")
