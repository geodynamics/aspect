#!/usr/bin/python

# This script can be used to plot errors for the
# instantaneous annulus benchmark. See the doc/
# folder for the output and a description of the
# benchmark

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

refinements = ["0","1","2","3","4","5","6","7"]
models = ["k_1","k_2","k_4","k_8"]
labels = ["k = 1","k = 2", "k = 4","k = 8"]
errors = ["u_L2","p_L2"]
ylabels = [r"$\|\boldsymbol{u} - \boldsymbol{u}_h\|_{L_2}$",r"$\|p - p_h\|_{L_2}$",r"$\|\rho - \rho_h\|_{L_2}$"]
markers=['o','X','P','v','s','D','<','>','^','+','x']
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

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

def setup_figure():
    figsize=(13,4)
    plt.rc('text', usetex=True)
    plt.rcParams['text.latex.preamble']='\\usepackage{relsize} \\usepackage{amsmath}'
    plt.rc('font', family='sanserif', size="16")
    plt.figure(figsize=figsize)
    return None

def plot_error_over_resolution(statistics, timestep, output_file):
    setup_figure()
    error_values = {}

    for model in models:
        error_values[model] = {}
        for error in errors:
            error_values[model][error] = []
            for refinement in refinements:
                error_values[model][error].append(statistics[refinement][model].iloc[timestep][error])

    scale_factors = [1.0,0.8]
    x = np.linspace(7e-3,1.0,100)
    y = x
    y2 = x*x
    y3 = x*x*x

    for i_error in range(2):
        ax = plt.subplot(1,2,i_error+1)
        ax.set_ybound(1e-7,0.5)
        ax.plot(x,scale_factors[i_error] * y, label='$h$', linestyle='--', color='grey')
        ax.plot(x,scale_factors[i_error] * y2, label='$h^2$', linestyle='-.', color='grey')
        ax.plot(x,scale_factors[i_error] * y3, label='$h^3$', linestyle=':', color='grey')

        ax.set_ylabel(ylabels[i_error])

        for i in range(len(models)):
            ax.loglog(h,error_values[models[i]][errors[i_error]], marker=markers[i], label=labels[i], color=colors[i])

        plt.xlabel("h")

    ax.legend(ncol=2, bbox_to_anchor=(1.05, 0.5), loc='upper left', borderaxespad=0.)
    plt.savefig(output_file, bbox_inches='tight',dpi=200)
    return None

statistics = {}

for refinement in refinements:
    statistics[refinement] = {}
    for model in models:
        statistics[refinement][model] = read_statistics("output-refinement_" + refinement + "_" + model + "/statistics")

plot_error_over_resolution(statistics, -1, "errors_annulus.svg")
