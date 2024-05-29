#!/usr/bin/env python3

""" This script reads the statistics from the output of the benchmark models and
writes them to a text file. It also plots the results and saves them to a png file.
See README.md for a description of how to run the benchmark models.
"""

import sys
import os
import re

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# Change current path to the directory of this script
os.chdir(os.path.dirname(__file__))

# Add the ASPECT root directory to the path so we can import from the aspect_data module
sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..'))
from contrib.python.scripts.aspect_data import read_statistics

# Define the cases and refinements that will appear in the output
cases = ["1a", "1b", "1c", "2a", "2b"]
refinements = [3, 4, 5, 6, 7]

# Read the reference statistics from file
statistics = pd.read_csv("reference_statistics.txt", sep='\s+', engine='python')

# Add columns for computed values
reference_model_name = "ASPECT"

def fill_statistics(case, refinement, statistics_path):
    """ Read the statistics of one model into the Pandas Dataframe."""
    print (statistics_path)
    model_statistics = read_statistics(statistics_path)
    
    Nu_ref = statistics.loc[(statistics['Case'] == case) & \
                            (statistics['#Cells'] == reference_model_name),'Nu'].to_numpy()[0]
    Nu = model_statistics.iloc[-1]["Outward heat flux density for boundary with id 3 (\"top\") (W/m)"]
    Nu_g = model_statistics.iloc[-1]["Outward heat flux (gradient) through boundary with indicator 3 (\"top\") (W)"]
    # The Nu_g output is not normalized by model width, do it here.
    # This only matters for case 2b, because that is the only one not on a unit square.
    model_width = 2.5 if case == "2b" else 1.0
    Nu_g /= model_width
    epsilon = np.abs(Nu - Nu_ref) / Nu_ref
    epsilon_g = np.abs(Nu_g - Nu_ref) / Nu_ref

    # Add the statistics to the dataframe
    statistics.loc[len(statistics)] = [case, \
                                        2**refinement, \
                                        1/2**refinement, \
                                        Nu, \
                                        epsilon, \
                                        Nu_g, \
                                        epsilon_g, \
                                        model_statistics.iloc[-1]["RMS velocity (m/s)"]]

def extrapolate_property(statistics, case, property):
    """ Extrapolate a property from the statistics of the finest mesh. """

    # Find the finest mesh
    # Unexpectedly case 1b resolution 32 and 128 have almost the
    # same value for Nu, and 64 deviates. This confuses the extrapolation.
    # Therefore we skip the second finest mesh (64) for extrapolating case 1b Nu.
    # This is likely a problem with the benchmark model output, not with the extrapolation.
    finest_index = -1
    coarser_index = finest_index - 2 if (case == "1b" and property == "Nu") else finest_index - 1
    coarsest_index = finest_index - 3 if (case == "1b" and property == "Nu") else finest_index - 2

    finest_property = statistics.loc[(statistics['Case'] == case) & \
                                    (statistics['#Cells'] == 2**refinements[finest_index]),property].to_numpy()[0]
    coarser_property = statistics.loc[(statistics['Case'] == case) & \
                                    (statistics['#Cells'] == 2**refinements[coarser_index]),property].to_numpy()[0]
    coarsest_property = statistics.loc[(statistics['Case'] == case) & \
                                    (statistics['#Cells'] == 2**refinements[coarsest_index]),property].to_numpy()[0]

    # Compute the convergence rate
    convergence_rate = np.log(np.abs(coarser_property - coarsest_property) / np.abs(finest_property - coarser_property)) / np.log(2**refinements[finest_index]/2**refinements[coarser_index])
    extrapolated_property = finest_property + (finest_property - coarser_property) / (2**(convergence_rate)-1)

    # Store the extrapolated property
    statistics.loc[(statistics['Case'] == case) & \
                   (statistics['#Cells'] == reference_model_name), property] = extrapolated_property

# Read the statistics from the output of the benchmark models
for case in cases:
    for refinement in refinements:
        case_directory = "output-case" + case + "_ref" + str(refinement)
        statistics_path = os.path.join(os.path.dirname(__file__), case_directory, "statistics")
        fill_statistics(case, refinement, statistics_path)

# Compute extrapolated values
for case in cases:
    for property in ["Nu", "Nu_g", "Vrms"]:
        extrapolate_property(statistics, case, property)

# Print the results to screen and to a text file
statistics.sort_values(by=['Case', '#Cells'], inplace=True)
print (statistics)
with open('statistics.txt', 'w') as out_file:
    statistics.to_string(out_file, header=True, index=False, float_format='%10.8g')

# Uncomment to update the reference_statistics.txt
# with open('reference_statistics.txt', 'w') as out_file:
#     statistics[(statistics['h'] == 0.0)].to_string(out_file, header=True, index=False, float_format='%10.8g')

# Uncomment to print a selection of the results into a LaTeX table
# This requires the Jinja2 package to be installed.
# statistics.to_latex('statistics.tex', \
#                     index=False, \
#                     float_format='%10.8g', \
#                     columns=['Case', '#Cells', 'Nu', 'Epsilon', 'Nu_g', 'Epsilon_g', 'Vrms'])

# Plot the results
figsize=(12,6)
fig, ax = plt.subplots(ncols=2, figsize=figsize)
prop={'size':12}
plt.rc('font', **prop)
cmap = plt.get_cmap("tab10")

# Plot model Nu values
results = statistics[(statistics['h'] != 0.0)]
for case, marker, color in zip(cases, Line2D.markers, cmap.colors[1:]):
    results[(results['Case'] == case)].plot(x='h', y='Nu_g', ax=ax[0], style='--', marker=marker, label="_Case " + case, color=color)

for case, marker, color in zip(cases, Line2D.markers, cmap.colors[1:]):
    results[(results['Case'] == case)].plot(x='h', y='Nu', ax=ax[0], style='-', marker=marker, label="Case " + case, color=color)

# Plot reference benchmark values
x = np.linspace(1/2**7,1/2**3,100)
xticks = [1.0/2**3, 1.0/2**4, 1.0/2**5, 1.0/2**6, 1/2**7]
xlabels = [str(xticks[0]), str(xticks[1]), str(xticks[2]), str(xticks[3]), str(xticks[4])]
linestyles = ['--', '-.', 'solid', ':', (-3,(7,10))]

for case, linestyle in zip(cases, linestyles):
    ref_index = statistics.index[(statistics['Case'] == case) & (statistics['#Cells'] == "Reference")].to_flat_index()[0]
    ref_Nu = statistics.loc[ref_index,'Nu']
    ref = ax[0].plot(x,ref_Nu * np.ones_like(x), label=case, linestyle=linestyle, linewidth=2, color='grey')

# Prettify image
ax[0].set_xscale('log')
ax[0].set_xticks(ticks=xticks, labels=xlabels)
ax[0].set_xticks(ticks=[], labels=[], minor=True)
ax[0].set_ylabel('Nu')
ax[0].set_xlabel('h')
ax[0].legend(ncols=2, title='Benchmarks / Reference values', bbox_to_anchor=(0.71,0.76))
ax[0].set_ybound([4.5,25])

# Plot reference benchmark relative errors
for case, marker, color in zip(cases, Line2D.markers, cmap.colors[1:]):
    results[(results['Case'] == case)].plot(x='h', y='Epsilon_g', ax=ax[1], style='--', marker=marker, label="_Case " + case, color=color)

for case, marker, color in zip(cases, Line2D.markers, cmap.colors[1:]):
    results[(results['Case'] == case)].plot(x='h', y='Epsilon', ax=ax[1], style='-', marker=marker, label="_Case " + case, color=color)

# Plot theoretical convergence rates
x = np.linspace(1/2**7,1/2**3,100)
y = 10 * x
y2 = 5 * x*x
y3 = 0.1 * x*x*x

ax[1].plot(x,y, label='$h$', linestyle='--', color='grey')
ax[1].plot(x,y2, label='$h^2$', linestyle='-.', color='grey')
ax[1].plot(x,y3, label='$h^3$', linestyle=':', color='grey')

# Prettify image
ax[1].set_xscale('log')
ax[1].set_xlabel('h')
ax[1].set_yscale('log')
ax[1].set_ylabel('$\epsilon_{Nu}$')
ax[1].set_ybound([1e-8,1])
ax[1].set_xticks(ticks=xticks, labels=xlabels)
ax[1].set_xticks(ticks=[], labels=[], minor=True)
ax[1].legend()

# Uncomment to show the plot instead of saving it to file
# plt.show()

plt.savefig("blankenbach.pdf", bbox_inches='tight',dpi=200)
plt.savefig("blankenbach.png", bbox_inches='tight',dpi=200)
