#!/usr/bin/env python3

""" This script reads the statistics from the output of the benchmark models and
writes them to a text file. It also plots the results and saves them to a png file.
See README.md for a description of how to run the benchmark models.
"""

import sys
import os
import re
import csv
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
formulations = ["ALA", "TALA"]
refinements = [3, 4, 5, 6, 7]
dissipation_numbers = [0.25, 0.5, 1.0]
rayleigh_numbers = ["1e4", "1e5"]

# Determine which model to use as reference for the convergence rates
# Other possible values are "KingCU" "KingVT" "KingUM"
# By default use the extrapolated value from ASPECT
reference_model_name = "ASPECT"

statistics = pd.read_csv("reference_statistics.txt", sep='\s+', engine='python')

def fill_statistics(formulation, Di, Ra, refinement, statistics_path):
    """ Read the statistics of one model into the Pandas Dataframe."""
    print (statistics_path)
    model_statistics = read_statistics(statistics_path)

    Nu_ref = statistics.loc[(statistics['Formulation'] == formulation) & \
                            (statistics['Di'] == Di) & \
                            (statistics['Ra'] == float(Ra)) & \
                            (statistics['Resolution'] == reference_model_name),'Nu'].to_numpy()[0]
    Nu = model_statistics.iloc[-1]["Outward heat flux through boundary with indicator 3 (\"top\") (W)"]
    Nu_g = model_statistics.iloc[-1]["Outward heat flux (gradient) through boundary with indicator 3 (\"top\") (W)"]

    # Adjust by adiabatic surface temperature
    T_avg = model_statistics.iloc[-1]["Average temperature (K)"] - 0.091
    # Adjust sign to match the conventions of the other benchmarks
    W = -model_statistics.iloc[-1]["Total adiabatic heating rate (W)"]

    epsilon = np.abs(Nu - Nu_ref) / Nu_ref
    epsilon_g = np.abs(Nu_g - Nu_ref) / Nu_ref

    statistics.loc[len(statistics)] = [formulation, \
                                       Di, \
                                       float(Ra), \
                                       2**refinement, \
                                       1/2**refinement, \
                                       Nu, \
                                       epsilon, \
                                       Nu_g, \
                                       epsilon_g, \
                                       model_statistics.iloc[-1]["RMS velocity (m/s)"], \
                                       T_avg, \
                                       model_statistics.iloc[-1]["Total shear heating rate (W)"], \
                                       W]
    
def extrapolate_property(statistics, formulation, Di, Ra, property):
    """ Extrapolate a property from the statistics of the finest mesh. """
    # Find the finest mesh
    finest_property = statistics.loc[(statistics['Formulation'] == formulation) & \
                                    (statistics['Di'] == Di) & \
                                    (statistics['Ra'] == float(Ra)) & \
                                    (statistics['Resolution'] == 2**refinements[-1]),property].to_numpy()[0]
    
    coarser_property = statistics.loc[(statistics['Formulation'] == formulation) & \
                                    (statistics['Di'] == Di) & \
                                    (statistics['Ra'] == float(Ra)) & \
                                    (statistics['Resolution'] == 2**refinements[-2]),property].to_numpy()[0]
    coarsest_property = statistics.loc[(statistics['Formulation'] == formulation) & \
                                    (statistics['Di'] == Di) & \
                                    (statistics['Ra'] == float(Ra)) & \
                                    (statistics['Resolution'] == 2**refinements[-3]),property].to_numpy()[0]

    # Compute the convergence rate
    convergence_rate = np.log(np.abs(coarser_property - coarsest_property) / np.abs(finest_property - coarser_property)) / np.log(2)

    extrapolated_property = finest_property + (finest_property - coarser_property) / (2**(convergence_rate)-1)

    # Store the extrapolated property
    statistics.loc[(statistics['Formulation'] == formulation) & \
                                    (statistics['Di'] == Di) & \
                                    (statistics['Ra'] == float(Ra)) & \
                                    (statistics['Resolution'] == reference_model_name), property] = extrapolated_property

for formulation in formulations:
    for Di in dissipation_numbers:
        for Rayleigh in rayleigh_numbers:
            for refinement in refinements:
                output_directory = "Di" + str(Di) + "_Ra" + str(Rayleigh) + "_refinement" + str(refinement) + "_" + formulation.lower()
                fill_statistics(formulation, Di, Rayleigh, refinement, output_directory + "/statistics")

# Compute extrapolated values
for formulation in formulations:
    for Di in dissipation_numbers:
        for Ra in rayleigh_numbers:
            for property in ["Nu", "Nu_g", "<T>", "Phi", "W", "Vrms"]:
                extrapolate_property(statistics, formulation, Di, Ra, property)

# Print the results to screen and to a text file
statistics.sort_values(by=['Formulation', 'Di', 'Ra', 'Resolution'], inplace=True)
print (statistics)
with open('statistics.txt', 'w') as out_file:
    statistics.to_string(out_file, header=True, index=False, float_format='%10.8g')

# Uncomment to update the reference_statistics.txt
# with open('reference_statistics.txt', 'w') as out_file:
#     statistics[(statistics['h'] == 0.0)].to_string(out_file, header=True, index=False, float_format='%10.8g')

# Uncomment to print a selection of the results into a LaTeX table
# This requires the Jinja2 package to be installed.
# statistics.loc[(statistics['Formulation'] == "ALA") & \
#                 (statistics['Ra'] == 1e5)].to_latex('statistics.tex', \
#                                                     index=False, \
#                                                     float_format='%10.8g', \
#                                                     columns=['Di', 'Resolution', 'Nu', 'Epsilon', 'Nu_g', 'Epsilon_g', 'Vrms', '<T>', 'Phi', 'W'])

# Plot the results
figsize=(12,6)
fig, ax = plt.subplots(ncols=2, figsize=figsize)
prop={'size':12}
plt.rc('font', **prop)
cmap = plt.get_cmap("tab10")

# Plot model Nu values
results =  statistics.loc[(statistics['Formulation'] == "ALA") & \
                           (statistics['Ra'] == 1e5) & \
                           (statistics['h'] != 0.0)]

for Di, marker, color in zip(dissipation_numbers, Line2D.markers, cmap.colors[1:]):
    results[(results['Di'] == Di)].plot(x='h', y='Nu_g', ax=ax[0], style='--', label="_Di " + str(Di), marker=marker, color=color)

for Di, marker, color in zip(dissipation_numbers, Line2D.markers, cmap.colors[1:]):
    results[(results['Di'] == Di)].plot(x='h', y='Nu', ax=ax[0], style='-', marker=marker, label="Di " + str(Di), color=color)

# Plot reference benchmark values
x = np.linspace(1/2**7,1/2**3,100)
xticks = [1.0/2**3, 1.0/2**4, 1.0/2**5, 1.0/2**6, 1.0/2**7]
xlabels = [str(xticks[0]), str(xticks[1]), str(xticks[2]), str(xticks[3]), str(xticks[4])]
linestyles = ['--', '-.', 'solid', ':', (-3,(7,10))]

reference_results = statistics.loc[(statistics['Formulation'] == "ALA") & \
                           (statistics['Ra'] == 1e5) & \
                           (statistics['Resolution'] == "KingVT")]

for Di, linestyle in zip(dissipation_numbers, linestyles):
    ref_Nu = reference_results.loc[(statistics['Di'] == Di),'Nu'].to_numpy()[0]
    ref = ax[0].plot(x,ref_Nu * np.ones_like(x), label='Ref Di ' + str(Di), linestyle=linestyle, linewidth=2, color='grey')

# Prettify image
ax[0].set_xscale('log')
ax[0].set_xticks(ticks=xticks, labels=xlabels)
ax[0].set_xticks(ticks=[], labels=[], minor=True)
ax[0].set_ylabel('Nu')
ax[0].set_xlabel('h')
ax[0].legend(ncols=2, title='Benchmarks / Reference values', loc='center left', bbox_to_anchor=(0.0,0.28))

# Plot reference benchmark relative errors
for Di, marker, color in zip(dissipation_numbers, Line2D.markers, cmap.colors[1:]):
    results[(results['Di'] == Di)].plot(x='h', y='Epsilon_g', ax=ax[1], style='--', label="_Di " + str(Di), marker=marker, color=color)

for Di, marker, color in zip(dissipation_numbers, Line2D.markers, cmap.colors[1:]):
    results[(results['Di'] == Di)].plot(x='h', y='Epsilon', ax=ax[1], style='-', marker=marker, label="_Di " + str(Di), color=color)

# Plot theoretical convergence rates
x = np.linspace(1/2**7,1/2**3,100)
y = 10 * x
y2 = 5 * x*x
y3 = 500.0 * x*x*x*x*x

ax[1].plot(x,y, label='$h$', linestyle='--', color='grey')
ax[1].plot(x,y2, label='$h^2$', linestyle='-.', color='grey')
ax[1].plot(x,y3, label='$h^5$', linestyle=':', color='grey')

# Prettify image
ax[1].set_xscale('log')
ax[1].set_xlabel('h')
ax[1].set_yscale('log')
ax[1].set_ylabel('$\epsilon_{Nu}$')
# ax[1].set_ybound([1e-8,1])
ax[1].set_xticks(ticks=xticks, labels=xlabels)
ax[1].set_xticks(ticks=[], labels=[], minor=True)
ax[1].legend()

# Uncomment to show the plot instead of saving it to file
# plt.show()

plt.savefig("king.pdf", bbox_inches='tight',dpi=200)
plt.savefig("king.png", bbox_inches='tight',dpi=200)
