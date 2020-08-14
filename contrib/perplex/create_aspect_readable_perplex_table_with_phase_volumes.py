from __future__ import absolute_import
from __future__ import print_function

import argparse
import os.path
from os import rename, remove
import sys
from subprocess import Popen, PIPE, STDOUT

"""
This standalone python file creates a merged P-T table from perplex output,
containing both the thermodynamic properties of the assemblage and
the volume fractions of each phase.

Before running this file, you should have run vertex on your build file
over the P-T range of interest at an appropriate resolution.

The following should be set in perplex_options.dat
proportions               vol
sample_on_grid            F

After this has been done, call this python script from the command line.
The following arguments are accepted by the file:
--werami_path path [path ...]
--project project name [project name ...]
--outfile output file
--n_pressures N_PRESSURES
--n_temperatures N_TEMPERATURES
[--pressure_range LOW_PRESSURE_LIMIT HIGH_PRESSURE_LIMIT]
[--temperature_range LOW_TEMPERATURE_LIMIT HIGH_TEMPERATURE_LIMIT]

The last two arguments are optional
(defaulting to the limits of the vertex-calculated grid),
but it is strongly recommended that you enter the limits
which are suitable for your problem.

In situations where a solvus is active
(such that a single solution is stable as two or more separate phases)
this code adds the components together.

All "/" and "-" characters in phase names are converted to underscores.
"""


def create_property_table(werami_path, project_name, outfile, n_pressures, n_temperatures, pressure_range=None, temperature_range=None):
    '''
    This function uses PerpleX's werami software to output a table file containing the following material properties.
    2 - Density (kg/m3)
    4 - Expansivity (1/K, for volume)
    19 - Heat Capacity (J/K/kg)
    13 - P-wave velocity (Vp, km/s)
    14 - S-wave velocity (Vs, km/s)
    18 - Enthalpy (J/kg)

    And another table file containing all of the phase modes (volumes or masses)
    '''

    print('Working on creating {0}x{1} P-T table file using werami. Please wait.\n'.format(n_pressures, n_temperatures))

    try:
        str2 = 'y\n{0} {1}\n{2} {3}\n'.format(pressure_range[0]/1.e5, pressure_range[1]/1.e5,
                                              temperature_range[0], temperature_range[1])
    except:
        print('Keeping P-T range the same as the original project range.\n')
        str2 = 'n\n'

    stdin=('{0:s}\n2\n'
           '2\nn\n'
           '4\nn\n'
           '19\nn\n'
           '13\nn\n'
           '14\nn\n'
           '18\nn\n'
           '0\n'
           '{1:s}'
           '{2:d} {3:d}\n'
           '0'.format(project_name, str2, n_pressures, n_temperatures))


    p = Popen(werami_path, stdout=PIPE, stdin=PIPE, stderr=STDOUT, text=True)
    stdout = p.communicate(input=stdin)[0]
    print(stdout)
    out = [s for s in stdout.split('\n') if "Output has been written to the" in s][0].split()[-1]
    rename(out, "properties.tmp")

    # Do the same for the phase modes (volumes)
    stdin=('{0:s}\n2\n'
           '25\nn\n'
           '{1:s}'
           '{2:d} {3:d}\n'
           '0'.format(project_name, str2, n_pressures, n_temperatures))

    p = Popen(werami_path, stdout=PIPE, stdin=PIPE, stderr=STDOUT, text=True)
    stdout = p.communicate(input=stdin)[0]
    print(stdout)
    out = [s for s in stdout.split('\n') if "Output has been written to the" in s][0].split()[-1]
    rename(out, "volumes.tmp")

    write_combined_output_file(outfile)
    remove("properties.tmp")
    remove("volumes.tmp")

def write_combined_output_file(outfile):
    with open(args.outfile[0], 'w') as outfile:
        with open('properties.tmp', 'r') as fp:
            i = 0
            while i < 11:
                outfile.write(fp.readline())
                i += 1

            with open('volumes.tmp', 'r') as fv:
                for i in range(11):
                    linev = fv.readline()

                # Now we get the numbers of columns
                linep = fp.readline()
                linev = fv.readline()

                np = int(linep.strip().split()[0])

                # Now we get the column labels
                linep = fp.readline()
                linev = fv.readline().replace('/', '_').replace('-', '_')
                all_phases = linev.strip().split()[2:]
                unique_phases = list(set(all_phases))

                nc_unique = len(unique_phases)
                outfile.write(str(np + nc_unique)+'\n')

                phase_indices = [[i for i, x in enumerate(all_phases) if x == phase]
                                 for phase in unique_phases]

                labels = linep.strip().split()
                labels.extend(["vol_fraction_"+phase for phase in unique_phases])

                outfile.write(' '.join(labels)+'\n')


                # Now we go through the rest of the lines
                linep = fp.readline()
                linev = fv.readline()
                while linep:
                    strw = linep.strip().split()
                    vs = list(map(float,
                                  linev.replace('NaN', "0").strip().split()[2:]))

                    fractions = ['{0:.6f}'.format(sum([vs[i] for i in indices])/100.) for indices in phase_indices]


                    strw += fractions
                    outfile.write(' '.join(strw)+'\n')

                    linep = fp.readline()
                    linev = fv.readline()


parser = argparse.ArgumentParser(description='Call werami to create a burnman-readable tab file.')

parser.add_argument('--werami_path', metavar='path', type=str, nargs='+', required=True,
                    help='The path to werami')
parser.add_argument('--project', metavar='project name', type=str, nargs='+', required=True,
                    help='The name of the project file (without the suffix)')
parser.add_argument('--outfile', metavar='output file', type=str, nargs='+', required=True,
                    help='The name of the output table file')
parser.add_argument('--n_pressures', type=int, nargs=1, required=True,
                    help='The number of pressure steps in the grid')
parser.add_argument('--n_temperatures', type=int, nargs=1, required=True,
                    help='The number of pressure steps in the grid')
parser.add_argument('--pressure_range', type=float, nargs=2,
                    help='Minimum and maximum values of pressure (Pa; optional)')
parser.add_argument('--temperature_range', type=float, nargs=2,
                    help='Minimum and maximum values of temperature (K; optional)')

args = parser.parse_args()
if not hasattr(args, 'pressure_range'):
    args.pressure_range = None
if not hasattr(args, 'temperature_range'):
    args.temperature_range = None

create_property_table(args.werami_path, args.project[0], args.outfile[0], args.n_pressures[0], args.n_temperatures[0], args.pressure_range, args.temperature_range)
