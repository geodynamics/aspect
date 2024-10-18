#!/usr/bin/env python3
# -*- coding: utf-8 -*-

""" This script reformats the given .prm files to move particle parameters
to the correct subsection.
"""

import sys
import os
import re
import argparse

__author__ = 'The authors of the ASPECT code'
__copyright__ = 'Copyright 2024, ASPECT'
__license__ = 'GNU GPL 2 or later'

# Add the ASPECT root directory to the path so we can import from the aspect_data module
sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..', '..'))
import python.scripts.aspect_input as aspect



def move_particle_parameters_to_own_subsection(parameters):
    """ Move the particle parameters to their own subsection. """

    # Collect existing parameters and delete old entries
    particle_params = dict({})
    if "Postprocess" in parameters:
        if "Particles" in parameters["Postprocess"]["value"]:
            particle_params = parameters["Postprocess"]["value"]["Particles"]
            del parameters["Postprocess"]["value"]["Particles"]

    # If there are any particle parameters
    if "value" in particle_params:
        # Make global section if necessary
        if not "Particles" in parameters:
            parameters["Particles"] = {"comment": "", "value" : dict({}), "type": "subsection"}

        # Move parameters one by one, throw an error if it already exists
        for parameter in particle_params["value"]:
            if not parameter in parameters["Particles"]["value"]:
                parameters["Particles"]["value"][parameter] = particle_params["value"][parameter]
            else:
                assert False, "Error: Duplicated parameter found with name: " + str(parameter)

    return parameters



def move_number_of_particles_to_correct_subsection(parameters):
    """ Move the number of particles to the correct subsection. """

    # Find the particle parameters and move
    if "Particles" in parameters:
        if "Number of particles" in parameters["Particles"]["value"]:
            parameter = parameters["Particles"]["value"]["Number of particles"]
            del parameters["Particles"]["value"]["Number of particles"]

            # We need to move the parameter, create the subsection if necessary
            if not "Generator" in parameters["Particles"]["value"]:
                parameters["Particles"]["value"]["Generator"] = {"comment": "", "value" : dict({}), "type": "subsection"}

            # Figure out if a generator is manually selected, if so move the parameter into their subsection
            if "Particle generator name" in parameters["Particles"]["value"]:
                generator = parameters["Particles"]["value"]["Particle generator name"]["value"]

                # if one of the generators was selected that use 'Number of particles' move it
                if generator == "uniform box":
                    if not "Uniform box" in parameters["Particles"]["value"]["Generator"]["value"]:
                        parameters["Particles"]["value"]["Generator"]["value"]["Uniform box"] = {"comment": "", "value" : dict({}), "type": "subsection"}
                    parameters["Particles"]["value"]["Generator"]["value"]["Uniform box"]["value"]["Number of particles"] = parameter

                elif generator == "uniform radial":
                    if not "Uniform radial" in parameters["Particles"]["value"]["Generator"]["value"]:
                        parameters["Particles"]["value"]["Generator"]["value"]["Uniform radial"] = {"comment": "", "value" : dict({}), "type": "subsection"}
                    parameters["Particles"]["value"]["Generator"]["value"]["Uniform radial"]["value"]["Number of particles"] = parameter

                elif generator == "random uniform":
                    if not "Random uniform" in parameters["Particles"]["value"]["Generator"]["value"]:
                        parameters["Particles"]["value"]["Generator"]["value"]["Random uniform"] = {"comment": "", "value" : dict({}), "type": "subsection"}
                    parameters["Particles"]["value"]["Generator"]["value"]["Random uniform"]["value"]["Number of particles"] = parameter

                elif generator == "probability density function":
                    if not "Probability density function" in parameters["Particles"]["value"]["Generator"]["value"]:
                        parameters["Particles"]["value"]["Generator"]["value"]["Probability density function"] = {"comment": "", "value" : dict({}), "type": "subsection"}
                    parameters["Particles"]["value"]["Generator"]["value"]["Probability density function"]["value"]["Number of particles"] = parameter

                # the parameter was not used by other generators. silently delete it
            else:
                # No generator was manually selected, move the parameter into the default generator subsection
                if not "Random uniform" in parameters["Particles"]["value"]["Generator"]["value"]:
                    parameters["Particles"]["value"]["Generator"]["value"]["Random uniform"] = {"comment": "", "value" : dict({}), "type": "subsection"}

                parameters["Particles"]["value"]["Generator"]["value"]["Random uniform"]["value"]["Number of particles"] = parameter

    # If 'random uniform' is selected or defaulted, but there is no 'random uniform' subsection, move the 'probability density function' subsection
    # this is necessary, because for a while the 'probability density function' section was used for the 'random uniform' generator as well
    for particle_section in ["Particles","Particles 2"]:
        if particle_section in parameters:
            if ("Particle generator name" in parameters[particle_section]["value"] and parameters[particle_section]["value"]["Particle generator name"]["value"] == "random uniform") \
                or not "Particle generator name" in parameters[particle_section]["value"]:
                if "Generator" in parameters["Particles"]["value"]:
                    # if the parameter does not already exist in random uniform
                    if not "Random uniform" in parameters[particle_section]["value"]["Generator"]["value"] \
                    and "Probability density function" in parameters[particle_section]["value"]["Generator"]["value"]:
                        subsection = parameters[particle_section]["value"]["Generator"]["value"]["Probability density function"]
                        parameters[particle_section]["value"]["Generator"]["value"]["Random uniform"] = subsection
                        del parameters[particle_section]["value"]["Generator"]["value"]["Probability density function"]
        
    return parameters

def move_particle_postprocess_parameters_back(parameters):
    """ Move the particle postprocessor parameters back to the postprocess subsection. """

    parameters_to_move = ["Time between data output", \
                          "Data output format", \
                          "Number of grouped files", \
                          "Write in background thread", \
                          "Temporary output location", \
                          "Exclude output properties"]

    # Find the particle parameters and move
    if "Particles" in parameters:
        for param in parameters_to_move:
            if param in parameters["Particles"]["value"]:
                # Create a new particles subsection in postprocess if necessary
                if "Particles" not in parameters["Postprocess"]["value"]:
                    parameters["Postprocess"]["value"]["Particles"] = {"comment": "", "value" : dict({}), "type": "subsection"}

                parameters["Postprocess"]["value"]["Particles"]["value"][param] = parameters["Particles"]["value"][param]
                del parameters["Particles"]["value"][param]

    return parameters

def remove_update_ghost_particles_parameter(parameters):
    """ Remove the parameter 'Update ghost particles'. """

    # Find the "Update ghost particles" parameter and remove
    if "Particles" in parameters:
        if "Update ghost particles" in parameters["Particles"]["value"]:
            del parameters["Particles"]["value"]["Update ghost particles"]

    return parameters

def main(input_file, output_file):
    """Echo the input arguments to standard output"""
    parameters = aspect.read_parameter_file(input_file)

    parameters = move_particle_parameters_to_own_subsection(parameters)
    parameters = move_number_of_particles_to_correct_subsection(parameters)
    parameters = move_particle_postprocess_parameters_back(parameters)
    parameters = remove_update_ghost_particles_parameter(parameters)

    aspect.write_parameter_file(parameters, output_file)



if __name__ == '__main__':
    parser = argparse.ArgumentParser(
                    prog='ASPECT .prm file reformatter',
                    description='Reformats ASPECT .prm files to follow our general formatting guidelines. See the documentation of this script for details.')
    parser.add_argument('input_file', type=str, help='The .prm file to reformat')
    parser.add_argument('output_file', type=str, help='The .prm file to write the reformatted file to')
    args = parser.parse_args()

    sys.exit(main(args.input_file, args.output_file))
