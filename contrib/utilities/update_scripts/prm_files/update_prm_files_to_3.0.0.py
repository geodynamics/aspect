#!/usr/bin/env python3
# -*- coding: utf-8 -*-

""" This script reformats the given .prm files to follow our
general formatting guidelines. These are:

- indent by two spaces for each subsection level, and indent comments as well
  as content
- ensure an empty line before a subsection or a comment, unless it is a
  subsection/comment following directly another subsection, or it is at the
  start of the file
- combine content of subsections that are duplicated
- if values are duplicated, merge them and use the last one in the file (as
  ASPECT does when it reads the file)
- retain as much user formatting in comments and parameter values as possible,
  e.g. spaces for padding values to align between adjacent lines. This is not
  always perfectly possible.
- retain broken lines (`\`) in values of parameters and comments, remove them
  from subsection or parameter names.
"""

import sys
import os
import re
import argparse

__author__ = 'The authors of the ASPECT code'
__copyright__ = 'Copyright 2023, ASPECT'
__license__ = 'GNU GPL 2 or later'

# Add the ASPECT root directory to the path so we can import from the aspect_data module
sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..', '..'))
import python.scripts.aspect_input as aspect



def line_contains_pattern_to_remove(line):
    """ Check if the given line contains the given pattern. """

    patterns_to_remove = ["# The parameters below this comment were created by the update script",
                          "# as replacement for the old 'Model settings' subsection. They can be",
                          "# safely merged with any existing subsections with the same name."]
    for pattern in patterns_to_remove:
        if re.match(pattern, line):
            return True

    return False



def reformat (parameters):
    for entry in parameters:
        if parameters[entry]["comment"] != "":
            comment_lines = parameters[entry]["comment"].split("\n")
            formatted_comment = ""
            for comment_line in comment_lines:
                if line_contains_pattern_to_remove(comment_line):
                    pass
                else:
                    if formatted_comment != "":
                        formatted_comment += "\n"
                    formatted_comment += comment_line
            parameters[entry]["comment"] = formatted_comment

            if isinstance(parameters[entry]["value"], dict):
                parameters[entry]["value"] = reformat(parameters[entry]["value"])

    return parameters



def move_particle_parameters_to_own_subsection(parameters):
    """ Move the particle parameters to their own subsection. """

    # Find the particle parameters and remove from the main dictionary
    particles = {}
    if "Postprocess" in parameters:
        if "Particles" in parameters["Postprocess"]["value"]:
            particles = parameters["Postprocess"]["value"]["Particles"]
            del parameters["Postprocess"]["value"]["Particles"]

    # Add the particle parameters to their own subsection before the postprocess subsection
    new_parameters = {}
    for key,value in parameters.items():
        if key == "Postprocess" and particles != {}:
            new_parameters["Particles"] = particles
        new_parameters[key] = value

    return new_parameters



def move_number_of_particles_to_correct_subsection(parameters):
    """ Move the number of particles to the correct subsection. """

    # Find the particle parameters and remove
    number_of_particles = None
    if "Particles" in parameters:
        if "Number of particles" in parameters["Particles"]["value"]:
            number_of_particles = parameters["Particles"]["value"]["Number of particles"]
            del parameters["Particles"]["value"]["Number of particles"]

    # Add the number of particles to the correct subsection
    if number_of_particles is not None:
            # which generator is used?
            particle_generator = aspect.get_parameter_value_from_subsection(parameters,"Particle generator name", ["Particles"])
            if particle_generator == None:
                particle_generator = "random uniform"

            # Find the correct subsection
            particle_generator_subsection = ""
            if particle_generator == "random uniform" or particle_generator == "probability density function":
                particle_generator_subsection = "Probability density function"
            if particle_generator == "uniform box":
                particle_generator_subsection = "Uniform box"
            if particle_generator == "uniform radial":
                particle_generator_subsection = "Uniform radial"
            
            # All other generators ignore this parameter, only add if necessary
            if particle_generator_subsection != "":
                aspect.set_parameter_value(parameters, "Number of particles", number_of_particles, ["Particles", "Generator", particle_generator_subsection])

    return parameters



def main(input_file, output_file):
    """Echo the input arguments to standard output"""
    parameters = aspect.read_parameter_file(input_file)

    parameters = reformat(parameters)
    parameters = move_particle_parameters_to_own_subsection(parameters)
    parameters = move_number_of_particles_to_correct_subsection(parameters)

    aspect.write_parameter_file(parameters, output_file)



if __name__ == '__main__':
    parser = argparse.ArgumentParser(
                    prog='ASPECT .prm file reformatter',
                    description='Reformats ASPECT .prm files to follow our general formatting guidelines. See the documentation of this script for details.')
    parser.add_argument('input_file', type=str, help='The .prm file to reformat')
    parser.add_argument('output_file', type=str, help='The .prm file to write the reformatted file to')
    args = parser.parse_args()

    sys.exit(main(args.input_file, args.output_file))
