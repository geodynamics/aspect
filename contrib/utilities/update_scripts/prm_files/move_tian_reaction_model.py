#!/usr/bin/env python3
# -*- coding: utf-8 -*-

""" This script reformats the given .prm files to move the Tian 2019 reaction
mode parameters to the correct subsection.
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



def create_tian2019_reaction_subsection(parameters):
    """ Move the Tian 2019 reaction parameters to their own subsection. """


    parameters_to_move = ["Maximum weight percent water in sediment", \
                          "Maximum weight percent water in MORB", \
                          "Maximum weight percent water in gabbro", \
                          "Maximum weight percent water in peridotite"]

    # Collect existing parameters and delete old entries
    reactive_fluid_params = dict({})
    if "Reactive Fluid Transport Model" in parameters["Material model"]["value"]:
        reactive_fluid_params = parameters["Material model"]["value"]["Reactive Fluid Transport Model"]
        for param in parameters_to_move:
            if "Tian 2019 model" not in reactive_fluid_params["value"]:
                reactive_fluid_params["value"]["Tian 2019 model"] = {"comment": "", "value" : dict({}), "type": "subsection"}
            
            if param in reactive_fluid_params["value"]:
                reactive_fluid_params["value"]["Tian 2019 model"]["value"][param] = reactive_fluid_params["value"][param]
                del reactive_fluid_params["value"][param]

    return parameters



def main(input_file, output_file):
    """Echo the input arguments to standard output"""
    parameters = aspect.read_parameter_file(input_file)

    parameters = create_tian2019_reaction_subsection(parameters)

    aspect.write_parameter_file(parameters, output_file)



if __name__ == '__main__':
    parser = argparse.ArgumentParser(
                    prog='ASPECT .prm file reformatter',
                    description='Reformats ASPECT .prm files to follow our general formatting guidelines. See the documentation of this script for details.')
    parser.add_argument('input_file', type=str, help='The .prm file to reformat')
    parser.add_argument('output_file', type=str, help='The .prm file to write the reformatted file to')
    args = parser.parse_args()

    sys.exit(main(args.input_file, args.output_file))
