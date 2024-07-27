#!/usr/bin/env python3
# -*- coding: utf-8 -*-

""" This script changes 'Use limiter for discontinuous temperature/composition solution'
to 'limiter for discontinuous temperature/composition solution' in order to allow users
to specify different limiters for different fields
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



def reformat_limiter_settings(parameters):
    """change 'Use limiter for discontinuous temperature/composition solution' 
    to 'limiter for discontinuous temperature/composition solution' in order to allow users
    to specify different limiters for different fields
    """

    old_entries = ["Use limiter for discontinuous temperature solution",
                   "Use limiter for discontinuous composition solution"]
    new_entries = ["Limiter for discontinuous temperature solution",
                   "Limiter for discontinuous composition solution"]

    # Go to the subsection
    if "Discretization" in parameters:
        if "Stabilization parameters" in parameters["Discretization"]["value"]:
            subsection = parameters["Discretization"]["value"]["Stabilization parameters"]["value"]

            for i in range(0,2):
                if old_entries[i] in subsection:

                    # Get the parameter value. Notice that the value may be followed by some
                    # extra comments, in which case we need to split them.
                    comments = subsection[old_entries[i]]["comment"]
                    values_and_extra_comments = subsection[old_entries[i]]["value"].split('#', 1)
                    values = values_and_extra_comments[0].rstrip()
                    extra_comments = '' if len(values_and_extra_comments) == 1 else " #" + values_and_extra_comments[1]

                    # Delete the old entry
                    del subsection[old_entries[i]]

                    new_values = ''
                    if i == 0:
                        # For temperature field, the parameter value is a boolean
                        new_values = "bound preserving" if values == "true" else "none"

                    else:
                        # For compositional fields, the parameter value may be a boolean or 
                        # a list of booleans.
                        split_values = values.split(',')
                        for i in range(0, len(split_values)):
                            if split_values[i] == "true":
                                new_values += "bound preserving"
                            else:
                                new_values += "none"

                            if i < len(split_values)-1:
                                new_values += ", "
                        
                    # append the new values with the extra comments
                    new_values += extra_comments
                    subsection[new_entries[i]] = {"comment"         : comments,
                                                  "value"           : new_values,
                                                  "type"            : "parameter",
                                                  "alignment spaces": 1}

    return parameters


        
def main(input_file, output_file):
    """Echo the input arguments to standard output"""
    parameters = aspect.read_parameter_file(input_file)
    parameters = reformat_limiter_settings(parameters)
    aspect.write_parameter_file(parameters, output_file)



if __name__ == '__main__':
    parser = argparse.ArgumentParser(
                    prog='ASPECT .prm file reformatter',
                    description='Reformats ASPECT .prm files to follow our general formatting guidelines. See the documentation of this script for details.')
    parser.add_argument('input_file', type=str, help='The .prm file to reformat')
    parser.add_argument('output_file', type=str, help='The .prm file to write the reformatted file to')
    args = parser.parse_args()

    sys.exit(main(args.input_file, args.output_file))
