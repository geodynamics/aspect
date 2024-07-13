#!/usr/bin/env python3
# -*- coding: utf-8 -*-

""" Routines to read and write input files (.prm files) of ASPECT """

import os

__author__ = 'The authors of the ASPECT code'
__copyright__ = 'Copyright 2023, ASPECT'
__license__ = 'GNU GPL 2 or later'


def get_parameter_value(parameters, name):
    """ Given a dictionary of parameters with a structure as
    the one created by read_parameter_file(), return the value
    of the parameter with the given name. Returns None if the
    parameter is not found.
    """

    if name in parameters:
        return parameters[name]["value"]
    else:
        for entry in parameters:
            if parameters[entry]["type"] == "subsection":
                value = get_parameter_value(parameters[entry]["value"], name)
                if value != None:
                    return value

    return None



def set_parameter_value(parameters, name, value):
    """ Given a dictionary of parameters with a structure as
    the one created by read_parameter_file(), set the value
    of the parameter with the given name to the given value.
    Returns 0 if the parameter was found and set, 1 otherwise.
    """

    if name in parameters:
        parameters[name]["value"] = value
        return 0
    else:
        for entry in parameters:
            if parameters[entry]["type"] == "subsection":
                if set_parameter_value(parameters[entry]["value"], name, value) == 0:
                    return 0

    return 1


def split_parameter_line(line):
    """ Read a 'set parameter' line and extract the name, value, and format. """
    
    equal_index = line.find('=')
    # Determine number of additional spaces left of equal sign.
    # We need to store these separately, because there might be two
    # lines setting the same parameter, but different number of spaces. They
    # should still map to the same parameter.
    alignment_spaces = len(line[:equal_index]) - len(line[:equal_index].rstrip())

    # skip initial word "set" in parameter name
    words_left_of_equal_sign = line[:equal_index].replace("\\\n","").split()
    param_name = words_left_of_equal_sign[1:]
    param_name = ' '.join(param_name)

    # strip spaces at end of value string, keep all inline comments
    param_value = line[equal_index+1:].rstrip()
    # if there is one or more spaces at the start of the string, remove one space.
    # keep additional spaces, because they are part of the formatting.
    # this way we can add one space when writing the file back and keep the formatting
    if param_value.startswith(' ') and len(param_value) > 1:
        param_value = param_value[1:]

    return param_name, param_value, alignment_spaces



def read_value_or_subsection(input_file, parameters):
    """ Read a value or a subsection from a parameter file into a parameter dictionary.
    """

    # Keep track of the comment lines to
    # add them to the next parameter or subsection
    accumulated_comment = ""

    for line in input_file:
        line = line.strip()

        # If line ends with \, read next line and append including \n
        if line != "":
            while line[-1] == "\\":
                line += "\n" + next(input_file).rstrip()

        words = line.replace("\\\n","").split()

        # Attach empty lines if we are inside a comment
        if line == "" or len(words) == 0:
            if accumulated_comment != "":
                accumulated_comment += "\n"

        # Store comments in accumulated_comment
        elif line[0] == "#":
            if accumulated_comment != "":
                accumulated_comment += "\n"
            accumulated_comment += line.replace("\\\n","")

        # If we encounter a 'subsection' line, store the subsection name
        # and recursively call this function to read the subsection
        elif words[0] == "subsection":
            subsection_name = ' '.join(words[1:]).replace("\\","")

            if subsection_name not in parameters:
                parameters[subsection_name] = {"comment": "", "value" : dict({}), "type": "subsection"}

            if parameters[subsection_name]["comment"] != "":
                parameters[subsection_name]["comment"] += "\n"

            parameters[subsection_name]["comment"] += accumulated_comment
            accumulated_comment = ""

            parameters[subsection_name]["value"].update(read_value_or_subsection(input_file, parameters[subsection_name]["value"]))

        # If we encounter a 'set' line, store the parameter name and value
        elif words[0] == 'set':
            name, value, alignment_spaces = split_parameter_line(line)

            parameters[name] = {"comment": accumulated_comment, "value": value, "alignment spaces": alignment_spaces, "type": "parameter"}
            accumulated_comment = ""

        elif words[0] == "include":
            name = ' '.join(words[1:])
            value = ''
            alignment_spaces = 0
            parameters[name] = {"comment": accumulated_comment, "value": value, "alignment spaces": alignment_spaces, "type": "include"}
            accumulated_comment = ""

        elif words[0] == "end":
            # sometimes there are comments at the end of files or subsection. store them in their own parameter
            if accumulated_comment != "":
                parameters["postcomment"] = {"comment": accumulated_comment, "value": "", "alignment spaces": 0, "type": "comment"}
                accumulated_comment = ""

            return parameters

        else:
            raise RuntimeError("Unrecognized first keyword in line: " + line)

    # sometimes there are comments at the end of files or subsection. store them in their own parameter
    if accumulated_comment != "":
        parameters["postcomment"] = {"comment": accumulated_comment, "value": "", "alignment spaces": 0, "type": "comment"}
        accumulated_comment = ""

    return parameters



def read_parameter_file(file):
    """ Read parameter file, and return a dictionary with all values. 
    
    The returned dictionary contains as keys the name of the subsection,
     include statement or parameter, and as value a dictionary with
     the following keys:
    - 'comment': the comment lines before the parameter or subsection
    - 'value': the value of the parameter as a string, or a dictionary
      with the same structure as the one returned by this function
      for a subsection. Empty for an include statement (the content
      of an include statement is its name).
    - 'alignment spaces': the number of spaces between name and equal sign
      in a parameter line. This is used to keep the formatting
      of the parameter file when writing it back.
    - 'type': either 'parameter', 'include', 'subsection', or 'comment'
      depending on the type of the entry.
    """

    parameters = dict({})
    
    with open(file, 'r') as f:
        read_value_or_subsection(f, parameters)
   
    return parameters



def write_comment_if_existent(comment, prepend_spaces, file):
    """ Print the comment lines, take care to use the correct indentation. 
    """
    if comment != "":
        for line in comment.split("\n"):
            if line != "":
                file.write(" " * prepend_spaces)
            file.write(line + "\n")


def write_value_or_subsection(parameters, prepend_spaces, file):
    """ Write a parameters dictionary recursively into a given file.
    prepend_spaces tracks the current indentation level.
    """

    # keep track of whether this is the first entry
    # this is important for some formatting
    first_entry = True

    # Print all parameter/include entries first
    for entry in parameters:
        if parameters[entry]["type"] == "parameter":
            # Add empty line before comments
            if first_entry == False and parameters[entry]["comment"] != "":
                file.write("\n")
            write_comment_if_existent(parameters[entry]["comment"], prepend_spaces, file)

            # Only add a space after equal sign if the value is not empty
            space_if_value = " " if parameters[entry]["value"] != "" else ""
            file.write(" " * prepend_spaces + "set " + entry + " " * parameters[entry]["alignment spaces"] + "=" \
                        + space_if_value + parameters[entry]["value"] + "\n")
            first_entry = False

        elif parameters[entry]["type"] == "include":
            # Add empty line before include
            if first_entry == False:
                file.write("\n")
            write_comment_if_existent(parameters[entry]["comment"], prepend_spaces, file)

            # Write the parameter
            file.write(" " * prepend_spaces + "include " + entry + "\n")
            # Add an empty line after include, except for the end of the file
            if entry != list(parameters)[-1]:
                file.write("\n")

            first_entry = False
    
    # Then print all subsections and end of subsection/file comments
    for entry in parameters:
        if parameters[entry]["type"] == "subsection":
            if first_entry == False:
                file.write("\n")
            write_comment_if_existent(parameters[entry]["comment"], prepend_spaces, file)

            file.write(" " * prepend_spaces + "subsection " + entry + "\n")
            prepend_spaces += 2
            # Recursively write the next subsection level
            write_value_or_subsection(parameters[entry]["value"], prepend_spaces, file)
            prepend_spaces -= 2
            file.write(" " * prepend_spaces + "end\n")
            first_entry = False

        elif parameters[entry]["type"] == "comment":
            if first_entry == False and parameters[entry]["comment"] != "":
                file.write("\n")
            write_comment_if_existent(parameters[entry]["comment"], prepend_spaces, file)
            first_entry = False

    return 0



def write_parameter_file(parameters, file):
    """ Given a dictionary of parameters with a structure as
    the one created by read_parameter_file(), write the dictionary
    to a file.
    """

    prepend_spaces = 0
    with open(file,'w') as f:
        write_value_or_subsection(parameters, prepend_spaces, f)
    return 0
