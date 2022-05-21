#!/usr/bin/python3

import json
from sys import argv, exit, stderr

def handle_subsection(data, cur_path):
    keys = list(data.keys())
    keys.sort()
    # Since there is no end to a markdown heading besides another markdown
    # heading of equal or greater size, we need to first do parameters and
    # aliases, and then do subsections otherwise parameters get nested
    # incorrectly
    for key in keys:
        path_str = key
        true_name = key.replace("_20", " ")
        
        if len(cur_path) != 0:
            path_str = ":".join(cur_path) + ":"  + key
        if "value" in data[key]:
            # this is a parameter
            print("(parameters:" + path_str + ")=")
            print("### __Parameter name:__ " + true_name)
            print("**Default value:** " +  data[key]["default_value"] + " \n")
            print("**Pattern:** " +  data[key]["pattern_description"] + " \n")
            print("**Documentation:** " + data[key]["documentation"] + " \n")
        elif "alias" in data[key]:
            # This is an alias for a parameter
            print("(parameters:" + path_str + ")=")
            aliased_name = data[key]["alias"]
            alias_path_str = ":".join(cur_path) + ":" + aliased_name.replace(" ", "_20") 
            print("### __Parameter name__: " +  true_name)
            print("**Alias:** [" + aliased_name + "](parameters:" + alias_path_str + ")\n")
            print("**Deprecation Status:** " + data[key]["deprecation_status"] + "\n")

    for key in keys:
        # Specific 20heats is a key in this json that isn't a parameter or subsection
        path_str = key
        if len(cur_path) != 0:
            path_str = ":".join(cur_path) + ":"  + key

        if (not "value" in data[key]) and (not "alias" in data[key]):
            # This is a subsection
            print("(parameters:" + path_str + ")=")
            new_path = cur_path + [key]
            section_name = path_str.replace("_20", " ").replace(":", "/")
            print("## **Parameters in section** " + section_name)
            handle_subsection(data[key], new_path)

   
if __name__  == "__main__":
    data = {}
    input_file = open(argv[1], "r")
    data = json.load(input_file)
    print("# Parameters\n\n")
    print("## **Parameters in section** \n\n")
    handle_subsection(data, [])
    exit(0)

