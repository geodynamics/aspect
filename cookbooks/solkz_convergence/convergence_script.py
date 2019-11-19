#!/usr/bin/env python

import os
import json
import subprocess


aspect_binary = "./aspect"
postprocess_binary = "aspect-postprocess"


default_prm_template = "./solkz_Q2_Q1.template.prm"
prm_dir = "prm_files"

directory_template = "::OUTPUT_DIR::"
global_ref_template = "::GLOBAL_REFINEMENT::"
adaptive_ref_template = "::ADAPTIVE_REFINEMENT::"

output_dir_fmt = "output_g{global_ref}_a{adaptive_ref}"


coarse_glob = "interpolated_values_{ref_level}_c_*.dat"
refined_glob = "interpolated_values_{ref_level}_r_*.dat"


def generate_prm(source_template, global_ref, adaptive_ref=0):
    # Generate filenames and ensure directory structure
    os.makedirs(prm_dir, exist_ok=True)
    stem = os.path.basename(source_template).split('.')[0]
    postfix = "_g{global_ref}_a{adaptive_ref}".format(
        global_ref=global_ref, adaptive_ref=adaptive_ref)
    output_directory = "output"+postfix
    dest_name = stem+postfix+".prm"
    dest_filename = os.path.join(prm_dir, dest_name)
    with open(source_template, "r") as template_file:
        with open(dest_filename, "w") as dest_file:
            for line in template_file:
                dest_file.write(
                    line.replace(directory_template, output_directory)
                    .replace(global_ref_template, str(global_ref))
                    .replace(adaptive_ref_template, str(adaptive_ref))
                )
    return ((global_ref, adaptive_ref), dest_filename, output_directory)


def generate_comparison_json(run_directories, json_dest):
    comparisions = []
    for a, b in zip(run_directories[:-1], run_directories[1:]):
        a_ref = a[0][0]
        a_dir = a[2]
        b_ref = b[0][0]
        b_dir = b[2]
        assert a_ref + 1 == b_ref
        entry = {
            "label": "$2^{"+str(b_ref)+"}",
            "A": [os.path.join(a_dir, refined_glob.format(ref_level=a_ref))],
            "B": [os.path.join(b_dir, coarse_glob.format(ref_level=b_ref))],
        }
        comparisions.append(entry)
    with open(json_dest, "w") as f:
        json.dump(comparisions, f)


if __name__ == "__main__":
    comparision_json = "./convergence_comparison.json"
    files = []
    for r in range(3, 7):
        files.append(generate_prm(default_prm_template, r))
    generate_comparison_json(files, comparision_json)
    for l in files:
        subprocess.call([aspect_binary, l[1]])
    subprocess.call([postprocess_binary, "convergence", comparision_json,
                     "-o", "comparison.csv"])
