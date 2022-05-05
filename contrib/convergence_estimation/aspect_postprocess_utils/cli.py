import sys
import json
import glob
import argparse

from aspect_postprocess_utils.data import box_point_hash
from aspect_postprocess_utils.readers import AspectLoadASCII
from aspect_postprocess_utils.diff import StateConvergence, AspectStateDiff

# Top level
arg_parser = argparse.ArgumentParser(description="Python postprocessing utilities for Aspect")
subparsers = arg_parser.add_subparsers(dest='subparser_name')

# Convergence args

convergence_args = subparsers.add_parser('convergence',
    description="""Compute the convergence rate from "state quadrature" output"""
)
convergence_args.add_argument(
    'input_file', type=argparse.FileType(mode='r'),
    help="JSON file with input data"
)
convergence_args.add_argument(
    '-o', '--outfile', type=argparse.FileType(mode='w'),
    help="Output file"
)

def cli_convergence(args):
    flist = json.load(args.input_file)
    args.input_file.close()

    convergence = StateConvergence()
    for item in flist:
        assert isinstance(item['label'], str)
        assert isinstance(item['A'], list)
        assert isinstance(item['B'], list)

        label = item['label']

        a_files = [f for g in item['A'] for f in glob.glob(g)]
        state_a = AspectLoadASCII.load_files(a_files)

        b_files = [f for g in item['B'] for f in glob.glob(g)]
        state_b = AspectLoadASCII.load_files(b_files)

        # Obtain consistent order
        state_a.hash_order(box_point_hash)
        state_b.hash_order(box_point_hash)

        # Add diff
        convergence.add_diff(AspectStateDiff(state_a, state_b), label = label)

    if args.outfile is None:
        convergence.csv_convergence(sys.stdout)
    else:
        convergence.csv_convergence(args.outfile)


def base_cli():
    args = arg_parser.parse_args()

    if args.subparser_name == "convergence":
        cli_convergence(args)
