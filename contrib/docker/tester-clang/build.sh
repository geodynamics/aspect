#!/bin/bash

# This script generates a docker image for the ASPECT tester.

docker build -t geodynamics/aspect-tester:noble-dealii-master-clang .
