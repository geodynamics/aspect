#!/bin/bash

# This script generates a docker image for the ASPECT tester.
# It requires a docker installation on the local machine, and the ability to
# communicate with the docker daemon without root user privileges (see the docker
# webpage for an explanation).

docker build -t geodynamics/aspect-tester:latest .

docker build -t geodynamics/aspect-tester:o1 -f Dockerfile.optimized .
