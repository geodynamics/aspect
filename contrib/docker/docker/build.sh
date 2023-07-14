#!/bin/bash

# This script generates a docker image from the latest aspect development version
# for the linux/amd64 and linux/arm64 architectures.
# It requires a docker installation on the local machine, and the ability to
# communicate with the docker daemon without root user privileges (see the docker
# webpage for an explanation).
# Note: This container is build from the developer version on Github, it does not use
# the local ASPECT folder. Therefore local changes are not included in the container.

# This builds the image for amd64 and arm64 platforms and pushed to dockerhub
#docker buildx build --no-cache -o type=registry --platform linux/amd64,linux/arm64 -t geodynamics/aspect .

# This builds the image locally for the current platform
docker build --no-cache -t geodynamics/aspect .
