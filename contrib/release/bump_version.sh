#!/bin/bash

# replaces the version information by the given string

if [[ $# != 1 ]]
then
  echo "Usage: ./bump_version.sh <versionstring>"
  exit
fi

VER=$1

echo "changing version info to '$VER':"

echo "VERSION ..."
echo "$VER" >../../VERSION

git add ../../VERSION
