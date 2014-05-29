#!/bin/bash

# replaces the version information in several places by the given string

if [[ $# != 1 ]]
then
  echo "Usage: ./bump_version.sh <versionstring>"
  exit
fi

VER=$1

echo "changing version info to '$VER':"

echo "VERSION ..."
echo "$VER" >../VERSION

echo "manual.tex ..."
sed -i "s/^.* %VERSION-INFO%\$/$VER %VERSION-INFO%/g" manual/manual.tex 

echo "global.h ..."
sed -i "s|version \(.*\)\\\\n\" //VERSION-INFO|version $VER\\\\n\" //VERSION-INFO|g" ../include/aspect/global.h

git add ../VERSION manual/manual.tex ../include/aspect/global.h
