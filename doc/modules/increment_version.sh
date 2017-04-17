#!/bin/bash

# execute with old and new version number. Example:
# ./increment_version.sh 1.3 1.4.0

OLDVERSION=$1
NEWVERSION=$2

if [[ $# != 2 ]]
then
  echo "Usage: ./increment_version.sh <oldversion> <newversion>"
  exit
fi

output=to-$NEWVERSION.h

if [ -e $output ]
then
  echo "ERROR: File 'to-$NEWVERSION' already exists."
  exit
fi

if [ ! -e to-$OLDVERSION.h ]
then
  echo "ERROR: File 'to-$OLDVERSION' doesn't exist."
  exit
fi


echo "Going from $OLDVERSION to $NEWVERSION"

# remove old changes.h if it exists
if [ -e changes.h ]
then
  rm changes.h
fi

# create new changes.h from bits in changes/
bash create_changes.sh
git rm changes/*

# move changes.h to to-$NEWVERSION
echo "$output ..."
cp changes_header.template $output
sed -i "s/NEWVERSION/$NEWVERSION/g" $output
sed -i "s/OLDVERSION/$OLDVERSION/g" $output
sed -n '/<ol>/,$p' changes.h >> $output

git add $output
rm changes.h

# create new current_changes_header:
echo "modules/current_changes_header ..."
cp current_changes_header.template current_changes_header
sed -i "s/VERSION/$NEWVERSION/g" current_changes_header
git add current_changes_header

echo "done."
