#!/bin/bash

# execute with old and new version number. Example:
# ./increment_version.sh 1.0 1.1

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

# move changes.h to to-$NEWVERSION
echo "$output ..."
cp changes_header.template $output
sed -i "s/NEWVERSION/$NEWVERSION/g" $output
sed -i "s/OLDVERSION/$OLDVERSION/g" $output
sed -n '/<ol>/,$p' changes.h >> $output

git add $output

# create new changes.h:
echo "modules/changes.h ..."
cp current_changes_header.template changes.h
sed -i "s/VERSION/$NEWVERSION/g" changes.h
git add changes.h

echo "done."
