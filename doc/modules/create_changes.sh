#!/bin/bash

# create changes.h from the bits in changes/ and the headers

if [[ $# != 0 ]]
then
  echo "Usage: ./create_changes.sh"
  exit
fi

output=changes.h

echo "Creating changes.h"

# copy header and information for every bit into changes.h
cp current_changes_header $output

for file in `ls -r changes`; do
  echo " *" >> $output
  sed -e '1s/^/<li> / ;s/^/ * /' changes/$file >> $output
done

cat current_changes_footer >> $output

echo "done."
