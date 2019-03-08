#!/bin/bash

# Purpose: Update the copyright year of every file based on the last
#          modification recorded in the git logs

# This is based on the copyright script as part of deal.II

if test ! -d source -o ! -d include -o ! -d cookbooks ; then
  echo "*** This script must be run from the top-level directory of ASPECT."
  exit
fi


files="
  $(find . -name CMakeLists.txt)
  $(find . | egrep '\.(cmake|cc|h|in)$')
"


for i in $files ; do
  # get the last year this file was modified from the git log.
  # we don't want to see patches that just updated the copyright
  # year, so output the dates and log messages of the last 3
  # commits, throw away all that mention both the words
  # "update" and "copyright", and take the year of the first
  # message that remains
  #
  # (it should be enough to look at the last 2 messages since
  # ideally no two successive commits should have updated the
  # copyright year. let's err on the safe side and take the last
  # 3 commits.)
  last_year=`git log -n 3 --date=short --format="format:%cd %s" $i | \
             egrep -i -v "update.*copyright|copyright.*update" | \
             head -n 1 | \
             perl -p -e 's/^(\d\d\d\d)-.*/\1/g;'`

  # get the first year this file was modified from the actual
  # file. this may predate the git log if the file was copied
  # from elsewhere
  first_year=`cat $i | egrep 'Copyright \(C\) [0-9]{4}' | \
              perl -p -e "s/.*Copyright \(C\) (\d{4}).*/\1/g;"`

  # print a status message. we really only have to update
  # the copyright year if the first and last year are
  # different
  echo "Processing $i: ${first_year} - ${last_year}"
  if test ! "${first_year}" = "${last_year}" ; then
    perl -pi -e "s/(Copyright \(C\) \d{4})( - \d{4})?(, \d{4}( - \d{4})?)*/\1 - ${last_year}/g;" $i
  fi
done

