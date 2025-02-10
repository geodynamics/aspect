# Compatibility of plugins

We strive to maintain compatibility for user written plugins with new versions
of ASPECT for as long as possible. However,
occasionally we have to restructure interface classes to improve
ASPECT further. This is in particular true for new
major versions. In order to allow running old plugins with newer
ASPECT versions, we provide scripts that can
automatically update existing plugins to the new syntax. Executing
`contrib/utilities/update_source_files.sh` with one or more plugin files as arguments will
create a backup of the old file (named `old_filename.bak`), and replace the
existing file with a version that should work with the current
ASPECT version. Using this script would look like
this:

``` ksh
bash contrib/utilities/update_source_files.sh cookbooks/finite_strain/finite_strain.cc
```

:::{note}
Not all text replacements are unique, and the structure of plugin files allows for constructs
the script can not properly parse. Thus, it is important that you check your updated plugin file for
errors. That being said, all plugin files in the main ASPECT repository are updated successfully
using this script.
:::
