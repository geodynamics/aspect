# Compatibility of input files with newer ASPECT versions

We strive to maintain compatibility for options in input files as long as
possible. However, occasionally we have to reorder, rename, or remove options
from parameter files to improve ASPECT further.
This is especially true for new major versions. In order to allow running old
parameter files with newer ASPECT versions we
provide scripts that can automatically update existing parameter files to the
new syntax. Executing `doc/update_prm_files.sh` with one or more parameter
files as arguments will create a backup of the old parameter file (named
`old_filename.bak`), and replace the existing file with a version that should
work with the current ASPECT version. Using
this script would look like this:

``` ksh
bash contrib/utilities/update_prm_files.sh cookbooks/convection_box.prm
```

:::{note}
Not all text replacements are unique, and the structure of input files allows for constructions
the script can not properly parse. Also we can not guarantee to preserve the structure and
position of comments, as it is not always clear to which part of the input file they refer. Thus, it
is important that you check your updated input file for errors. That being said, all input files in
the main ASPECT repository are updated successfully using this script.
:::
