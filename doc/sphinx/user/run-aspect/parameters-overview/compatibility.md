# Compatibility of input files with newer ASPECT versions

#### Compatibility of input files with newer ASPECT versions

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

<div class="center">

</div>
