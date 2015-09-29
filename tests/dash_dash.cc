#include <aspect/simulator.h>

// create a function that is run upon loading the plugin.  as discussed in the
// corresponding .prm file, this function simply calls ASPECT again, and then
// terminates the original instance of ASPECT
int f()
{
  // call ASPECT with "--" and pipe an existing input file into it.
  //
  // notes:
  // - the box_origin.prm file contains no explicit output directory
  //   (nor would that be helpful here), so we also pipe the name
  //   for a temporary output directory into ASPECT
  // - if this directory does not exist, this would trigger
  //   ASPECT's 'directory appears not to exist' warning.
  //   consequently, remove the directory if it existed before
  //   and re-create it as an empty directory
  system ("cd output-dash_dash ; "
          "(cat " ASPECT_SOURCE_DIR "/tests/box_origin.prm "
          " ; "
          " echo 'set Output directory = output.tmp' "
          " ; "
          " rm -rf output.tmp ; mkdir output.tmp "
          ") "
          "| ../../aspect -- ");
  exit (0);
}

// run this function by initializing a global variable by it
int i = f();
