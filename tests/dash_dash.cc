/*
  Copyright (C) 2022 by the authors of the ASPECT code.

  This file is part of ASPECT.

  ASPECT is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2, or (at your option)
  any later version.

  ASPECT is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with ASPECT; see the file LICENSE.  If not see
  <http://www.gnu.org/licenses/>.
*/

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
