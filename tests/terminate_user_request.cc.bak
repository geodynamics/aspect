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

#include <stdlib.h>

// Create a function that is run upon loading the plugin. We delete the
// restart.mesh file that is used as a trigger for terminating
// ASPECT. Otherwise we might have a left over file that triggers termination
// immediately.
int f()
{
  system ("rm -f output-terminate_user_request/restart.mesh");
  return 42;
}

// run this function by initializing a global variable by it
int i = f();
