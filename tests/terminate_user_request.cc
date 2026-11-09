// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2022-2025 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------

#include <stdlib.h>

// Create a function that is run upon loading the plugin. We delete the
// file that is used as a trigger for terminating
// ASPECT. Otherwise we might have a left over file that triggers termination
// immediately.
int f()
{
  system ("rm -f output-terminate_user_request/restart/last_good_checkpoint.txt");
  return 42;
}

// run this function by initializing a global variable by it
int i = f();
