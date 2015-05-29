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
