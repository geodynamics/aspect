(part:dev_manual:chap:developing_for_the_GWB:sec:creating_new_plugins)=
Creating new plugins
====================

It is possible that the user wants to implement a different distribution of properties (temperature or composition) than the ones available in the GWB. This can be achieved by creating a new plugin with the steps outlined subsequently. We will discuss two ways to write a plugin: one is to base it on an existing plugin, and another is to write the plugin from scratch.

##  Modify an existing plugin

Identify an existing implementation that nearly does what you want and copy the contents into another file, say `myplugin.cc`, in the same directory. It is also recommended that you make a copy of the corresponding header file into `myplugin.h`. This allows other plugins to use the declarations in `myplugin`. Then, rename the class in both the header and the source files to `Myplugin` (the class names are written in the Camel case by GWB convention) and modify the header guard accordingly.

As an example, assume that you want to implement faults in your world that have compositional value that varies according to a hat function. In this case, the closest existing functionality is the `smooth` fault composition that varies following a hyperbolic tangent function. You can then rename the class to say Hat and modify the implementation of the relevant function, i.e., `get_composition()`, such that it now returns a compositional value based on the hat function.

## Write a plugin from scratch

Create your plugin file that contains the declarations and implementations of all the members of your class, which must be derived from the relevant `Interface` class. For example, if you are creating a new temperature computation for a `continental_plate_model`, then you would include the function declarations of `features/continental_plate_models/fault_models/interface.h`. This implies that you need to have the functions `declare_entries`, `parse_entries`, and `get_temperature` and their definitions for your plugin. Similar to above, it is recommended that you split the function declarations into a header file, say `myplugin.h`.

-----
Register the plugin at the bottom of `myplugin.cc` in `WB_REGISTER_*` that instantiates the plugin, documents it, and makes it available to the parameter file handlers. You can do this by using `WB_REGISTER_*`(name of the class, name of the plugin used in the input file). Finally, compile the plugin, using `cmake .` and then `make` in the build directory.
