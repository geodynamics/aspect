(part:dev_manual:chap:developing_for_the_GWB:sec:creating_new_plugin_systems)=
Creating new plugin system
==========================

Suppose you want to create a new feature or tectonic unit in your world besides the available features, i.e., subducting plate, fault, mantle layer etc. You can do this by creating a new plugin system that includes the possible distributions of `composition`, `grains`, and `temperature` within that feature. This is quite an advanced change to GWB, and it would be useful to first contact the developers to discuss the implementation before doing it yourself.

Below are the following broad steps you would need to follow to implement a new plugin system:

1. Write a class of your feature, say `myfeature`, that is derived from the `interface` class of the GWB features. You can do this by first copying the source (`.cc`) and the header (`.h`) files of an existing feature.  Then, rename all the instances of the class to `myfeature`, modify the name of the feature in the constructor, add the feature in WB_REGISTER_FEATURE, and modify the header guards.
2. Create a folder of your feature inside `source/world_builder/features`. Next, add the subfolders for all the properties available in GWB, i.e., `composition`, `grains`, and `temperature` within `myfeature`.
3. Add `interface` class for each of the properties into the file `interface.cc` within their separate subfolders. You can copy the existing `interface.cc` in the available tectonic features. Do the same for the header files, `interface.h` in `include/world_builder/features`.
4. Change the namespace and header guard in the `interface` class to reflect `myfeature`.
5. Write a class that describes for all the properties how you want their respective distribution. Make sure that this class is derived from the interface class. You can either copy the available implementations of these properties in the existing features or write your own following the steps in {ref}`part:dev_manual:chap:developing_for_the_GWB:sec:creating_new_plugins`.