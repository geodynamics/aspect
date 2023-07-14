# Plugins and additional files

In case you need other files (like shared libraries) to run your cookbook, you
have to create a new folder in the [cookbooks](https://github.com/geodynamics/aspect/tree/main/cookbooks) directory that is named after
your cookbook (with words divided by underscores). {ref}`sec:extending:write-a-plugin`
explains how to add a `CMakeLists.txt` file to that directory so that your
plugin can be compiled easily (see the bullet point starting with "Put
the `my_plugin.cc` file into a directory of your choice..."). Note that
after you have copied and renamed the [doc/plugin-CMakeLists.txt](https://github.com/geodynamics/aspect/blob/main/doc/plugin-CMakeLists.txt)
file, you have to modify it in the following way: in the command
`SET(TARGET "my_plugin")`, replace `"my_plugin"` by the name you want your
shared library to have (usually the name of the cookbook), and in
`ADD_LIBRARY(${TARGET} SHARED source_1.cc source_2.cc)`, replace
`source_1.cc source_2.cc` by the name of your .cc file.
