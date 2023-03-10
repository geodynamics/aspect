(sec:extending:write-a-plugin)=
# How to write a plugin

Before discussing what each kind of plugin actually has to implement (see the
next subsection), let us briefly go over what you actually have to do when
implementing a new plugin. Essentially, the following steps are all you need
to do:

-   Create a file, say `my_plugin.cc` that contains the declaration of the
    class you want to implement. This class must be derived from one of the
    `Interface` classes we will discuss below. The file also needs to contain
    the implementation of all member functions of your class.

    As discussed above, it is possible &ndash; but not necessary &ndash; to
    split this file into two: a header file, say `my_plugin.h`, and the
    `my_plugin.cc` file (or, if you prefer, into multiple source files). We do
    this for all the existing plugins in ASPECT
    so that the documentation of these plugins shows up in the
    doxygen-generated documentation. However, for your own plugins, there is
    typically no need for this split. The only occasion where this would be
    useful is if some plugin actually makes use of a different plugin (e.g.,
    the implementation of a gravity model of your own may want to query some
    specifics of a geometry model you also implemented); in that case the
    *using* plugin needs to be able to see the declaration of the class of the
    *used* plugin, and for this you will need to put the declaration of the
    latter into a header file.

-   At the bottom of the `my_plugin.cc` file, put a statement that
    instantiates the plugin, documents it, and makes it available to the
    parameter file handlers by registering it. This is always done using one
    of the `ASPECT_REGISTER_*` macros that will be discussed in the next
    subsections; take a look at how they are used in the existing plugins in
    the ASPECT source files.

-   You need to compile the file. There are two ways by which this can be
    achieved:

    -   Put the `my_plugin.cc` into one of the
        ASPECT source directories and call `cmake .`
        followed by `make` to ensure that it actually gets compiled. This
        approach has the advantage that you do not need to worry much about
        how the file actually gets compiled. On the other hand, every time you
        modify the file, calling `make` requires not only compiling this one
        file, but also linking the ASPECT binary.
        Furthermore, when you upgrade from one version of
        ASPECT to another, you need to remember to
        copy the `my_plugin.cc` file.

    -   Put the `my_plugin.cc` file into a directory of your choice and
        compile it into a shared library, which will be loaded by ASPECT.
        While one could compile the `.cc` file with a command like

             # NOTE: Do not do this, but use the cmake command below!
             g++ -I/path/to/aspect/headers -I/path/to/deal.II/headers \
                 -fPIC -shared my_plugin.cc -o my_plugin.so

        on Linux, but the command may be different on other systems.
        Having to remember all of these
        pieces is a hassle, and a much easier way is in fact to set up a
        small CMake project for this. To this end, simply copy the file
        [doc/plugin-CMakeLists.txt](https://github.com/geodynamics/aspect/blob/main/doc/plugin-CMakeLists.txt)
        to the directory where you have your
        plugin source files and rename it to `CMakeLists.txt`.

        You can then just run the commands

            cmake -DAspect_DIR=/path/to/aspect/build/ .
            make

        and it should compile your plugin files into a shared library
        `libmy_plugin.so`. A concrete example of this process is discussed in
        {ref}`sec:benchmark-run`. Of course, you may want to
        choose different names for the source files `source_1.cc`, `source_2.cc`
        or the name of the plugin `my_plugin` inside the `CMakeLists.txt`.

        In essence, what the few lines inside the `CMakeLists.txt` do, is that
        they find an
        ASPECT installation (i.e., the directory where
        you configured and compiled it as discussed in
        {ref}`cha:installation`) in either the directory
        explicitly specified in the `Aspect_DIR` variable passed to `cmake`, the
        shell environment variable `ASPECT_DIR`, or just one directory up. It then
        sets up compiler paths and similar, and the following lines simply define
        the name of a plugin, list the source files for it, and define everything
        that's necessary to compile them into a shared library. Calling
        `make` on the command line then simply compiles everything.

        Now you
        only need to tell ASPECT to load this
        shared library at startup so that the plugin becomes available at run
        time and can be selected from the input parameter file. This is done
        using the `Additional shared libraries` parameter in the input file,
        see {ref}`parameters:Additional_20shared_20libraries`:

            set Additional shared libraries = ./libmy_plugin.so

        This approach to build your own shared library has the upside that you can
        keep all files that define new plugins in your own directories where
        you also run the simulations, also making it easier to keep around
        your plugins as you upgrade your ASPECT
        installation. On the other hand, compiling the file into a shared
        library is a bit more that you need to do yourself. Nevertheless, this
        is the preferred approach.

        Recently, ASPECT learned to build debug and optimized builds from a single
        build directory as described in {ref}`sec:run-aspect:debug-mode`. Since then,
        plugins will be compiled into `libmy_plugin.debug.so`, `libmy_plugin.release.so`,
        or both depending on how ASPECT was configured. You do not need to modify the
        `Additional shared libraries` to point to one or the other library, because ASPECT will automatically load the corresponding file if you just specify
        `./libmy_plugin.so`.


:::{note}
Complex projects built on ASPECT often require plugins of more than just one kind. For
example, they may have plugins for the geometry, the material model, and for postprocessing.
In such cases, you can either define multiple shared libraries by repeating the calls to `PROJECT`,
`ADD_LIBRARY` and `ASPECT_SETUP_PLUGIN` for each shared library in your `CMakeLists.txt` file
above, or you can just compile all of your source files into a single shared library. In the latter
case, you only need to list a single library in your input file, but each plugin will still be selectable
in the various sections of your input file as long as each of your classes has a corresponding
`ASPECT_REGISTER_*` statement somewhere in the file where you have its definition. An even
simpler approach is to just put everything into a single file – there is no requirement that different
plugins are in separate files, though this is often convenient from a code organization point of
view.
:::

:::{note}
If you choose to compile your plugins into a shared library yourself, you will need to
recompile them every time you upgrade your ASPECT installation since we do not guarantee
that the ASPECT application binary interface (ABI) will remain stable, even if it may not be
necessary to actually change anything in the *implementation* of your plugin
:::
