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
        file, but also link ASPECT.
        Furthermore, when you upgrade from one version of
        ASPECT to another, you need to remember to
        copy the `my_plugin.cc` file.

    -   Put the `my_plugin.cc` file into a directory of your choice and
        compile it into a shared library yourself. This may be as easy as
        calling

             # NOTE: do not do this, but use the cmake command below!
             g++ -I/path/to/aspect/headers -I/path/to/deal.II/headers \
                 -fPIC -shared my_plugin.cc -o my_plugin.so

        on Linux, but the command may be different on other systems. Now you
        only need to tell ASPECT to load this
        shared library at startup so that the plugin becomes available at run
        time and can be selected from the input parameter file. This is done
        using the `Additional shared libraries` parameter in the input file,
        see {ref}`sec:3.1`. This approach has the upside that you can
        keep all files that define new plugins in your own directories where
        you also run the simulations, also making it easier to keep around
        your plugins as you upgrade your ASPECT
        installation. On the other hand, compiling the file into a shared
        library is a bit more that you need to do yourself. Nevertheless, this
        is the preferred approach.

        In practice, the compiler line above can become tedious because it
        includes paths to the ASPECT and
        deal.II header files, but possibly also other
        things such as Trilinos headers, etc. Having to remember all of these
        pieces is a hassle, and a much easier way is in fact to set up a
        mini-CMake project for this. To this end, simply copy the file
        [doc/plugin-CMakeLists.txt] to the directory where you have your
        plugin source files and rename it to `CMakeLists.txt`.

    You can then just run the commands

         cmake -DAspect_DIR=/path/to/aspect/build/ .
         make

    and it should compile your plugin files into a shared library
    `my_plugin.so`. A concrete example of this process is discussed in
    {ref}`sec:benchmark-run`. Of course, you may want to
    choose different names for the source files `source_1.cc`, `source_2.cc`
    or the name of the plugin `my_plugin`.

    In essence, what these few lines do is that they find an
    ASPECT installation (i.e., the directory where
    you configured and compiled it, which may be the same directory as where
    you keep your sources, or a different one, as discussed in
    {ref}`sec:installation`) in either the directory
    explicitly specified in the `Aspect_DIR` variable passed to `cmake`, the
    shell environment variable `ASPECT_DIR`, or just one directory up. It then
    sets up compiler paths and similar, and the following lines simply define
    the name of a plugin, list the source files for it, and define everything
    that's necessary to compile them into a shared library. Calling
    `make` on the command line then simply compiles everything.
