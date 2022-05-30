(sec:benchmark-run)=
# Running benchmarks that require code

Some of the benchmarks require plugins like custom material models, boundary
conditions, or postprocessors. To not pollute ASPECT with all these
purpose-built plugins, they are kept separate from the more generic plugins in
the normal source tree. Instead, the benchmarks have all the necessary code in
`.cc` files in the benchmark directories. Those are then compiled into a shared
library that will be used by ASPECT if it is referenced in a `.prm`
file. Let's take the SolCx benchmark as an example (see Section {ref}`sec:benchmarks:solcx`).
The directory contains:

-   `solcx.cc` -- the code file containing a material model
    "SolCxMaterial" and a postprocessor "SolCxPostprocessor",

-   `solcx.prm` -- the parameter file referencing these plugins,

-   `CMakeLists.txt` -- a cmake configuration that allows you to
    compile `solcx.cc`.

To run this benchmark you need to follow the general outline of
steps discussed in {ref}`sec:extending:write-a-plugin`. For the current case, this
amounts to the following:

1.   Move into the directory of that particular benchmark:
     ```{code-block} console
     $ cd benchmarks/solcx
     ```

2.   Set up the project:
     ```{code-block} console
     $ cmake .
     ```
     By default, `cmake` will look for the ASPECT binary and other
     information in a number of directories relative to the current one.
     If it is unable to pick up where ASPECT was built and installed, you can
     specify this directory explicitly this using `-D
     Aspect\_DIR=$<$...$>$` as an additional flag to `cmake`, where
     `$<$...$>$` is the path to the build directory.

3.   Build the library:
     ```{code-block} console
     $ make
     ```
     This will generate the file `libsolcx.so`.

Finally, you can run ASPECT with `solcx.prm`:
```{code-block} console
 $ ../../aspect solcx.prm
```
where again you may have to use the appropriate path to get to the ASPECT
executable. You will need to run ASPECT from the current directory because
`solcx.prm` refers to the plugin as `./libsolcx.so`, i.e., in
the current directory.
