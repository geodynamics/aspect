(part:user_manual:chap:installation:sec:other_installation_methods)=
Other installation methods
=======

Choose how you would like to install: as part of a...

::::::{tab-set}

:::::{tab-item} C++ program

For this case, there are two options. Either compile the library and link it, or directly compile the library source files within your project. Compiling the Geodynamic World Builder separately and linking it through the C++ wrapper is the recommended option.

::::{tab-set}
:::{tab-item} Compile the library and link it
First, follow the instructions from {ref}`part:user_manual:chap:installation:sec:stand-alone-install`.
The library can be found in the `build/lib` directory with the name `libWorldBuilder.a`.
Link your project with this file.

The only file you need to include in your code is `world_builder/wrapper_cpp.h`. Initialize a variable of type `wrapper_cpp::WorldBuilderWrapper` with a valid std::string pointing to the World Builder file. Then implement in the correct locations the temperature and composition querying functions. It is also possible to directly use the World class in the World Builder, by including `world_builder/world.h`.

:::

:::{tab-item} Directly compile the library source files within your project
(implemented-in-aspect)=

This is the way it is implemented in ASPECT.
Although it is possible to just copy the library source and head directory files, it is recommended to clone the whole World Builder project into a directory in your project (possibly as a git submodule).
This way the World Builder can be easily updated. Add the World Builder files to cmake/make/compiler.

The only file you need to include in your code is `world_builder/world.h`.
Initialize a variable of type `WorldBuilder::World` with a valid std::string pointing to the world builder file.
Then implement in the correct locations the temperature and composition querying functions.

:::
::::

:::::

:::::{tab-item} C program
First, follow the instructions from {ref}`part:user_manual:chap:installation:sec:stand-alone-install`. The library can be found in the `build/lib/` directory with the name `libWorldBuilder.a`. Link your project with this file.

The only file you need to include in your code is `world_builder/wrapper_c.h`.
Create a void variable which is a pointer to a pointer and set it so NULL (e.g. `void **ptr_ptr_world`), and create a pointer to a c-string (e.g. `char *world_builder_file`).
Pass these variables to the `create_world` function.
This function will create the World Builder world.
Now this pointer can be used to call the temperature and composition querying functions.

```{important}
When done with the World Builder, call the `release_world` function.
This will clean up the memory used by the World Builder.
```

:::::

:::::{tab-item} Fortran program
First, follow the instructions from {ref}`part:user_manual:chap:installation:sec:stand-alone-install`.
The library can be found in the `build/mod/` directory with the name `worldbuilder.mod`.
Link your project with this file.
Include the module into your project. The only thing you need to care for when creating the world is to provide a file name which ends with `//C_NULL_CHAR`.
Then call the `create_world` function with the variable cworld and the file name variable.
The Fortran module takes care of the world pointer internally.
When the World Builder world is created, the temperature and composition functions can be called at will.

```{important}
When done with the World Builder, call the `release_world` function.
This will clean up the memory used by the World Builder.
```

To be more clear, we show here an example Fortran program using the Geodynamic World Builder.

```{code-block} console
---
caption: Example of how to interface the World Builder with Fortran code
---
program test
use WorldBuilder

IMPLICIT NONE

  ! Declare the types which will be needed.
  REAL*8 :: temperature,x=120e3,y=500e3,z=0,depth=0,gravity = 10
  INTEGER :: composition_number = 3
  REAL*8 :: composition
  character(len=256) :: file_name = "path/to/world_builder_file"//C_NULL_CHAR

  ! Show how to call the functions.
  CALL create_world(cworld, file_name)

  write(*, *) '2d temperature:'
  CALL temperature_2d(cworld,x,z,depth,gravity,temperature)
  write(*, *) 'temperature in Fortran = ', temperature

  write(*, *) '3d temperature:'
  CALL temperature_3d(cworld,x,y,z,depth,gravity,temperature)
  write(*, *) 'temperature in Fortran = ', temperature

    write(*, *) '2d composition:'
  CALL composition_2d(cworld,x,z,depth,composition_number,composition)
  write(*, *) 'composition in Fortran = ', composition

  write(*, *) '3d composition:'
  CALL composition_3d(cworld,x,y,z,depth,composition_number,composition)
  write(*, *) 'composition in Fortran = ', composition

  CALL release_world(cworld)
END program
```

A more extensive example for how to link Fortran code with the World Builder can be found in the example directory.
This includes instructions on how to compile the example.

:::::

:::::{tab-item} Python program
You only need to include the module called gwb.
To be more clear, we show here a python example using the Geodynamic World Builder.
```{code-block} python
---
caption: Example how to interface the World Builder with Python
---
from gwb import WorldBuilderWrapper

filename = "../../tests/data/continental_plate.wb"
world_builder = WorldBuilderWrapper(filename);

print("2d temperature:")
print("Temp. = ", world_builder.temperature_2d(120.0e3,500.0e3,0,10));
print("3d temperature:")
print("Temp.  = ", world_builder.temperature_3d(120.0e3,500.0e3,0,0,10));
print("2d composition:")
print("Comp. = ", world_builder.composition_2d(120.0e3,500.0e3,0,3));
print("3d composition:")
print("Comp. = ", world_builder.composition_3d(120.0e3,500.0e3,.0e3,0e3,3));
```

:::::

:::::{tab-item} ASPECT
The Geodynamic World Builder is installed by default in ASPECT starting with the release 2.2.0. If you want to use a newer version of the Geodynamic World Builder in ASPECT, you can point ASPECT to the location of the World Builder directory.

```{important}
The Geodynamic World Builder is implemented in ASPECT version 2.1.0 as a submodule using a method similar to what is described under the tab [`cpp program` &#8594; `Directly compile the library source files within your project`](implemented-in-aspect), where ASPECT can be directly compiled with the World Builder source files.
The only thing you have to do is make sure that the submodule is actually loaded.

If the World Builder submodule is initialized, ASPECT's cmake configuration will automatically find it and use it. When cloning the ASPECT repository, add the `-recursive` flag (e.g., `git clone -recursive git@github.com:geodynamics/aspect.git`) to automatically initialize the git submodule.
If you already cloned ASPECT, use the command `git submodule update -init -recursive`.
```

When ASPECT has been successfully compiled with the World Builder, set the ASPECT global input parameter `World builder file` to the World Builder file location.
:::::
::::::
