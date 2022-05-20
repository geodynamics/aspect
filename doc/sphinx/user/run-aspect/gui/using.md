# Using ASPECT-GUI

When configuring ASPECT after executing the
above steps, it should automatically pick up the location of the
parameter-GUI, and will create a new script named `aspect-gui` within the
build folder. If this does not happen, it is possible to hand over the
location of the parameter-GUI as a cmake variable during configuration (e.g.
`cmake -D PARAMETER_GUI_EXECUTABLE=path_to_your_executable`).

The aspect-gui script can be executed from any folder either with no argument
or with one argument that contains the path to an existing parameter file. The
script will run ASPECT with the given (or an
empty) parameter file, generate a database of existing input parameters, and
open the parameter-GUI program with this database. The resulting window looks
similar to Fig.&nbsp;`[5]`. If an existing parameter file was given, all
parameter fields are pre-filled with the values set in the file instead of the
default values. In the program's main window you can change parameters
as necessary, and then save the file as a new parameter file. This parameter
file can then be used to start an ASPECT model
as usual. Note that it is possible to prepare and execute the parameter file
with different versions of ASPECT, e.g. if you
prepare parameter files on a local machine, and execute the model on a remote
compute cluster. Note however that if the two
ASPECT versions contain different default values or
parameter names have changed, this can lead to unexpected model behavior or
even unusable parameter files.

```{figure-md} fig:aspect-gui
<img src="aspect-gui.png" title="fig:" id="fig:aspect-gui" style="width:40.0%" alt="The parameter GUI lists all available parameter options, and allows to change and save them into a new parameter file. Input fields know about the type of the variable and will display useful options to change them (e.g. drop-down menus, file dialogs, text fields)." />

The parameter GUI lists all available parameter options, and allows to change and save them into a new parameter file. Input fields know about the type of the variable and will display useful options to change them (e.g. drop-down menus, file dialogs, text fields).
```

.
