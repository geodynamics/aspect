# Installing parameter-GUI

#### Installing parameter-GUI

The deal.II parameter-GUI program can be
downloaded at <https://github.com/dealii/parameter_gui>, and is compiled using
the `cmake` program just like ASPECT itself.
The program has no dependencies except for the Qt development libraries that
should be available as packages for most Linux distributions and can also be
obtained for all major operating systems at
<https://www.qt.io/download-open-source/>.

Example steps for installing the parameter-gui could look as follows:

1.  Download the program from <https://github.com/dealii/parameter_gui>.

2.  Prepare a Makefile by running `cmake .` in the source folder.

3.  Compile the program by running `make`.

4.  Make sure to set the environment variable `PARAMETER_GUI_DIR` to the
    directory that contains the parameter-GUI executable (optional). This will
    allow ASPECT to automatically enable the
    GUI during configuration.

##### Installing on macOS

On a mac machine with recent macOS Sierra 10.12.4, Qt development libraries of
version 4.x.x at the libraries' official website
<https://www.qt.io/download-open-source/> may fail to install. Alternatively,
you can install `qt4` through Homebrew (also see instruction here
<https://github.com/cartr/homebrew-qt4>)

    brew tap cartr/qt4
    brew tap-pin cartr/qt4
    brew install qt@4

or install it through Mac Ports (<https://www.macports.org/>)

    sudo port install qt4-mac

Then you can follow the Linux user instructions provided previously to
download and install dealii parameter-GUI. Before running `cmake .`, you may
need to either pass the path of `qt4` and specify the value of variable
`QT_LIBRARIES` to the directory that contains the libraries of `qt4` or add
those information into your `.bash_profile`. For example, for installation
through Mac Ports, you can set the following into your `.bash_profile`

    export PATH="$PATH:/opt/local/libexec/qt4"
    export QT_LIBRARIES="/opt/local/libexec/qt4"
