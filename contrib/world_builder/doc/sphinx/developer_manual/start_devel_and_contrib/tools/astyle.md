(part:dev_manual:chap:start_developing_and_contribute:sec:tools:subsec:astyle)=
Astyle
=====

Astyle (Artistic Style) is a tool that automates the formatting of code. Executing astyle will automatically indent new code to the standard that the Geodynamic WorldBuilder uses, and should be done before each pull request to ensure that the indenter tests do not fail. 

GWB requires the use of astyle v2.04, which can be downloaded from sourceforge [here](https://sourceforge.net/projects/astyle/files/astyle/astyle%202.04/) for mac, linux or windows. 

An example one line command for how to build and install the correct version of astyle on linux is shown below, followed by a more general description on building astyle and linking it to the GWB.

`mkdir astyle && cd astyle && wget 'https://sourceforge.net/projects/astyle/files/astyle/astyle 2.04/astyle_2.04_linux.tar.gz' && tar -zxvf astyle_2.04_linux.tar.gz && cd astyle/build/gcc && make && sudo make install`

After the download is complete, simply unzip the files in your desired directory, and compile the astyle executable. To do this, run `make` in the `$ASTYLE_DIR/build/$COMPILER`. Here `$ASTYLE_DIR` is the path to the unzipped astyle directory, and `$COMPILER` is the directory within `$ASTYLE_DIR/build/` that corresponds to the compiler you will use to generate astyle (gcc, clang, intel, etc.). This should generate the astyle executable in `$ASTYLE_DIR/build/$COMPILER/bin/`. For the GWB to see this executable, you now need to add this to your `$PATH`. For mac, open your `~/.bash_profile` and add the following line:

`export PATH=$ASTYLE_DIR/build/$COMPILER/bin/:$PATH`

Now, if you navigate to `$WORLDBUILDER_SOURCE_DIR/build` and run `make indent`, astyle will run and automatically format all files within `$WORLDBUILDER_SOURCE_DIR`.