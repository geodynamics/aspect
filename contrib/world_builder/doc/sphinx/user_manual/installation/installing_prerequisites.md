(part:user_manual:chap:installation:sec:installing_prerequisites)=
Installing the prerequisites
============================

When installing the World Builder on a system, make sure CMake is installed.
If you want a Fortran wrapper, then also make sure a Fortran compiler is installed (GFortran is the preferred option, but other Fortran compilers should work as well).
If you want a Python interface, make sure you have Python 3 installed. Some World Builder apps can make use of MPI parallelization, although the main library generally makes no use of it. If you want to use MPI parallelization in the apps make sure that a MPI library *including its development header files* is installed on your system.

There are many ways to install the prerequisites of the World Builder per operating system.
For each system, we show some options which we know work.
If it doesn't work for your system, please let us know through GitHub issues (<https://github.com/GeodynamicWorldBuilder/WorldBuilder/issues>) and/or make a pull request to describe how to fix the issue for your system.

::::::{tab-set}
:::::{tab-item} Linux (Debian based)
1. Ensure that all modules are up-to-date by running `sudo apt update` and `sudo apt upgrade`.
2. Run in a terminal `sudo apt install [module]` to install all the required modules. These include `cmake`, `gcc`, `g++`, and `libopenmpi-dev`. Use `which [module]` to locate their installation path.
3. Set the environment variables using `export CC=gcc; export CXX=g++;`.
4. (Optional) For a Fortran wrapper, run `sudo apt install gFortran`.
5. (Optional) For a Python wrapper, run `sudo apt install swig python3-setuptools`.
6. (Optional) For MPI parallelization, run `sudo apt install libopenmpi-dev`.

:::::

:::::{tab-item} OSX

Choose one of the options. "No Fortran compiler" is recommended for normal use:
::::{tab-set}

:::{tab-item} No Fortran compiler

The world builder can create a Fortran wrapper to be able to link to Fortran programs.
If you don't explicitly need a Fortran wrapper to be built, you can use the following instructions:
1. If not already installed, install Homebrew. Run in a terminal: `/usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"`
2. Run in a terminal `brew install cmake`
3. For a Python wrapper, also run in a terminal `brew install swig`

Note that if you already had a Fortran compiler installed and CMake can find it, the Fortran wrapper will still be built.
:::

:::{tab-item} Before Xcode 10 with Fortran compiler
1. If not already installed, install Homebrew. Run in a terminal: `/usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"`
2. Run in a terminal `brew install cmake`
3. Run in a terminal `brew install gcc || true && brew link –overwrite gcc`
4. Run in a terminal `export LDFLAGS="-I/usr/local/opt/llvm/lib"`
5. Run in a terminal `export CPPFLAGS="-I/usr/local/opt/llvm/include"`
6. Run in a terminal `export CC=/usr/local/opt/llvm/bin/clang`
7. Run in a terminal `export CXX=/usr/local/opt/llvm/bin/clang++`
8. Run in a terminal `export FC=gFortran"`
9. For a Python wrapper, also run in a terminal `brew install swig`
:::

:::{tab-item} From Xcode 10 with Fortran compiler
You will need to check where clang is installed.
It should be installed in a folder called /Applications/X-code.X.X.app/ or something similar, where the X represents a number.
For Xcode 10 for example, use the following:
1. If not already installed, install Homebrew. Run in a terminal: `/usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"`
2. Run in a terminal `brew install cmake`
3. Run in a terminal `brew install gcc || true && brew link –overwrite gcc`
4. Run in a terminal `export LDFLAGS="-L/Applications/Xcode-10.0.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/lib"`
5. Run in a terminal `export CPPFLAGS="-I/Applications/Xcode-10.0.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include"`
6. Run in a terminal `export CC=/Applications/Xcode-10.0.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/clang`
7. Run in a terminal `export CXX=/Applications/Xcode-10.0.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/clang++`
8. Run in a terminal `export FC=gFortran"`
9. For a Python wrapper, also run in a terminal `brew install swig`
:::

::::
:::::

:::::{tab-item} Windows

There are three main ways to install on Windows.
The recommended way is to use Linux subsystems for Windows (see <https://docs.microsoft.com/en-us/windows/wsl/install-win10>).
In this case, start Linux in the Windows terminal and follow the Linux installation description.

If you want to have a native installation we recommend using  Visual Studio to compile the world builder. The only problem that we are aware of is using Fortran with Visual Studio. The problem here is that you need to install a Fortran compiler somehow, so if you know how to do that, please contribute.

::::{tab-set}

:::{tab-item} Visual Studio
1. If not already installed, install Chocolatey (<https://chocolatey.org>). In a PowerShell, you can install it with the following command (in one line): `Set-ExecutionPolicy Bypass -Scope Process -Force; iex ((New-Object System.Net.WebClient).DownloadString(’https://chocolatey.org/install.ps1’))`
2. Run in a terminal `choco install cmake`
3. For a Python wrapper, run in a terminal `choco install python`
4. For a Python wrapper, run in a terminal `choco install swig` 
:::
::::
:::::
::::::
