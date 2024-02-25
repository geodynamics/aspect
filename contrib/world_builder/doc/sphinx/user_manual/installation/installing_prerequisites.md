(part:user_manual:chap:installation:sec:installing_prerequisites)=
Installing the prerequisites
============================

When installing the World Builder on a system, make sure CMake is installed.
If you want a Fortran wrapper, then also make sure a Fortran compiler is installed (GFortran is the preferred option, but other Fortran compilers should work as well).
If you want a Python interface, make sure you have Python 3 installed.

There are many ways to install the prerequisites of the World Builder per operating system.
For each system, we show some options which we know work.
If it doesn't work for your system, please let us know through GitHub issues (<https://github.com/GeodynamicWorldBuilder/WorldBuilder/issues>) and/or make a pull request to describe how to fix the issue for your system.

::::::{tab-set}
:::::{tab-item} Linux (Debian based)

1. Run in a terminal `sudo apt install cmake`
2. For a Fortran wrapper, also run in a terminal `sudo apt install gFortran`
3. For a Python wrapper, also run in a terminal `sudo apt install swig python3-setuptools`

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

If you want to have a native installation, the two main options are using MinGW or Visual Studio.
In both cases, it might be possible to install both the Fortran wrapper and the Python wrapper, but we have not gotten it to work on our tester setup.
Currently we know that with MinGW you can create a successful Fortran wrapper, and with Visual Studio you can create a successful Python wrapper.
The problem with Python in MinGW is not entirely clear, but it seems that only Visual Studio compilers are supported.
So it may or may not be able to find the GWB Python module when it is compiled and installed.
The problem with Fortran with Visual Studio is that you need to install a Fortran compiler somehow, so if you know how to do that, please contribute.

::::{tab-set}

:::{tab-item} MinGW
1. If not already installed, install Chocolatey (<https://chocolatey.org>). In a PowerShell, you can install it with the following command (in one line): `Set-ExecutionPolicy Bypass -Scope Process -Force; iex ((New-Object System.Net.WebClient).DownloadString(’https://chocolatey.org/install.ps1’))`
2. Run in a terminal `choco install msys2`
3. Open a mingw64 terminal
4. Run in a mingw64 terminal `pacman –noconfirm -Syu`
5. Run in a mingw64 terminal `pacman -S mingw-w64-x86_64-toolchain`
6. Run in a mingw64 terminal `pacman –noconfirm -S dos2unix`
:::

:::{tab-item} Visual Studio
1. If not already installed, install Chocolatey (<https://chocolatey.org>). In a PowerShell, you can install it with the following command (in one line): `Set-ExecutionPolicy Bypass -Scope Process -Force; iex ((New-Object System.Net.WebClient).DownloadString(’https://chocolatey.org/install.ps1’))`
2. Run in a terminal `choco install cmake`
3. For a Python wrapper, run in a terminal `choco install python`
4. For a Python wrapper, run in a terminal `choco install swig` 
:::

::::
:::::
::::::
