(cha:installation)=
# Installation

There are three distinct ways to install ASPECT - compilation from
source, installing a virtual machine, and using a Docker container -
each providing distinct advantages and disadvantages. In this section we
describe all three options and start with a summary of their properties to
guide users to an informed decision about the best option for their purpose.

```{table} Features of the different installation options of ASPECT.
:name: tab:install

|             Feature             | Compile & Install |      Virtual Machine       | Docker Container |
|:-------------------------------:|:-----------------:|:--------------------------:|:----------------:|
|         Speed overhead          |        0%         |            30%             |    0-5%    |
|          Disk overhead          |     0&nbsp;GB     |         1&nbsp;GB          |   200&nbsp;MB    |
|       Knowledge required        |       Much        |        Very Little         |      Little      |
|    Root privileges required     |        No         | No (installed VM software) |    Partially     |
| Embedded in native environment  |        Yes        |             No             |    Partially     |
|          MacOS support          |        Yes        |            Yes             |       Yes        |
|         Windows support         |        No         |            Yes             |       Yes        |
|      Local parallelization      |        Yes        |            Yes             |       Yes        |
| Massively parallel computations |        Yes        |             No             |        No        |
|        Modifying ASPECT         |     Possible      |          Possible          |     Possible     |
|    Configuring dependencies     |     Possible      |             No             |        No        |
```

The available options can be best presented in form of typical use cases:

1.  Virtual Machine (ASPECT beginner and
    tutorial participant): The virtual machine image provides a fully prepared
    user environment that contains installations of
    ASPECT, all required libraries, and visualization
    software on top of a full Linux environment. This way beginning users and
    tutorial participants can work in a unified environment, thus minimizing
    installation time and technical problems. Due to the overhead of
    virtualizing a full operating system this installation typically needs
    more space, and is approximately 30&nbsp;% slower than a native
    installation. Additionally working in a virtual machine
    &lsquo;feels' differently from working in your usual desktop
    environment. The virtual machine can be run on all host operating systems
    that can run a virtualization software like VirtualBox (e.g. Linux, Apple
    MacOS, Microsoft Windows).

2.  Docker Container (advanced user with no need to configure/change the
    underlying libraries, possibly changing parts of
    ASPECT): Docker containers are lightweight
    packages that only encapsulate the minimal dependencies to run an
    application like ASPECT on top of the host
    operating system. They allow easy installation and usage of
    ASPECT in a unified environment, while relying on
    the user's operating system to provide visualization software and
    model input data. When compared to the virtual machine it is simple to
    exchange files between the host operating system and the docker container,
    and it provides the benefit to work in the desktop environment you are
    used to. They have very little overhead in terms of memory and speed
    compared to virtual machines, and allow for reproducible computations. The
    container is set up with a standard ASPECT
    installation, but this can be modified by advanced users (source code
    development within the container is possible).

3.  Compile & Install (advanced users and developers with the need to
    reconfigure underlying libraries or running massively parallel models):
    The most advanced option is to compile and install
    ASPECT from source. This allows maximal control
    over the underlying libraries like Trilinos
    and deal.II, as well as easy modifications
    to ASPECT by recompiling a modified source
    directory. Our installation instructions cover most Linux and MacOS
    operating systems, but incompatibilities on individual systems can always
    occur and make the installation more cumbersome. If you are planning to
    run massively parallel computations on a compute cluster this is likely
    your only option. Since clusters usually have a very individual setup, it
    is always a good idea to ask IT support staff for help when installing
    ASPECT, to avoid hard to reproduce setup
    problems, and performance penalties.






:::{toctree}
docker-container/index.md
virtual-machine/index.md
local-installation/index.md
:::
