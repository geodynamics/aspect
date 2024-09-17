
# Developing ASPECT within a container

The workflow described previously does not include advice on how to modify
ASPECT inside the container. We recommend a slightly
different workflow for advanced users that want to modify parts of
ASPECT. The ASPECT
docker container itself is built on top of a deal.II
container that contains all dependencies for
compiling ASPECT. Therefore it is possible to
run the deal.II container, mount an ASPECT
source directory from your host system and compile it inside of the container.
An example workflow could look as follows (assuming you navigated in a
terminal into the modified ASPECT source
folder):

``` ksh
docker pull geodynamics/aspect:latest
docker run -it -v "$(pwd):/home/dealii/aspect:ro" geodynamics/aspect:latest
```

Inside of the container you now find a read-only
ASPECT directory that contains your modified source
code. You can compile and run a model inside the container, e.g. in the
following way:

``` ksh
mkdir aspect-build
cd aspect-build
cmake -DCMAKE_BUILD_TYPE=Debug $HOME/aspect
make -j4
./aspect $HOME/aspect/cookbooks/shell_simple_2d/shell_simple_2d.prm
```

To avoid repeated recompilations of the ASPECT
source folder we recommend to reuse the so prepared container instead of
starting new containers based on the deal.II
image. This can be achieved by the following commands outside of the
container:

``` ksh
docker ps -a # Find the name of the running / recently closed container in the output
docker restart CONTAINER_NAME
docker attach CONTAINER_NAME
```

For more information on the differences between using images and containers,
and how to attach additional terminals to a running container, we refer to the
docker documentation (e.g.
<https://docs.docker.com/engine/getstarted/step_two/>).
