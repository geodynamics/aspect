
# Developing ASPECT within a container


The above given workflow does not include advice on how to modify
ASPECT inside the container. We recommend a slightly
different workflow for advanced users that want to modify parts of
ASPECT. The ASPECT
docker container itself is build on top of a <span
class="smallcaps">deal.II</span> container that contains all dependencies for
compiling ASPECT. Therefore it is possible to
run the deal.II container, mount an ASPECT
source directory from your host system and compile it inside of the container.
An example workflow could look as following (assuming you navigated in a
terminal into the modified ASPECT source
folder):

``` ksh
docker pull tjhei/dealii:v9.2.0-full-v9.2.0-r2-gcc5
docker run -it -v "$(pwd):/home/dealii/aspect:ro" \
  tjhei/dealii:v9.2.0-full-v9.2.0-r2-gcc5 bash
```

Inside of the container you now find a read-only
ASPECT directory that contains your modified source
code. You can compile and run a model inside the container, e.g. in the
following way:

``` ksh
mkdir aspect-build
cd aspect-build
cmake -DCMAKE_BUILD_TYPE=Debug -DDEAL_II_DIR=$HOME/deal.II-install $HOME/aspect
./aspect $HOME/aspect/cookbooks/shell_simple_2d/shell_simple_2d.prm
```

To avoid repeated recompilations of the ASPECT
source folder we recommend to reuse the so prepared container instead of
starting new containers based on the DEAL.II
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
