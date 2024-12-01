
# Installing Docker and downloading the ASPECT image

Docker is a lightweight virtualization software that allows to ship
applications with all their dependencies in a simple way. It is outside of the
scope of this manual to explain all possible applications of Docker, and we
refer to the introduction (<https://www.docker.com/what-docker>) and
installation and quickstart guides (<https://www.docker.com/products/docker>)
on the Docker website for more detailed descriptions of how to set up and use
the docker engine. More importantly Docker provides a marketplace for
exchanging prepared docker images (called Docker Hub). After setting up the
docker engine, downloading a precompiled ASPECT
image from Docker Hub is as simple as typing in a terminal:

``` ksh
docker pull geodynamics/aspect
```

The precompiled version of ASPECT downloaded through the command above
will be from a recent commit on the main (development) branch of ASPECT.

A list of Docker containers with precompiled versions of ASPECT from earlier
releases and download instructions can also be found on Docker Hub
((<https://hub.docker.com/r/geodynamics/aspect/tags>)).

For example, the command to download the container for ASPECT v.2.5.0 is
``` ksh
docker pull geodynamics/aspect:v2.5.0
```

Note that the transfer size of the compressed image containing
ASPECT and all its dependencies is a few GB.
When extracted the image requires less than 10 GB of disk space.

Running the image is as simple as typing
``` ksh
docker run -it geodynamics/aspect
```
and then executing the ASPECT binary using
``` ksh
./aspect
```

:::{note}
If you are running Docker under Windows, we strongly recommend running the
image by using the command line (open "Command prompt" or "Terminal") instead
of launching through the Docker graphical interface.
:::
