
# Installing Docker and downloading the ASPECT image

#### Installing Docker and downloading the ASPECT image

Docker is a lightweight virtualization software that allows to ship
applications with all their dependencies in a simple way. It is outside of the
scope of this manual to explain all possible applications of Docker, and we
refer to the introduction (<https://www.docker.com/what-docker>) and
installation and quickstart guides (<https://www.docker.com/products/docker>)
on the Docker website for more detailed descriptions of how to set up and use
the docker engine. More importantly Docker provides a marketplace for
exchanging prepared docker images (called Docker Hub). After setting up the
docker engine downloading a precompiled ASPECT
image from Docker Hub is as simple as typing in a terminal:

``` ksh
docker pull geodynamics/aspect
```

Note that the transfer size of the compressed image containing <span
ASPECT and all its dependencies is about 900&nbsp;MB.
When extracted the image requires about 3.2&nbsp;GB of disk space.
