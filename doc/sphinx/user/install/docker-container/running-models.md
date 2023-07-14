
# Running ASPECT models

Although it is possible to use the downloaded
ASPECT docker image in a number of different ways, we
recommend the following workflow:

1.  Create your ASPECT input file in a folder
    of your choice (also containing any input data that is required
    by your model) and navigate in a terminal into that directory.

2.  Run the docker image and mount the current directory as a read-only volume
    into the docker container.[^footnote1] This is accomplished by specifying the -v
    flag followed by the absolute path on the host machine, colon, absolute
    path within the docker container, colon, and specifying read-only
    permissions as in the example below.

    Make sure your parameter file specifies a model output directory *other*
    than the input directory, e.g., `/home/dealii/aspect/model_output`. When
    you have started the container run the aspect model inside the container.
    Note that there are two ASPECT versions in the container:
    `aspect` and `aspect-release`. For a discussion of the different
    versions see {ref}`sec:run-aspect:debug-mode`, in
    essence: You should run `aspect` first to check your model for errors,
    then run `aspect-release` for a faster model run.

    To sum up, the steps you will want to execute are:

    ``` ksh
    docker run -it -v "$(pwd):/home/dealii/aspect/model_input:ro" \
      geodynamics/aspect:latest bash
    ```

    Within the container, simply run your model by executing:

    ``` ksh
    aspect model_input/your_input_file.prm
    ```

3.  After the model has finished (or during the model run if you want to check
    intermediate results), copy the model output out of the container into your
    current directory. For this you need to find the name or ID of the docker
    container by running `docker ps -a` in a separate terminal first. Look for
    the most recently started container to identify your current
    ASPECT container.

    Commands that copy the model output to the current directory could be:

    ``` ksh
    docker ps -a # Find the name of the running / recently closed container in the output
    docker cp CONTAINER_NAME:/home/dealii/aspect/model_output .
    ```

4.  The output data is saved inside your container even after the computation
    finishes and even when you stop the container. After you have copied the
    data out of the container, you should therefore delete the container to
    avoid duplication of output data. Even after deleting you will always be
    able to start a new container from the downloaded image following step 2.
    Deleting the finished container can be achieved by the
    `docker container prune` command that removes any container that is not
    longer running.

    :::{note}
    If you own other finished containers that you want to keep, use `docker container
    rm CONTAINER_NAME` to only remove the container named `CONTAINER_NAME`.
    :::

    To remove all finished containers use the following command:

    ``` ksh
    docker container prune
    ```

    Alternatively only remove a particular container:

    ``` ksh
    docker container rm CONTAINER_NAME
    ```

You are all set. Repeat steps 1-4 of this process as necessary when updating
your model parameters.

[^footnote1]: Note that it is possible to mount a directory as writeable
into the container. However, this often causes file permission conflicts
between your host operating system and the container. Therefore, we
recommend this slightly more cumbersome, but also more reliable workflow.
