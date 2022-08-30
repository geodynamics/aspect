# Running ASPECT models

The internal setup of the virtual machine is similar to the Docker container
discussed above, except that it contains a full-featured desktop environment.
Also note that the user name is `ubuntu`, not `dealii` as in the Docker
container. Again there are multiple ways to use the virtual machine, but we
recommend the following workflow:

1.  Create your ASPECT input file in the shared
    folder and start the virtual machine.

2.  Navigate in a terminal to your model directory.

3.  Run your model using the provided ASPECT
    executable:

    ``` ksh
    ~/aspect/aspect your_input_file.prm
    ```

4.  The model output should automatically appear on your host machine in the
    shared directory.

5.  After you have verified that your model setup is correct, you might want
    to consider recompiling ASPECT in release
    mode to increase the speed of the computation. See {ref}`sec:run-aspect:debug-mode`
    for a discussion of debug and release mode.

6.  Visualize your model output either inside of the virtual machine (ParaView
    and VisIt are pre-installed), or outside on your host system.

You are all set. Repeat steps 1-6 of this process as necessary when updating
your model parameters.
