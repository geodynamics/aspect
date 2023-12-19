Author: Lorraine J. Hwang, Ian Rose, Juliane Dannberg and the ASPECT development community

# Introduction to ASPECT

This notebook is based on tutorials by J. Dannberg that provide a basic introduction to ASPECT.
This notebook demonstrates the onset of convection and the Nusselt-Rayleigh number relationship.

To run, copy the contents of this directory to your workspace.


## Running using ASPECT Jupyter Notebooks tool

The current version is verified to run within the ASPECT Jupyter Notebooks tool which can be launched from the CIG website:

https://geodynamics.org/resources/aspectnotebook


## Running on your desktop

If you run in your local JupyterLab environment, note the following dependencies:

 * Python 3.8.5
 * IPython 7.19.0
 * Jupyterlab 3.2.1
 * Jupyter widget extension 1.0.0 (see Appendix B)
 * matplotlib 3.3.2
 * numpy 1.19.2
 * tables
 * ipympl
 * scipy
 * glob
 
 
## Running using Docker
The current version of the docker container supports an older version of ASPECT.  The model_input directory contains two files with `_old` in its name which will work with this container.

Use at your own risk.

Here are the old instructions

A. Install Docker

B. Download the docker image
> docker pull tjhei/aspect-jupyter

1. Start Docker
2. In a terminal window, navigate to directory where your model input is located:
>	cd my_directory/ModelInput
3. Load the Jupyter extension for widgets:
>jupyter nbextension enable widgetsnbextension --py --sys-prefix
4. Spin up Timo’s container with Jupyter notebooks and map your current working directory as input: 
>	docker run -it -v "$(pwd):/home/dealii/aspect/model_input:ro" -d -p 8888:8888 --name tmpnb-aspect-jupyter tjhei/aspect-jupyter start-notebook.sh --NotebookApp.token=' ' 
5. In a browser, type:
>	http://localhost:8888/tree?
6. Upload this notebook and the image files.

Remember you are now running in a Docker container and NOT your desktop.

**Common problems:**

	unexpected error address already in use.
	change your local host address 8887:8888
	change this to 8887 in step 5

You have to remove (or rename) that container to be able to reuse that name.
Easiest solution (somewhat heavy handed) stop your container, prune, and check.

> docker stop CONTAINER_ID

> docker container prune 

> docker ps -a
 
## Packages install
 
The heat flux slider requires the installation of two additional packages.

**Adding Jupyter widgets**

You should have added this in Step 3 above:
> jupyter nbextension enable widgetsnbextension --py --sys-prefix

**Installing tables**
> conda install tables
 
 