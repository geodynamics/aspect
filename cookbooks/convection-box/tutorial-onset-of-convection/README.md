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



## Packages install

The heat flux slider requires the installation of two additional packages.

**Adding Jupyter widgets**

You should have added this in Step 3 above:
> jupyter nbextension enable widgetsnbextension --py --sys-prefix

**Installing tables**
> conda install tables
