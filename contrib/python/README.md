Python contributions for ASPECT 
==============================

This folder contains python scripts, functions, and notebooks that ASPECT users can find useful.

Python packages needed:
- numpy
- pandas
- xmltodict
- requests
- matplotlib
- pyvista
- scipy
- jupyter

1- Install Python with anaconda: 
https://www.anaconda.com/distribution/

2- The provided YAML configuration file <environment.yml> can be used to create the <py_aspect> environment and install all the packages listed above:
>> conda env create -f environment.yml

3- Activate the new environment:
>> conda activate py_aspect

Instead of using the pre-configured environment, users can:

    i) Create a new python environment with conda:
    >> conda create --name py_aspect python=3.10
    Users can change the name of the environment (i.e., py_aspect) and the python version (i.e., 3.10).

    ii) Activate conda environment:
    >> conda activate py_aspect

    iii)- Install required packages:
    >> conda install numpy
    >> conda install pandas
    >> conda install xmltodict
    >> conda install requests
    >> conda install matplotlib
    >> conda install scipy
    >> conda install -c conda-forge pyvista
    >> conda install jupyter

To use jupyter notebook:
    1- Create an ipykernel from the <py_aspect> environment:
    >> conda install ipykernel

    2- Install the kernel for the <py_aspect> environment:
    ipython kernel install --user --name=py_aspect

    3- Start the server for jupyter notebook:
    >> jupyter notebook

    OR Start the server for jupyter lab:
    >> jupyter lab

    It will open a page on your default web browser. Double click on the file .ipynb and you are all set!
