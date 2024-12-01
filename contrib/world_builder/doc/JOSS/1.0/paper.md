---
title: 'The Geodynamic World Builder: A planetary structure creator for the geosciences'
tags:
  - C++
  - CPP
  - C
  - Fortran
  - Python
  - Geodynamics
  - Earth science
  - Tectonics
  - Seismology
authors:
  - name: Menno R. T. Fraters
    orcid: 0000-0003-0035-7723
    affiliation: "1" # (Multiple affiliations must be quoted)
  - name: Magali I. Billen
    orcid: 0000-0002-7316-1791
    affiliation: "2"
  - name: Rene Gassmöller
    orcid: 0000-0001-7098-8198
    affiliation: "3,4"
  - name: Arushi Saxena
    orcid: 0000-0003-2930-3497
    affiliation: "2,5"
  - name: Timo Heister
    orcid: 0000-0002-8137-3903
    affiliation: "5"
  - name: Haoyuan Li
    orcid: 0000-0003-0676-9884
    affiliation: "2"
  - name: Daniel Douglas
    orcid: 0000-0002-7871-018X
    affiliation: "6"
  - name: Juliane Dannberg
    orcid: 0000-0003-0357-7115
    affiliation: "3,4"
  - name: Wolfgang Bangerth
    orcid: 0000-0003-2311-9402
    affiliation: "7"
  - name: Yijun Wang
    orcid: 0000-0002-7637-3239
    affiliation: "8"
affiliations:
 - name: GFZ German Research Centre for Geosciences, Potsdam, Germany
   index: 1
 - name: UC Davis, USA
   index: 2
 - name: GEOMAR Helmholtz-Zentrum für Ozeanforschung, Kiel, Germany
   index: 3
 - name: Universiy of Florida, USA
   index: 4
 - name: Clemson University, USA
   index: 5
 - name: New Mexico Tech, USA
   index: 6
 - name: Colorado State University, USA
   index: 7
 - name: University of Oslo, Norway
   index: 8
date: 11 March 2024
bibliography: paper.bib
---

# Summary

Many Earth science applications require the discretization, parameterization, and/or visualization of complex geological features in a 3D geometry in global or regional settings. A prime example are geodynamic models that have to make assumptions about the Earth's thermal and chemical structure and the geometry of different features such as plates, subducted slabs, or mantle plumes. This structure is needed in instantaneous models, as model initial conditions, or to test different hypotheses and compare model predictions to observations. Other examples are the creation of an Earth velocity structure for seismic forward modeling and hypothesis-testing, or the visualization of tectonic features in structural geology.

The Geodynamic World Builder (GWB) has been designed to make the creation of complex parameterized models significantly easier. It can also be used to design simple synthetic models and it supports the use of several types of datasets to set up models. Besides setting up initial conditions for geodynamic models, the GWB can also visualize complex 3D geologic, tectonic, and geodynamic settings.

# Statement of need

Today's computational resources, infrastructure, and numerical methods allow for the creation of complex numerical models that closely resemble specific locations on the Earth, using 3D geometries and high resolutions. However, the related increase in complexity has also made setting up these more detailed regional or global models exceedingly difficult, especially in three dimensions. Furthermore, investigating the model dynamics often requires testing different scenarios involving variations in model geometry, thermal, or chemical structure, or other model assumptions. Although studies with such complex models have been published, the practical realization of these model setups often have one or many of the following disadvantages:

1. The configuration is not human-readable.
2. The software is not easily modifiable and extendable.
3. The model setup is not portable to other computing systems or reproducible in other software frameworks.
4. The model setup is not shareable with other users.

These issues lead to a number of problems with the reproducibility and reliability of modeling studies, which threaten to undermine the predictive power and usefulness of modeling results, and highlight the need for an easy, efficient, and robust way to describe model setups. The GWB has been designed to address these challenges, by creating human readable, parameterized, portable, reproducible, and shareable geodynamic model setups. Critically, the GWB comes with its own programs to visualize the constructed model through applications like Paraview. Creating the models requires no programming knowledge. Therefore, the GWB can be easily used to visualize tectonic and geodynamic settings for publications, teaching, and public outreach. 

![A workflow diagram for how a world builder file can be used to create and visualize a geodynamic model.](workflow_diagram.png)

The GWB has been used in several published studies to model global fault patterns, plumes, and plate dynamics [@Saxena_Dannberg_etal_2023; @Gea_Negredo_etal_2023; @Sandiford_Craig_timothy_2023; and @vanderWiel_Hinsbergen_etal_2024]. Other tools to solve this problem have emerged at around the same time as the first GWB release [@Fraters_Thieulot_etal_2019]. Examples include GeomIO [@Bauville_Baumann_2019; @Spang_Baumann_2022], which uses an approach based on vector graphics; Easy (https://easyinit.readthedocs.io/), which uses a more generic function-based approach; UWGeodynamics [@Beucher_Moresi_etal_2019], which is specifically designed for Underworld [@Moresi_Dufour_2002]; and GemPy [@Varga_Schaaf_2019; @Schaaf_Varga_2021], which is designed for structural modeling. The GWB was designed to be a more general planetary structure creator, using the methods shown below.

# Methods

To address the challenges outlined in the previous section, the Geodynamic World Builder implements specific code and world parameterization principles.

## GWB Code Principles
The GWB's software architecture is built around the following principles:

1. A single text-based input file.
2. Code, language and platform independence:
    1. Support for **Linux**, **macOS** and **Windows**;
    2. Official interfaces for **C++**, **C**, **FORTRAN** and **Python**.
3. Safe parallel execution.
4. Readable and extensible software modules.
5. Strict version numbering to ensure reproducible results.

These principles are implemented in an object-oriented C++ code with interfaces to other programming languages. All parts a user might want to modify are implemented as plugin systems using interface classes that decouple individual modules and allow the user to easily extend the code with new features. In addition, the GWB includes an extensive automated test suite with benchmarks, integration, and unit tests with high code coverage, memory checking, automatic code indentation, and a spell checker to keep the GWB in a healthy state.

## GWB World Parameterization Principles

The GWB's world parameterization principles are built around the idea that a complex model region can be split into individual tectonic features. These tectonic features can be parameterized by defining their location and geometry in terms of points, lines, or areas in a map view. For example, a continental plate can be represented as an area on a map, and the GWB user defines this area. A fault is a linear feature on a map, so the user can define the fault trace as a line at the surface. Users can also provide additional information for a feature, such as a spatially variable thickness or dip angle. The GWB then uses these parameters to create the 3D geometry of the feature, defining its volume. Furthermore, users can attach one or many models to those volumes to define additional properties such as thermal or chemical structure. These can be very simple models, such as a uniform temperature distribution; or follow a more complex distribution, such as a half space cooling model, or a McKenzie model [@McKenzie_1970], or a mass conserving slab temperature model [@Billen_Fraters_AGU_2023]. 

All these tectonic features are bundled in a single input file in standard JSON format, which is human readable, writeable, and editable. The main idea behind this design of the GWB is that users can easily create, modify, and visualize complex parameterized geodynamic or tectonic settings.

## Example 
Below we show an example input file for a Cartesian model that contains a single feature, namely a subducting plate.

```json
{
  "version": "1.0",
  "coordinate system": {"model":"cartesian"},
  "features":
  [
    { 
      "model": "subducting plate",  "name": "Slab",  "dip point": [0,0],
      "coordinates": [[1500e3,1000e3],[1600e3,350e3],[1500e3,0]],
      "segments": [{"length":300e3, "thickness": [100e3], "angle": [0,60]}],
      "temperature models": [{"model": "plate model", "plate velocity": 0.02}],
    }
  ]
}
```

A more complicated example (only requiring 85 lines, and can be found [here](https://github.com/GeodynamicWorldBuilder/WorldBuilder/blob/GWB-v1.0.0/doc/sphinx/_static/gwb_input_files/BST_19_spherical_models.wb)) featuring a spherical geometry, a spatially variable subducting plate, continental plate, oceanic plate and plume can be seen in Fig 2.

![A schematic example of what can be built with 85 lines of a GWB input file formatted in the same way as in the example input file shown above. This includes a slab with variable dip and thickness along strike and down dip, subducting under an oceanic plate on the right side of the ridge, as well as a passive continental margin with variable thickness, and a mantle plume beneath the ridge. The temperature profile of the continent is linear, the oceanic plates are defined by a half-space cooling model, the slab temperature is defined by a mass conserving temperature model and the plume adds heat based on a Gaussian around the center.\label{fig:example}](../../sphinx/_static/images/user_manual/basic_starter_tutorial/BST_19.png)


# Acknowledgements

We would like to acknowledge all other contributors to the project, especially Lorraine Hwang, Rebecca Fildes, and John Naliboff for their advice and support for this project throughout the years. We would also like to acknowledge NSF for their funding and support through grants EAR-1620618, OCE-1948902, EAR-0949446, EAR-1550901, EAR-1925677, and EAR-2149126.

M. Fraters also acknowledges the support of the Department of Earth and Planetary Sciences at UC Davis and the Department of Geological Sciences at the University of Florida, where much of the research presented was completed while he was a postdoctoral scholar. 

# References


