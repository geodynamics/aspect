# ASPECT

```{admonition} Community Project
:class: information

ASPECT is a community software project. A list of ASPECT developers and contributors is available {ref}`here <sec:authors>`. Contributions to software or documentation by every user are welcome and encouraged. See [here](https://github.com/geodynamics/aspect/blob/main/CONTRIBUTING.md) for how to contribute.
```

```{image} _static/images/aspect_logo.png
:alt: ASPECT Logo
:width: 80%
:align: center
```

## About ASPECT

 ASPECT is a code to simulate problems in thermal convection. Its primary focus is on the simulation of processes in the Earth's mantle, but its design is more general than that. The primary aims developing ASPECT are:

* **Usability and extensibility**: Simulating mantle convection is a difficult problem characterized not only by complicated and nonlinear material models but, more generally, by a lack of understanding which parts of a much more complicated model are really necessary to simulate the defining features of the problem. This uncertainty requires a code that is easy to extend by users to support the community in determining what the essential features of convection in the Earth's mantle are.
* **Modern numerical methods**: We build ASPECT on numerical methods that are at the forefront of research in all areas -- adaptive mesh refinement, linear and nonlinear solvers, stabilization of transport-dominated processes. This implies complexity in our algorithms, but also guarantees highly accurate solutions while remaining efficient in the number of unknowns and with CPU and memory resources.
* **Parallelism**: Many convection processes of interest are characterized by small features in large domains -- for example, mantle plumes of a few tens of kilometers diameter in a mantle almost 3,000 km deep. Such problems require hundreds or thousands of processors to work together. ASPECT is designed from the start to support this level of parallelism.
* **Building on others' work**: Building a code that satisfies above criteria from scratch would likely require several 100,000 lines of code. This is outside what any one group can achieve on academic time scales. Fortunately, most of the functionality we need is already available in the form of widely used, actively maintained, and well tested and documented libraries. Thus, ASPECT builds immediately on top of the deal.II library for everything that has to do with finite elements, geometries, meshes, etc.; and, through deal.II on Trilinos for parallel linear algebra and on p4est for parallel mesh handling.
* **Community**: We believe that a large project like ASPECT can only be successful as a community project. Every contribution is welcome and we want to help you so we can improve ASPECT together.

ASPECT is published under the GNU GPL v2 or newer license.

## Table of Contents
```{toctree}
---
maxdepth: 1
---
user/index.md
parameters/index.md
user/developer_documentation.md
user/authors.md
references.md
```
