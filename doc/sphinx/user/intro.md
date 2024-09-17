(cha:intro)=
# Introduction

ASPECT - short for Advanced Solver for Planetary Evolution, Convection, and Tectonics - is a code intended to solve the equations that describe thermally driven convection with a focus on doing so in the context of convection in the Earth's mantle. The scope of this software has grown to include other planetary bodies and other processes including tectonics, lithospheric deformation, two-phase flow, and core convection. This is reflected in renaming ASPECT
from "Advanced Solver for Problems in Earth's ConvecTion" to its current name in 2023.

It is developed by computational scientists all over the world based on the following principles:

-   *Usability and extensibility:* Simulating mantle convection is a difficult problem characterized not only by complicated and nonlinear material models but, more generally, by a lack of understanding of which parts of a much more complicated model are really necessary to simulate the defining features of the problem.
To name just a few examples:

    -   Mantle convection is often solved in a spherical shell geometry, but the Earth is not a sphere - its true shape on the longest length scales is dominated by polar oblateness, but deviations from spherical shape relevant to convection patterns may go down to the length scales of mountain belts, mid-ocean ridges or subduction trenches.
    Furthermore, processes outside the mantle like crustal depression during glaciations can change the geometry as well.

    -   Rocks in the mantle flow on long time scales, but on shorter time scales they behave more like a visco-elasto-plastic material as they break and as their crystalline structure heals again.
    The mathematical models discussed in {ref}`cha:methods` can therefore only be approximations.

    -   If pressures are low and temperatures high enough, rocks melt, leading to all sorts of new and interesting behavior.

    This uncertainty in what problem one actually wants to solve requires a code that is easy to extend by users to support the community in determining what are the essential features of convection in the Earth's mantle.
    Achieving this goal also opens up possibilities outside the original scope, such as the simulation of convection in exoplanets or the icy satellites of the gas giant planets in our solar system.

-   *Modern numerical methods:* We build ASPECT on numerical methods that are at the forefront of research in all areas - adaptive mesh refinement, linear and nonlinear solvers, and stabilization of transport-dominated processes.
This implies complexity in our algorithms, but also guarantees highly accurate solutions while remaining efficient in the number of unknowns and with CPU and memory resources.

-   *Parallelism:* Many convection processes of interest are characterized by small features in large domains - for example, mantle plumes of a few tens of kilometers diameter in a mantle almost 3,000 km deep.
Such problems can not be solved on a single computer but require dozens or hundreds of processors to work together.
ASPECT is designed from the start to support this level of parallelism.

-   *Building on others' work:* Building a code that satisfies the above criteria from scratch would likely require several 100,000 lines of code.
This is outside what any one group can achieve on academic time scales.
Fortunately, most of the functionality we need is already available in the form of widely used, actively maintained, and well tested and documented libraries, and we leverage these to make ASPECT a much smaller and easier to understand system.
Specifically, ASPECT builds immediately on top of the deal.II library (see <https://www.dealii.org/>) for everything that has to do with finite elements, geometries, meshes, etc.; and, through deal.II on Trilinos (see <http://trilinos.org/>) for parallel linear algebra and on p4est (see <http://www.p4est.org/>) for parallel mesh handling.

-   *Community:* We believe that a large project like ASPECT can only be successful as a community project.
Every contribution is welcome and we want to help you so we can improve ASPECT together.

Combining all of these aspects into one code makes for an interesting challenge.
We hope to have achieved our goal of providing a useful tool to the geodynamics community and beyond!

:::{note}
ASPECT is a community project. As such, we encourage contributions from the community to improve this code over time.
Natural candidates for such contributions are implementations of new plugins as discussed in {ref}`sec:extending:plugin-types` since they are typically self-contained and do not require much knowledge of the details of the remaining code.
Obviously, however, we also encourage contributions to the core functionality in any form! If you have something that might be of general interest, please contact us.
:::

:::{note}
ASPECT will only solve problems relevant to the community if we get feedback from the community on things that are missing or necessary for what you want to do.
Let us know by personal email to the developers, or open a topic on our forum hosted at <https://community.geodynamics.org/c/aspect>!
:::

## Referencing ASPECT

As with all scientific work, funding agencies have a reasonable expectation that if we ask for continued funding for this work, we need to demonstrate relevance.
In addition, many have contributed to the development of ASPECT and deserve credit for their work.
To this end, we ask that you cite the appropriate references if you publish results that were obtained to some part using ASPECT.
For what exactly to cite and suggestions for acknowledgments, please see **<https://aspect.geodynamics.org/cite.html>**.

Also see {cite:t}`aspectmanual,aspect-doi-v1.5.0,aspect-doi-v2.0.0,aspect-doi-v2.0.1,kronbichler:etal:2012,heister:etal:2017`.

## Acknowledgments

The development of ASPECT has been funded through a variety of grants to the authors.
Most immediately, it has been supported through the Computational Infrastructure in Geodynamics (CIG), initially by the CIG-I grant (National Science Foundation Award No. EAR-0426271, via The California Institute of Technology) and later by the CIG-II and CIG-III grants (National Science Foundation Awards No. EAR-0949446, EAR-1550901, and EAR-2149126 via the University of California Davis).
In addition, the libraries upon which ASPECT builds heavily have been supported through many other grants that are equally gratefully acknowledged.

Please acknowledge CIG as follows:

:::{important}
ASPECT is hosted by the Computational Infrastructure for Geodynamics (CIG) which is supported by the National Science Foundation award EAR-2149126.
:::

The ASPECT community as a whole, and a number of the primary developers in particular, owe great thanks to Louise Kellogg who, when she was the head of CIG, was a strong supporter of the ASPECT project.
Louise loved how collaborative the ASPECT development model was, and how many people contributed.
Louise passed away far too early in 2019, but her support lives on in the spirit of this project.
