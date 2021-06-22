[![Codacy Badge](https://api.codacy.com/project/badge/Grade/a047af8dae6e498b8900d0ccdd2b726b)](https://www.codacy.com/app/MFraters/WorldBuilder?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=GeodynamicWorldBuilder/WorldBuilder&amp;utm_campaign=Badge_Grade)
[![AppVeyor Build status](https://ci.appveyor.com/api/projects/status/8amaw31qwwlo33vs?svg=true)](https://ci.appveyor.com/project/MFraters/worldbuilder)
[![actions](https://github.com/GeodynamicWorldBuilder/WorldBuilder/actions/workflows/test.yml/badge.svg?branch=master)](https://github.com/GeodynamicWorldBuilder/WorldBuilder/actions?query=branch%3Amaster)
[![Documentation](https://codedocs.xyz/GeodynamicWorldBuilder/WorldBuilder.svg)](https://codedocs.xyz/GeodynamicWorldBuilder/WorldBuilder/index.html)
[![Coverage Status](https://coveralls.io/repos/github/GeodynamicWorldBuilder/WorldBuilder/badge.svg?branch=master)](https://coveralls.io/github/GeodynamicWorldBuilder/WorldBuilder?branch=master)
[![codecov](https://codecov.io/gh/GeodynamicWorldBuilder/WorldBuilder/branch/master/graph/badge.svg)](https://codecov.io/gh/GeodynamicWorldBuilder/WorldBuilder)
[![Mastodon Follow](https://img.shields.io/mastodon/follow/106136314313793382?domain=https%3A%2F%2Fsocial.mfraters.net&style=social)](https://social.mfraters.net/@world_builder)

# The Geodynamic World Builder (GWB)
## What is the Geodynamic World Builder?
The Geodynamic World Builder(GWB) is an open source code library intended to set up initial conditions for computational geodynamic models and/or visualize complex 3d teconic setting, in both Cartesian and Spherical geometries. The inputs for the JSON-style parameter file are not mathematical, but rather a structured nested list describing tectonic features, e.g. a continental, an oceanic or a subducting plate. Each of these tectonic features can be assigned a specific temperature profile (e.g. plate model) or composition label (e.g. uniform). For each point in space, the GWB can return the composition and/or temperature. It is written in C++, but can be used in almost any language through its C and Fortran wrappers. Various examples of 2D and 3D subduction settings are presented.

For more information see the the [GWB site](https://geodynamicworldbuilder.github.io/), see the automatically generated extensive [online User Manual](https://gwb.mfraters.net/manual/manual.pdf) or the automatically generated [code documentation](https://codedocs.xyz/GeodynamicWorldBuilder/WorldBuilder/index.html).

## What do I do when I have a question or want to request a feature?
If you have a question about the code and you can not find the answer easily in the documentation or the question is already raised as an [issue](https://github.com/GeodynamicWorldBuilder/WorldBuilder/issues), please let us know by opening an [issue on github](https://github.com/GeodynamicWorldBuilder/WorldBuilder/issues/new). This is also the preferred method for asking for new features for GWB.

## I found a mistake in the documentation or code, what should I do?
Please do not keep it to yourself and let us know. Others will also profit from mistakes being found and fixed. Even if it is just a typo in the documentation feel free to raise an issue or, even better, make a pull request to fix the issue.

## How to cite?
The developers of the Geodynamic World Builder request that you cite the following publication:

Fraters, M., Thieulot, C., van den Berg, A., and Spakman, W.: The Geodynamic World Builder: a solution for complex initial conditions in numerical modelling, Solid Earth, [https://doi.org/10.5194/se-10-1785-2019](https://doi.org/10.5194/se-10-1785-2019), 2019.

And cite the specific version of the software used. Version 0.4.0 can be cited as:

Menno Fraters and others. 2021, June 22. The Geodynamic World Builder v0.4.0. Zenodo. [https://doi.org/10.5281/zenodo.5014808](https://doi.org/10.5281/zenodo.5014808).

## How can I follow the progress of this project
There are multiple ways in which you can follow this project:
 1. Watch the repository on github. This will give you an update of what happening in the repository. This happens automatically and you can set it up to notify you for different kind of events. 
 2. Follow the [world builder Mastodon account](https://social.mfraters.net/@world_builder). This is a manually updated feed which updates when there are new release or major new features merged in the main branch. More general new related to the project may also be posted.
 3. Subscribe the the Mastodon RSS feed: https://social.mfraters.net/@world_builder.rss. This will show exactly the same information as the mastodon account, but you can use any RSS reader.
 4. Visit the [world builder website](https://geodynamicworldbuilder.github.io/). Besides all kind of useful information and links, it also contains a RSS feed viewer the world builder Mastodon account.
