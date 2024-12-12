# Release Tasklist

## Leading up to a release
- Send out an email about problems or outstanding patches
- Go through the list of TODOs in the source code and see what can be done
- Go through the list of issues marked as bugs (https://github.com/geodynamics/aspect/issues?q=is%3Aissue+is%3Aopen+label%3Abug) and see which have to be fixed
- Go through the list of open pull requests and decide which ones have to go into the release, postpone all others
- Check that the used deal.II version for the Docker container in [contrib/docker/docker/Dockerfile](https://github.com/geodynamics/aspect/blob/main/contrib/docker/docker/Dockerfile) and in the manual is appropriate for the release
- Check that [README.md](https://github.com/geodynamics/aspect/blob/main/README.md) and https://aspect.geodynamics.org/ is up-to-date
and the links are working
- Run (and be patient), if any cookbooks/benchmarks fail, find a fix or open an issue as broken:

  ```
  cd benchmarks && make -f check.mk BUILD=$BUILDDIR -j4
  cd cookbooks && make -f check.mk BUILD=$BUILDDIR -j4
  ```

- Find and fix doxygen errors:

  ```
  git checkout -b pre-release-tasks
  find include -name "*h" -print | xargs -n 1 $DEALSRCDIR/contrib/utilities/checkdoxygen.py
  git commit -a -m "doxygen fixes"
  ```

- Build doxygen and manual, check for missing labels and for warnings, fix if possible:

  ```
  cd doc
  make aspect.tag
  cd sphinx
  make html
  cd ../..
  ```

- Check and fix doxygen documentation. some of the changes of the script will destroy intentional indentation. Go through the list of changes manually and decide which ones to include.
  ```
  find . -name "*.h" -not -wholename "*/doc/modules/*" -not -wholename "*/contrib/world_builder/*" -print | while read file;do $DEALSRCDIR/contrib/utilities/wrapcomments.py $file >temp;mv temp $file;done
  git add -p
  git checkout .
  ```

- Fix formatting, copyright years:

  ```
  ./contrib/utilities/indent
  ./contrib/release/update_copyright.sh
  git commit -a -m "doxygen formatting, update copyright years"
  ```

- Make sure all CI workflows on the main branch pass: https://github.com/geodynamics/aspect/actions?query=branch%3Amain

- Create a pull request with the pre release tasks

## Create a release pull-request
- determine new version roughly following semantic versioning: http://semver.org/
  - format is X.Y.Z for a release, X.Y.Z-pre for the dev version or X.Y.Z-rcW for release candidates
  - backwards incompatible changes require incrementing X, adding features incrementing Y

  setup version numbers:

  ```
  export OLDVER=3.0.0
  export VER=3.1.0
  export VERSHORT=3.1
  export NEXTVER=3.2.0-pre
  # Make sure DEALSRCDIR is set to correct path
  export DEALSRCDIR=$DEAL_II_DIR
  ```

- create branch for main PR to update changes.h in doc/modules:

  ```
  git checkout -b post-release-$VER
  cd doc/modules; rm -f changes/*~; ./increment_version.sh $OLDVER $VER;cd ../..
  cd contrib/release;./bump_version.sh $NEXTVER; cd ..
  git commit -m "release task: update version and changes.h"
  ```

- create a branch, bump version (note, make sure the PR above is included):

  ```
  git checkout post-release-$VER && \
  git checkout -b aspect-$VERSHORT && \
  cd contrib/release && ./bump_version.sh $VER-rc0 && cd ../.. && \
  git commit -m "release task: update version info"
  ```

- compile aspect, make sure you have a symlink in the main directory for the next step
  - make sure the WorldBuilder is using the included version

- update parameters and documentation:

  ```
  cd doc && ./update_parameters.sh && cd sphinx && make html && cd ../.. && \
  git add doc/sphinx/parameters && \
  git commit -m "release task: update manual"
  ```

- Tag a release candidate (RC):

  ```
  export TAG=$VER-rc1
  cd contrib/release && ./bump_version.sh $TAG && cd ../.. && \
  git commit -m "version $TAG" && \
  git tag -s v$TAG -m "version $TAG"
  ```

- Tag the release:

  ```
  export TAG=$VER
  cd contrib/release && ./bump_version.sh $TAG && cd ../.. && \
  git commit -m "version $TAG" && \
  git tag -s v$TAG -m "version $TAG"
  ```

- create a tar file:
  ```
  export PREFIX=aspect-$TAG && rm -rf $PREFIX.tar.gz && \
  git archive --format=tar.gz --prefix=$PREFIX/ HEAD >$PREFIX.tar.gz
  ```

- build pdf doc (temporary by building html one page and print to pdf until
  we fix the sphinx pdf):
  ```
  cd doc/sphinx && make singlehtml && cd ../..
  firefox ./doc/sphinx/_build/singlehtml/index.html
  # print to pdf
  ```

- final testing by extracting tarball, compiling, and running:

  ```
  tar xf $PREFIX.tar.gz
  cd $PREFIX
  docker run --rm -it -v `pwd`:/home/dealii/aspect \
        tjhei/dealii:v9.2.0-full-v9.2.0-r2-gcc5 /bin/bash

  mkdir build; cd build
  cmake -G "Ninja" -D ASPECT_RUN_ALL_TESTS=ON ~/aspect
  ninja
  ctest -j 8 -V
  ```

- make public (branch and tag):

  ```
  git push upstream aspect-$VERSHORT
  git push upstream v$TAG
  ```

- sign:

  ```
  gpg --detach-sign --armor aspect-$TAG.tar.gz
  gpg --detach-sign --armor aspect-manual-$TAG.pdf
  sha1sum aspect-$TAG.tar.gz aspect-manual-$TAG.pdf >sha1sum-$TAG.txt
  ```

- create a release on github, upload .tar.gz
- update website (www branch):
  - header.include: add link to changes
  - index.html: add news entry
  - cite.html: change to current version (2x)
- create zenodo release for source code:
  - https://zenodo.org/deposit?page=1&size=20
  - title: ASPECT v2.0.0
  - license: GPL 2
  - check CIG comments: https://github.com/geodynamics/best_practices/blob/master/ZenodoBestPractices.md
  - add to "Computational Infrastructure for Geodynamics" community
  - update Zenodo button on main readme (see badge button on the right of zenodo page)
  - doc/sphinx/references.bib: add new zenodo entry
- add to github release:
    - Zenodo button
    - [![pdf manual](https://img.shields.io/badge/get-PDF-green.svg)](https://doi.org/10.6084/m9.figshare.4865333)
    - [![online manual](https://img.shields.io/badge/online-manual-red)](https://aspect-documentation.readthedocs.io/en/v2.5.0/)
- create figshare DOI for manual (just upload a new version as the same entry)
  - update doc/sphinx/references.bib entry
- update doc/sphinx/references.bib with src and manual doi
- update aspect.geodynamics.org/cite.html and citing.html in www repo:
  - add new version in citing.html, search for "<option"
  - doc/make_cite_html.py:
    - add new version, update doc/zenodo dois
  - run aspect/doc/ python3 make_cite_html.py add to www
- update http://geodynamics.org/cig/software/aspect/:
  - update current release number
  - create entry for the new release
  - update the list of contributors

- update the spack installation package with the latest tarball,
  see https://github.com/spack/spack/pull/13830 for an example:
      spack checksum aspect

- announce on
  - cig-all@geodynamics.org
  - https://community.geodynamics.org/c/aspect
  - dealii@googlegroups.com

## List of prior release notes

Announcement for 3.0.0 (Nov 6, 2024)
-----------------------------------------
We are pleased to announce the release of ASPECT 3.0.0. ASPECT is the Advanced
Solver for Planetary Evolution, Convection, and Tectonics. It uses modern
numerical methods such as adaptive mesh refinement, multigrid solvers, and
a modular software design to provide a fast, flexible, and extensible mantle
convection solver. ASPECT is available from

                   https://aspect.geodynamics.org/

and the release is available from

        https://geodynamics.org/resources/aspect

and

        https://github.com/geodynamics/aspect/releases/tag/v3.0.0

Among others this release includes the following significant changes:

- ASPECT has been renamed from "Advanced Solver for Problems in Earth's
  ConvecTion" to "Advanced Solver for Planetary Evolution, Convection, and
  Tectonics" to reflect that the scope of ASPECT has grown beyond mantle
  convection.
  (Timo Heister on behalf of all maintainers)

- ASPECT now includes version 1.0.0 of the Geodynamic World Builder and no
  longer supports World Builder versions older than 0.5.0.
  (Menno Fraters and other contributors)

- ASPECT can now be coupled to the landscape evolution code FastScape to deform
  the surface through erosion and sediment deposition. Solution variables can
  now also be output on the surface mesh of the model domain.
  (Derek Neuharth, Anne Glerum)

- ASPECT can now compute crystal-preferred orientation of mineral fabrics.
  DREX-like calculations are used to compute anisotropy tensors and
  distributions of mineral orientations.
  (Menno Fraters, Xiaochuan Tian)

- A sea level postprocessor for glacial isostatic adjustment modeling has been
  added. It computes the sea level based on the free surface topography, ocean
  basin, ice thickness, and perturbed gravitational potential.
  (Maaike Weerdesteijn, John Naliboff)

- There is now a new material model that is designed to advect fluids and
  compute fluid release and absorption based on different models for fluid-rock
  interaction. New melt-rock interactions have been added.
  (Daniel Douglas, Juliane Dannberg, Grant Block, John Naliboff)

- ASPECT now requires deal.II 9.5 or newer. ASPECT is also compatible with
  deal.II 9.6, including new features and performance improvements.
  (Rene Gassmoeller, Timo Heister)

- ASPECT now by default builds a debug and an optimized (release) version of
  the executable in the same build directory.
  (Rene Gassmoeller, Wolfgang Bangerth, Timo Heister)

- ASPECT now has a Visual Studio Code extension, which provides syntax
  highlighting and auto-completion for input parameter files. The old
  Parameter GUI has been removed as it was no longer maintained.
  (Zhikui Guo, Timo Heister, Rene Gassmoeller)

- ASPECT now supports compositional fields with different discretizations
  (continuous or discontinuous) and different polynomial degrees in the same
  model. Compositional fields can now be solved using a different list of
  assemblers for each field, effectively allowing to add additional terms to
  each advection equation.
  (Timo Heister, Juliane Dannberg)

- ASPECT now outputs the physical units of quantities into .pvtu files.
  (Wolfgang Bangerth)

- The geometric multigrid (GMG) solver described in Clevenger and Heister,
  2021, has become ASPECT's new default Stokes solver. The previous algebraic
  multigrid (AMG) option is still available.
  (Conrad Clevenger, Jiaqi Zhang, Timo Heister, Rene Gassmoeller)

- ASPECT now utilizes solvers for ordinary differential equations from the
  SUNDIALS ARKODE library for grain-size evolution and other purposes.
  (Juliane Dannberg, Wolfgang Bangerth, Bob Myhill, Rene Gassmoeller)

- All ASPECT plugin classes and plugin systems are now derived from common base
  classes. This unifies class interfaces across plugin systems and allows for
  removal of duplicate code and documentation.
  (Wolfgang Bangerth)

- The particle subsystem has been overhauled. Most particle parameters have
  moved. Multiple particle systems (with different properties) can be active
  in the same model. Existing input files can be updated with the update
  scripts.
  (Rene Gassmoeller, Menno Fraters, Timo Heister)

- 15 new cookbooks and benchmark cases have been added.
  (Many authors, see link below)

- Many deprecated input options and source code functions have been removed.
  Many bugs and inconsistencies have been fixed.
  (Many authors, see link below).

A complete list of all changes and their authors can be found at
  https://aspect.geodynamics.org/doc/doxygen/changes_between_2_85_80_and_3_80_80.html

We are thankful for all feature and model contributions, code reviews,
forum posts, bug reports, and general help provided by members of our
community. Your contributions have helped make ASPECT what is it today.

Wolfgang Bangerth, Juliane Dannberg, Menno Fraters, Rene Gassmoeller,
Anne Glerum, Timo Heister, Bob Myhill, John Naliboff, Cedric Thieulot,
and many other contributors.


Announcement for 2.5.0 (July 8, 2023)
-----------------------------------------
We are pleased to announce the release of ASPECT 2.5.0. ASPECT is the Advanced
Solver for Problems in Earth's ConvecTion. It uses modern numerical methods such
as adaptive mesh refinement, multigrid solvers, and a modular software design to
provide a fast, flexible, and extensible mantle convection solver. ASPECT is
available from

                   https://aspect.geodynamics.org/

and the release is available from

        https://geodynamics.org/resources/aspect

and

        https://github.com/geodynamics/aspect/releases/tag/v2.5.0

Among others this release includes the following significant changes:

- ASPECT now includes version 0.5.0 of the Geodynamic World Builder.
  (Menno Fraters and other contributors)

- ASPECT's manual has been converted from LaTeX to Markdown to be hosted as a
  website at https://aspect-documentation.readthedocs.io.
  (Chris Mills, Mack Gregory, Timo Heister, Wolfgang Bangerth, Rene
  Gassmoeller, and many others)

- New: ASPECT now requires deal.II 9.4 or newer.
  (Rene Gassmoeller, Timo Heister)

- ASPECT now supports a DebugRelease build type that creates a debug build and
  a release build of ASPECT at the same time.  It can be enabled by setting the
  CMake option CMAKE_BUILD_TYPE to DebugRelease or by typing "make debugrelease".
  (Timo Heister)

- ASPECT now has a CMake option ASPECT_INSTALL_EXAMPLES that allows building
  and install all cookbooks and benchmarks. ASPECT now additionally installs
  the data/ directory. Both changes are helpful for installations that are used
  for teaching and tutorials.
  (Rene Gassmoeller)

- Changed: ASPECT now releases the memory used for storing initial conditions
  and the Geodynamic World Builder after model initialization unless an
  owning pointer to these objects is kept. This reduces the memory footprint
  for models initialized from large data files.
  (Wolfgang Bangerth)

- Added: Various helper functions to distinguish phase transitions for
  different compositions and compositional fields of different types.
  (Bob Myhill)

- Added: The 'adiabatic' initial temperature plugin can now use a spatially
  variable top boundary layer thickness read from a data file or specified as a
  function in the input file. Additionally, the boundary layer temperature can
  now also be computed following the plate cooling model instead of the
  half-space cooling model.
  (Daniel Douglas, John Naliboff, Juliane Dannberg, Rene Gassmoeller)

- New: ASPECT now supports tangential velocity boundary conditions with GMG for
  more geometries, such as 2D and 3D chunks.
  (Timo Heister, Haoyuan Li, Jiaqi Zhang)

- New: Phase transitions can now be deactivated outside a given temperature
  range specified by upper and lower temperature limits for each phase
  transition. This allows implementing complex phase diagrams with transitions
  that intersect in pressure-temperature space.
  (Haoyuan Li)

- New: There is now a postprocessor that outputs the total volume of the
  computational domain. This can be helpful for models using mesh deformation.
  (Anne Glerum)

- New: Added a particle property 'grain size' that tracks grain size evolution
  on particles using the 'grain size' material model.
  (Juliane Dannberg, Rene Gassmoeller)

- Fixed: Many bugs, see link below for a complete list.
  (Many authors. Thank you!).

A complete list of all changes and their authors can be found at
  https://aspect.geodynamics.org/doc/doxygen/changes_between_2_84_80_and_2_85_80.html

Wolfgang Bangerth, Juliane Dannberg, Menno Fraters, Rene Gassmoeller,
Anne Glerum, Timo Heister, Bob Myhill, John Naliboff,
and many other contributors.


Announcement for 2.4.0 (July 25, 2022)
-----------------------------------------
We are pleased to announce the release of ASPECT 2.4.0. ASPECT is the Advanced
Solver for Problems in Earth's ConvecTion. It uses modern numerical methods such
as adaptive mesh refinement, multigrid solvers, and a modular software design to
provide a fast, flexible, and extensible mantle convection solver. ASPECT is
available from

                   https://aspect.geodynamics.org/

and the release is available from

        https://geodynamics.org/resources/aspect

and

        https://github.com/geodynamics/aspect/releases/tag/v2.4.0

Among others this release includes the following significant changes:

- New: ASPECT now requires deal.II 9.3.0 or newer, and cmake 3.1.0 or newer.
  (Timo Heister)

- New: The matrix-free GMG Stokes solver now works for problems with
  free-surface boundaries and elasticity.
  (Jiaqi Zhang, Anne Glerum, Timo Heister, John Naliboff)

- New: The matrix-free GMG Stokes preconditioner is now implemented for the
  Newton solver.
  (Timo Heister, Menno Fraters, Jiaqi Zhang)

- New: Visualization postprocessors now record the physical units of the
  quantity they compute, and this information is also output into visualization
  files with a sufficiently new version of deal.II.
  (Wolfgang Bangerth)

- New: Where possible, when using large data tables as input (e.g., for initial
  conditions specified as tables), these data are now stored only once on each
  node in memory areas that is accessible by all MPI processes on that node.
  (Wolfgang Bangerth)

- New: There is now a new material model for melting in the lowermost mantle.
  It can be used to reproduce the results of Dannberg et al. (2021).
  (Juliane Dannberg)

- New: The geoid postprocessor can now handle a deforming mesh, in addition to the
  already existing option from the dynamic topography postprocessor output.
  (Maaike Weerdesteijn, Rene Gassmoeller, Jacky Austermann)

- New: There is now a 'static' option for the temperature field that is set-up
  similarly to the 'static' option for compositional fields. This allows the
  temperature field to be constant over time so you can still advect and build
  up elastic stresses.
  (Rebecca Fildes, Magali Billen)

- Changed: The least squares particle interpolation plugins now provide a bound
  preserving slope limiter that respects local bounds on each cell.
  (Mack Gregory, Gerry Puckett, Rene Gassmoeller)

- New: Add an advection field method that advects a compositional field
  according to Darcy's Law.
  (Daniel Douglas)

- New: The material model 'dynamic_friction' has been integrated into a new
  rheology model friction_models that can be used together with the
  visco_plastic material model.
  (Esther Heckenbach)

- New: ASPECT now has a ThermodynamicTableLookup equation of state plugin,
  which allows material models to read in one or more Perple_X or HeFESTo table
  files.
  (Bob Myhill)

- Changed: The initial composition model called 'ascii data' can now read in 3d
  ascii datasets into a 2d model and slice the dataset in a user controlled
  plane. This allows it to make high-resolution 2d models of problems that use
  observational data (such as seismic tomography models).
  (Juliane Dannberg, Rene Gassmoeller)

- New: Added a new postprocessor which computes the parameter "Mobility"
  following Lourenco et al., 2020.
  (Elodie Kendall, Rene Gassmoeller, Anne Glerum and Bob Myhill)

- Improved: Particle operations have been significantly accelerated, in
  particular in combination with a recent deal.II version (9.4.0 or newer).
  (Rene Gassmoeller)

- New: Add a benchmark for load induced flexure with options for specifying
  sediment and rock material infilling the flexural moat.
  (Daniel Douglas)

- New: ASPECT now has a cookbook that uses the gravity postprocessor to
  compute gravity generated by S40RTS-based mantle density variations.
  (Cedric Thieulot)

- New: ASPECT now has a cookbook that shows how velocities can be prescribed
  at positions specified by an ASCII input file.
  (Bob Myhill)

- New: There is now a cookbook of kinematically driven oceanic subduction in 2D
  with isoviscous materials and without temperature effects. The cookbook model
  setup is based on Quinquis (2014).
  (Anne Glerum)

- New: There is now a cookbook that visualizes the phase diagram from results
  of a model run. This includes examples from the Visco-Plastic and Steinberger
  material model.
  (Haoyuan Li and Magali Billen)

- New: There is now a cookbook that reproduces convection models with a phase
  function from Christensen and Yuen, 1985.
  (Juliane Dannberg)

- Fixed: Many bugs, see link below for a complete list.
  (Many authors. Thank you!).

A complete list of all changes and their authors can be found at
  https://aspect.geodynamics.org/doc/doxygen/changes_between_2_83_80_and_2_84_80.html

Wolfgang Bangerth, Juliane Dannberg, Menno Fraters, Rene Gassmoeller,
Anne Glerum, Timo Heister, Bob Myhill, John Naliboff,
and many other contributors.


Announcement for 2.3.0 (June 30, 2021)
-----------------------------------------
We are pleased to announce the release of ASPECT 2.3.0. ASPECT is the Advanced
Solver for Problems in Earth's ConvecTion. It uses modern numerical methods such
as adaptive mesh refinement, multigrid solvers, and a modular software design to
provide a fast, flexible, and extensible mantle convection solver. ASPECT is
available from

                   https://aspect.geodynamics.org/

and the release is available from

        https://geodynamics.org/cig/software/aspect/

and

        https://github.com/geodynamics/aspect/releases/tag/v2.3.0

Among others this release includes the following significant changes:

- New: ASPECT now requires deal.II 9.2.0 or newer.
  (Timo Heister)

- New: ASPECT has a new, reproducible logo.
  (Rene Gassmoeller, Juliane Dannberg)

- New: Mesh deformation now also works in combination with particles. Instead
  of the end of the timestep, particles are now advected before solving
  the compositional field advection equations. In iterative advection
  schemes, the particle location is restored before each iteration.
  (Anne Glerum, Rene Gassmoeller, Robert Citron)

- New: ASPECT now supports the creation of visualization postprocessors that
  only output data on the surface of a model. An example is the "surface
  stress" visualization postprocessor.
  (Wolfgang Bangerth)

- New: A new class TimeStepping::Manager to control time stepping with a plugin
  architecture has been added. The architecture allows to repeat time steps
  if the time step length changes significantly.
  (Timo Heister)

- New: A mesh refinement plugin that allows to set regions of minimum and
  maximum refinement level between isosurfaces of solution variables.
  (Menno Fraters and Haoyuan Li)

- New: There is a new nullspace removal option 'net surface rotation',
  which removes the net rotation of the surface.
  (Rene Gassmoeller)

- New: Particle advection can now be used in combination with the repetition of
  timesteps. Before each repetition the particles are restored to their previous
  position.
  (Anne Glerum)

- New: There is a new property in the depth average postprocessor that averages
  the mass of a compositional field (rather than its volume).
  (Juliane Dannberg)

- New: The Drucker Prager rheology module now has an option to include a
  plastic damper, which acts to stabilize the plasticity formulation. At
  sufficient resolutions for a given plastic damper viscosity, the plastic
  shear band characteristics will be resolution independent.
  (John Naliboff and Cedric Thieulot)

- New: ASPECT can now compute viscosity values depending on the values
  of phase functions for an arbitrary number of phases.
  (Haoyuan Li, 2020/08/06)

- New: Added calculation for temperature-dependent strain healing in the strain
  dependent rheology module.
  (Erin Heilman)

- New: Added new rheology module, which computes the
  temperature dependent Frank Kamenetskii viscosity approximation.
  (Erin Heilman)

- New: ASPECT now includes a CompositeViscoPlastic rheology module. This
  rheology is an isostress composite of diffusion, dislocation and Peierls
  creep rheologies and optionally includes a damped Drucker-Prager plastic
  element. The rheology module for Peierls creep includes a formulation to
  compute the exact Peierls viscosity, using an internal Newton-Raphson
  iterative scheme.
  (Bob Myhill, John Naliboff and Magali Billen)

- New: There is a new visualization postprocessor 'principal stress',
  which outputs the principal stress values and directions at every point in the
  model.
  (Rene Gassmoeller)

- New: Added the functionality to compute averages in user defined depth layers
  (e.g. lithosphere, asthenosphere, transition zone, lower mantle) to the depth
  average postprocessor and the lateral averaging plugin.
  (Rene Gassmoeller)

- New: The 'spherical shell' geometry model now supports periodic boundary
  conditions in polar angle direction for a 2D quarter shell (90 degree opening
  angle).
  (Kiran Chotalia, Timo Heister, Rene Gassmoeller)

- New: A new particle interpolator based on quadratic least squares has been
  added.
  (Mack Gregory, Gerry Puckett)

- New: There is now a mesh deformation plugin "diffusion" that can be used to
  diffuse surface topography in box geometry models.
  (Anne Glerum)

- Bug fixes to: Steinberger and Calderwood viscosity profile, particle
  generation, viscous strain weakening, incompressible equation of state,
  pressure sign convention, Neumann heat flow boundaries with the Newton
  solver, viscosity on the adiabat for extended Boussinesq approximation
  models, and many more.
  (many authors)

A complete list of all changes and their contributing authors can be found at
  https://aspect.geodynamics.org/doc/doxygen/changes_between_2_82_80_and_2_83_80.html

Wolfgang Bangerth, Juliane Dannberg, Menno Fraters, Rene Gassmoeller,
Anne Glerum, Timo Heister, Bob Myhill, John Naliboff,
and many other contributors.


Announcement for 2.2.0 (June 30, 2020)
-----------------------------------------
We are pleased to announce the release of ASPECT 2.2.0. ASPECT is the Advanced
Solver for Problems in Earth's ConvecTion. It uses modern numerical methods such
as adaptive mesh refinement, multigrid, and a modular software design to
provide a fast, flexible, and extensible mantle convection solver. ASPECT is
available from

                   https://aspect.geodynamics.org/

and the release is available from

        https://geodynamics.org/cig/software/aspect/

and

        https://github.com/geodynamics/aspect/releases/tag/v2.2.0

This release includes the following significant changes:

- New: There is a new matrix-free Stokes solver which uses geometric multigrid.
  This method is significantly faster than the default algebraic multigrid
  preconditioner and uses less memory. Free surface and melt transport are not
  yet implemented.
  (Thomas C. Clevenger, Timo Heister)

- New: There is now a new approximation for the compressible convection
  models that is called 'projected density field'.
  (Rene Gassmoeller, Juliane Dannberg, Timo Heister, Wolfgang Bangerth)

- Changed: The Geodynamic World Builder has been updated to version 0.3.0.
  (Menno Fraters)

- Changed: ASPECT now requires deal.II version 9.0.0 or newer.
  (Timo Heister, Rene Gassmoeller)

- New: There is a new, alternative stabilization method for the advection equation 
  called SUPG.
  (Thomas C. Clevenger, Rene Gassmoeller, Timo Heister, Ryan Grove)

- Changed: The entropy viscosity method for stabilizing the advection equations
  was substantially improved leading to less artificial diffusion in particular
  close to boundaries.
  (Rene Gassmoeller)

- New: The 'visco plastic' material model now has an option to simulate 
  viscoelastic-plastic deformation. The 'viscoelastic plastic' material
  model has been superseded and removed.
  (John Naliboff, Dan Sandiford)

- New: The "Free surface" functionality has been generalized and is now part of
  "Mesh deformation". This change is incompatible to old parameter files that
  used the free surface. 
  (Rene Gassmoeller, Anne Glerum, Derek Neuharth, Marine Lasbleis)

- New benchmarks: entropy equation, viscoelastic cantilever, buoyancy-driven
  viscoelastic plate stress, advection in annulus, slab detachment benchmark,
  several advection benchmarks, rigid shear, polydiapirs, surface loading.
  (Wolfgang Bangerth, Fiona Clerc, Juliane Dannberg, Daniel Douglas, Rene
  Gassmoeller, Timo Heister, Garrett Ito, Harsha Lokavarapu, John Naliboff,
  Elbridge G. Puckett, Cedric Thieulot)

- Incompatibility: The option to use PETSc for linear algebra has been removed
  until further notice.
  (Timo Heister)

- New: If the user has the libdap libraries installed then input data can be
  pulled from the server instead of a local file.
  (Kodi Neumiller, Sarah Stamps, Emmanuel Njinju, James Gallagher)

- New: Implement the "no Advection, single Stokes" and 
  "single Advection, iterated Newton Stokes" solver schemes. 
  (Timo Heister, Anne Glerum)

- New: The chunk geometry model can now incorporate initial
  topography from an ascii data file.
  (Anne Glerum)

- New: The 'depth average' postprocessor now additionally computes the laterally
  averaged density of vertical mass flux for each depth slice in the model.
  (Rene Gassmoeller)

- Changed: The gravity point values postprocessor has been significantly extended.
  (Ludovic Jeanniot, Cedric Thieulot)

- New: There is now a general class
  `MaterialModel::Utilities::PhaseFunction` that can be used to model
  phase transitions using a smooth phase function.
  (Rene Gassmoeller, John Naliboff, Haoyuan Li)

- New: ASPECT now includes a thermodynamically self-consistent compressible
  material model, that implements the Modified Tait equation of state that is
  described in Holland and Powell, 2011.
  (Bob Myhill)

- New: The material models can now outsource the computation of the viscosity
  into a separate rheology model.
  (Rene Gassmoeller)

- New: ASPECT now includes initial temperature and initial composition plugins
  that use ASCII data files to define the initial temperature or composition
  at a series of layer boundaries.
  (Sophie Coulson, Anne Glerum, Bob Myhill)

- New: Extended spherical shell geometry model to include custom mesh schemes.
  (Ludovic Jeanniot, Marie Kajan, Wolfgang Bangerth)

- New: There is a new termination criterion that cancels the model run
  when a steady state average temperature is reached.
  (Rene Gassmoeller, Juliane Dannberg, Eva Bredow)

- Bug fixes to : parallel hdf5 output, chunk geometry model, initial
  topography modules, gplates boundary velocity plugin.
  (many authors)

A complete list of changes and their contributing authors can be found at
  https://aspect.geodynamics.org/doc/doxygen/changes_between_2_81_80_and_2_82_80.html

Wolfgang Bangerth, Juliane Dannberg, Rene Gassmoeller, Timo Heister,
Jacqueline Austermann, Menno Fraters, Anne Glerum, John Naliboff,
and many other contributors.


Announcement for 2.1.0 (April 29, 2019)
-----------------------------------------
We are pleased to announce the release of ASPECT 2.1.0. ASPECT is the Advanced
Solver for Problems in Earth's ConvecTion. It uses modern numerical methods such
as adaptive mesh refinement, multigrid, and a modular software design to
provide a fast, flexible, and extensible mantle convection solver. ASPECT is
available from

                   https://aspect.geodynamics.org/

and the release is available from

        https://geodynamics.org/cig/software/aspect/

and

        https://github.com/geodynamics/aspect/releases/tag/v2.1.0

This release includes the following significant changes:
- New: ASPECT has a new plugin system that allows it to prescribe a fixed
  heat flux (instead of prescribing the temperature) at the model boundaries.
- New: Compositional fields can optionally be advected with the melt velocity.
- New: There is now a visualization postprocessor that outputs the compaction
  length, the characteristic length scale of melt transport.
- New: ASPECT can optionally use the Geodynamic World Builder
  (https://github.com/GeodynamicWorldBuilder/WorldBuilder/) to create complex
  initial conditions for temperature and composition.
- New: ASPECT can now read in a depth-dependent vs to density conversion file, which
  can be used with the included tomography model plugins.
- New: ASPECT can now read in a depth-dependent initial temperature from file.
- New: The 'ascii data' and 'function' boundary velocity plugins now allow
  velocities to be specified along spherical (up, east, north) unit vectors.
- New: Added a visualization plugin that directly outputs the strain rate tensor.
- New: ASPECT can now call PerpleX to calculate material properties, phase
  amounts and compositions on-the-fly. This model is provided as a
  proof-of-concept; more efficient procedures are required for production runs.
- New: ASPECT now outputs a dynamically generated URL based on used features to
  ask people to cite appropriate papers.
- New: ASPECT has two visualization postprocessors which calculate and output
  the grain lag angle and the infinite strain axis (ISA) rotation timescale,
  respectively. These two quantities can be used to calculate the grain
  orientation lag parameter of Kaminski and Ribe (G3, 2002).
- Improved: The artificial diffusion term that is added in the entropy
  viscosity method to the temperature and composition equations is now computed
  as the maximum of the physical diffusion and entropy viscosity instead of the
  sum. This reduces numerical diffusion for the temperature field.
- New: Compositional fields can now be prescribed to a value that is computed
  in the material model as an additional output at every time step.
- Changed: The heat flux through boundary cells is now computed using the
  consistent boundary flux method as described in Gresho, et al. (1987), which
  is much more accurate than the previously used method.
- New: ASPECT can now calculate gravity anomalies in addition to the geoid.
- New: ASPECT now outputs a file named original.prm in the output directory
  with the exact content of the parameter it got started with.
- New: Added basic support for a volume-of-fluid interface tracking advection
  method in 2D incompressible box models. The VoF method is an efficient method
  to track a distinct compositional field without artificial diffusion.
- New: There is now an option to output visualization data as higher order
  polynomials. This is an improvement in accuracy and requires less disk space
  than the 'Interpolate output' option that was available before. However the
  new output can only be read by ParaView version 5.5 and newer and is
  therefore disabled by default.
- New: Several new benchmark cases were added.
- Many other fixes and smaller improvements.

A complete list of changes and their contributing authors can be found at
  https://aspect.geodynamics.org/doc/doxygen/changes_between_2_80_80_and_2_81_80.html

Wolfgang Bangerth, Juliane Dannberg, Rene Gassmoeller, Timo Heister,
Jacqueline Austermann, Menno Fraters, Anne Glerum, John Naliboff,
and many other contributors.


Announcement for 2.0.1 (June 21, 2018)
-----------------------------------------
We are pleased to announce the release of ASPECT 2.0.1. ASPECT is the Advanced
Solver for Problems in Earth's ConvecTion. It uses modern numerical methods such
as adaptive mesh refinement, multigrid, and a modular software design to
provide a fast, flexible, and extensible mantle convection solver. ASPECT is
available from

                   https://aspect.geodynamics.org/

and the release is available from

        https://github.com/geodynamics/aspect/releases/tag/v2.0.1

This release is a bugfix release for 2.0.0 and includes the following fixes:
- Fixed: The 'compositional heating' heating plugin had a parameter
  'Use compositional field for heat production averaging' that was used
  inconsistently with its description. Its first entry did not correspond
  to the background field, but to the first compositional field, and the
  last value was ignored. This is fixed now, the first entry is used for
  the background field, and all following values determine whether to include
  the corresponding compositional fields.
- Fixed: The 'depth dependent' material model did not properly initialize
  the material model it uses as a base model. This caused crashes if the
  base model requires an initialization (such as the 'steinberger' material
  model). This is fixed now by properly initializing the base model.
- Fixed: The advection assembler for DG elements was not thread-safe,
  which led to wrong results or crashes if a discontinuous temperature
  or composition discretization was combined with multithreading.
- Disabled: The particle functionality was not tested when combined with a
  free surface boundary, and this combination is currently not supported. This
  limitation is now made clear by failing for such setups with a descriptive
  error message.

Wolfgang Bangerth, Juliane Dannberg, Rene Gassmoeller, Timo Heister,
Jacqueline Austermann, Menno Fraters, Anne Glerum, John Naliboff,
and many other contributors.


Announcement for 2.0.0 (May 10, 2018)
-----------------------------------------
We are pleased to announce the release of ASPECT 2.0.0. ASPECT is the Advanced
Solver for Problems in Earth's ConvecTion. It uses modern numerical methods such
as adaptive mesh refinement, multigrid, and a modular software design to
provide a fast, flexible, and extensible mantle convection solver. ASPECT is
available from

                   https://aspect.geodynamics.org/

and the release is available from

        https://github.com/geodynamics/aspect/releases/tag/v2.0.0

This release includes the following changes:
- New: Newton solver and defect correction Picard iterations for nonlinear
  problems (for the Stokes system)
- Melt solver: Overhaul leading to improved performance and stability, better
  integration with other plugins
- New: Material model with grain size evolution
- New: Boundary temperature plugin with evolving core-mantle boundary
  temperature based on the heat flux through the core-mantle boundary
- New: ASPECT can now compute the geoid in 3D spherical shell geometry
- New: Operator splitting for reactions between compositional fields
- New: Added a PREM gravity profile
- Improved: Significantly reduced memory consumption in models that use many
  compositional fields
- Improved: A large number of performance improvements for preconditioners,
  assembly, seismic tomography initial conditions, and lateral averaging
- Improved: More flexibility for boundary and initial conditions, different
  plugins can be combined
- Improved: The dynamic topography postprocessor now uses the consistent
  boundary flux method for computing surface stresses, which is significantly
  more accurate
- New: Additional RHS force terms in the Stokes system can be added
- New particle interpolators: nearest neighbor, bilinear least squares,
  harmonic average
- New: Graphical user interface for the creation and modification of input
  parameter files
- Many other fixes and small improvements.
- Rework: Updated parameter and section names to make them more consistent and
  easier to understand. A script for updating parameter and source files is
  provided with the release.

A complete list of changes can be found at
  https://aspect.geodynamics.org/doc/doxygen/changes_between_1_85_80_and_2_80_80.html

Wolfgang Bangerth, Juliane Dannberg, Rene Gassmoeller, Timo Heister,
Jacqueline Austermann, Menno Fraters, Anne Glerum, John Naliboff,
and many other contributors.


Announcement for 1.5.0 (February 28, 2017)
-----------------------------------------
We are pleased to announce the release of ASPECT 1.5.0. ASPECT is the Advanced
Solver for Problems in Earth's ConvecTion. It uses modern numerical methods such
as adaptive mesh refinement, multigrid, and a modular software design to
provide a fast, flexible, and extensible mantle convection solver. ASPECT is
available from

                   https://aspect.dealii.org/

and the release is available from

        https://github.com/geodynamics/aspect/releases/tag/v1.5.0

This release includes the following changes:

- New: Choice between different formulations for the governing equations
  including Boussinesq and anelastic liquid approximation.
- New: Melt transport (two-phase flow).
- Particles: new generators, ghost exchange, performance improvements,
  interpolation to fields.
- New: Nondimensional material model for incompressible (using the Boussinesq
  approximation) and compressible computations (with ALA or TALA) for
  nondimensionalized problems. This can be used for benchmark problems like
  Blankenbach, King, etc..
- New: Optional DG method for temperature/composition.
- Adiabatic conditions: rework, now includes a reference density profile.
- Free surface: overhaul.
- New cookbooks: continental extension, finite strain, BurnMan interface,
  active tracers.
- New benchmarks: TanGurnis, Blankenbach, King.
- New: viscoplastic material model.
- Material model interface cleanup.
- Assembly performance improvements.
- New: memory statistics postprocessor.
- New: initial topography plugins.
- Many other fixes and small improvements.

A complete list of changes can be found at
  https://aspect.dealii.org/doc/doxygen/changes_between_1_84_80_and_1_85_80.html
and the release is also available from
  https://github.com/geodynamics/aspect/releases/tag/v1.5.0

Wolfgang Bangerth, Juliane Dannberg, Rene Gassmoeller, Timo Heister, and many
other contributors.


Announcement for 1.4.0 (May 15, 2016)
-----------------------------------------
We are pleased to announce the release of ASPECT 1.4.0. ASPECT is the Advanced
Solver for Problems in Earth's ConvecTion. It uses modern numerical methods such
as adaptive mesh refinement, multigrid, and a modular software design to
provide a fast, flexible, and extensible mantle convection solver. ASPECT is
available from

                   https://aspect.dealii.org/

and the release is available from

        https://github.com/geodynamics/aspect/releases/tag/v1.4.0

This release includes the following changes:
- Complete overhaul of the particle architecture.
- Modularization of the assembly process.
- Large improvements to serial and parallel performance of particle handling.
- Support for traction boundary conditions.
- Support for time-dependent mesh refinement criteria.
- Flexible support for modifying the underlying FEM variables for the PDE.
- Optional DG support for temperature and compositional fields.
- More robust parallel file I/O.
- Support for tangential mesh velocity boundaries in free surface
  computations.
- Support for anisotropic viscosity.
- Various fixes to free surface computations: checkpointing, better
  stabilization, crash fixes
- New Chunk Geometry model.
- New material property averaging options.
- Complete rewrite of the heating model infrastructure.
- Several new cookbooks.
- Several new postprocessors.
- Support for signals in various locations that allows plugins to
  inspect/manipulate things inside the core application.
- Several new mesh refinement plugins.
- Improved spherical interpolation of data used in the GPlates plugin. 
- Many other fixes and small improvements.

A complete list of changes can be found at
  https://aspect.dealii.org/doc/doxygen/changes_between_1_83_and_1_84_80.html
and the release is also available from
  https://github.com/geodynamics/aspect/releases/tag/v1.4.0

Wolfgang Bangerth, Timo Heister, and many other contributors.





Announcement for 1.3 (May 18, 2015)
-----------------------------------------
We are pleased to announce the release of ASPECT 1.3. ASPECT is the Advanced
Solver for Problems in Earth's ConvecTion. It uses modern numerical methods such
as adaptive mesh refinement, multigrid, and a modular software design to
provide a fast, flexible, and extensible mantle convection solver. ASPECT is
available from

                   https://aspect.dealii.org/

and the release is available from

        https://github.com/geodynamics/aspect/releases/tag/v1.3

This release includes the following changes:
- New: Averaging of material properties between the quadrature points of a
  cell. This greatly increases the stability of solutions in simulations with
  spatially varying coefficients, and also greatly accelerates the solution,
  at times up to a factor of ten.
- Corrections to the entropy stabilization scheme for compositional fields.
- Fixed and extended: Removal of rigid body translations and rotations when
  the simulation has a nullspace.
- New: VTU visualization output can now be grouped into an arbitrary number of
  files per time step.
- Various fixes to the nonlinear solver residual computation.
- New visualization postprocessors that can output the shear stress and full
  stress tensors.
- Fixes to the latent heat formulation.
- New 'ascii data' plugins for boundary and initial conditions.
- New mass flux statistics postprocessor.
- Many other fixes and small improvements.

A complete list of changes can be found at
  https://aspect.dealii.org/doc/doxygen/changes_between_1_82_and_1_83.html
and the release is also available from
  https://github.com/geodynamics/aspect/releases/tag/v1.3

Wolfgang Bangerth, Timo Heister, and many other contributors.




Announcement for 1.2 (January 25rd, 2015)
-----------------------------------------
We are pleased to announce the release of ASPECT 1.2. ASPECT is the Advanced
Solver for Problems in Earth's ConvecTion. It uses modern numerical methods such
as adaptive mesh refinement, multigrid, and a modular software design to
provide a fast, flexible and extensible mantle convection solver. ASPECT is
available from

                   http://aspect.dealii.org/

This release includes the following changes:
- methods to have truly concentric shells without distorted cells
- linear solver improvements which give up to 2x speedup in some cases
- boundary names like "top" instead of numbers are now supported in input files
- new example: free surface computation with a crust as a stagnant lid
- new benchmarks: Davies et al. and Burstedde et al.
- new initial condition: S40RTS perturbation based on shear wave data
- fixes to minimum/maximum refinement plugins
- better error messages when linear solvers fail
- many other fixes and small improvements (direct solver, file output,
  checkpointing, etc.)

A complete list of changes can be found at
  http://aspect.dealii.org/doc/doxygen/changes_between_1_81_and_1_82.html

Wolfgang Bangerth, Timo Heister, and many other contributors.



Announcement for 1.1 (June 1st, 2014)
---------------------
We are pleased to announce the release of ASPECT 1.1. ASPECT is the Advanced
Solver for Problems in Earth ConvecTion. It uses modern numerical methods such
as adaptive mesh refinement, multigrid, and a modular software design to
provide a fast, flexible and extensible mantle convection solver. ASPECT is
available from

                   http://aspect.dealii.org/

This release combines the incredible number of changes made during the
ASPECT hackathon in May 2014, as well as other changes made since the
last release. In particular, it includes the following changes:
- fixes to Stokes solver tolerance (users might need to change settings!)
- free surface computations
- optional direct solver support
- new material models (multi-component, simple compressible, ...)
- finer mesh refinement control using plugins
- updates to dynamic topography postprocessors
- updates to gplates boundary conditions
- compositional fields can now react with each other
- new cookbooks: reaction between fields, free surfaces
- improved PETSc support
- boundary conditions can now be time dependent 
- support for radiogenic heating models

Wolfgang Bangerth, Timo Heister, and many other contributors.


Announcement for 1.0:
---------------------
We are pleased to announce the release of ASPECT 1.0. ASPECT is the Advanced
Solver for Problems in Earth ConvecTion. It uses modern numerical methods such
as adaptive mesh refinement, multigrid, and a modular software design to
provide a fast, flexible and extensible mantle convection solver. ASPECT is
available from

                   http://aspect.dealii.org/

This release contains a large number of changes, incorporating about one year
of work of the main authors and many contributors:
- a lot of new documentation
- new examples (2d/3d shells, GPlates, ...)
- option to remove rotational/translational modes
- big performance improvements
- compositional fields: reactions, boundary conditions
- support for periodic meshes
- output of dynamic topography
- incorporation of latent heat
- experimental PETSc support

Wolfgang Bangerth, Timo Heister, and many other contributors.





Announcement for 0.3:
---------------------
We are pleased to announce the release of ASPECT 0.3. ASPECT is the Advanced
Solver for Problems in Earth ConvecTion. It uses modern numerical methods such
as adaptive mesh refinement, multigrid, and a modular software design to
provide a fast, flexible and extensible mantle convection solver. ASPECT is
available from

                   http://aspect.dealii.org/

This is primarily a bug fix release. It fixes a bug in the computation of adiabatic initial temperature conditions when specifying inner or outer boundary layers (thanks to Juliane Dannberg for finding and fixing it). The doxygen-generated documentation now also automatically cross-links to deal.II classes wherever they are referenced in ASPECT.

Wolfgang Bangerth, Timo Heister.




Announcement for 0.2:
---------------------
We are pleased to announce the release of ASPECT 0.2. ASPECT is the Advanced
Solver for Problems in Earth ConvecTion. It uses modern numerical methods such
as adaptive mesh refinement, multigrid, and a modular software design to
provide a fast, flexible and extensible mantle convection solver. ASPECT is
available from

                    http://aspect.dealii.org/

This release adds a significant number of new features, including:
   - support for active and passive "compositional" fields (courtesy of
     Juliane Dannberg)
   - more flexibility to output only some variables, drastically reducing
     the amount of data written
   - support for user-defined mesh refinement criteria
   - support for GPlates-generated velocity boundary conditions (courtesy
     of Rene Gassmoeller)
   - support for passive tracer particles
   - provision of an "introspection" module as part of the source code
     to make it easier and less error-prone to write extension code
In addition, the manual has been significantly expanded, with many new
cookbooks, and we have fixed a number of bugs. A list of changes is
available here:

  http://aspect.dealii.org/doc/doxygen/changes_between_0_81_and_0_82.html

Wolfgang Bangerth, Timo Heister, and many other contributors.


