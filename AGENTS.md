# ASPECT: The Advanced Solver for Planetary Evolution, Convection, and Tectonics.

A parallel, extensible finite element code.

## Overview
- ASPECT is a C++ CMake project
- unit tests are in /unit_tests/ and are implemented using the Catch2 framework; they are run via ``aspect --test`` (also wired as a ctest target)
- tests are in /tests/, can be run using ``ctest``. By default only a small subset runs; enable the full suite with ``make setup_tests``. Also see /doc/sphinx/user/extending/testing/running-tests.md
- instructions for code formatting, use ``./contrib/utilities/indent`` at /CONTRIBUTING.md

## Layout
- ``source/`` and ``include/aspect/`` — paired implementation and headers
- ``tests/`` — integration / regression tests
- ``unit_tests/`` — unit tests
- ``cookbooks/`` — example programs / tutorials
- ``benchmarks/`` — similar, but for benchmarks
- ``doc/sphinx/`` — user documentation
- ``contrib/`` — vendored third-party code; do not casually edit

## Plugins
New physics is usually added as a plugin under matching directories in ``source/`` and ``include/aspect/``. See ``doc/sphinx/user/extending/write-a-plugin.md`` and ``doc/sphinx/user/extending/idea-of-plugins.md``.
