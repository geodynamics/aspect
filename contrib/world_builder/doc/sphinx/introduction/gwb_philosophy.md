(part:introduction:chap:GWB_philosophy)=
GWB philosophy
==============

Setting up complex 2D models of geodynamic settings has in the past been hard but feasible. With the advance and the widespread availability of ever more powerful computers to the geodynamic community, detailed regional 3D problems are now within computational range. Although successful attempts have been published, such model setups often have one or many of the following issues:

1. Code is not readable (even hard for their developers)
2. Code/initial conditions is/are not modifiable (even hard for their developers)
3. Code is not extendable
4. Code is not portable or reproducible in other codes
5. Code is not shareable which makes everyone reinvent the wheel

**The Geodynamic World Builder (GWB) is designed to solve these problems.**


To solve these problems, the GWB is built with both a code and a user philosophy in mind. The code philosophy is designed to solve the extendibility, portability and shareability issues and the user philosophy is designed to solve the readability and modifiability issues.

## GWB Code Philosophy

The code philosophy is built around the following points:

1. A single text-based input file
2. Code, language and platform independent
   1. Supports **Linux**, **OSX** and **Windows**
   2. Can interface with **C++**, **C**, **FORTRAN** and **Python**
3. Up-to-date manual and code documentation
4. Safe to use in parallel codes
5. Readable and extensible (modular) codes
6. Strict version numbering to ensure reproducible results

Following these points will help to create a clean, portable, extendable code with reproducible results. This is of course not everything needed to reach such results. For example, having integration and unit tests with high code coverage and automatic code indentation are important to keep the GWB in a healthy state.

## GWB User Philosophy

The user philosophy is built around the following points:

1. Tectonic features can be parameterized by lines and area
2. These features implicitly define a volume 
3. To which a model can be assigned describing
   1. temperature
   2. composition (a label for a material)
   3. Crystal Preferred Orientation
   4. etc.
4. Parameterized by a human readable JSON file

The main idea behind these points is to design the GWB so that users can easily create complex parameterized initial conditions for their geodynamic setting. How this works should become more clear when reading the rest of the manual.
