(sec:benchmarks:solubility)=
# 1D test for water solubility

The implementation of the two-phase flow equations in ASPECT does not only describe
the flow of a silicate melt through porous rocks, but can also be used to model the
flow of other fluids, such as water. Water in the Earth's mantle can either be bound
in minerals, or it can be present as a fluid phase that can migrate relative to the
solid rock. The fraction of this free water then represents the porosity of the
material. Therefore, the terms describing melting and freezing of a silicate melt
in the equations would instead describe the release and reabsorption of water
(or other volatiles) into the rock, which are governed by the water (volatile) solubility.

This test shows how a model for water solubility can be implemented in a material
model in ASPECT, and demonstrates the mass of water is conserved as water is released,
migrates and is being reabsorbed. Specifically, in this test material with a water
content of 1\% flows into the model from the bottom. The model has three layers with
different water solubility, specifically, and infinite solubility above 30 km depth
and below 60 km depth, and a zero solubility in between. When the upwelling material
reaches the boundary at 60 km depth, water is being released and can move relative to
the solid as a free fluid phase. Due to its lower density, it moves upwards faster than
the solid, until it reaches 30 km depth where it is reabsorbed into the solid. In
steady state, the water content should be the same in the top and bottom layer, and in
the middle layer it should be lower, proportional to the ratio between solid and melt
velocity.

Specifically, the material properties are chosen in such a way that the free water
moved twice as fast as the solid, so the porosity in the middle layer should be half
of what flows in from the bottom, or 0.5\%. This is illustrated in the visualizations
below (Figure {numref}`fig:solubility`).

```{figure-md} fig:solubility
<img src="solubility.*" alt="Evolution" />

Evolution of the bound water content (top) and the free water (bottom). Initially, all water is bound. Since the middle layer has a zero solubility, water is released and starts migrating upward with respect to the solid, initially leading to a higher bound water content when it reaches the top layer. As the model reaches steady state, the (bound) water content in the top and bottom layer are equal (1\%) and the water content in the middle layer is 0.5\%.
```
Note that this example is specifically for water, but the general concept can be applied
to a fluid of an arbitrary composition.

The benchmark can be run using the parameter files in
[benchmarks/solubility/](https://github.com/geodynamics/aspect/blob/main/benchmarks/solubility).
The material model is implemented in [benchmarks/solubility/plugin/solubility.cc](https://github.com/geodynamics/aspect/blob/main/benchmarks/solubility/plugin/solubility.cc). Consequently, this code needs to be compiled into a shared library before you can run the test.


``` {literalinclude} ../solubility.prm
```
