(part:user_manual:chap:concepts:sec:features)=
Features
================

Tectonic features are a core concept in the World Builder. The name "tectonic feature" is often shorted to just "feature" for simplicity. A tectonic feature in the GWB represents a tectonic/geodynamic unit. This can be a continental plate, an oceanic plate, a subducting plate, a fault, etc. 

```{note}
The term plate is used here in a loose definition. A more correct term would probably be continental lithosphere, oceanic plate lithosphere, etc. If you have strong opinion on this, feel free to contribute to [the open issue on GitHub](https://github.com/GeodynamicWorldBuilder/WorldBuilder/issues/210).
```

Conceptually, the GWB tries to stay as close as possible to how a geologist would think about the problem. Ideally, GWB allows a user without any or with little programming and modeling experience to take a map of a region, identify the tectonic units and create a World Builder file with little to no help from a computational geodynamicist. This lets the user, for example, make interesting 3D renderings of their region for visualization purposes or as a starting point for a discussion about how to model it.

In the next few sections, we will go into more detail of how this works. We will start with showing how features can be used to create 3D setups and then go into the three main types of features: {ref}`part:user_manual:chap:concepts:sec:area_features`, {ref}`part:user_manual:chap:concepts:sec:line_features` and {ref}`part:user_manual:chap:concepts:sec:point_features`. 

