```{tags}
category:cookbook
```

(sec:cookbooks:phase-transition-kinetics)=

# Using Non-Equilibrium Thermodynamics to Drive Phase Transformations

_Buchanan Kerswell and Rene Gassmöller contributed this cookbook._

You can find the input parameter file for this model at [https://github.com/geodynamics/aspect/blob/main/cookbooks/phase_transition_kinetics/simple-subduction.prm](https://www.github.com/geodynamics/aspect/blob/main/cookbooks/phase_transition_kinetics/simple-subduction.prm)

## Overview

This cookbook demonstrates how to model phase transitions in a simple subduction model using reaction kinetics driven by non-equilibrium thermodynamics. This [`PhaseTransitionKinetics`](https://github.com/geodynamics/aspect/blob/main/cookbooks/phase_transition_kinetics/phase-transition-kinetics.h) material model differs from existing material models in `ASPECT` that apply the [`PhaseFunction`](https://github.com/geodynamics/aspect/blob/4a0743e738e65c3c8b371b4e8579e304f855ec0d/include/aspect/material_model/utilities.h#L756) class (e.g., see the [Convection in a 2d box with a phase transition](https://aspect-documentation.readthedocs.io/en/latest/user/cookbooks/cookbooks/christensen_yuen_phase_function/doc/christensen_yuen_phase_function.html) cookbook), which assumes phase transitions proceed over a finite depth (or pressure) interval under hydrostatic conditions after {cite:t}`christensen:yuen:1985`. As an alternative approach, the `PhaseTransitionKinetics` material model uses operator splitting to govern phase transitions under non-hydrostatic conditions---in this cookbook we focus on the olivine &rarr; wadsleyite reaction. In this simple subduction example, cold material flows from the top boundary of a 300 $\times$ 200 km 2D box at a fixed temperature and velocity. The top part of the model consists of pure olivine, which becomes metastable as it flows downwards into the wadsleyite stability field. In this cookbook we ignore the reverse wadsleyite &rarr olivine reaction because the mantle flow is unidirectional.

Rather than explicitly prescribing a univariant reaction boundary as in the `PhaseFunction` class, the `PhaseTransitionKinetics` material model formulates the phase transition as a function of the excess Gibbs free energy relative to an adiabatic reference state:

```{math}
:label: eqn:driving-force

\Delta G = \Delta G_0 + \hat{P} \Delta V - \hat{T} \Delta S
```

where $\Delta G$ is excess Gibbs free energy driving the phase transition, $\Delta G_0 = \Delta G_b - \Delta G_a$ is the difference in Gibbs free energy between the reacting phases $a$ and $b$, $\Delta V = V_b - V_a$ and $\Delta S = S_b - S_a$ are the differences in molar volume and entropy of the phases, respectively, $\hat{P} = P - \bar{P}$ is the nonadiabatic ("dynamic") pressure, and $\hat{T} = T - \bar{T}$ is the nonadiabatic temperature. The thermodynamic variables $\Delta G_0$, $\Delta V$, and $\Delta S$ are evaluated along an isentropic adiabat ($\bar{P}$, $\bar{T}$) using the thermodynamic data and equations of state from {cite:t}`stixrude:lithgow-bertelloni:2011`. These thermodynamic data are stored in the table [`ol-wad-driving-force-profile.txt`](https://github.com/geodynamics/aspect/blob/main/cookbooks/phase_transition_kinetics/thermodynamic-driving-force-profile-olivine-wadsleyite.txt), which is referenced using the [`AsciiDataProfile`](https://github.com/geodynamics/aspect/blob/a947bac5a92025ed0f12bf4da70757e78a531b25/include/aspect/structured_data.h#L717) class. Figure {numref}`fig:ascii-data-profile` shows contours of the thermodynamic driving force $\Delta G$ as a function of depth and the degree of local dynamic pressure $\hat{P}$ and heating $\hat{T}$.

```{figure-md} fig:ascii-data-profile
<img src="thermodynamic-driving-force-olivine-wadsleyite.*" style="width:80.0%" />

Thermodynamic driving force for the reaction olivine &rarr; wadsleyite as a function of depth. The pressure- and temperature-dependent terms in Equation {math:numref}`eqn:driving-force` are shown in (a) and (b), respectively, with the other term held fixed. The $\Delta G = 0$ contour indicates thermodynamic equilibrium. The olivine &rarr; wadsleyite reaction begins when $\Delta G$ is negative (purple contours) and no reaction occurs where $\Delta G$ is positive (orange contours). At a fixed depth, positive dynamic pressure $\hat{P}$ tends to drive the olivine &rarr; wadsleyite reaction, while excess local heating (positive $\hat{T}$) tends to hinder it. $\Delta G_0$, $\Delta V$, and $\Delta S$ for olivine and wadsleyite were computed with the mineral physics toolkit `BurnMan` {cite}`myhill:etal:2023` using thermodynamic data and equations of state from {cite:t}`stixrude:lithgow-bertelloni:2011`.
```

The reaction rate is then calculated as:

```{math}
:label: eqn:reaction-rate

\frac{dX}{dt} =
  \begin{cases}
    0, & \text{if } \Delta G > 0 \text{ or } X \geq 1 \\
    Q \, \Delta G \, (1 - X), & \text{if } \Delta G < 0 \text{ and } X < 1
  \end{cases}
```

where $Q$ is a kinetic prefactor constant (units: mol/J/s) and $X$ is the mass fraction of phase $b$ (the product phase). The simple subduction models below show a metastable olivine layer forming with a thickness that varies from $\leq$ 10–15 km for $Q$ = 1e$^{-9}$ mol/J/s (Figure {numref}`fig:simple-subduction-Q9`) to $\leq$ 35–40 km for $Q$ = 1e$^{-10}$ mol/J/s (Figure {numref}`fig:simple-subduction-Q10`). Virtually no metastable olivine layer forms for $Q \leq$ 1e$^{-8}$ because the reaction rate $dX/dt$ is too high relative to the advection timestep in the simulation (i.e., the phases are always in equilibrium).

```{figure-md} fig:simple-subduction-Q9
<img src="simple-subduction-Q9-evolution.*" style="width:80.0%" />

Evolution of a simple subduction model with an olivine &rarr; wadsleyite phase transition and a kinetic prefactor $Q$ of 1e$^{-9}$. Cold material flows from the top boundary at a fixed temperature and velocity (left column) and olivine transforms to wadsleyite (center column). The phase transition is governed by a reaction rate driven by non-equilibrium thermodynamics (Equation {math:numref}`eqn:reaction-rate`). A high-density metastable wadsleyite layer and a low-density metastable olivine layer form above and below the equilibrium reaction boundary, respectively (right column). The morphology and sharpness of the phase transition layer develop dynamically as the reaction rate responds to changes in non-hydrostatic ("dynamic") pressure and local heating.
```

```{figure-md} fig:simple-subduction-Q10
<img src="simple-subduction-Q10-evolution.*" style="width:80.0%" />

Evolution of a simple subduction model with an olivine &rarr; wadsleyite phase transition and a kinetic prefactor $Q$ of 1e$^{-10}$. All other parameters are the same as in Figure {numref}`fig:simple-subduction-Q9`.
```

## Assumptions and limitations

The `PhaseTransitionKinetics` material model assumes that reactions begin at univariant PT boundaries and proceed until the metastable phase is entirely absent (i.e., $\Delta G \leq$ 0 until $X$ = either 0 or 1). This implies that there are no divariant PT fields where the reacting phases (or mineral assemblages) remain in equilibrium by adjusting their volume fractions or bulk compositions. The `PhaseTransitionKinetics` material model also does not account for latent heat of reaction, which would affect the reaction rate through the $\hat{T}$ term in Equation {math:numref}`eqn:driving-force` {cite}`{e.g.,}schmeling:etal:1999`. Thus, the `PhaseTransitionKinetics` material model is not intended to be used with other material models that use lookup tables (e.g., the [`Steinberger`](https://aspect-documentation.readthedocs.io/en/latest/user/cookbooks/cookbooks/steinberger/doc/steinberger.html) material model) or the `PhaseFunction` class (e.g., the [`Latent heat`](https://aspect-documentation.readthedocs.io/en/latest/user/cookbooks/cookbooks/latent-heat/doc/latent-heat.html#latent-heat-benchmark) material model). Rather, the utility of the `PhaseTransitionKinetics` material model is to explore feedbacks between mantle dynamics and reaction kinetics, the morphology and PT-dependence of metastable mantle boundary layers, and their implications for and seismic structure and mantle flow.

The `PhaseTransitionKinetics` material model also ignores the dependence of $\alpha$, $\beta$, and $\text{Cp}$ on the reaction rate (i.e., the "metamorphic" terms in equations 7, 9, and 10 of {cite:t}`stixrude:lithgow-bertelloni:2022`). Instead, it computes $\alpha$, $\beta$, and $\text{Cp}$ as the arithmetic average using the mass fractions of olivine and wadsleyite (i.e., the "isomorphic" terms in equations 7, 9, and 10 of {cite:t}`stixrude:lithgow-bertelloni:2022`). Our model therefore underestimates the effective bulk $\alpha$, $\beta$, and $\text{Cp}$.

Moreover, Equations {math:numref}`eqn:driving-force` and {math:numref}`eqn:reaction-rate` define a simple linear model that relates macro-scale thermodynamics to reaction rates. One could also consider adding other terms to the right-hand-side of Equation {math:numref}`eqn:reaction-rate` that scale micro-scale reaction mechanisms to macro-scale reaction rates, such as Avrami-type equations {cite}`avrami:1939,devaux:etal:2000,guest:etal:2003,guest:etal:2004`. These alternative formulations could introduce strong feedbacks between mantle flow and reaction rates, especially if the reaction rate equation is stress-dependent {cite}`liu:etal:1998,wheeler:2018`. Testing different reaction rate models is as simple as swapping out Equation {math:numref}`eqn:reaction-rate` in `PhaseTransitionKinetics` for another expression, thereby extending the functionality of the `PhaseTransitionKinetics` to include other models of reaction kinetics.

## The input parameter file

The [`simple-subduction.prm`](https://github.com/geodynamics/aspect/blob/main/cookbooks/phase_transition_kinetics/simple-subduction.prm) file sets up the simple subduction model shown in Figures {numref}`fig:simple-subduction-Q9` and {numref}`fig:simple-subduction-Q10`. From top to bottom, these subsections include:

A set of global parameters for the solver and other settings

```
# This cookbook illustrates the effect of a thermodynamic driving force on
# the olivine --> wadsleyite phase transition. It features a simplified
# subducting slab, induced by cold material entering the model with a
# prescribed temperature. Due to pressure and a phase transition the
# material compresses and the velocity changes to conserve mass.

set Dimension                              = 2
set Use years instead of seconds           = true
set End time                               = 75e6
set Output directory                       = results/simple_subduction_Q9
set Adiabatic surface temperature          = 1706
set Surface pressure                       = 1e10
set Nonlinear solver scheme                = iterated Advection and Stokes
set Nonlinear solver tolerance             = 1e-4
set Use operator splitting                 = true
set Max nonlinear iterations               = 100
set CFL number                             = 0.7
set Resume computation                     = auto
set Additional shared libraries            = $ASPECT_SOURCE_DIR/cookbooks/phase_transition_kinetics/libphase-transition-kinetics.so

subsection Formulation
  set Mass conservation    = projected density field
  set Temperature equation = real density
end

subsection Solver parameters
  subsection Stokes solver parameters
    set Linear solver tolerance             = 1e-6
    set GMRES solver restart length         = 200
    set Number of cheap Stokes solver steps = 5000
  end

  subsection Operator splitting parameters
    set Reaction time step                 = 5e-4
  end
end

subsection Discretization
  set Use locally conservative discretization      = true
  set Use discontinuous composition discretization = true
  set Use discontinuous temperature discretization = true
end
```

The 2D box geometry and initial mesh resolution

```
subsection Geometry model
  set Model name = box

  subsection Box
    set X extent      = 300e3
    set Y extent      = 200e3
    set X repetitions = 3
    set Y repetitions = 2
  end
end

subsection Mesh refinement
  set Initial adaptive refinement        = 0
  set Initial global refinement          = 6
  set Time steps between mesh refinement = 0
end
```

A constant boundary velocity that flows from top to bottom

```
subsection Boundary velocity model
  set Prescribed velocity boundary indicators = top:function, right x:function, left x:function, bottom x:function

  subsection Function
    set Function expression = 0.0; -1.5e-2
    set Variable names      = x, y
  end
end
```

A compositional boundary condition that ensures olivine and wadsleyite flow into the model at appropriate depths

```
subsection Boundary composition model
  set Fixed composition boundary indicators = top, right, bottom, left
  set List of model names = function

  subsection Function
    set Function expression = if(y<110e3, 1, 0); 3300
    set Variable names      = x, y
  end
end
```

A temperature field that models a cool lithospheric slab entering the model at the top boundary with a -500 K temperature difference from the adiabatic temperature

```
subsection Boundary temperature model
  set List of model names                   = initial temperature
  set Fixed temperature boundary indicators = top, left
end

subsection Initial temperature model
  set List of model names = adiabatic, function

  subsection Function
    set Function expression = -500 * exp(-((y-200e3)*(y-200e3)+(x-150e3)*(x-150e3))/(2*35e3*35e3))
    set Variable names      = x, y
  end

  subsection Adiabatic
    subsection Function
      set Function expression = 0; 3500
    end
  end
end

subsection Heating model
  set List of model names =  adiabatic heating, shear heating
end
```

An initial lithostatic pressure at the top boundary equal to approximately 300 km of mantle rock with a constant gravitational acceleration of 10 m/s$^2$

```
subsection Boundary traction model
  set Prescribed traction boundary indicators = right y:initial lithostatic pressure, left y:initial lithostatic pressure, bottom y:initial lithostatic pressure

  subsection Initial lithostatic pressure
    set Representative point = 300e3, 0
  end
end

subsection Gravity model
  set Model name = function

  subsection Function
    set Function expression = 0.0; -10.0
  end
end
```

The adiabatic conditions are computed from the initial composition model, which defines a two-layered mantle with olivine at depths above 110 km (from the base of the model), and wadsleyite below. The composition field `X_field` stores the amount of wadsleyite and is updated by the reaction rate $dX/dt$ with a separate "reaction" timestep that is smaller than the advection timestep (see {ref}`sec:benchmark:operator-splitting`). We use the projected density approximation {cite}`gassmoller:etal:2020` formulation of the continuity equation, which requires a prescribed `density_field`

```
subsection Adiabatic conditions model
  subsection Compute profile
    set Number of points = 10000
  end
end

subsection Initial composition model
  set List of model names = function

  subsection Function
    set Function expression = if(y<110e3, 1, 0); 3300
    set Variable names      = x, y
  end
end

subsection Compositional fields
  set Number of fields            = 2
  set Names of fields             = X_field, density_field
  set Compositional field methods = field, prescribed field
end
```

The `PhaseTransitionKinetic` material model that requires a data file that contains columns for pressure, temperature, $\rho_a$, $\rho_b$, $\alpha_a$, $\alpha_b$, $\text{Cp}_a$, $\text{Cp}_b$, $\Delta G$, $\Delta V$, and $\Delta S$ for a reaction of interest. The data file can optionally store seismic wave velocities and their temperature-derivatives. The thermodynamic data profile could represent a reaction between polymorphic phases as demonstrated in this cookbook, or it could represent a reaction between two fixed mineral assemblages. The kinetic prefactor $Q$ scales the reaction rate according to Equation {math:numref}`eqn:reaction-rate`

```
subsection Material model
  set Model name = phase transition kinetics

  subsection Phase transition kinetics
    set Kinetic prefactor Q        = 1e-9

    set Viscosity                  = 1e21
    set Minimum viscosity          = 1e19
    set Maximum viscosity          = 1e24
    set Thermal viscosity exponent = 0
    set Transition depths          = 110e3
    set Viscosity prefactors       = 1, 1

    set Data directory = $ASPECT_SOURCE_DIR/cookbooks/phase_transition_kinetics/
    set Data file name = thermodynamic-driving-force-profile-olivine-wadsleyite.txt
  end
end
```

Finally, the model output settings define how to write the solution for visualization

```
subsection Postprocess
  set List of postprocessors = visualization, velocity statistics, temperature statistics, pressure statistics, material statistics, depth average

  subsection Visualization
    set Output format                 = vtu
    set Time between graphical output = 1e6
    set Number of grouped files       = 0
    set List of output variables      = material properties, adiabat, nonadiabatic temperature, nonadiabatic pressure, stress, shear stress, strain rate, stress second invariant, named additional outputs

    subsection Material properties
      set List of material properties = density, thermal expansivity, compressibility, specific heat, viscosity
    end

    subsection Principal stress
      set Use deviatoric stress = true
    end
  end

  subsection Depth average
    set Time between graphical output = 1e6
    set Number of zones               = 50
  end
end

subsection Checkpointing
  set Time between checkpoint  = 1800
end
```
