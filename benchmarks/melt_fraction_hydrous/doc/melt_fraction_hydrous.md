```{tags}
category:benchmark
feature:2d
feature:cartesian
feature:melt
```
(sec:benchmarks:melt-fraction-hydrous)=
# Melt fraction hydrous benchmark
This benchmark validates the hydrous extension of the melt fraction
visualization postprocessor. The implementation extends the existing
anhydrous melting parameterization by accounting for the effect of
water on peridotite melting following {cite:t}`ball:etal:2022`, building on
the melting formulation of {cite:t}`katz:etal:2003`.

The domain is a 2d Cartesian box that is not intended to represent any
geological domain, but rather a numerical representation of pressure-temperature
space. The horizontal coordinate is used to set temperatures ranging from 1100 °C to 2000 °C,
whereas the vertical coordinate spans 0–8 GPa of pressure.
This way this benchmark evaluates the computed melt fraction and
compares the numerical results with published references.


## Example setup
The relevant parameters for the hydrous melting formulation are in the parameter file:
```{literalinclude} melt_fraction_hydrous.prm
```
For the validation comparisons, the Python scripts use the parameter values corresponding to each reference case (Katz et al. (2003) or Ball et al. (2022)). The corresponding ASPECT simulations use the same values to ensure a fair direct comparison.

## Hydrous formulation
In Katz et al. (2003) dry formulation, melting prior to clinopyroxene exhaustion is expressed in terms of the normalized temperature
```{math}
T'_{\mathrm{dry}} = \frac{T-T_{\mathrm{solidus}}}{T_{\mathrm{lherzolite\ liquidus}}-T_{\mathrm{solidus}}},
```

from which the melt fraction is calculated as
```{math}
F = \left(T'_{\mathrm{dry}}\right)^{\beta_1}.
```

The hydrous extension accounts for the depression of the solidus caused by water. The temperature depression is described as
```{math}
\Delta T_{\mathrm{H_2O}} = K X_{\mathrm{H_2O}}^{\gamma},
```

where $X_{\mathrm{H_{2}O}}$ is the water concentration in the melt, and the corresponding normalized temperature and melt fraction become
```{math}
:label: eq:hydrousT'
T'_{\mathrm{hydrous}} = \frac{T-\left(T_{\mathrm{solidus}}-\Delta T_{\mathrm{H_2O}}\right)}{T_{\mathrm{lherzolite\ liquidus}}-T_{\mathrm{solidus}}}\quad\text{and}\quad
F = \left(T'_{\mathrm{hydrous}}\right)^{\beta_1}.
```

On the other hand, the melt water concentration depends on the melt fraction according to
```{math}
:label: eq:water-concentration

X_{\mathrm{H_2O}} = \frac{X_{\mathrm{H_2O}}^{\mathrm{bulk}}}{
D_{\mathrm{H_2O}} + F\left(1-D_{\mathrm{H_2O}}\right)}.
```

As a consequence, this hydrous formulation introduces an **implicit dependence**: the melt fraction controls the water concentration in the melt, while the water concentration modifies the solidus depression and therefore the melt fraction. This requires an **iterative solution**, as discussed below.

The formulation of Katz et al. (2003) for the post-clinopyroxene-exhaustion regime is also followed in our implementation. The reaction coefficient and the melt fraction at which cpx is exhausted vary with pressure as
```{math}
:label: eq:rcpx
R_{\mathrm{cpx}}(P) = r_0 + r_1 P \quad\text{and}\quad
F_{\mathrm{cpx-out}} =
\frac{M_{\mathrm{cpx}}}{R_{\mathrm{cpx}}(P)}.
```

Once this happens, the parameterization switches to the post-cpx branch:
```{math}
F = F_{\mathrm{cpx-out}} + \left(1-F_{\mathrm{cpx-out}}\right)
\left(\frac{T-T_{\mathrm{cpx-out}}}{T_{\mathrm{liquidus}}-T_{\mathrm{cpx-out}}}\right)^{\beta_2},
```

where $T_{\mathrm{cpx-out}}$ is the clinopyroxene exhaustion temperature. Setting a different value for the parameter $\beta_2$ is a modification introduced by Ball et al. (2022) to allow for different progression of melting in the cpx branch compared to the opx branch.

By adding the term $K X_{\mathrm{H_2O}}^{\gamma}$ on the fraction numerator, the expression takes its hydrous form.

Therefore, this hydrous formulation introduces the additional parameters `Bulk H2O ppm`, `D H2O`, `K H2O`, `gamma H2O`, `beta2`, `k1 H2O`, `k2 H2O`, and `lambda H2O`. They describe, respectively, the bulk water content, water partitioning coefficient, the magnitude and exponent of the water-induced solidus depression, and the water saturation limit coefficient and exponent.

### Iterative algorithm:
Briefly explained, the implicit melt fraction equation is solved in two steps. First, a bracket containing the current melt fraction is identified by searching the interval from $F=0$ to $F_{\mathrm{cpx-out}}$ (or from $F_{\mathrm{cpx-out}}$ to $F=1$ after clinopyroxene exhaustion), progressively reducing the search increment from $10^{-1}$ to $10^{-9}$. Once a sign-changing bracket is found, this post processor uses the Boost TOMS Algorithm 748 root finder {cite:p}`alefeld:etal:1995` to find $F$. In contrast, the original implementation of Ball et al. (2022) uses Brent's method {cite:p}`brent:1973`.


### Mantle compositions
Ball et al. (2022) defines three mantle configurations: **Primitive**, **50/50**, and **Depleted**. They use the same values of $B_1$ and $\beta_{2}$, but differ in the clinopyroxene mass fraction $M_{\mathrm{cpx}}$, the coefficients $r_0$ and $r_1$ of the reaction coefficient $R_{\mathrm{cpx}}(P)$, and the bulk water content.

The respective values are:

| Configuration | $M_{\mathrm{cpx}}$ | $r_0$ | $r_1$ | Bulk H$_2$O |
|---|---:|---:|---:|---:|
| Primitive | 0.17127832 | 0.99133219 | -0.12357699 | 280 ppm |
| 50/50 | 0.17022480 | 1.09774656 | -0.14365651 | 200 ppm |
| Depleted | 0.16783327 | 1.24715140 | -0.17266209 | 100 ppm |


## Validation
### Katz et al. (2003)
The implementation is first validated against the original parameterization of {cite:t}`katz:etal:2003` by reproducing their Figure 4, where melt fraction is compared as a function of temperature at a constant pressure of 1 GPa for several bulk water contents (0, 200, 500, 1000, and 3000 ppm).

Comparisons between two different resolutions are presented to visualize the convergence of the ASPECT implementation toward the original Katz et al. (2003).

```{figure} Katz_1GPa_comparison_Res7.png
:alt: Comparison between Katz et al. (2003) and ASPECT (resolution 7) melt fractions at 1 GPa.
:width: 80%

Comparison between Katz et al. (2003) and ASPECT (resolution 7) melt fractions at 1 GPa.
```

```{figure} Katz_1GPa_comparison_Res9.png
:alt: Comparison between Katz et al. (2003) and ASPECT (resolution 9) melt fractions at 1 GPa.
:width: 80%

Comparison between Katz et al. (2003) and ASPECT (resolution 9) melt fractions at 1 GPa.
```

The figures were generated with
{download}`plot_katz.py <plot_katz.py>`.

### Ball et al. (2022)
The implementation is further validated against the code of {cite:t}`ball:etal:2022` using the three mantle configurations described above.
It is worth mentioning that for each ASPECT temperature value, the Ball et al. (2022) formulation is evaluated exactly at the same pressure and temperature, allowing a direct point-by-point comparison.

**First**, with Fig 4 of Katz et al. (2003) as inspiration, melt fraction is again compared as a function of temperature at a constant pressure of 1 GPa.

The **second** validation compares the implementations over the pressure-temperature domain computing their difference in melt fraction, defined as
```{math}
\Delta F = F_{\mathrm{ASPECT}} - F_{\mathrm{Ball}},
```

The differences observed are negligible, with $\Delta F$ remaining very close to zero throughout the tested range, indicating that the ASPECT implementation reproduces the Ball et al. (2022) formulation with only minor numerical differences (probably due to our output accuracy) for the evaluated conditions. The plots were practically identical in all mantle configurations so only one of them is shown as an example. These comparisons refer to the original hydrous implementation, before the water saturation extension described below.

```{figure} Ball_1GPa_comparison.png
:alt: Comparison between Ball et al. (2022) and ASPECT melt fractions at 1 GPa.
:width: 80%

1st validation: Comparison between Ball et al. (2022) implementation and ASPECT at 1 GPa.
```

The figure was generated with
{download}`plot_ball.py <plot_ball.py>`.

```{figure} Ball_primitive_PT_comparison.png
:alt: Difference between ASPECT's and Ball's melt fractions.
:width: 80%

2nd validation: Difference between ASPECT and Ball et al. (2022) melt fractions.
```

The figures were generated with
{download}`plot_ball_PT.py <plot_ball_PT.py>`.


### Extension to when $R_{\mathrm{cpx}} \leq 0$
As mentioned above, the reaction coefficient is defined as $R_{\mathrm{cpx}} = r_0 + r_1 P$.

Since $r_{\mathrm{1}}$ is negative for the mantle compositions considered in Ball et al. (2022) paper (it specifies experiments show F_{\mathrm{cpx-out}} increases with pressure), $R_{\mathrm{cpx}}$ decreases with pressure and can reach zero when it gets sufficiently high.

By equation {eq}`eq:rcpx`, the original expression becomes undefined or does not make physical sense when $R_{\mathrm{cpx}}=0$ and $R_{\mathrm{cpx}}<0$ respectively. This could occur for example for the depleted mantle at approximately 7.22 GPa.

This ASPECT implementation therefore sets
```{math}
F_{\mathrm{cpx-out}} = 1 \qquad \text{when} \qquad R_{\mathrm{cpx}} \leq 0.
```

In the code, assuming $M_{\mathrm{cpx}}geq 0$, $F_{\mathrm{cpx-out}}$ is constrained to the physical range $0 \leq F_{\mathrm{cpx-out}} \leq 1$, meaning that once the expression would exceed 1, it is kept at 1, and once it drops below 0, it is kept at 0. Thus, this extension is consistent with the limiting behaviour as $R_{\mathrm{cpx}}$ decreases and avoids the artificial jump. It was also verified that the extension maintains the original model behaviour when $R_{\mathrm{cpx}}>0$.

To verify that this extension does not introduce an artificial discontinuity in melt fraction, $F$ is evaluated around the pressure at which $R_{\mathrm{cpx}}=0$. The four selected temperatures (approximately 1700, 1750, 1800, and 1900 °C) span conditions from no melting to high melt fractions. Only the example of the depleted mantle is shown here although the same output arises with the 50/50 setting (due to its parametrization, the primitive mantle is not affected by this issue).

```{figure} extension_Rcpx0_depleted.png
:alt: Melt fraction across the pressure at which R_cpx equals zero for the depleted mantle configuration.
:width: 80%

Melt fraction across $R_{\mathrm{cpx}}=0$ for depleted mantle at approximately 1700, 1750, 1800, and 1900 °C. No artificial discontinuity in $F$ is observed at the transition.
```
The figure was generated with
{download}`plot_extension_Rcpx0.py <plot_extension_Rcpx0.py>`.

## Extension of water saturation correction to post clinopyroxene exhaustion.
### Motivation
In the previous hydrous implementation, as in {cite:t}`ball:etal:2022`, the water-saturation correction was applied only when the original calculation used the clinopyroxene regime. The branch used when the clinopyroxene was exhausted did not call this correction. Since there is no reason in Katz et al. (2003) for this difference, this extension adds this call. This is therefore a new extension of the saturation correction, not an unchanged reproduction of Ball et al. (2022).

The sections below develop the mathematical foundation of the implementation and the numerical behaviour observed in the ASPECT benchmark.

### Mathematical formulation
#### Saturation-limited temperature depression
The water concentration in melt as a function of melt fraction is {eq}`eq:water-concentration`.

The saturation limit is
```{math}
:label: eq:water_saturation

X_{\mathrm{H_2O}}^{\mathrm{sat}}=k_{1}P^{\lambda}+k_{2}P,
```

and then, the induced temperature depression is
```{math}
\Delta T_{\mathrm{sat}}=K_{\mathrm{H_2O}}
\left(X_{\mathrm{H_2O}}^{\mathrm{sat}}\right)^{\gamma}.
```

If $X_{\mathrm{H_2O}}(F)\le X_{\mathrm{H_2O}}^{\mathrm{sat}}$, the function returns the original melt fraction $F$ unchanged.

#### Selection of the corrected melting regime
The water saturation correction potentially reduces the original melt fraction, implying that a melt fraction initially calculated in the post-cpx-exhaustion regime may fall into the cpx regime after correction.

To determine which of the regimes, and hence which set of equations, to use, we first calculate a *candidate melt fraction* $F_1$ using the pre-cpx-exhaustion equation {eq}`eq:hydrousT'` when water saturation is exceeded:

```{math}
F_1 = \left(\frac{T + \Delta T_{\mathrm{sat}} - T_{\mathrm{solidus}}}{T_{\mathrm{lherz,liquidus}} - T_{\mathrm{solidus}}}\right)^{\beta_1}.
```

For $T+\Delta T_{\mathrm{sat}}\leq T_{\mathrm{solidus}}$, we set $F_1=0$.

We then compare $F_1$ with $F_{\mathrm{cpx\_out}}$:

- If $F_1\leq F_{\mathrm{cpx\_out}}$, the candidate belongs to the `cpx` regime and is returned directly.
- If $F_1 > F_{\mathrm{cpx\_out}}$, the candidate melt fraction exceeds the validity range of the `cpx` equation, and we must instead calculate the corrected melt fraction using the `opx` equation:

```{math}
F_2=F_{\mathrm{cpx\_out}}+(1-F_{\mathrm{cpx\_out}})
\left(\frac{T+\Delta T_{\mathrm{sat}}-T_{\mathrm{cpx\_out}}}{T_{\mathrm{liquidus}}-T_{\mathrm{cpx\_out}}}\right)^{\beta_2}
```

**Why is this regime check sufficient?**

The temperature at clinopyroxene exhaustion is defined as:

```{math}
:label: eq:T_cpx_out

T_{\mathrm{cpx\_out}}=T_{\mathrm{solidus}}+F_{\mathrm{cpx\_out}}^{1/\beta_1}(T_{\mathrm{lherz,liquidus}}-T_{\mathrm{solidus}}).
```

Since $\beta_1>0$ and $T_{\mathrm{lherz,liquidus}}>T_{\mathrm{solidus}}$, the `cpx` melt fraction increases with effective temperature. Therefore, following equation {eq}`eq:T_cpx_out`:
```{math}
F_1>F_{\mathrm{cpx\_out}}\iff T+\Delta T_{\mathrm{sat}}>T_{\mathrm{cpx\_out}}
```

Thus, whenever we select the `opx` equation, its temperature-dependent term is positive, assuming $\beta_2>0$ and $T_{\mathrm{liquidus}}>T_{\mathrm{cpx\_out}}$.

For $0\leq F_{\mathrm{cpx\_out}}<1$, this guarantees $F_2>F_{\mathrm{cpx\_out}}$.

The corrected result therefore remains in the `opx` regime and no further checks are needed.

(sec:continuity_cpx_out)=
#### Continuity at clinopyroxene exhaustion

At the boundary $F_1=F_{\mathrm{cpx\_out}}$, rearranging the `cpx` equation gives $T+\Delta T_{\mathrm{sat}}=T_{\mathrm{cpx\_out}}$. Substituting this into the `opx` equation gives
```{math}
F_2=F_{\mathrm{cpx\_out}}+(1-F_{\mathrm{cpx\_out}})
\underbrace{\left(\frac{0}{T_{\mathrm{liquidus}}-T_{\mathrm{cpx\_out}}}\right)^{\beta_2}}_{=0}=F_{\mathrm{cpx\_out}}.
```

So the two saturation-limited expressions match continuously at clinopyroxene exhaustion (assuming positive denominators and strictly positive exponents).

#### Bounding values and limiting cases
With $0\le F_{\mathrm{cpx\_out}}\le1$, a valid original solution $0\le F\le1$, positive temperature denominators, $\beta_1,\beta_2>0$, and a nonnegative water-induced temperature depression ($\Delta T_{\mathrm{sat}}$), the corrected melt fraction is also in $[0,1]$:

The lower bound follows directly: the `cpx` candidate is either zero or a nonnegative power and in `opx`, $F_2$ is the sum of nonnegative terms. For the upper bound, the saturation-limited water concentration must be lower than the original water concentration, so $\Delta T_{\mathrm{sat}}\le\Delta T_{\mathrm{H_2O}}(F)$, and then, if the original solution is $F$, since the melt fraction increases with increasing temperature, the correction can not increase it: $F_{\mathrm{corrected}}\le F\le1$.

If $F_{\mathrm{cpx\_out}}=1$, the regime must always be `cpx`, there is no `opx` interval to correct. If $F_{\mathrm{cpx\_out}}=0$, then $T_{\mathrm{cpx\_out}}=T_{\mathrm{solidus}}$ and positive melting falls totally in the `opx` regime. The root-bracketing logic for these limiting cases was inspected separately and no change specific to this extension was necessary.

### Validation against the previous implementation
The extended implementation was evaluated using the existing `melt_fraction_postprocessor_hydrous` test. The "Previous $F$" reference corresponds to the version of the hydrous melt fraction postprocessor prior to this water saturation extension. Representative results are shown below:

| Case | Pressure (GPa) | Temperature (K) | Previous $F$ | New $F$ |
|---|---:|---:|---:|---:|
| No water oversaturation | 0.25 | 1457.53 | 0.0871563 | 0.0871563 |
| Oversaturation, remaining in `cpx` | ~0 | 1457.53 | 0.108299 | 0.108299 |
| Oversaturation, remaining in `opx` | ~0 | 1513.78 | 0.197329 | 0.188145 |
| Correction crossing from `opx` to `cpx` | ~0 | 1485.65 | 0.174103 | 0.157759 |

At 0.25 GPa and 1457.53 K, the calculated water content is below the saturation limit, and both implementations agree on melt fraction value. At the surface and 1457.53 K, the original corrected result lies in `cpx` and is unchanged by this extension. At the surface and 1513.78 K, the new implementation reduces the fraction while remaining in `opx`. At 1485.65 K, the previous fraction lies above $F_{\mathrm{cpx\_out}}\simeq0.172776$, but the new corrected fraction lies below it. So everything matches our expectations.

#### Global comparison
The validation script {download}`validate_water_saturation.py <validate_water_saturation.py>` compares the Gnuplot output of the original ASPECT implementation, preserved as {download}`original_hydrous_melt_fraction.gnuplot <../original_hydrous_melt_fraction.gnuplot>`, with the output of the extended implementation test. It uses the primitive-mantle parameters of this benchmark and generates the comparisons presented below:

| Metric | Result |
|---|---:|
| Gnuplot data rows compared | 4096 |
| Unchanged rows | 4054 |
| Rows with decreased $F$ | 42 |
| Rows with increased $F$ | 0 |
| Maximum reduction in $F$ | 0.016344 |
| Changed rows originally classified as `opx` | 42 |
| Changed rows with calculated oversaturation | 42 |
| Changed rows crossing `opx` to `cpx` | 2 |

All changes are confined to the surface ($y=8$), where the fifth Gnuplot column shows $P\simeq0$ Pa. Because Gnuplot output may list the same mesh vertex multiple times, the row counts should not be interpreted as distinct mesh points.

#### Plots
With the parameter values and pressure resolution used in this benchmark, differences between the current and previous version of the postprocessor are observed **only** at the surface ($P\simeq0$). This is explained by the pressure dependence of the water saturation limit. However, the correction itself is not restricted to zero pressure. With the current parameters, it may also affect melt fractions at small, nonzero pressures that are not sampled by the benchmark. Different parameter choices could extend the pressure range over which the correction becomes relevant.

At the surface, the water saturation limit is approximately zero, and therefore the saturation-limited temperature depression is also approximately zero:
```{math}
X_{\mathrm{H_2O}}^{\mathrm{sat}}\simeq0 \quad\Rightarrow\quad \Delta T_{\mathrm{sat}}\simeq0.
```

Consequently, considering the {download}`melt_fraction_hydrous.prm <melt_fraction_hydrous.prm>` values, the melt fraction in the extended implementation changes from `cpx` to `opx` at approximately $T_{\mathrm{cpx-out}}=1493.57\ \mathrm{K}$. In contrast, without the extension, melting is allowed in the `opx` regime at lower temperatures. This explains why differences already appear at 1485.65 K: the original solution is in `opx`, whereas the extended solution is still in `cpx`.

The vertical markers show the regime boundaries of the final melt fractions, after applying the water-saturation correction where implemented. At 1482.97 K, the original implementation switches discontinuously from saturation-corrected `cpx` to uncorrected `opx`, and at 1493.57 K, the extended implementation switches continuously between the two saturation-corrected regimes.

At higher temperatures, both implementations eventually reach complete melting $(F=1)$, and their results coincide again. The last sampled temperature showing differences is 2048.15 K, close to the peridotite liquidus at zero pressure (2053.15 K).

The comparison plots {numref}`fig:surface_melt_fraction_difference` and {numref}`fig:surface_melt_fraction_comparison` illustrate these differences.

The benchmark's temperature spacing (approximately 28 K) does not allow us to examine continuity immediately at the transition from `cpx` to `opx`, or determine the exact temperatures at which differences between implementations begin and disappear. Continuity at the transition is instead demonstrated analytically in {ref}`sec:continuity_cpx_out`.

```{figure} surface_melt_fraction_difference.png
:alt: Difference of melt fractions before and after extending the water saturation treatment at the surface.
:width: 80%
:name: fig:surface_melt_fraction_difference

Difference of the previous and extended hydrous melt fraction
implementations at the surface.
```

```{figure} surface_melt_fraction_comparison.png
:alt: Comparison of melt fractions before and after extending the water saturation treatment at the surface.
:width: 80%
:name: fig:surface_melt_fraction_comparison

Comparison of the previous and extended hydrous melt fraction
implementations at the surface.
```
