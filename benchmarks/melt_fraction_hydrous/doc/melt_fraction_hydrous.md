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
water on peridotite melting following Ball et al. (2022), building on
the melting formulation of Katz et al. (2003).

The domain is a 2d Cartesian box that is not intended to represent any
geological domain, but rather a numerical representation of pressure-temperature
space. The horizontal coordinate is used to set temperatures ranging from 1100 °C to 2000 °C,
whereas the vertical coordinate spans 0–8 GPa of pressure.
This way this benchmark evaluates the computed melt fraction and
compares the numerical results with published references.


## Example setup
The relevant parameters for the hydrous melting formulation are in the parameter file:
```{literalinclude} melt_fraction_hydrous.prm.
```
However, the posterior validation files (.py and figures) shown here make use of specific parameter configurations, although all derived from the same setup, in order to make fair comparisons.

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

and the corresponding normalized temperature and melt fraction become
```{math}
T'_{\mathrm{hydrous}} = \frac{T-\left(T_{\mathrm{solidus}}-\Delta T_{\mathrm{H_2O}}\right)}{T_{\mathrm{lherzolite\ liquidus}}-T_{\mathrm{solidus}}}\quad\text{and}\quad
F = \left(T'_{\mathrm{hydrous}}\right)^{\beta_1}.
```

On the other hand, the melt water concentration depends on the melt fraction according to
```{math}
X_{\mathrm{H_2O}} = \frac{X_{\mathrm{H_2O}}^{\mathrm{bulk}}}{
D_{\mathrm{H_2O}} + F\left(1-D_{\mathrm{H_2O}}\right)}.
```

As a consequence, this hydrous formulation introduces an **implicit dependence**: the melt fraction controls the water concentration in the melt, while the water concentration modifies the solidus depression and therefore the melt fraction. This requires an **iterative solution**, as discussed below.

Katz et al. (2003) formulation after clinopyroxene-exhaustion is kept in this approach. The reaction coefficient and melt fraction at which cpx is exhausted vary with pressure as
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

where $T_{\mathrm{cpx-out}}$ is the clinopyroxene exhaustion temperature. The parameter `beta 2` as exponent is a novelty Ball et al. (2022) introduces with respect to Katz at al. (2003).

By adding the term $K X_{\mathrm{H_2O}}^{\gamma}$ on the fraction numerator, the expression takes its hydrous form.

Therefore, this hydrous formulation introduces the additional parameters `Bulk H2O ppm`, `D H2O`, `K H2O`, `gamma H2O`, `k1 H2O`, `k2 H2O`, and `lambda H2O`. They describe, respectively, the bulk water content, water partitioning coefficient, the magnitude and exponent of the water-induced solidus depression, and the water saturation limit coefficient and exponent.

### Iterative algorithm:
Briefly explained, the implicit melt fraction equation is solved in two steps. First, a bracket containing the root is identified by searching the interval from $F=0$ to $F_{\mathrm{cpx-out}}$ (or from $F_{\mathrm{cpx-out}}$ to $F=1$ after clinopyroxene exhaustion), progressively reducing the search increment from $10^{-1}$ to $10^{-9}$. Once a sign-changing bracket is found, this post processor uses the Boost TOMS Algorithm 748 root finder to find F. In contrast, the original Ball et al. (2022) implementation uses Brent's method.


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
The implementation is first validated against the original Katz et al. (2003) parameterization by reproducing their Figure 4, where melt fraction is compared as a function of temperature at a constant pressure of 1 GPa for several bulk water contents (0, 200, 500, 1000, and 3000 ppm).

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
The implementation is further validated against the Ball et al. (2022) code using the three mantle configurations described above.
It is worth mentioning that for each ASPECT temperature value, the Ball et al. (2022) formulation is evaluated exactly at the same pressure and temperature, allowing a direct point-by-point comparison.

**First**, with Fig 4 of Katz et al. (2003) as inspiration, melt fraction is again compared as a function of temperature at a constant pressure of 1 GPa.

The **second** validation compares the implementations over the pressure-temperature domain computing their difference in melt fraction, defined as
```{math}
\Delta F = F_{\mathrm{ASPECT}} - F_{\mathrm{Ball}},
```

The differences observed are negligible, with $\Delta F$ remaining very close to zero throughout the tested range, indicating that the ASPECT implementation reproduces the Ball et al. (2022) formulation with only minor numerical differences (probably due to our output accuracy) for the evaluated conditions. The plots were practically identical in all mantle configurations so only one of them is shown as an example.

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

Since $r_{\mathrm{1}}$ is negative for the mantle compositions considered in Ball et al. (2022) paper (it specifies experiments show F_{\mathrm{cpx_out}} increases with pressure), $R_{\mathrm{cpx}}$ decreases with pressure and can reach zero when it gets sufficiently high.

By equation {eq}`eq:rcpx`, the original expression becomes undefined or does not make physical sense when $R_{\mathrm{cpx}}=0$ and $R_{\mathrm{cpx}}<0$ respectively. This could occur for example for the depleted mantle at approximately 7.22 GPa.

This ASPECT implementation therefore sets
```{math}
F_{\mathrm{cpx-out}} = 1 \qquad \text{when} \qquad R_{\mathrm{cpx}} \leq 0.
```

In the code, $F_{\mathrm{cpx-out}}$ is constrained to the physical range $0 \leq F_{\mathrm{cpx-out}} \leq 1$, meaning that once the expression would exceed 1, it is kept at 1, and once it drops below 0, it is kept at 0. Thus, this extension is consistent with the limiting behaviour as $R_{\mathrm{cpx}}$ decreases and avoids the artificial jump. It was also verified that the extension mantains the original model behaviour when $R_{\mathrm{cpx}}>0$.

To verify that this extension does not introduce an artificial discontinuity in melt fraction, $F$ is evaluated around the pressure at which $R_{\mathrm{cpx}}=0$. The four selected temperatures (approximately 1700, 1750, 1800, and 1900 °C) span conditions from no melting to high melt fractions. Only the example of the depleted mantle is shown here although the same output arises with the 50/50 setting (due to its parametrization, the primitive mantle is not affected by this issue).

```{figure} extension_Rcpx0_depleted.png
:alt: Melt fraction across the pressure at which R_cpx equals zero for the depleted mantle configuration.
:width: 80%

Melt fraction across $R_{\mathrm{cpx}}=0$ for depleted mantle at approximately 1700, 1750, 1800, and 1900 °C. No artificial discontinuity in $F$ is observed at the transition.
```
The figure was generated with
{download}`plot_extension_Rcpx0.py <plot_extension_Rcpx0.py>`.


## References
The aforementioned:
- Katz, R. F., Spiegelman, M., & Langmuir, C. H. (2003). A new parameterization of hydrous mantle melting. *Geochemistry, Geophysics, Geosystems*, 4(9), 1073. https://doi.org/10.1029/2002GC000433

- Ball, P. W., Duvernay, T., & Davies, D. R. (2022). A coupled geochemical-geodynamic approach for predicting mantle melting in space and time. *Geochemistry, Geophysics, Geosystems*, 23(4), e2022GC010421. https://doi.org/10.1029/2022GC010421
