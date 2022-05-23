#### Reproducing rheology of Morency and Doin, 2004

*This section was contributed by Jonathan Perry-Houts*

Modeling interactions between the upper mantle and the lithosphere can be
difficult because of the dynamic range of temperatures and pressures involved.
Many simple material models will assign very high viscosities at low
temperature thermal boundary layers. The pseudo-brittle rheology described in
(Morency and Doin 2004) was developed to limit the strength of lithosphere at
low temperature. The effective viscosity can be described as the harmonic mean
of two non-Newtonian rheologies:
$$v_{\text{eff}} = \left(\frac{1}{v_{\text{eff}}^v}+\frac{1}{v_{\text{eff}}^p}\right)^{-1}$$
where $$\begin{aligned}
  v_{\text{eff}}^v = B \left(\frac{\dot{\epsilon}}{\dot{\epsilon}_\text{ref}}\right)^{-1+1/n_v}
  \exp\left(\frac{E_a +V_a \rho_m g z}{n_v R T}\right),
  \\
  v_{\text{eff}}^p = (\tau_0 + \gamma \rho_m g z) \left( \frac{\dot{\epsilon}^{-1+1/n_p}}
  {\dot{\epsilon}_\text{ref}^{1/n_p}} \right),\end{aligned}$$ where $B$ is a
scaling constant; $\dot{\epsilon}$ is defined as the quadratic sum of the
second invariant of the strain rate tensor and a minimum strain rate,
$\dot{\epsilon}_0$; $\dot{\epsilon}_\text{ref}$ is a reference strain rate;
$n_v$, and $n_p$ are stress exponents; $E_a$ is the activation energy; $V_a$
is the activation volume; $\rho_m$ is the mantle density; $R$ is the gas
constant; $T$ is temperature; $\tau_0$ is the cohesive strength of rocks at
the surface; $\gamma$ is a coefficient of yield stress increase with depth;
and $z$ is depth.

By limiting the strength of the lithosphere at low temperature, this rheology
allows one to more realistically model processes like lithospheric
delamination and foundering in the presence of weak crustal layers. A similar
model setup to the one described in (Morency and Doin 2004) can be reproduced
with the files in the directory [cookbooks/morency_doin_2004]. In
particular, the following sections of the input file are important to
reproduce the setup:

``` prmfile
```

<figure>
<embed src="cookbooks/morency_doin_2004/doc/morency_doin_2004_fig1.pdf" id="fig:md-1" /><figcaption aria-hidden="true"><em>Approximate reproduction of figure 1 from <span class="citation" data-cites="MD04">(Morency and Doin 2004)</span> using the &#x2018;morency doin&#x2019; material model with almost all default parameters. Note the low-viscosity Moho, enabled by the low activation energy of the crustal component.</em></figcaption>
</figure>

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-MD04" class="csl-entry">

Morency, C., and M.-P. Doin. 2004. "Numerical Simulations of the Mantle
Lithosphere Delamination." *Journal of Geophysical Research: Solid Earth
(1978&ndash;2012)* 109 (B3).

</div>

</div>

  [cookbooks/morency_doin_2004]: cookbooks/morency_doin_2004
