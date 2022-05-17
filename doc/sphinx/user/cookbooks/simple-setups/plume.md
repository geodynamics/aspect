#### Plume in a 2D chunk

*This section was contributed by Cedric Thieulot and Paul Bremner.*

This cookbook is inspired by Kellogg & King (1997) (Kellogg and King 1997) but
is not an attempt at reproducing the results of their publication. Their study
was entirely carried out in dimensionless form. Here, however, we choose to
build a similar experiment based on Earth dimensions and material properties.
Furthermore, in this cookbook we run strictly in 2D, whereas the original
paper makes use of axisymmetry to restrict the problem to two velocity degrees
of freedom (Kiefer and Hager 1992).

The two-dimensional domain is a section of an annulus, i.e. a 2D chunk with
$3\pi/8 \leq \phi \leq \pi/2$. Free-slip boundary conditions are imposed on
all boundaries. The inner and outer radii are
$R\textsubscript{inner}=3480~\si{\km}$ and
$R\textsubscript{outer}=6371~\si{\km}$, respectively.

Temperature boundary conditions are $T=T\textsubscript{surf}$ at the surface
(outer boundary), insulating on the sides, and $T=T\textsubscript{cmb}$ at the
inner boundary except along a patch $\phi > 7\pi/16$ where
$T=T\textsubscript{patch}=T\textsubscript{surf}+\Delta T$. Note that in the
original publication the portion of the inner boundary that is not the patch
is also insulating. By default, cannot accommodate two different boundary
condition types on one boundary. So, we will prescribe plausible temperatures
along the whole boundary, instead: $T\textsubscript{surf}=0\si{\celsius}$,
$\Delta T=3000\si{\celsius}$, and $T\textsubscript{cmb}=2750\si{\celsius}$
(see for instance Steinberger & Calderwood (Steinberger and Calderwood 2006)).

The domain contains a single fluid described by the &lsquo;visco
plastic&rsquo; material model, with thermal expansion
$\alpha=3\times 10^{-5}~\si{\per\kelvin}$, heat capacity
$C_p=1250~\si{\joule\per\kg\per\kelvin}$, thermal diffusivity
$\kappa=5.5\times 10^{-7}~\si{\square\meter\per\second}$ (corresponding to the
thermal conductivity value
$k=\kappa \rho_0 C_p=2.25~\si{\watt\per\meter\per\kelvin}$), reference
temperature $T_0=1023~\si{\kelvin}$, reference density
$\rho_0=3250~\si{\kg\per\cubic\meter}$, and viscosity
$\eta_0=1.25\times 10^{23}$. Gravity is constant throughout the mantle with
$g=9.81~\si{\meter\per\square\second}$. With these parameters we find that the
dimensionless Rayleigh number is $$Ra 
=\frac{\rho_0 g \alpha \Delta T (R\textsubscript{outer}-R\textsubscript{inner})^3}{\kappa \eta_0}
=10^6.$$ In Kellogg & King (1997), the dimensionless viscosity is
temperature-dependent and defined as: $$\eta'(T')
=\eta_0' \exp\left[ \frac{E}{R \Delta T} \left( \frac{1}{T'+T_0'} -\frac{1}{1+T_0'}  \right)   \right],$$
where $\eta_0'=1$ is the dimensionless viscosity, $E$ is the activation
energy, and $R$ is the universal gas constant. The dimensionless temperature
$T_0'$ is the surface temperature $T_{surf}$ divided by the temperature drop
across the shell $\Delta T$. Assuming that the authors used the common
relationship
$$T'=\frac{T-T\textsubscript{surf}}{T\textsubscript{patch}-T\textsubscript{surf}}= \frac{T-T\textsubscript{surf}}{\Delta T}.$$
then multiplying the equation above by $\eta_0$, it follows that: $$\eta(T)
=\eta_0 \exp\left[ \frac{E}{R} \left( \frac{1}{T-T\textsubscript{surf}  + T\textsubscript{surf}} 
-\frac{1}{T\textsubscript{patch}-T\textsubscript{surf} + T\textsubscript{surf}}  \right)   \right]
=\eta_0 \exp\left[ \frac{E}{R} \left( \frac{1}{T} -\frac{1}{T\textsubscript{patch}}\right) \right]$$
so that $\eta(T_{patch})=\eta_0$.

Kellogg & King (1997) investigated three cases:
$E/R \Delta T = \{0,0.25328,3\}$. Setting
$R=8.31~\si{\joule\per\kelvin\per\mol}$ and $\Delta T=3000\si{\kelvin}$, the
activation energy becomes $E=\{ 0, 6317.6 , 74829.6 \}\si{\joule\per\mole}$,
lower than typical values (above 200&nbsp;kJ, see for example Karato & Wu
(Karato and Wu 1993)).

The viscosity expression can be written as $$\eta(T)  
= \frac12 \underbrace{2 \eta_0  \exp\left( -\frac{E}{R T\textsubscript{patch}} \right) }_{A^{-1}} \exp \frac{E}{R T}$$
which is effectively a diffusion creep-type viscosity. We find that
$A^{-1} = \{ 2.5\times 10^{22}, 1.98\times 10^{22} , 1.6\times 10^{21} \}\si{\pascal\second}$,
or $A = \{ 4\times 10^{-24}, 5.05\times 10^{-24},  6.26\times 10^{-23}\}$. As
in Kellogg & King (1997) the viscosity is limited to
$\eta\textsubscript{max}=1000\eta_0$.

We run the model to steady state for the reason given in Redmond & King
(Redmond and King 2004) (which is similar to Kellogg & King (Kellogg and King
1997), only applied to Mars): &ldquo;We use steady state calculations so that
we can separate time-dependent from parameter-dependent effects. Once again,
it is unlikely that any Mars-sized, or larger, planetary body is in steady
state. These calculations mainly serve as a guide, allowing us to determine
the relationship between surface observations and internal parameters.&rdquo;

Steady-state fields are shown in Fig.&nbsp;[9][]. We find that the velocity
fields are similar between the isoviscous and weakly temperature-dependent
cases, where the convection cell occupies most of the domain. In contrast, the
strongly temperature-dependent experiment showcases a large viscosity zone a
few hundred kilometers thick (the thermal lithosphere), forcing the convection
cell below it.i In other words, the thick lid insulates the convection cell,
raising the temperature within it.

<img src="cookbooks/plume_2D_chunk/doc/vel1.*" title="fig:" id="fig:plume-diff-creep" style="width:5cm" alt="Plume in a 2D chunk. Columns from left to right: isoviscous case, weakly temperature dependent case, and strongly temperature-dependent case. Rows from top to bottom: Velocity, temperature, and viscosity field at steady state. Angular opening of \pi/8." />
<img src="cookbooks/plume_2D_chunk/doc/vel2.*" title="fig:" id="fig:plume-diff-creep" style="width:5cm" alt="Plume in a 2D chunk. Columns from left to right: isoviscous case, weakly temperature dependent case, and strongly temperature-dependent case. Rows from top to bottom: Velocity, temperature, and viscosity field at steady state. Angular opening of \pi/8." />
<img src="cookbooks/plume_2D_chunk/doc/vel3.*" title="fig:" id="fig:plume-diff-creep" style="width:5cm" alt="Plume in a 2D chunk. Columns from left to right: isoviscous case, weakly temperature dependent case, and strongly temperature-dependent case. Rows from top to bottom: Velocity, temperature, and viscosity field at steady state. Angular opening of \pi/8." />  
<img src="cookbooks/plume_2D_chunk/doc/T1.*" title="fig:" id="fig:plume-diff-creep" style="width:5cm" alt="Plume in a 2D chunk. Columns from left to right: isoviscous case, weakly temperature dependent case, and strongly temperature-dependent case. Rows from top to bottom: Velocity, temperature, and viscosity field at steady state. Angular opening of \pi/8." />
<img src="cookbooks/plume_2D_chunk/doc/T2.*" title="fig:" id="fig:plume-diff-creep" style="width:5cm" alt="Plume in a 2D chunk. Columns from left to right: isoviscous case, weakly temperature dependent case, and strongly temperature-dependent case. Rows from top to bottom: Velocity, temperature, and viscosity field at steady state. Angular opening of \pi/8." />
<img src="cookbooks/plume_2D_chunk/doc/T3.*" title="fig:" id="fig:plume-diff-creep" style="width:5cm" alt="Plume in a 2D chunk. Columns from left to right: isoviscous case, weakly temperature dependent case, and strongly temperature-dependent case. Rows from top to bottom: Velocity, temperature, and viscosity field at steady state. Angular opening of \pi/8." />  
<img src="cookbooks/plume_2D_chunk/doc/eta1.*" title="fig:" id="fig:plume-diff-creep" style="width:5cm" alt="Plume in a 2D chunk. Columns from left to right: isoviscous case, weakly temperature dependent case, and strongly temperature-dependent case. Rows from top to bottom: Velocity, temperature, and viscosity field at steady state. Angular opening of \pi/8." />
<img src="cookbooks/plume_2D_chunk/doc/eta2.*" title="fig:" id="fig:plume-diff-creep" style="width:5cm" alt="Plume in a 2D chunk. Columns from left to right: isoviscous case, weakly temperature dependent case, and strongly temperature-dependent case. Rows from top to bottom: Velocity, temperature, and viscosity field at steady state. Angular opening of \pi/8." />
<img src="cookbooks/plume_2D_chunk/doc/eta3.*" title="fig:" id="fig:plume-diff-creep" style="width:5cm" alt="Plume in a 2D chunk. Columns from left to right: isoviscous case, weakly temperature dependent case, and strongly temperature-dependent case. Rows from top to bottom: Velocity, temperature, and viscosity field at steady state. Angular opening of \pi/8." />

Obtaining a steady state is contingent on the narrow angular opening. We find
that simply increasing the angular opening from $\pi/8$ to $\pi/4$ yields only
a statistical steady state, as multiple downwellings occur near the side but
the system never stabilizes (see Fig.&nbsp;[12][]). Also, decreasing $\eta_0$
by a factor 10 would yield $Ra=10^7$. In this case, too, a statistical steady
state is reached (not shown here).

<img src="cookbooks/plume_2D_chunk/doc/exp1_22.*" title="fig:" id="fig:plume-angular-opening" style="width:5cm" alt="Plume in a 2D chunk: Temperature at the end of the run. From left to right: Angular opening of \pi/8, \pi/4 and \pi/2. The first two have reached a steady state while the third one has not." />
<img src="cookbooks/plume_2D_chunk/doc/exp1_45.*" title="fig:" id="fig:plume-angular-opening" style="width:5cm" alt="Plume in a 2D chunk: Temperature at the end of the run. From left to right: Angular opening of \pi/8, \pi/4 and \pi/2. The first two have reached a steady state while the third one has not." />
<img src="cookbooks/plume_2D_chunk/doc/exp1_90.*" title="fig:" id="fig:plume-angular-opening" style="width:5cm" alt="Plume in a 2D chunk: Temperature at the end of the run. From left to right: Angular opening of \pi/8, \pi/4 and \pi/2. The first two have reached a steady state while the third one has not." />

Fig.&nbsp;[13][] shows the time evolution of the root mean square velocity as
a function of time. As mentioned above, no active planet is at steady state so
the time on the horizontal axis is not really meaningful. Also, it is easy to
show that the path to steady state (if at all attained) is vastly influenced
by the initial temperature field. We find that the average velocities are
smaller in the models with larger activation energy, which is reasonable since
on average viscosities increase with an increase in activation energy $E$.
Also, although it looks like the system reaches steady state for the $\pi/2$
opening angle it ultimately proves unstable and becomes chaotic, akin to the
figure on the cover of this manual. Please check the corresponding video
<https://youtu.be/FN2BBmbiA8E> to see the system for the entire duration of
the simulation.

<figure>
<img src="cookbooks/plume_2D_chunk/doc/vrms.*" id="fig:plume-diff-creep-vrms" style="width:10cm" alt="Plume in a 2D chunk: Root mean square velocity for each experiment." /><figcaption aria-hidden="true"><em>Plume in a 2D chunk: Root mean square velocity for each experiment.</em></figcaption>
</figure>

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-KW93" class="csl-entry">

Karato, S. I., and P. Wu. 1993. &ldquo;Rheology of the Upper Mantle: A
Synthesis.&rdquo; *Science* 260: 771&ndash;78.

</div>

<div id="ref-keki97" class="csl-entry">

Kellogg, Louise H, and Scott D King. 1997. &ldquo;The Effect of Temperature
Dependent Viscosity on the Structure of New Plumes in the Mantle: Results of a
Finite Element Model in a Spherical, Axisymmetric Shell.&rdquo;
*Earth&nbsp;Planet.&nbsp;Sci.&nbsp;Lett.* 148 (1-2): 13&ndash;26.
<https://doi.org/10.1016/S0012-821X(97)00025-3>.

</div>

<div id="ref-kiha92" class="csl-entry">

Kiefer, W. S., and B. Hager. 1992. &ldquo;<span class="nocase">Geoid anomalies
and dynamic topography from convection in cylindrical geometry: applications
to mantle plumes on Earth and Venus</span>.&rdquo; *Geophy.&nbsp;J.&nbsp;Int.*
108: 198&ndash;214.

</div>

<div id="ref-reki04" class="csl-entry">

Redmond, Hannah L, and Scott D King. 2004. &ldquo;A Numerical Study of a
Mantle Plume Beneath the Tharsis Rise: Reconciling Dynamic Uplift and
Lithospheric Support Models.&rdquo; *Journal of Geophysical Research: Planets*
109 (E9). <https://doi.org/10.1029/2003JE002228>.

</div>

<div id="ref-stca06" class="csl-entry">

Steinberger, B., and A. R. Calderwood. 2006. &ldquo;<span
class="nocase">Models of large-scale viscous flow in the Earth&rsquo;s mantle
with constraints from mineral physics and surface observations</span>.&rdquo;
*Geophy.&nbsp;J.&nbsp;Int.* 167: 1461&ndash;81.
<https://doi.org/10.1111/j.1365-246X.2006.03131.x>.

</div>

</div>

  [9]: #fig:plume-diff-creep
  [12]: #fig:plume-angular-opening
  [13]: #fig:plume-diff-creep-vrms
