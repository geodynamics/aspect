#### Polydiapirism

*This section was contributed by Cedric Thieulot.*

Diapirs are a type of geologic intrusion in which a more mobile and ductily
deformable material (e.g., salt) is emplaced into (brittle) overlying rocks.
As salt domes are capable of trapping petroleum and natural gas these
structures have been extensively studied (Jackson and Hudec 2017).

We consider in this experiment the three-layer viscous Rayleigh-Taylor
instability proposed by Weinberg and Schmeling (Weinberg and Schmeling 1992)
and we focus in what follows on the case II of that publication. The domain is
a 2D Cartesian box of size $2.24~\si{m} \times 1~\si{m}$. Gravity is
Earth-like ($9.81~\si{\meter\per\square\second}$). Boundary conditions are
free-slip on the sides and top and no-slip at the bottom. All three layers are
initially horizontal. The top layer (fluid 1) has a thickness of
$0.75~\si{m}$, a viscosity $\eta_1=100~\si{\pascal\second}$ and a density
$\rho_1=100~\si{\kg\per\cubic\meter}$. The middle layer (fluid 2) has a
thickness $0.125~\si{\meter}$ with $\rho_2=90~\si{\kg\per\cubic\meter}$ and
$\eta_2=1~\si{\pascal\second}$. The bottom layer (fluid 3) has a thickness
$0.125~\si{\meter}$ with $\rho_3=89~\si{\kg\per\cubic\meter}$ and
$\eta_3=1~\si{\pascal\second}$. The two interfaces between the layers are
perturbed by a random noise of amplitude $\pm 0.001~\si{m}$. Since fluid 3 is
lighter than fluid 2 and fluid 2 is lighter than fluid 1, both interfaces are
unstable. We observe that interface 2-3 deforms first, produces domes which
are subsequently incorporated in the domes being generated at the interface
1-2, as shown in Figure&nbsp;[4][]. The root mean square velocity
(Figure&nbsp;[5][]) shows two slopes in the early stages ($t<50~\si{\second}$)
corresponding to the two different growth rates of the interfaces, as
explained by linear stability analysis (Weinberg and Schmeling 1992; Ramberg
1981).

<img src="diapirs0000.png" title="fig:" id="fig:polydiapirs_density" style="width:48.0%" alt="Figure" />
<img src="diapirs0005.png" title="fig:" id="fig:polydiapirs_density" style="width:48.0%" alt="Figure" />
<img src="diapirs0010.png" title="fig:" id="fig:polydiapirs_density" style="width:48.0%" alt="Figure" />
<img src="diapirs0015.png" title="fig:" id="fig:polydiapirs_density" style="width:48.0%" alt="Figure" />

```{figure-md} fig:polydiapirs_vrms
<img src="vrms.svg" style="width:75.0%" />

Polydiapirism benchmark: Root mean square velocity as a function of time
```

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-jahu17" class="csl-entry">

Jackson, M., and M. Hudec. 2017. *Salt Tectonics: Principles and Practice*.
Cambridge: Cambridge University Press.
<https://doi.org/10.1017/9781139003988>.

</div>

<div id="ref-ramb81" class="csl-entry">

Ramberg, H. 1981. *Gravity, Deformation, and the Earth&rsquo;s Crust: In
Theory, Experiments and Geological Application*. Academic Press, London,
214pp.

</div>

<div id="ref-wesc92" class="csl-entry">

Weinberg, R. F., and H. Schmeling. 1992. &ldquo;<span
class="nocase">Polydiapirs: multiwavelength gravity structures</span>.&rdquo;
*Journal of Structural Geology* 14 (4): 425&ndash;36.

</div>

</div>

  [4]: #fig:polydiapirs_density
  [5]: #fig:polydiapirs_vrms
