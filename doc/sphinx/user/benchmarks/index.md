(cha:benchmarks)=
# Benchmarks

Benchmarks are used to verify that a solver solves the problem correctly,
i.e., to *verify* correctness of a code.[3] Over the past decades, the
geodynamics community has come up with a large number of benchmarks. Depending
on the goals of their original inventors, they describe stationary problems in
which only the solution of the flow problem is of interest (but the flow may
be compressible or incompressible, with constant or variable viscosity, etc),
or they may actually model time-dependent processes. Some of them have
solutions that are analytically known and can be compared with, while for
others, there are only sets of numbers that are approximately known. We have
implemented a number of them in <span class="smallcaps">ASPECT</span> to
convince ourselves (and our users) that <span class="smallcaps">ASPECT</span>
indeed works as intended and advertised. Some of these benchmarks are
discussed below. Numerical results for several of these benchmarks are also
presented in a number of papers (such as (Kronbichler, Heister, and Bangerth
2012; Heister et al. 2017; Tosi et al. 2015; Fraters et al. 2019)) in much
more detail than shown here.

Before going on with showing these benchmarks, let us mention that the data
shown below (and in the papers mentioned above) reflect the state of <span
class="smallcaps">ASPECT</span> at a particular time. On the other hand, <span
class="smallcaps">ASPECT</span> has become more accurate and faster over time,
for example by implementing better stabilization schemes for the advection
equations and improving assembly and solver times. We occasionally update
sections of the manual, but when reading through the sections on individual
benchmarks below, it is worthwhile keeping in mind that <span
class="smallcaps">ASPECT</span> may yield different (and often better) results
than the one shown.

:::{toctree}
benchmarks/onset-of-convection/doc/onset-of-convection.md
cookbooks/van-keken/doc/van-keken.md
cookbooks/van-keken-vof/doc/van-keken-vof.md
cookbooks/bunge_et_al_mantle_convection/doc/bunge_et_al_mantle_convection.md
benchmarks/rayleigh_taylor_instability/doc/rayleigh_taylor_instability.md
benchmarks/polydiapirs/doc/polydiapirs.md
benchmarks/sinking_block/doc/sinking_block.md
benchmarks/solcx/doc/solcx.md
benchmarks/solkz/doc/solkz.md
benchmarks/inclusion/doc/inclusion.md
benchmarks/burstedde/doc/burstedde.md
benchmarks/slab_detacahment/doc/slab_detacahment.md
benchmarks/hollow_sphere/doc/hollow_sphere.md
benchmarks/annulus/doc/annulus.md
cookbooks/stokes/doc/stokes.md
benchmarks/viscosity_grooves/doc/viscosity_grooves.md
cookbooks/latent-heat/doc/latent-heat.md
benchmarks/davies_et_al/doc/davies_et_al.md
benchmarks/crameri_et_al/doc/crameri_et_al.md
benchmarks/solitary_wave/doc/solitary_wave.md
benchmarks/operator_splitting/doc/operator_splitting.md
benchmarks/tosi_et_al_2015_gcubed/doc/tosi_et_al_2015_gcubed.md
benchmarks/layeredflow/doc/layeredflow.md
benchmarks/doneahuerta/doc/doneahuerta.md
benchmarks/advection/doc/advection.md
benchmarks/yamauchi_takei_2016_anelasticity/doc/yamauchi_takei_2016_anelasticity.md
benchmarks/gravity_thin_shell/doc/gravity_thin_shell.md
benchmarks/gravity_thick_shell/doc/gravity_thick_shell.md
benchmarks/gravity_mantle/doc/gravity_mantle.md
benchmarks/buiter_et_al_2016_jsg/doc/buiter_et_al_2016_jsg.md
:::
