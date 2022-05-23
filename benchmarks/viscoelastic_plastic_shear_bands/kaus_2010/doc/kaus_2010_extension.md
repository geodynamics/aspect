(sec:benchmarks:the-kaus_2010-brick)=
# The Kaus (2010) Brick

*This section was contributed by Marcel Saaro, Cedric Thieulot and John Naliboff.*

This setup is based on from {cite:t}`kaus10` with the only exception for not using brittle strain weakening. The prm file is located in `/benchmarks/viscoelastic\_plastic\_shear\_bands`. The domain is a Cartesian box of size $(L_x \times L_y)=(\SI{40}{km} \times \SI{10}{km})$.
The setup is shown in Fig. {numref}`fig:kaus_brick`. 
The viscous inclusion of size $(\SI{800}{m} \times \SI{400}{m})$ is centered at the bottom of the domain. Its viscosity is $\eta_i=\SI{1e20}{\pascal\second}$.

The brick material is characterised by an (elasto-)visco-plastic
rheology following a Drucker-Prager yield criterion with a cohesion $c=\SI{40}{\mega\pascal}$ and angle of 
friction $\phi=\SI{30}{\degree}$. The elastic shear modulus is set to $G=\SI{50e10}{\pascal}$.

Both materials (brick and inclusion) have a density of $\rho=\SI{2700}{\kg\per\cubic\meter}$. The gravity is pointing downwards with $g=\SI{10}{\meter\per\square\second}$.
The flow is assumed to be isothermal and incompressible.

The boundary conditions are free slip at the bottom and sides and stress free at the top (free surface). The horizontal component of the velocity is prescribed on the sides such that it results in a background strain rate of $\dot{\varepsilon}=\SI{2e-15}{\per\second}$. By reversing the velocity directions on both sides the brick can be put in compression or extension. 

The mesh is composed of $100\times25$ repetitions, i.e. a resolution of $\SI{400}{\m}$ in the absence of further refinement. The model runs for $\SI{20}{\kilo year}$.

```{figure-md} fig:kaus_brick
<img src="./kaus_2010_brick_setup.*" alt="kaus_brick"  width="100%"/>

Setup for the Kaus (2010) brick in compression.
```

Figure {numref}`kaus_brick_result_time_evolution` shows the effective strain rate field at 3 different times as obtained with the python script included in the corresponding folder. Elastic stresses build up until the yield criterion is reached and plastic deformation is activated, ultimately generating two conjugated shear bands rooted on the inclusion. The plots are produced with the script `kaus_2010_time_evolution.py`.

```{figure-md} fig:kaus_brick_result_time_evolution
<img src="./kaus_2010_brick_result_time_evolution.*" alt="kaus_brick"  width="100%"/>

Time evolution of the shear bands in the extensional case of the Kaus (2010) brick.
```

The profile at the height $\frac{5}{16}$ of the brick is displayed in Fig. {numref}`fig:kaus_brick_result_profile`. It is generated with the script `kaus_2010_analysis.py`.

```{figure-md} fig:kaus_brick_result_profile
<img src="./kaus_2010_brick_result_strain_rate_profile.*" alt="kaus_brick"  width="100%"/>

Strain rate profile at height $\frac{5}{16}$ of the Kaus (2010) brick.
```
