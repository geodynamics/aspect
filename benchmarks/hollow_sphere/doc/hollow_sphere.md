# The hollow sphere benchmark

This benchmark is based on {cite:t}`thieulot:2017` in which an analytical
solution to the isoviscous incompressible Stokes equations is derived in a
spherical shell geometry. The velocity and pressure fields are as follows:
```{math}
\begin{aligned}
v_r(r,\theta)      &=& g(r) \cos \theta, \\
v_\theta(r,\theta) &=& f(r) \sin \theta, \\
v_\phi(r,\theta)   &=& f(r) \sin \theta, \\
p(r,\theta)        &=& h(r) \cos \theta ,
\end{aligned}
```
where
```{math}
\begin{aligned}
f(r) &= &\frac{\alpha}{r^2} + \beta r, \\
g(r) &=& -\frac{2}{r^2} \left(  \alpha \ln r + \frac{\beta}{3}  r^3  + \gamma \right),   \\
h(r) &= &\frac{2\mu_0}{r} g(r),\end{aligned}
```
with
```{math}
\begin{aligned}
\alpha&=&-\gamma \frac{R_2^3-R_1^3}{R_2^3 \ln R_1 - R_1^3 \ln R_2}, \\
\beta &=& -3\gamma \frac{\ln R_2 - \ln R_1  }{R_1^3 \ln R_2 - R_2^3 \ln R_1}.
\end{aligned}
```
These two parameters are chosen so that $v_r(R_1)=v_r(R_2)=0$, i.e. the
velocity is tangential to both inner and outer surfaces. The gravity vector is
radial and of unit length, while the density is given by:
```{math}
\rho(r,\theta)=  \left(   \frac{\alpha}{r^4}  (8 \ln r -6) +  \frac{8\beta}{3r}  +8 \frac{\gamma}{r^4}  \right) \cos\theta.
```
We set $R_1=0.5$, $R_2=1$ and $\gamma=-1$. The pressure is zero on both
surfaces so that the surface pressure normalization is used. The boundary
conditions that are used are simply the analytical velocity prescribed on both
boundaries. The velocity and pressure fields are shown in {numref}`fig:hollow-sphere`.

{numref}`fig:hollow-sphere-errors` shows the velocity and pressure errors in the $L_2$-norm as a
function of the mesh size $h$ (taken in this case as the radial extent of the
elements). As expected we recover a third-order convergence rate for the
velocity and a second-order convergence rate for the pressure.


```{figure-md} fig:hollow-sphere
<img src="hollow-sphere.png" />

 Velocity and pressure fields for the hollow sphere benchmark.
```

```{figure-md} fig:hollow-sphere-errors
<img src="errors_hollowsphere.svg" />

Velocity and pressure errors in the $L_2$-norm as a function of the mesh size.
```
