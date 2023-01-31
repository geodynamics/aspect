(sec:methods:approximate-equations:tala)=
# The truncated anelastic liquid approximation (TALA)

The *truncated anelastic liquid approximation (TALA)* further simplifies the ALA by assuming that the variation of the density due to pressure variations is small, i.e., that
```{math}
\begin{aligned}
  \rho(p,T) \approx  \bar\rho (1 - \bar \alpha T').
\end{aligned}
```
This does not mean that the density is not pressure dependent - it will, for example, continue to be
depth dependent because the hydrostatic pressure grows with depth.
It simply means that the deviations from the reference pressure are assumed to be so small that they do not matter in describing the density.
Because the pressure variation $p'$ is induced by the flow field (the static component pressure is already taken care of by the hydrostatic pressure), this assumption in essence means that we assume the flow to be very slow, even beyond the earlier assumption that we can neglect inertial terms when deriving {math:numref}`eq:stokes-1`-{math:numref}`eq:stokes-2`.

This further assumption then transforms{math:numref}`eq:stokes-ALA-1`-{math:numref}`eq:stokes-ALA-2` into the following equations:
```{math}
\begin{aligned}
  -\nabla \cdot \left[2\eta \left(\varepsilon(\mathbf u)
                                  - \frac{1}{3}(\nabla \cdot \mathbf u)\mathbf 1\right)
                \right] + \nabla p' &=  -\bar \alpha \bar\rho T' \mathbf g
                & \qquad  & \textrm{in $\Omega$},  \\
  \nabla \cdot (\bar\rho \mathbf u) &= 0  & \qquad  & \textrm{in $\Omega$}.
\end{aligned}
```
The energy equation is the same as in the ALA case.
