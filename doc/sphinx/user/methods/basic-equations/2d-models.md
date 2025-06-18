(sec:methods:2d-models)=
# Two-dimensional models

ASPECT allows solving both two- and three-dimensional models via a parameter in the input files, see also {ref}`sec:run-aspect:2d-vs-3d`.
At the same time, the world is unambiguously three-dimensional.
This raises the question what exactly we mean when we say that we want to solve two-dimensional problems.

The notion we adopt here - in agreement with that chosen by many other codes - is to think of two-dimensional models in the following way: We assume that the domain we want to solve on is a two-dimensional cross section (parameterized by $x$ and $y$ coordinates) that extends infinitely far in both negative and positive $z$ direction.
Further, we assume that the velocity is zero in $z$ direction and that all variables have no variation in $z$ direction.
As a consequence, we ought to really think of these two-dimensional models as three-dimensional ones in which the $z$ component of the velocity is zero and so are all $z$ derivatives.

If one adopts this point of view, the Stokes equations {math:numref}`eq:stokes-1`-{math:numref}`eq:stokes-2` naturally simplify in a way that allows us to reduce the $3+1$ equations to only $2+1$, but it makes clear that the correct description of the compressible strain rate is still $\varepsilon(\mathbf u) - \frac{1}{3}(\nabla \cdot \mathbf u)\mathbf 1$, rather than using a factor of $\frac{1}{2}$ for the second term.
(A derivation of why the compressible strain rate tensor has this form can be found in {cite:t}`schubert:etal:2001`, sec. 6.5.)

Because we think of these 2d models as just a slice through an infinite cylinder along the $z$-axis, we should still see what happens if we express the strain rate as a $3\times 3$ tensor. Under the assumptions above, we have
```{math}
\begin{aligned}
  \varepsilon(\mathbf u) = \begin{pmatrix} \tfrac{\partial u_x}{\partial x} & \tfrac 12 \tfrac{\partial u_x}{\partial y} + \tfrac 12 \tfrac{\partial u_y}{\partial x} & 0 \\
    \tfrac 12 \tfrac{\partial u_x}{\partial y} + \tfrac 12 \tfrac{\partial u_y}{\partial x} & \tfrac{\partial u_y}{\partial y} & 0 \\
    0 & 0 & 0 \end{pmatrix}
\end{aligned}
```
with the expected zeros in the last row and column: the last row contains derivatives of $u_z=0$, whereas the last column contains $z$-derivatives of all velocities (which are also zero under our assumption). Because the strain rate tensor is zero for all $z$-related entries, our assumption about what 2d models represent is often called "[plane strain](https://en.wikipedia.org/wiki/Infinitesimal_strain_theory#Plane_strain)" in mechanics, indicating that nonzero strain components are all in the $x$-$y$ plane. Interestingly, though, for compressible materials, the full compressible strain rate tensor then reads
```{math}
\begin{aligned}
  \varepsilon(\mathbf u) - \frac{1}{3}(\nabla \cdot \mathbf u)\mathbf 1 = \begin{pmatrix} \tfrac 23 \tfrac{\partial u_x}{\partial x}  - \tfrac 13 \tfrac{\partial u_y}{\partial y} & \tfrac 12 \tfrac{\partial u_x}{\partial y} + \tfrac 12 \tfrac{\partial u_y}{\partial x} & 0 \\
    \tfrac 12 \tfrac{\partial u_x}{\partial y} +  \tfrac 12 \tfrac{\partial u_y}{\partial x}  & \tfrac 23 \tfrac{\partial u_y}{\partial y} - \tfrac 13 \tfrac{\partial u_x}{\partial x}  & 0 \\
  0 & 0 & - \tfrac 13 \tfrac{\partial u_y}{\partial y} - \tfrac 13 \tfrac{\partial u_x}{\partial x} \end{pmatrix}.
\end{aligned}
```
The entry in the $(3,3)$ position of this tensor may be surprising.
It disappears, however, when taking the (three-dimensional) divergence of the stress, as is done in {math:numref}`eq:stokes-1`, because the divergence applies the $z$ derivative to all elements of the last row - and the assumption above was that all $z$ derivatives are zero; consequently, whatever lives in the third row of the strain rate tensor does not matter.

Two other comments are in order:

* Above, we have described how 2d models are implemented in ASPECT, namely using the "plane strain" approximation. There are, however, other approximations. Specifically, we could have assumed, for example, that where changes in temperature or pressure occur, the resulting change in volume is partitioned in a 2:1 ratio between the $x-y$-plane and $z$-dimension. In this approximation, one would obtain nonzero velocities $u_z$ that are *proportional* to the distance from the $x$-$y$-plane in which we simulate the material. This is related to the "[plane stress](https://en.wikipedia.org/wiki/Plane_stress)" approximation (which, strictly speaking, has $\sigma_{iz}=0$). We specifically do not adopt this perspective.

* It is interesting to consider what happens when we think of the plane strain model if the material is incompressible. Incompressibility of the conceptual 3d object implies that $\nabla \cdot \vec u = \frac{\partial u_x}{\partial x} + \frac{\partial u_y}{\partial y} + \frac{\partial u_z}{\partial z}=0$. Under our assumptions, the last term is zero, and as a consequence, in our two-dimensional cross-section of the infinite three-dimensional domain, we also have $\nabla_{2d} \cdot \vec u_2d} = \frac{\partial u_x}{\partial x} + \frac{\partial u_y}{\partial y}=0$; i.e., the material is also incompressible in the 2d cross section. Perhaps equally importantly, if a material expands its three-dimensional volume by a factor of $\alpha$, under the plane strain assumption, it must also expand its two-dimensional (cross-sectional) area by a factor of $\alpha$. A second consequence is that if you consider the full strain tensor for incompressible materials, i.e., $\varepsilon(\mathbf u) - \frac{1}{3}(\nabla \cdot \mathbf u)\mathbf 1$, then under the plane strain assumption, it does not actually matter whether one uses $\frac{1}{3}$ or $\frac{1}[2}$ (as one might naively do for two-dimensional models): $\nabla \cdot \mathbf u$ is zero, and so the factor we multiply that term by is unimportant. Of course, if the material *is* compressible, then $\nabla \cdot \mathbf u$ is in general not zero, and there *is* a difference between the two choices for the factor. In those cases, it is important to be clear what notion we adopt for the "2d problem"; the assumptions we have laid out above then imply that even for 2d problems, the choice of $\frac{1}{3}$ is correct and should for example be used when computing a deviatoric strain rate for any purpose.


# Comments on the final set of equations

ASPECT solves the equations given in the previous sections in essentially the form stated.
In particular, the form given in {math:numref}`eq:stokes-1` implies that the pressure $p$ we compute is in fact the *total pressure*, i.e., the
sum of hydrostatic pressure and dynamic pressure (however, see {ref}`sec:methods:static-v-dynamic` for more information on this, as well as the extensive discussion of this issue in {cite:t}`kronbichler:etal:2012`).
Consequently, it allows the direct use of this pressure when looking up pressure dependent material parameters.
