# Stokes solver



Problems in ASPECT require solving the Stokes
problem below for velocity $u$ and pressure
$p$

```{math}
:label: eq:strong-form
-\nabla \cdot (2 \eta(x)\varepsilon(u))+\nabla p &= f \quad \text{ in } \Omega \\
-\nabla \cdot u&=0 \quad\text{ in } \Omega
```

equipped with boundary conditions.


Here, $x$ lies in our
domain $\Omega$ and $\eta(x)$ is our
viscosity. $\varepsilon(u)=\frac{1}{2}(\nabla u+
\nabla u^T)$ is known as the strain rate tensor.

In this document, we discuss the case
of Dirichlet boundary conditions
$u=g_D$ on $\Gamma_D$ and Neumann
$n \cdot (pI-2\varepsilon(u))=g_N$ on
$\Gamma_N$. However,
ASPECT also supports
other boundary conditions.


To enable the
computation of a numerical solution,
and to loosen the regularity
restrictions imposed on the velocity
space by the strong formulation,
we must first derive the weak form of
the Stokes equations. This derivation
can be found in step-22 of deal.II.

To derive the weak form,
let $u \in V_g:=\{\phi \in H^1(\Omega)^d: \phi_{\Gamma_D}=g_d\}$,
$v \in V_0=\{\phi \in H^1(\Omega)^d: \phi_{\Gamma_D}=0\}$. We call $V_g$ the solution
space and $V_0$ the test space.

Multiplying the first equation in the strong form by $v \in V_0$ and integrating by parts:

```{math}
:label: eq:stokes-weak-form-1
\begin{aligned}
(-\nabla \cdot (2 \eta(x)\varepsilon(u)),v)+(\nabla p,v)&=
(2\eta(x)\varepsilon(u),\nabla v)-(n \otimes v,2\eta(x)\varepsilon(u))_{d\Omega}\\
&-(p,\nabla \cdot v)+(n \cdot v,p)_{d\Omega}\\
&=(2\eta(x)\varepsilon(u),\varepsilon(v))-(p,\nabla \cdot v)+(v,n \cdot (pI-2\varepsilon(u))_{\Gamma_N})\\
&=(f,v),
\end{aligned}
```

where we use that $(n \otimes v,2\eta(x)\varepsilon(u))_{\Gamma_D}+(n\cdot v,p)_{\Gamma_D}=0$
due to the definition of $V_0$, and the observation that
$(2\eta(x)\varepsilon(u),\nabla v)=(2\eta(x)\varepsilon(u),\varepsilon(v))$.

Of course, the $\nabla \cdot u =0$ must also be
incorporated into the
weak formulation. Multiplying by $q \in Q=L^2$,
we get $(q,\nabla \cdot u)=0$.

Combining these two results,
we have the weak form

```{math}
:label: eq:stokes-weak-form-2
\begin{aligned}
(\varepsilon(v),2\varepsilon(u))-(\nabla \cdot v,p)
-(q,\nabla \cdot u)&=(v,f)-(v,g_N)_{\Gamma_N}
\end{aligned}
```
where we apply the Neumann boundary
condition to get the $(v,g_N)_{\Gamma_N}$ term
from $(v,n \cdot(pI-2\varepsilon(u)_{\Gamma_N}))$.

We proceed now to the
discrete weak form. In order
to guarantee existence and uniqueness
of the discrete velocity and
pressure solution, we
search for $u_h \in Q_{k+1}^d$
and $p_h \in Q_k$, with $k \geq 1$. The use
$Q_{k+1}$ for velocity and  $Q_k$ for pressure for $k \geq 1$ satisfies
the inf-sup condition
and is typically known as the Taylor-Hood pair.
Other options are available in ASPECT,
including discontinuous pressure and Q1-Q1
stabilized.

This yields the linear system




```{math}
:label: eq:discrete-weak-form
\begin{aligned}
    Mx&=\begin{pmatrix}A & B^T \\
    B & 0
    \end{pmatrix}
    \begin{pmatrix} u\\
    p
    \end{pmatrix} =
    \begin{pmatrix}
    f \\
    0
    \end{pmatrix}.
\end{aligned}
```

Here, $A$ is the viscous stress matrix,
$B^T$ is the discrete gradient operator on the
pressure space, and $B$ is the divergence operator
on the velocity space.


We solve this linear system iteratively with
preconditioner

```{math}
:label: eq:preconditioning-1
    P^{-1}&=\begin{pmatrix}
    \widetilde{A}^{-1} & \widetilde{A}^{-1}B^TS^{-1} \\
    0      & -S^{-1}
    \end{pmatrix}
```

applied from the right.



The motivation for this is that

```{math}
    :label: eq:preconditioning-2
    MP^{-1}&= \begin{pmatrix}
    A\widetilde{A}^{-1} & 0 \\
    BA^{-1} & BA^{-1}B^TS^{-1}
    \end{pmatrix} .
```

Notice that $A\widetilde{A}^{-1}$ in the top left block and  $BA^{-1}B^TS^{-1}$ in the bottom right block. This $BA^{-1}B^T$ is often referred to as the negative Schur Complement.
The effectiveness of the preconditioner depends significantly on how well $S^{-1}$ approximates
the inverse of the negative Schur complement.


# Pressure Scaling

ASPECT uses dimensional units,
and consequently a pressure scaling factor
based on reference viscosity
and length scale is required in
order to ensure the iterative
solver enforces both
equations in the Stokes system. To
understand why, we need to
examine the iterative linear solver
{cite}`kronbichler:etal:2012`.



An iterative method such as GMRES will continue iteration until
```{math}
\lvert\lvert\begin{pmatrix}
    F_u-AU^{(k)}-B^TP^{(k)} \\
    -BU^{(k)}
\end{pmatrix}\rvert\rvert< \text{ tol }.
```




We must be aware that the $U$ block corresponding to velocity and
$P$ block corresponding to pressure will
have different units, and this generally leads
to one of the equations,
likely $\nabla \cdot u=0$, not being
enforced. This is because the
two equations have vastly different
numerical values, meaning
one will contribute to the residual
significantly more than the other.


For instance, assume the residual of the first block
has associated units $\frac{Pa}{m} m^{\text{dim}}$
where $Pa=\frac{kg}{ms^2}$
and the residual of the
second block has units $\frac{m^{\text{dim}}}{s}$.
Then the norm of the residual
vector has units
$m^{\text{dim}-1}\sqrt{(Pa)^2+(\frac{m}{s})^2}$.



To remedy this, we recall our Stokes system of the form
```{math}
:label: eq: pressure-scaling-1
\begin{aligned}
        -\nabla \cdot (2 \eta(x)
        \varepsilon(u))+\nabla p & = f\\
        \nabla \cdot u&=0.
\end{aligned}
```
We introduce the pressure scaling $\lambda:=\frac{\eta}{L}$ and scale $\nabla \cdot u$ as follows:

```{math}
:label: eq: pressure-scaling-2
\begin{aligned}
 \lambda \nabla \cdot u &=0
 \end{aligned},
```

where $L$ is determined by length scale and computed reference viscosity.

However, notice that this destroys the symmetry
we have between the $B$ and $B^T$
blocks of our Stokes system. To remedy this, let $\widehat{p}=\lambda^{-1}p$.
Then we can rewrite our system as

```{math}
:label: eq: pressure-scaling-3
\begin{aligned}
    -\nabla \cdot (2 \eta \varepsilon(u))+\nabla (\lambda \widehat{p}) & = f \quad&x \in \Omega\\
    \lambda \nabla \cdot u &=0 \quad&x \in \Omega
\end{aligned}
```

We immediately recover $p$ from
$\widehat{p}$ after solving.
