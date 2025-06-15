#Solver



Problems in ASPECT require solving the Stokes
problem below for velocity $u$ and pressure
$p$

```{math}
:label: eq:aligned
-\nabla \cdot (2 \eta(x)\varepsilon(u))+\nabla p &= f \quad \text{ in } \Omega \\
\nabla \cdot u&=0 \quad\text \in \Omega
```

with $u=g_D$ on $\Gamma_D$ and
$n \cdot (pI-2\varepsilon(u))=g_N$ on 
$\Gamma_N$. Here, $\eta(x)$ is our
viscosity and $\Omega$ is our domain. Of
course, $x \in \Omega$. Of course,
we must discretize the problem
and will not do so from the strong form 
of the Stokes problem provided above. This derivation
can be found in step-22 of deal.II.

To derive the weak form, 
let $u \in V_g:=\{\phi \in H^1(\Omega)^d: \phi_{\Gamma_D}=g_d\}$,
$v \in V_0=\{\phi \in H^1(\Omega)^d: \phi_{\Gamma_D}=0\}. We call $V_g$ the solution
space and $V_0$ the test space.

Multiplying the first equation in the strong form by $v \in V_0$ and integrating by parts:

```{math}
:label: eq:aligned
(-\nabla \cdot (2 \eta(x)\varepsilon(u)),v)+(\nabla p,v)&=
(2\eta(x)\varepsilon(u),\nabla v)-(n \otimes v,2\eta(x)\varepsilon(u))_{d\Omega}-(p,\nabla \cdot v)
+(n \cdot v,p)_{d\Omega}\\
&=(2\eta(x)\varepsilon(u),\varepsilon(v))-(p,\nabla \cdot v)+(v,n \cdot (pI-2\varepsilon(u))_{\Gamma_N})\\
&=(f,v),
```

where we use that $(n \otimes v,2\eta(x)\varepsilon(u))_{\Gamma_D}+(n\cdot v,p)_{\Gamma_D}=0$
due to the definition of $V_0$, and the observation that
$(2\eta(x)\varepsilon(u),\nabla v)=(2\eta(x)\varepsilon(u),\varepsilon(v))$.

Of course the $\nabla \cdot u =0$ must also be addressed. Multiplying by $q \in Q=L^2$,
we get $(q,\nabla \cdot u)=0$.

Combining these two results,
we have the weak form

```{math}
:label: eq:aligned
(\varepsilon(v),2\varepsilon(u))-(\nabla \cdot v,p)
-(q,\nabla \cdot u)=(v,f)-(v,g_N)_{\Gamma_N}
```
where we apply the Neumann boundary
condition to get the $(v,g_N)_{\Gamma_N}$ term.

Of course, we must discretize this weak form.
To do so, care has to be 
taken in order to guarantee
the existence of a unique solution.




```{math}
:label: eq:aligned
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
```

Here, $A$ is the viscous stress matrix,
$B^T$ comes from integrating 
grad on pressure space by parts after
testing with $v$, and $B$ comes from
integrating div $u$ by parts after testing
with $q$.


We have a right preconditioner of the form
```{math}
:label: eq:aligned
    P^{-1}&=\begin{pmatrix}
    \widetilde{A}^{-1} & \widetilde{A}^{-1}B^TS^{-1} \\
    0      & -S^{-1}
    \end{pmatrix}.
```


The motivation for this is that

```{math}
    :label: eq:aligned
    MP^{-1}&= \begin{pmatrix}
    A\widetilde{A}^{-1} & 0 \\
    BA^{-1} & BA^{-1}B^TS^{-1}
    \end{pmatrix} .
```

Notice we have $A\widetilde{A}^{-1}$ in the top left block and  $BA^{-1}B^TS^{-1}$ in the bottom right block - 
we want to choose $S$ such that $S^{-1} \approx (BA^{-1}B^T)^{-1}$. Note that 
this $BA^{-1}B^T$ is often referred to as the negative Schur Complement, so we want $S^{-1}$ to serve
as the negative Schur Complement inverse.

Pressure Scaling: ASPECT uses dimensionalized units
and a pressure scaling factor is required in
order to ensure the iterative
solver enforces both 
equations in the Stokes system. 
The pressure scaling factor
is computed
based on viscosity and other parameters.

An iterative method such as GMRES will continue iteration until
```{math}
\lvert\lvert\begin{pmatrix}
    F_u-AU^{(k)}-B^TP^{(k)} \\
    BU^{(k)}
\end{pmatrix}\rvert\rvert< \text{ tol }.
```

We must note that the $U$ block corresponding to velocity and 
$P$ block corresponding to pressure will 
have different units, and this generally leads 
to one of the equations,
likely $\text{div}\mathbf{u}=0$, not being
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
:label: eq:aligned
        -\nabla \cdot (2 \eta(x)
        \varepsilon(u))+\nabla p & = f\\
        \nabla \cdot u&=0.
```
We introduce the pressure scaling $\lambda:=\frac{\eta}{L}$ and scale $\nabla \cdot u$ as follows:

```{math}
:label: eq:aligned
 \lambda \nabla \cdot u =0,
```

where $L$ is determined by length scale and computed reference viscosity.

However, notice that this destroys the symmetry we have in our block Stokes matrix
with the $B$ block. To remedy this, let $\widehat{p}=\lambda^{-1}p$ \cite{kronbichler:etal:2012}.
Then we can rewrite our system as 

```{math}
:label: eq:aligned
    -\nabla \cdot (2 \eta \varepsilon(u))+\nabla (\lambda \widehat{p})  = f \quad&x \in \Omega\\
    \lambda \nabla \cdot u=0 \quad& x \in \Omega.
```