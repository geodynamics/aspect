(sec:cookbooks:CPO_induced_anisotropic_viscosity)=
# CPO induced anisotropic viscosity

*This section was contributed by Yijun Wang and Ágnes Király.*

This cookbook explains how to use the CPO-induced anisotropic viscosity material model to set up a shear box model.

## Introduction

Individual crystals of the mineral olivine reorganized their orientations into the crystal-preferred orientations (CPO) under deformation. The directional-dependence (anisotropy) in the viscous properties of olivine crystals suggests that the effective viscosity for olivine rocks/aggregates is different when deformations occur in different directions relative to the CPO. This CPO-induced anisotropic viscosity material model computes anisotropic viscosity based on the CPO evolution predicted by D-Rex ({cite}`fraters_billen_2021_cpo`; {cite}`kaminski2004`) and modifies the subsequent deformation.

Our constitutive equation for the relationship between the strain rate and stress using the anisotropic viscosity tensor is adapted from {cite:t}`signorelli:etal:2021`:

```{math}
:label: eqn:one
\dot{\varepsilon}_{ij} = \gamma J(\sigma_{ij})^{(n-1)} A_{ij} \sigma_{ij}
```
```

where $\gamma$ is the part of fluidity (the inverse of viscosity) which is temperature- and grain-size dependent:

```{math}
:label: eqn:two
\gamma=\gamma_0 exp \left(\frac{-Q}{RT} \right) /d^m
```

$\gamma_0=1.1\times 10^{5}$ is the isotropic fluidity, $Q=530 kJ/mol$ is the activation energy, $R=8.314 m^3 \cdot Pa \cdot K^{−1} \cdot mol^{−1}$ is the gas constant, $d=0.001 m$ is the grain size, and $m=0.73$ is the grain size exponent. These values for olivine are taken from rock experiments performed by {cite:t}`hansen:etal:2016` and {cite:t}`HK04`. $J(\sigma_{ij})$ is the equivalent yield stress, where $\sigma_{ij}$ is the deviatoric (anisotropic) stress computed using the tensorial and scalar component of the anisotropic viscosity:

```{math}
:label: eqn:three
J(\sigma_{ij})=(F(\sigma_{11} - \sigma_{22})^2+G(\sigma_{22} - \sigma_{33})^2+H(\sigma_{33} - \sigma_{11})^2+2L\sigma_{23}^2+2M\sigma_{13}^2+2N\sigma_{12}^2)^{1/2}$$
```

and $A_{ij}$ is the anisotropic tensor of fluidity in Kelvin notation:

```{math}
:label: eqn:four
A_{ij}=\frac{2}{3} \left[
\begin{matrix}
F+H & -F & -H & 0 & 0 & 0 \\
-F & G+F & -G & 0 & 0 & 0 \\
-H & -G & H+G & 0 & 0 & 0 \\
0 & 0 & 0 & L & 0 & 0 \\
0 & 0 & 0 & 0 & M & 0 \\
0 & 0 & 0 & 0 & 0 & N
\end{matrix} \right]
```

$$J(\sigma_{ij})$$ and $A_{ij}$ are computed using Hill coefficients $H, J, K, L, M,$ and $N$ {cite}`hill:1948`, which are obtained from regression analysis of a texture database constructed with olivine textures from laboratory experiments, shear box models, and subduction models (Kiraly et al., in rev.).

In this material model plugin, strain rate, density, temperature, and other parameters are taken as input to compute the viscosity, which is passed into the Stokes system to compute the stress. As a result, we adapt {math:numref}`eqn:one` to be:

```{math}
:label: eqn:five
\sigma_{ij} = \frac{1}{\gamma J(\sigma_{ij})^{(n-1)}} * inv(A_{ij}) * \dot\varepsilon_{ij}
```

Since the Hill coefficients are defined in the microscopic CPO reference frame, and parameters computed in ASPECT are in the macroscopic model reference frame, several reference frame conversions are needed. First, we need to rotate $\sigma_{ij}$ in {math:numref}`eqn:three` from the model reference frame to the CPO reference frame so that $J(\sigma_{ij})$ is in the CPO reference frame. This is achieved by constructing a matrix from the eigenvectors corresponding with the largest eigenvalues of the covariance matrix for the a-, b-, and c-axis of olivine textures and then we assign the nearest orthogonal matrix to be the rotation matrix R:

```{math}
:label: eqn:six
R = \left[
\begin{matrix}
\verb|max_eigenvector|\_{a1} & \verb|max_eigenvector|\_{b1} & \verb|max_eigenvector|\_{c1} \\
\verb|max_eigenvector|\_{a2} & \verb|max_eigenvector|\_{b2} & \verb|max_eigenvector|\_{c2} \\
\verb|max_eigenvector|\_{a3} & \verb|max_eigenvector|\_{b3} & \verb|max_eigenvector|\_{c3} \\
\end{matrix} \right]
```

We compute the rotation matrix R on the particles and further convert it to Euler angles for computation and memory efficiency. These properties need to be interpolated from particles to fields to be used in the material model (Fig.1). As a result, the anisotropic viscosity material model requires at least one particle in each cell so that all cells can have the texture parameters (Euler angles and eigenvalues) for constructing the rotation matrix R and compute the Hill coefficients. In the material model, the interpolated Euler angles are converted to the rotation matrix again. We use the same notation R to describe the rotation matrix used in the material model in the following paragraphs.

The inverse of the A tensor then needs to be rotated to the model reference frame. Since $inv(A)$ is the Kelvin notation of the rank-4 tensor, we apply the Kelvin notation representation of the R rotation matrix, $R_K$, on $inv(A_{ij})$:

```{math}
:label: eqn:seven
R_K = \left[
\begin{matrix}
R_{11}^2 & R_{12}^2 & R_{13}^2 & \sqrt2* R_{12}* R_{13} & \sqrt2* R_{11}* R_{13} & \sqrt2* R_{11}* R_{12} \\
R_{21}^2 & R_{22}^2 & R_{23}^2 & \sqrt2* R_{22}* R_{23} & \sqrt2* R_{21}* R_{23} & \sqrt2* R_{21}* R_{22} \\
R_{31}^2 & R_{32}^2 & R_{33}^2 & \sqrt2* R_{32}* R_{33} & \sqrt2* R_{31}* R_{33} & \sqrt2* R_{31}* R_{32} \\
\sqrt2* R_{21}* R_{31} & \sqrt2* R_{23}* R_{32} & \sqrt2* R_{23}* R_{33} & R_{22}* R_{33}+R_{23}* R_{32} & R_{21}* R_{33}+R_{23}* R_{31} & R_{21}* R_{32}+R_{22}* R_{31} \\
\sqrt2* R_{11}* R_{31} & \sqrt2* R_{12}* R_{32} & \sqrt2* R_{13}* R_{33} & R_{12}* R_{33}+R_{13}* R_{32} & R_{11}* R_{33}+R_{13}* R_{31} & R_{11}* R_{32}+R_{12}* R_{31} \\
\sqrt2* R_{11}* R_{21} & \sqrt2* R_{12}* R_{22} & \sqrt2* R_{13}* R_{23} & R_{12}* R_{23}+R_{13}* R_{22} & R_{11}* R_{23}+R_{13}* R_{32} & R_{11}* R_{22}+R_{12}* R_{21}
\end{matrix} \right]
```

The final equation involving all reference frame conversions is:

```{math}
:label: eqn:eight
\sigma_{ij} = \frac{1}{\gamma J(R\prime*\sigma_{ij}* R)^{(n-1)}} * (R_K * inv(A_{ij}) * R_K\prime) * \dot\varepsilon_i
```

We save $\frac{1}{\gamma J(R\prime*\sigma_{ij}* R)^{(n-1)}}$ as the material model viscosity output, which we call the scalar viscosity. The tensorial part of anisotropic viscosity, $R_K * inv(A_{ij}) * R_K\prime$ is called the stress-strain director and is stored in the additional outputs. Due to the dependence of the scalar viscosity on stress and the stress is determined by the scalar viscosity from the previous time step, the computation of the scalar viscosity is not stable and its value naturally fluctuates up and down. This fluctuation ultimately causes a numerical instability associated with the prediction of the anisotropic viscosity. Therefore, we damp the scalar viscosity using a non-linear Newtonian iteration, and in each iteration, only half of the change in the scalar viscosity is applied until the result is stable.



## Compiling requirement

Since the determinant of $A_{ij}$ is 0, $A_{ij}$ is a singular, non-invertible matrix. We find the MoorePenrose pseudoinverse of the matrix $A_{ij}$ as an approximate of the inverse of $A_{ij}$, using the SCALAPACK package in deal.ii. Thus it is required to link ASPECT with a deal.ii version with the scalapack package included and you can turn on the scalapack package when compiling deal.ii, for example, with `./candi.sh -j 8 --packages="p4est trilinos sundials scalapack dealii"`.

## Input

The usage of the AV material model is demonstrated with a 3d simple shear box model, where its dimension is $1 \times 1 \times 1 $ (non-dimensionalized). The shear strain rate is set to
$0.5$. The origin is the center of the box, and one Olivine particle with 1000 grains sits at the origin to track CPO developments for computation of anisotropic viscosity parameters.

Since the AV material model computes viscosity based on the evolving CPO stored on particles, several setup requirements must be met:

- **Particles per cell**: Each computational cell must contain at least one particle, to allow interpolation of the CPO particle property. This is achieved by setting (in the Particles subsection):

```{literalinclude} min_particles_per_cell.part.prm
```

- **CPO particle property**: The CPO particle property must be stored for use by the AV model. This requires enabling the particle and crystal preferred orientation postprocessors and the relevant subsuctions for them, including the CPO Bingham Average plugin, which calculates the Hill coefficients:

```{literalinclude} cpo_particle_property.part.prm
```

Note: These settings are similar to those used for simulations involving CPO alone. However, for the AV model, it is essential to set `Use rotation matrix = false` in the CPO Bingham Average subsection, so that the CPO is represented using Euler angles, as required.

- **Compositional fields**: The eigenvalues and Euler angles of the CPO tensor are stored in compositional fields. This requires the following input file section:

```{literalinclude} compositional_field.part.prm
```

In the `CPO induced Anisotropic Viscosity` material model subsection, all parameters have reasonable default values and do not need to be manually specified unless customization is needed.
