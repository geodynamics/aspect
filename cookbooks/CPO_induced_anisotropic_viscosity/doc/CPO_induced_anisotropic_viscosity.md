(sec:cookbooks:CPO_induced_anisotropic_viscosity)=
# CPO induced anisotropic viscosity

*This section was contributed by Yijun Wang and Ágnes Király.*

This cookbook explains how to use the CPO-induced anisotropic viscosity material model to set up a shear box model.

## Introduction
 
The directional-dependence (anisotropy) in the viscous properties of olivine suggests that the effective viscosity for olivine would be different when deformations occur in different directions relative to the CPO. This CPO-induced anisotropic viscosity material model computes AV based on the CPO evolution predicted by D-Rex (Fraters & Billen *et al.* in {cite:t}`fratersbillen:etal:2021`; Kaminski *et al.* in {cite:t}`kaminski:etal:2004`) and modifies the subsequent deformation. 

Our constitutive equation for the relationship between the strain rate and stress using the anisotropic viscosity tensor is adapted from Signorelli et al. (2021):
\varepsilon\ \dot_ij=\gamma\begin J(\sigma_ij\ )〗^((n-1) ) A_ij σ_ij,									(1)
where \gamma is the part of fluidity which is temperature- and grain-size dependent:
\gamma=\gamma_0\ exp(-Q/RT)/d^m,									(2)
\gamma_0=1.1*105 is the isotropic fluidity, Q=530 kJ/mol is the activation energy, R=8.314 m3⋅Pa⋅K−1⋅mol−1 is the gas constant, d=0.001 m is the grain size, and m=0.73 is the grain size exponent. These values are taken from Hansen et al. (2016) and Hirth & Kohlstedt (2003). J(\sigma_ij\ ) is the equivalent yield stress, where \sigma_ij is the deviatoric (anisotropic) stress computed using the tensorial and scalar component of the anisotropic viscosity from the previous time step stored in the prescribed field (Figure 1):
J(\sigma_ij\ )=(F(\sigma_11\ \ -\ \sigma_22\ )^2+G(\sigma_22\ \ -\ \sigma_33\ )^2+H(\sigma_33\ \ -\ \sigma_11\ )^2+2L\begin\sigma_23〗^2+2M〖σ_13〗^2+2N〖σ_12〗^2 )^(1\/2),	(3)
and Aij is the anisotropic tensor of fluidity in Kelvin notation: 
A_ij=\ \ 2/3\ [\begin{matrix}(\begin{matrix}(\begin{matrix}(F+H@-F)&\begin{matrix}(-F@G+F))&\begin{matrix}(\begin{matrix}(-H@-G)&\ \ \ \ \ \begin{matrix}(0@0))&\begin{matrix}(0&0@0&0)@\begin{matrix}(\begin{matrix}(-H\ \ \ \ @0\ \ )&\begin{matrix}(-G@0\ \ ))&\begin{matrix}(\begin{matrix}(H+G@0\ \ \ \ \ \ \ \ )&\begin{matrix}(0@L))&\begin{matrix}(0&0@0&0)@\begin{matrix}(\begin{matrix}(0@0)&\begin{matrix}(\ \ \ \ \ \ 0@\ \ \ \ \ 0))&\begin{matrix}(\begin{matrix}(0@0)&\begin{matrix}(\ \ \ \ \ \ \ \ 0@\ \ \ \ \ \ \ \ 0))&\begin{matrix}(M&0@0&N))],						(4)
J(\sigma_ij\ ) and Aij are computed using Hill coefficients H, J, K, L, M, and N (Hill, 1948), which are obtained from regression analysis of a texture database constructed with olivine textures from laboratory experiments, shear box models, and subduction models (Paper III). 
In the ASPECT material model plugin, strain rate, density, temperature, and other parameters are taken as input to compute the viscosity, which is passed into the Stokes system to compute the stress (Figure 1). As a result, we adapt eq. (1) to be:
\sigma_ij=1/(\gamma\begin J(\sigma_ij)〗^(n-1) )*inv(A_ij )*ε ̇_ij,								(5)
Since the determinant of A_ij is 0, A_ij is a singular, non-invertible matrix. We find the Moore-Penrose pseudoinverse of the matrix A_ij as an approximate of the inverse of A_ij.
Since the Hill coefficients are defined in the microscopic CPO reference frame, and parameters computed in ASPECT are in the macroscopic model reference frame, several reference frame conversions are needed. First, we need to rotate \sigma_ij in eq.(3) from the model reference frame to the CPO reference frame so that J(\sigma_ij\ ) is in the CPO reference frame. This is achieved by constructing a matrix from the eigenvectors corresponding with the largest eigenvalues of the covariance matrix for the a-, b-, and c-axis of olivine textures and then we assign the nearest orthogonal matrix to be the rotation matrix R:
R=[\begin{matrix}(max_eigenvector_a_1&max_eigenvector_b_1&max_eigenvector_c_1@max_eigenvector_a_2&max_eigenvector_b_2&max_eigenvector_c_2@max_eigenvector_a_3&max_eigenvector_b_3&max_eigenvector_c_3\ )],		(6)
We compute the rotation matrix R on the particles and further convert it to Euler angles for computation and memory efficiency. These properties need to be interpolated from particles to fields to be used in the material model (Fig.1). As a result, the anisotropic viscosity material model requires at least one particle in each cell so that all cells can have the texture parameters (Euler angles and eigenvalues) for constructing the rotation matrix R and compute the Hill coefficients. In the material model, the interpolated Euler angles are converted to the rotation matrix again. We use the same notation R to describe the rotation matrix used in the material model in the following paragraphs. 
The inverse of the A tensor then needs to be rotated to the model reference frame. Since inv(A) is the Kelvin notation of the rank-4 tensor, we apply the Kelvin notation representation of the R rotation matrix, R_K, on inv(A_ij):
R_K=[\begin{matrix}(R_11^2&R_12^2&R_13^2&\sqrt2\ast R_12\ast R_13&\sqrt2\ast R_11\ast R_13&\sqrt2\ast R_11\ast R_12@R_21^2&R_22^2&R_23^2&\sqrt2\ast R_22\ast R_23&\sqrt2\ast R_21\ast R_23&\sqrt2\ast R_21\ast R_22@R_31^2&R_32^2&R_33^2&\sqrt2\ast R_32\ast R_33&\sqrt2\ast R_31\ast R_33&\sqrt2\ast R_31\ast R_32@\sqrt2\ast R_21\ast R_31&\sqrt2\ast R_23\ast R_32&\sqrt2\ast R_23\ast R_33&R_22\ast R_33+R_23\ast R_32&R_21\ast R_33+R_23\ast R_31&R_21\ast R_32+R_22\ast R_31@\sqrt2\ast R_11\ast R_31&\sqrt2\ast R_12\ast R_32&\sqrt2\ast R_13\ast R_33&R_12\ast R_33+R_13\ast R_32&R_11\ast R_33+R_13\ast R_31&R_11\ast R_32+R_12\ast R_31@\sqrt2\ast R_11\ast R_21&\sqrt2\ast R_12\ast R_22&\sqrt2\ast R_13\ast R_23&R_12\ast R_23+R_13\ast R_22&R_11\ast R_23+R_13\ast R_32&R_11\ast R_22+R_12\ast R_21\ )]
The final equation involving all reference frame conversions is:
\sigma_ij=1/(\gamma\begin J(R^\prime\ast\sigma_ij\ast R)〗^(n-1) )*(R_K*inv(A_ij )*R_K')*(ε_i ) ̇,						(7)
We save 1/(\gamma\begin J(R^\prime\ast\sigma_ij\ast R)〗^(n-1) ) as the material model viscosity output, which we call the scalar viscosity. The tensorial part of anisotropic viscosity, R_K\ast inv(A_ij\ )\ast R_K\prime is called the stress-strain director and is stored in the compositional fields. Due to the dependence of the scalar viscosity on stress and the stress is determined by the scalar viscosity from the previous time step (Fig.1), the computation of the scalar viscosity is not stable and its value naturally fluctuates up and down. This fluctuation ultimately causes a numerical instability associated with the prediction of the anisotropic viscosity. Therefore, we damp the scalar viscosity using a non-linear Newtonian iteration, and in each iteration, only half of the change in the scalar viscosity is applied until the result is stable (Fig.1):




## Input

We prescribe simple shear in a 3d Cartesian box/cube
with dimensions of $1 \times 1 \times 1 $ $[m^3]$. The shear strain rate is set to
$\dot{\epsilon}_{xz} = -5\times 10^{-6} [s^{-1}]$. The particle with olivine/enstatite volume fraction of 0.7/0.3 is placed
right at the center of the cubic box so it stays stationary. The DRex implementation
keep tracks of rotations of crystal grains within the particle under macroscopic deformation.

```{figure-md} fig:model_setup_3D_box
<img src="model_setup_3D_box.png" style="width:80.0%" />

Description
```

The model computes how crystal grains rotate and align under simple shear and the pole figures visualize this rotation.


## Model results


## References:
- Fraters, M. R. T., & Billen, M. I. (2021).
On the implementation and usability of crystal preferred orientation evolution in geodynamic modeling. Geochemistry, Geophysics, Geosystems, 22(10), e2021GC009846.

- Kaminski, E., Ribe, N. M., & Browaeys, J. T. (2004).
D-Rex, a program for calculation of seismic anisotropy due to crystal lattice preferred orientation in the convective upper mantle. Geophysical Journal International, 158(2), 744–752.

- Karato, S., Jung, H., Katayama, I., & Skemer, P. (2008).
Geodynamic significance of seismic anisotropy of the upper mantle: New insights from laboratory studies. Annu. Rev. Earth Planet. Sci., 36, 59–95.
