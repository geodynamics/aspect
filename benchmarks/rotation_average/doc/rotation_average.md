```{tags}
category:benchmark
feature:2d
feature:cartesian
feature:particles
feature:quaternions
```

(sec:benchmarks:rotation_average)=
# Particle Interpolator for Rotations and Crystal Orientation
*This section was contributed by Theo Häußler*

This benchmark describes a particle interpolator for rotations represented by quaternions. The interpolator was written to interpolate an average crystal orientation with an arbitrary crystal symmetry group onto the model grid. The option `Symmetry group = triclinic`, i.e. no crystal symmetry, corresponds to a standard average for rotations (see {cite}`markley:etal:2007`).

## Definition of Crystal Orientation
Any rotation in 3D can be represented by a 3x3 orthogonal matrix that preserves orientation, formally the group $SO(3)$. For a crystal with certain symmetries there is operations, which leave the crystal orientation invariant. For example a 180°-rotation around $z$ is defined by the rotation, $\mathbf{P}_z = \mathrm{diag}(-1,-1,1)$, and would leave a crystal with monoclinic symmetry around $z$ invariant. Following {cite}`man:2022` Chapter 6, one can define the set of all symmetry operations as,
```{math}
G_{cr} = \{\mathbf{P}_1,\ \dots ,\ \mathbf{P}_{N_{cr}}\} \subset SO(3) \text{ ,}
```
with $\mathbf{P_1}=\mathbf{I}$ and $N_{cr}$ the number of crystal symmetries. Each crystal orientation then has $N_{cr}$ representations such that one crystal orientation is defined by the set,
```{math}
\mathbf{R}G_{cr} = \{ \mathbf{R}\mathbf{P}_i:P_i\in G_{cr} \} \text{ .}
```
Formally, the space of crystal orientation is then defined as a quotient set of the space of rotations,
```{math}
SO(3)/G_{cr} = \{\mathbf{R}G_{cr}:\mathbf{R}\in SO(3)\} \text{ .}
```
We can choose one specific representation of the quotient space, e.g. in terms of a range of euler angles. Then we can send all rotation into this range of values by applying the symmetry operations. In crystallography people say they send all rotations into a fundamental zone. A fundamental zone is not unique and can be centered around any point in the space of rotations.

To define an average between crystal orientations we need to define a distance between two objects in the space $SO(3)/G_{cr}$. If we have two orientations $\mathbf{R}_1G_{cr}$ and $\mathbf{R}_2G_{cr}$, we find a distance by picking the symmetry operations $\mathbf{P}_i$ of one of the crystal orientations such that it is closest to the other. We can also say we send one crystal orientation into the fundamental zone around the other orientation. Following ({cite}`man:2022` Chapter 6.4) we can write,
```{math}
d_{SO(3)/G_{cr}} (\mathbf{R}_1G_{cr}, \mathbf{R}_2G_{cr}) = \min_{\mathbf{P} \in G_{cr}} d(\mathbf{R}_1, \mathbf{R}_2\mathbf{P} ) \text{ ,}
```
here $d(\cdot, \cdot)$ is a metric in the space of rotations.

## Definition of an Average

An average in the space of rotations is defined as the rotation, which minimizes the squared distance to all samples (e.g. {cite}`moakher:2002`),
```{math}
:label: eqn:rotation_average
\bar{\mathbf{R}} =  \underset{\mathbf{R}\in SO(3)}{\rm argmin}\sum_n d^2(\mathbf{R}_i,\mathbf{R}) \text{ .}
```
The two definitions of the metric distance we consider are the euclidean distance (e.g. straight line through the ambient space),
```{math}
:label: eqn:rotation_euclidean_norm
d_{\rm{eucl.}}(\mathbf{R}_1,\mathbf{R}_2) = ||\mathbf{R}_1- \mathbf{R}_2||_F
```
and the geodesic or Riemannian distance (e.g. arc following the curvature of the space),
```{math}
:label: eqn:rotation_riemann_norm
d_{\rm{geod.}}(\mathbf{R}_1,\mathbf{R}_2) =  ||\log(\mathbf{R}_1^T\mathbf{R}_2)||_F \text{ .}
```
In both cases $|| \cdot ||_F$ denotes the Frobenius norm, $||\mathbf{A}||_F:= \sqrt{\mathrm{Tr}(\mathbf{A}^T \mathbf{A}) }$. For a general crystal symmetry, the minimization problem we have to solve becomes,
```{math}
:label: eqn:rotation_symmetry_average
\bar{\mathbf{R}} =  \underset{\mathbf{R}\in SO(3)}{\rm argmin}\sum_i d_{SO(3)/G_{cr}}^2(\mathbf{R}_i G_{cr},\mathbf{R}G_{cr}) = \underset{\mathbf{R}\in SO(3)}{\rm argmin}\sum_i \min_{\mathbf{P}_i \in G_{cr}} d^2(\mathbf{R}_i\mathbf{P}_i,\mathbf{R}) \text{ .}
```
This is a combinatoric optimization problem for finding the optimal fundamental zone where all rotations are closest to one another. We can also rewrite the problem as a minimization problem for the objective function,
```{math}
:label: eqn:objective_function
F(\mathbf{R}, \{\mathbf{P}_i\}) =  \sum_i d^2(\mathbf{R}_i\mathbf{P}_i,\mathbf{R}) \text{ .}
```
With the objective function we can compare different estimates for the average.

## Quaternions

To efficiently solve this minimization problem we represent rotations by unit quaternions, $\mathbf{q} \in \mathbb{H}$. Similar to complex numbers, quaternions can be seen as 4-dimensional numbers, $\mathbf{q} = w + x\mathbf{i} + y\mathbf{j} + z\mathbf{k}$. The quaternion algebra (e.g. $\mathbf{i}\mathbf{j} = \mathbf{k}$) defines the quaternion product $\circ$, which can be related to an active rotation of a vector $\mathbf{v}$. By representing the vector as a quaternion with no scalar part $\mathbf{v} = 0 + v_x\mathbf{i} + v_y\mathbf{j} + v_z\mathbf{k}$, we can write,
```{math}
\mathbf{v}' = \mathbf{R}\cdot \mathbf{v} = \mathbf{q}\circ \mathbf{v} \circ \mathbf{q}^{-1} \text{ .}
```
With some algebra we can write the right-hand side as a rotation matrix and find the function $\mathbf{R}(\mathbf{q})$ and the backwards transform $\mathbf{q}(\mathbf{R})$ (https://en.wikipedia.org/wiki/Quaternions_and_spatial_rotation). Physically, unit quaternions are an axis-angle representation of a rotation. The scalar part gives $\cos(\alpha/2)$ and the vector part gives $\sin(\alpha/2)\mathbf{n}$. Here $\alpha$ is the rotation angle and $\mathbf{n}$ is the normalized rotation axis. Each rotation can then be represented by a quaternion, $\mathbf{q}$ or $-\mathbf{q}$, because rotating by the negative angle around the negative axis is exactly the same. Therefore, quaternions are a double cover of the space of rotations, and we will work in the half-space where $q.w > 0$.
Just as for rotation matrices we can write a rotation by first $\mathbf{q}_1$ and then $\mathbf{q}_2$ as $\mathbf{q}' = \mathbf{q}_2 \circ \mathbf{q}_1$. If we rewrite the symmetry operations in terms of quaternions, $\mathbf{g}_i = \mathbf{q}(\mathbf{P}_i)$, one crystal orientation then becomes the set $\mathbf{q}G_{cr} = \{\mathbf{q}\circ\mathbf{g}_i:g_i \in G_{cr} \}$. For orthotropic symmetry the set of symmetry operations is (e.g. {cite}`kagan:1991`), $G_{\rm ortho.} = \{\pm 1, \pm \mathbf{i}, \pm \mathbf{j}, \pm \mathbf{k}\}$.

Quaternions also enable us to efficiently compute the Riemannian distance as {cite}`huynh:2009`,
```{math}
d_{\rm geod.}(\mathbf{q}_1, \mathbf{q_2}) = 2\arccos(|\mathbf{q}_1\cdot\mathbf{q}_2|) \text{ .}
```
With the riemannian metric and the symmetry operations we can define the objective function in terms of quaternions as,
```{math}
F(\bar{\mathbf{q}}, \{\mathbf{g}_i\}) = \frac{1}{\sum_j w_j}\sum_i w_i[\ 2\arccos(|\bar{\mathbf{q}}\cdot(\mathbf{q}_i\circ\mathbf{g}_i)|)\ ]^2 \text{ ,}
```
where $\cdot$ is the dot product for vectors in $\mathbb{R}^4$. This is a slight generalization, as $w_i$ are weights associated with each rotation (quaternion). In the future, a distance weighting could be implemented based on these.

Quaternions are also convenient as there exists a well defined average. {cite}`markley:etal:2007` show that the minimization defined in Equ. {math:numref}`eqn:rotation_average` for the euclidean norm Equ. {math:numref}`eqn:rotation_euclidean_norm` can be reformulated into,
```{math}
    :label: eqn:markley_average_maximization
    \bar{\mathbf{q}} = \underset{\mathbf{q}\in \mathbb{H}}{\rm argmax} \ \mathbf{q}^T\mathbf{K}\mathbf{q}
```
with the 4x4 traceless matrix,
```{math}
    \mathbf{K} = 4\sum_i w_i\mathbf{q}_i\mathbf{q}_i^T - \mathbf{I}_4\sum_i w_i \text{ .}
```
The eigenvector corresponding to the maximum eigenvalue of $\mathbf{K}$ solves the maximization problem Equ. {math:numref}`eqn:markley_average_maximization`. As the matrix is traceless the largest eigenvalue will always be positive. I will give the result of this procedure the name $\bar{\mathbf{q}}_{\rm mk07}$ and define it as,
```{math}
    \bar{\mathbf{q}}_{\rm mk07}(\{\mathbf{q}_i\}, \{w_i\}) =  \underset{\mathbf{q}\in \mathbb{H}}{\rm argmin} \sum_i w_i d^2_{\rm eucl.}(\mathbf{R}(\mathbf{q}_i), \mathbf{R}(\mathbf{q}))
```
**Note: This average is valid for general rotations** or crystals with triclinic symmetry (no symmetry objects except the identity). In the interpolator this function is called markley average.

## Algorithm

I propose the following algorithm to solve the optimization problem {math:numref}`eqn:rotation_symmetry_average`. Our mean value and the set of symmetry operations will be denoted as $\bar{\mathbf{q}}^n$ and $\{\mathbf{g}^n_i\}$ for each iteration $n$. $n$ is restricted to a maximum number of iterations set to 50 as I have not seen more than 20 iteration for 1000 crystal orientations.

  1. Pick an initial guess $\bar{\mathbf{q}}^0$ and compute the objective function $F^0(\bar{\mathbf{q}}^{0}, \{\mathbf{g}_i^0\})$ with $\mathbf{g}_i^0 = (1, \mathbf{0})$ (no symmetrization).
  2. Project all rotations into the fundamental zone around the guess.
     - For each quaternion $\mathbf{q}_i$, compute the geodesic distances $d_{\rm geod.}(\mathbf{q}_iG_{cr}, \bar{\mathbf{q}}^n)$ for all elements in the set $\mathbf{q}_i G_{cr}$.
     - Pick $\mathbf{q}_i' = \mathbf{q}_i\circ\mathbf{g}_i^n$ such that $d_{\rm geod.}(\mathbf{q}_i', \bar{\mathbf{q}}^n)$ is minimal. In other words, find the representation in $\mathbf{q}_i G_{cr}$ that is closest to the current mean value.
  3. Compute the markley average $\bar{q}^{n+1} = \bar{\mathbf{q}}_{\rm mk07}(\{\mathbf{q}_i'\})$.
  4. Compute the objective function $F^{n+1}(\bar{\mathbf{q}}^{n+1}, \{\mathbf{g}_i^n\})$.
  5. Check if the objective function got smaller, $F^n > F^{n+1}$ ?
     - true : repeat 2,3 and 4
     - false: return $\bar{\mathbf{q}}^{n+1}$

The choice of the initial guess for the combinatoric optimization matters for a small number of particles. Here we choose to guess the rotation with the maximum geodesic distance from the initial mean value. For a unimodal or antipodal distribution this is likely to be close to the actual mean. More information and justification for these choices will be given in a publication at a later point.

## Benchmark cases

These benchmarks consist of two models, each set in a unit box in 2D.
1. The first case is a shear box with a prescribed Stokes solution. See the parameter file [olivineA.prm](https://www.github.com/geodynamics/aspect/blob/main/cookbooks/crystal_preferred_orientation_olivine_fraters_billen_2021/olivineA.prm) and the documentation [crystal_preferred_orientation_olivine_fraters_billen_2021.md](https://www.github.com/geodynamics/aspect/blob/main/cookbooks/crystal_preferred_orientation_olivine_fraters_billen_2021/doc/crystal_preferred_orientation_olivine_fraters_billen_2021.md).
2. The second case is a convection box with an evolving temperature and velocity field. See the parameter file [convection-box](https://www.github.com/geodynamics/aspect/blob/main/cookbooks/convection-box/convection-box.prm) and the documentation [convection-box.md](https://www.github.com/geodynamics/aspect/blob/main/cookbooks/convection-box/doc/convection-box.md).

For each of them we track the crystal preferred orientation on particles and interpolate the output rotation to the grid. The changes done to the prm files mentioned above are as follows.

I added compositional fields and mapped particle properties to activate the particle interpolator.

```{literalinclude} compositional_field.part.prm
```

With the Postprocessor and the Particles we generate an average crystal orientation and specify the particle interpolator.

```{literalinclude} particle.part.prm
```

The full prm is located in [/benchmarks/rotation_average/convection_box_orthorhombic.prm](https://www.github.com/geodynamics/aspect/blob/main/benchmarks/rotation_average/convection_box_orthorhombic.prm).

To compare to the previously implemented way of averaging rotations I also included an adaptation of the convection box case with euler angles

```{literalinclude} euler_angle.part.prm
```

The prm changing the compositional fields and interpolator is located in [/benchmarks/rotation_average/convection_box_euler.prm](https://www.github.com/geodynamics/aspect/blob/main/benchmarks/rotation_average/convection_box_euler.prm).

The output presented below can be derived from quaternions by recovering the rotation angle $\alpha$ and rotation axis $\mathbf{n}$. They can then be used in the [Rodrigues formula](https://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula) to rotate the 100-axis (or a-axis) from the crystal frame $\mathbf{v}_a^{cr} = (1,0,0)^T$ into the laboratory frame $\mathbf{v}_a$,
```{math}
\begin{aligned}
    \alpha &= 2\mathrm{arccos}(q.w) \text{ ,} \\
    \mathbf{n} &= \frac{1}{\sin(\alpha/2)}(q.x,\ q.y,\ q.z)^T \text{ ,} \\
    \mathbf{v}_a &= \cos(\alpha) \mathbf{v}_a^{cr} + \sin(\alpha)(\mathbf{n}\times\mathbf{v}_a^{cr}) + (1-\cos\alpha) \mathbf{n}(\mathbf{n}\cdot\mathbf{v}_a^{cr}) \text{ .}
\end{aligned}
```
If you do not want to compute the rotation axis and angle directly, you can also use the first column of the active rotation matrix defined in `aspect::Utilities::Quaternions::quaternion_to_rotation_matrix` and retrieve the a-axis as,

```{math}
\mathbf{v}_a = (q_x^2 - q_y^2 - q_z^2 + q_w^2,\ 2(q_x q_y + q_z q_w),\ 2(q_x q_z - q_y q_w))^T \text{ .}
```

For euler angles the first row of the passive rotation matrix defined in `aspect::Utilities::zxz_euler_angles_to_rotation_matrix` gives back the a-axis.
```{math}
\mathbf{v}_a = (\cos(\phi_2)\cos(\phi_1) - \cos(\theta)\sin(\phi_1)\sin(\phi_2),\ -( \cos(\phi_2)\sin(\phi_1) + \cos(\theta)\cos(\phi_1)\sin(\phi_2)),\ -\sin(\phi_2)\sin(\theta) \ )^T
```

### Shear Box Benchmark

```{figure-md} fig:rotation_average_shearbox_t0
<img src="shearbox_comparison_t0.*" style="width:100.0%" />

Expected output of the shear box model at time zero. White arrows are the applied velocity boundary condition. Dark red arrows are the interpolated rotation axis. Light red is the interpolated 100-axis. Dark green is the rotation axis on the particles. Light green is the 100-axis on the particles.
```
At time zero one can nicely see how the symmetry average is able to interpolate rotations while taking into account that an orientation of $[\mathbf{v}_1, \mathbf{v}_2, \mathbf{v}_2]$ is equivalent to an orientation of $[-\mathbf{v}_1, -\mathbf{v}_2, \mathbf{v}_2]$.

```{figure-md} fig:rotation_average_shearbox_tlast
<img src="shearbox_comparison_tlast.*" style="width:100.0%" />

Expected output of the shear box model at time t=3.5. White arrows are the applied velocity boundary condition. Dark red arrows are the interpolated rotation axis. Light red is the interpolated 100-axis. Dark green is the rotation axis on the particles. Light green is the 100-axis on the particles.
```
At the last time step one can see that both averages align as the crystal orientation of the average crystal orientation of the 3 particles is close to each other and no rotation needs to be reprojected into the fundamental zone of another.

### Unit Convection Box Benchmark

This is the expected output for the 3 interpolation algorithms applied to the unit convection box.
We compare 3 methods to average a crystal orientation subject to orthorhombic symmetry.1. We output the rotation as Euler angles and use a cell average.
2. Using ageneral rotation average (triclinic average/markley average)
3. Using the new average for orthorhombic crystal symmetry.
By now the only implementation of a rotation average was the Euler angle average (possibly with linear least squares). The triclinic average is a well established method in the aviation and compute graphics industry and the orthorhombic average is a new experimental feature described in the section above.

```{figure-md} fig:rotation_average_t1
<img src="average_comparison.*" style="width:100.0%" />

Expected output of a 2d convection box model in its spin up phase. White arrows are the velocity. Black lines are the 100-axis tracked on particles. White lines is the average 100-axis interpolated to the model grid.
```
In the figure one can see the a-axis or 100-axis orientation of olivine. The black lines are the expected output tracked on particles and the white lines are the interpolated output on the fields. As we use a cell-wide average and discontinuous elements, there are multiple lines at each evaluation point.

In areas where crystal orientations are close to one another, the triclinic and orthorhombic averages align, whereas in more dynamic areas with spatial variations in the flow over the length of one grid cell, they tend to disagree.

```{figure-md} fig:rotation_average_tlast
<img src="average_comparison_tlast.*" style="width:100.0%" />

Expected output of a 2d convection box model in its steady state phase. White arrows are the velocity. Black lines are the 100-axis tracked on particles. White lines is the average 100-axis interpolated to the model grid.
```
A more thorough evaluation of the new algorithm to coupled anisotropic viscous flow will be given in a future paper.

## possible next steps

- implement a distance weighting
- enable using different symmetry groups for different minerals
- generalize how the quaternions are called instead of requiring that cpo bingham average is active and the properties are named cpo mineral x q.w/q.n
