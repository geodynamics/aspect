(parameters:Solver_20parameters)=
# Solver parameters


## **Subsection:** Solver parameters


(parameters:Solver_20parameters/Composition_20solver_20tolerance)=
### __Parameter name:__ Composition solver tolerance
**Default value:** 1e-12

**Pattern:** [Double 0...1 (inclusive)]

**Documentation:** The relative tolerance up to which the linear system for the composition system gets solved. See &lsquo;Stokes solver parameters/Linear solver tolerance&rsquo; for more details.

(parameters:Solver_20parameters/Temperature_20solver_20tolerance)=
### __Parameter name:__ Temperature solver tolerance
**Default value:** 1e-12

**Pattern:** [Double 0...1 (inclusive)]

**Documentation:** The relative tolerance up to which the linear system for the temperature system gets solved. See &lsquo;Stokes solver parameters/Linear solver tolerance&rsquo; for more details.

(parameters:Solver_20parameters/AMG_20parameters)=
## **Subsection:** Solver parameters / AMG parameters
(parameters:Solver_20parameters/AMG_20parameters/AMG_20aggregation_20threshold)=
### __Parameter name:__ AMG aggregation threshold
**Default value:** 0.001

**Pattern:** [Double 0...1 (inclusive)]

**Documentation:** This threshold tells the AMG setup how the coarsening should be performed. In the AMG used by ML, all points that strongly couple with the tentative coarse-level point form one aggregate. The term strong coupling is controlled by the variable aggregation\_threshold, meaning that all elements that are not smaller than aggregation\_threshold times the diagonal element do couple strongly. The default is strongly recommended. There are indications that for the Newton solver a different value might be better. For extensive benchmarking of various settings of the AMG parameters in this section for the Stokes problem and others, see https://github.com/geodynamics/aspect/pull/234.

(parameters:Solver_20parameters/AMG_20parameters/AMG_20output_20details)=
### __Parameter name:__ AMG output details
**Default value:** false

**Pattern:** [Bool]

**Documentation:** Turns on extra information on the AMG solver. Note that this will generate much more output.

(parameters:Solver_20parameters/AMG_20parameters/AMG_20smoother_20sweeps)=
### __Parameter name:__ AMG smoother sweeps
**Default value:** 2

**Pattern:** [Integer range 0...2147483647 (inclusive)]

**Documentation:** Determines how many sweeps of the smoother should be performed. When the flag elliptic is set to true, (which is true for ASPECT), the polynomial degree of the Chebyshev smoother is set to this value. The term sweeps refers to the number of matrix-vector products performed in the Chebyshev case. In the non-elliptic case, this parameter sets the number of SSOR relaxation sweeps for post-smoothing to be performed. The default is strongly recommended. There are indications that for the Newton solver a different value might be better. For extensive benchmarking of various settings of the AMG parameters in this section for the Stokes problem and others, see https://github.com/geodynamics/aspect/pull/234.

(parameters:Solver_20parameters/AMG_20parameters/AMG_20smoother_20type)=
### __Parameter name:__ AMG smoother type
**Default value:** Chebyshev

**Pattern:** [Selection Chebyshev|symmetric Gauss-Seidel ]

**Documentation:** This parameter sets the type of smoother for the AMG. The default is strongly recommended for any normal runs with ASPECT. There are some indications that the symmetric Gauss-Seidel might be better and more stable for the Newton solver. For extensive benchmarking of various settings of the AMG parameters in this section for the Stokes problem and others, see https://github.com/geodynamics/aspect/pull/234.

(parameters:Solver_20parameters/Advection_20solver_20parameters)=
## **Subsection:** Solver parameters / Advection solver parameters
(parameters:Solver_20parameters/Advection_20solver_20parameters/GMRES_20solver_20restart_20length)=
### __Parameter name:__ GMRES solver restart length
**Default value:** 50

**Pattern:** [Integer range 1...2147483647 (inclusive)]

**Documentation:** This is the number of iterations that define the GMRES solver restart length. Increasing this parameter makes the solver more robust and decreases the number of iterations. Be aware that increasing this number increases the memory usage of the advection solver, and makes individual iterations more expensive.

(parameters:Solver_20parameters/Diffusion_20solver_20parameters)=
## **Subsection:** Solver parameters / Diffusion solver parameters
(parameters:Solver_20parameters/Diffusion_20solver_20parameters/Diffusion_20length_20scale)=
### __Parameter name:__ Diffusion length scale
**Default value:** 1.e4

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)]

**Documentation:** Set a length scale for the diffusion of advection fields if the &ldquo;prescribed field with diffusion&rdquo; method is selected for a field. More precisely, this length scale represents the square root of the product of diffusivity and time in the diffusion equation, and controls the distance over which features are diffused. Units: \si{\meter}.

(parameters:Solver_20parameters/Matrix_20Free)=
## **Subsection:** Solver parameters / Matrix Free
(parameters:Solver_20parameters/Matrix_20Free/Execute_20solver_20timings)=
### __Parameter name:__ Execute solver timings
**Default value:** false

**Pattern:** [Bool]

**Documentation:** Executes different parts of the Stokes solver repeatedly and print timing information. This is for internal benchmarking purposes: It is useful if you want to see how the solver performs. Otherwise, you don&rsquo;t want to enable this, since it adds additional computational cost to get the timing information.

(parameters:Solver_20parameters/Matrix_20Free/Output_20details)=
### __Parameter name:__ Output details
**Default value:** false

**Pattern:** [Bool]

**Documentation:** Turns on extra information for the matrix free GMG solver to be printed.

(parameters:Solver_20parameters/Newton_20solver_20parameters)=
## **Subsection:** Solver parameters / Newton solver parameters
(parameters:Solver_20parameters/Newton_20solver_20parameters/Max_20Newton_20line_20search_20iterations)=
### __Parameter name:__ Max Newton line search iterations
**Default value:** 5

**Pattern:** [Integer range 0...2147483647 (inclusive)]

**Documentation:** The maximum number of line search iterations allowed. If the criterion is not reached after this number of iterations, we apply the scaled increment even though it does not satisfy the necessary criteria and simply continue with the next Newton iteration.

(parameters:Solver_20parameters/Newton_20solver_20parameters/Max_20pre_2dNewton_20nonlinear_20iterations)=
### __Parameter name:__ Max pre-Newton nonlinear iterations
**Default value:** 10

**Pattern:** [Integer range 0...2147483647 (inclusive)]

**Documentation:** If the &rsquo;Nonlinear Newton solver switch tolerance&rsquo; is reached before the maximal number of Picard iterations, then the solver switches to Newton solves anyway.

(parameters:Solver_20parameters/Newton_20solver_20parameters/Maximum_20linear_20Stokes_20solver_20tolerance)=
### __Parameter name:__ Maximum linear Stokes solver tolerance
**Default value:** 1e-2

**Pattern:** [Double 0...1 (inclusive)]

**Documentation:** The linear Stokes solver tolerance is dynamically chosen for the Newton solver, based on the Eisenstat Walker (1994) paper (https://doi.org/10.1137/0917003), equation 2.2. Because this value can become larger than one, we limit this value by this parameter.

(parameters:Solver_20parameters/Newton_20solver_20parameters/Nonlinear_20Newton_20solver_20switch_20tolerance)=
### __Parameter name:__ Nonlinear Newton solver switch tolerance
**Default value:** 1e-5

**Pattern:** [Double 0...1 (inclusive)]

**Documentation:** A relative tolerance with respect to the residual of the first iteration, up to which the nonlinear Picard solver will iterate, before changing to the Newton solver.

(parameters:Solver_20parameters/Newton_20solver_20parameters/SPD_20safety_20factor)=
### __Parameter name:__ SPD safety factor
**Default value:** 0.9

**Pattern:** [Double 0...1 (inclusive)]

**Documentation:** When stabilizing the Newton matrix, we can encounter situations where the coefficient inside the elliptic (top-left) block becomes negative or zero. This coefficient has the form $1+x$ where $x$ can sometimes be smaller than $-1$. In this case, the top-left block of the matrix is no longer positive definite, and both preconditioners and iterative solvers may fail. To prevent this, the stabilization computes an $\alpha$ so that $1+\alpha x$ is never negative and so that always $0\le \alpha \le 1$.  On the other hand, we also want to stay away from $1+\alpha x=0$, and so modify the choice of $\alpha$ by a factor $c$ between zero and one so that if $c<1$, we are assured that $1+\alpha x>0$, i.e., bounded away from zero. If $c=1$, we allow $1+\alpha x=0$, i.e., an unsafe situation. If $c=0$, then $\alpha$ is always set to zero which guarantees the desired property that $1+\alpha x=1>0$, but at the cost of a diminished convergence rate of the Newton method.

(parameters:Solver_20parameters/Newton_20solver_20parameters/Stabilization_20preconditioner)=
### __Parameter name:__ Stabilization preconditioner
**Default value:** SPD

**Pattern:** [Selection SPD|PD|symmetric|none ]

**Documentation:** This parameters allows for the stabilization of the preconditioner. If one derives the Newton method without any modifications, the matrix created for the preconditioning is not necessarily Symmetric Positive Definite. This is problematic (see {cite}`fraters:etal:2019`). When &lsquo;none&rsquo; is chosen, the preconditioner is not stabilized. The &lsquo;symmetric&rsquo; parameters symmetrizes the matrix, and &lsquo;PD&rsquo; makes the matrix Positive Definite. &lsquo;SPD&rsquo; is the full stabilization, where the matrix is guaranteed Symmetric Positive Definite.

(parameters:Solver_20parameters/Newton_20solver_20parameters/Stabilization_20velocity_20block)=
### __Parameter name:__ Stabilization velocity block
**Default value:** SPD

**Pattern:** [Selection SPD|PD|symmetric|none ]

**Documentation:** This parameters allows for the stabilization of the velocity block. If one derives the Newton method without any modifications, the matrix created for the velocity block is not necessarily Symmetric Positive Definite. This is problematic (see {cite}`fraters:etal:2019`). When &lsquo;none&rsquo; is chosen, the velocity block is not stabilized. The &lsquo;symmetric&rsquo; parameters symmetrizes the matrix, and &lsquo;PD&rsquo; makes the matrix Positive Definite. &lsquo;SPD&rsquo; is the full stabilization, where the matrix is guaranteed Symmetric Positive Definite.

(parameters:Solver_20parameters/Newton_20solver_20parameters/Use_20Eisenstat_20Walker_20method_20for_20Picard_20iterations)=
### __Parameter name:__ Use Eisenstat Walker method for Picard iterations
**Default value:** false

**Pattern:** [Bool]

**Documentation:** If set to true, the Picard iteration uses the Eisenstat Walker method to determine how accurately linear systems need to be solved. The Picard iteration is used, for example, in the first few iterations of the Newton method before the matrix is built including derivatives of the model, since the Picard iteration generally converges even from points where Newton&rsquo;s method does not.

Once derivatives are used in a Newton method, ASPECT always uses the Eisenstat Walker method.

(parameters:Solver_20parameters/Newton_20solver_20parameters/Use_20Newton_20failsafe)=
### __Parameter name:__ Use Newton failsafe
**Default value:** false

**Pattern:** [Bool]

**Documentation:** When this parameter is true and the linear solver fails, we try again, but now with SPD stabilization for both the preconditioner and the velocity block. The SPD stabilization will remain active until the next timestep, when the default values are restored.

(parameters:Solver_20parameters/Newton_20solver_20parameters/Use_20Newton_20residual_20scaling_20method)=
### __Parameter name:__ Use Newton residual scaling method
**Default value:** false

**Pattern:** [Bool]

**Documentation:** This method allows to slowly introduce the derivatives based on the improvement of the residual. If set to false, the scaling factor for the Newton derivatives is set to one immediately when switching on the Newton solver. When this is set to true, the derivatives are slowly introduced by the following equation: $\max(0.0, (1.0-(residual/switch\_initial\_residual)))$, where switch\_initial\_residual is the residual at the time when the Newton solver is switched on.

(parameters:Solver_20parameters/Operator_20splitting_20parameters)=
## **Subsection:** Solver parameters / Operator splitting parameters
(parameters:Solver_20parameters/Operator_20splitting_20parameters/Reaction_20solver_20relative_20tolerance)=
### __Parameter name:__ Reaction solver relative tolerance
**Default value:** 1e-6

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)]

**Documentation:** The relative solver tolerance used in the ARKode reaction solver. This tolerance is used to adaptively determine the reaction step size. For more details, see the ARKode documentation. This parameter is only used if the &lsquo;ARKode&rsquo; reaction solver type is used. Units: none.

(parameters:Solver_20parameters/Operator_20splitting_20parameters/Reaction_20solver_20type)=
### __Parameter name:__ Reaction solver type
**Default value:** ARKode

**Pattern:** [Selection ARKode|fixed step ]

**Documentation:** This parameter determines what solver will be used when the reactions are computed within the operator splitting scheme. For reactions where the reaction rate is a known, finite quantity, the appropriate choice is &lsquo;ARKode&rsquo;, which uses an ODE solver from SUNDIALs ARKode (adaptive-step additive Runge Kutta ODE solver methods) to compute the solution. ARKode will pick a reasonable step size based on the reaction rate and the given &lsquo;Reaction solver relative tolerance&rsquo;. However, in some cases we have instantaneous reactions, where we know the new value of a compositional field (and the reaction rate would be infinite), or reaction where we need to know or be able to control the step size we use to compute the reactions. In theses cases, it is appropriate to use the &lsquo;fixed step&rsquo; scheme, a method that a forward Euler scheme and a fixed number of steps given by the &lsquo;Reaction time step&rsquo; and &lsquo;Reaction time steps per advection step&rsquo; parameters.

(parameters:Solver_20parameters/Operator_20splitting_20parameters/Reaction_20time_20step)=
### __Parameter name:__ Reaction time step
**Default value:** 1000.0

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)]

**Documentation:** Set a time step size for computing reactions of compositional fields and the temperature field in case operator splitting is used. This is only used when the parameter &ldquo;Use operator splitting&rdquo; is set to true and when the &lsquo;fixed step&rsquo; reaction solver type is used. The reaction time step must be greater than 0. If you want to prescribe the reaction time step only as a relative value compared to the advection time step as opposed to as an absolute value, you should use the parameter &ldquo;Reaction time steps per advection step&rdquo; and set this parameter to the same (or larger) value as the &ldquo;Maximum time step&rdquo; (which is 5.69e+300 by default). Units: Years or seconds, depending on the &ldquo;Use years in output instead of seconds&rdquo; parameter.

(parameters:Solver_20parameters/Operator_20splitting_20parameters/Reaction_20time_20steps_20per_20advection_20step)=
### __Parameter name:__ Reaction time steps per advection step
**Default value:** 0

**Pattern:** [Integer range 0...2147483647 (inclusive)]

**Documentation:** The number of reaction time steps done within one advection time step in case operator splitting is used. This is only used if the parameter &ldquo;Use operator splitting&rdquo; is set to true and when the &lsquo;fixed step&rsquo; reaction solver type is used. If set to zero, this parameter is ignored. Otherwise, the reaction time step size is chosen according to this criterion and the &ldquo;Reaction time step&rdquo;, whichever yields the smaller time step. Units: none.

(parameters:Solver_20parameters/Stokes_20solver_20parameters)=
## **Subsection:** Solver parameters / Stokes solver parameters
(parameters:Solver_20parameters/Stokes_20solver_20parameters/Force_20nonsymmetric_20A_20block_20solver)=
### __Parameter name:__ Force nonsymmetric A block solver
**Default value:** false

**Pattern:** [Bool]

**Documentation:** This parameter determines whether to enforce a solver that supports nonsymmetric matrices when solving the inner $A$ block of the Stokes system. By default ASPECT recognizes cases where the A block is nonsymmetric automatically, and chooses an appropriate solver. However, if the inner A block solver does not converge, this parameter can be set to &rsquo;true&rsquo; to force the use of a solver that can handle nonsymmetric matrices.

(parameters:Solver_20parameters/Stokes_20solver_20parameters/GMRES_20solver_20restart_20length)=
### __Parameter name:__ GMRES solver restart length
**Default value:** 100

**Pattern:** [Integer range 1...2147483647 (inclusive)]

**Documentation:** This is the number of iterations that define the GMRES solver restart length. Increasing this parameter helps with convergence issues arising from high localized viscosity jumps in the domain. Be aware that increasing this number increases the memory usage of the Stokes solver, and makes individual Stokes iterations more expensive.

(parameters:Solver_20parameters/Stokes_20solver_20parameters/IDR_28s_29_20parameter)=
### __Parameter name:__ IDR(s) parameter
**Default value:** 2

**Pattern:** [Integer range 1...2147483647 (inclusive)]

**Documentation:** This is the sole parameter for the IDR(s) Krylov solver and will dictate the number of matrix-vector products and preconditioner applications per iteration (s+1) and the total number of temporary vectors required (5+3*s). For s=1, this method is analogous to BiCGStab. As s is increased this method is expected to converge to GMRES in terms of matrix-vector/preconditioner applications to solution.

(parameters:Solver_20parameters/Stokes_20solver_20parameters/Krylov_20method_20for_20cheap_20solver_20steps)=
### __Parameter name:__ Krylov method for cheap solver steps
**Default value:** GMRES

**Pattern:** [Selection GMRES|IDR(s) ]

**Documentation:** This is the Krylov method used to solve the Stokes system. Both options, GMRES and IDR(s), solve non-symmetric, indefinite systems. GMRES guarantees the residual will be reduced in each iteration while IDR(s) has no such property. On the other hand, the vector storage requirement for GMRES is dependent on the restart length and can be quite restrictive (since, for the matrix-free GMG solver, memory is dominated by these vectors) whereas IDR(s) has a short term recurrence. Note that the IDR(s) Krylov method is not available for the AMG solver since it is not a flexible method, i.e., it cannot handle a preconditioner which may change in each iteration (the AMG-based preconditioner contains a CG solve in the pressure space which may have different number of iterations each step).

(parameters:Solver_20parameters/Stokes_20solver_20parameters/Linear_20solver_20A_20block_20tolerance)=
### __Parameter name:__ Linear solver A block tolerance
**Default value:** 1e-2

**Pattern:** [Double 0...1 (inclusive)]

**Documentation:** A relative tolerance up to which the approximate inverse of the $A$ block of the Stokes system is computed. This approximate $A$ is used in the preconditioning used in the GMRES solver. The exact definition of this block preconditioner for the Stokes equation can be found in {cite}`kronbichler:etal:2012`.

(parameters:Solver_20parameters/Stokes_20solver_20parameters/Linear_20solver_20S_20block_20tolerance)=
### __Parameter name:__ Linear solver S block tolerance
**Default value:** 1e-6

**Pattern:** [Double 0...1 (inclusive)]

**Documentation:** A relative tolerance up to which the approximate inverse of the $S$ block (i.e., the Schur complement matrix $S = BA^{-1}B^{T}$) of the Stokes system is computed. This approximate inverse of the $S$ block is used in the preconditioning used in the GMRES solver. The exact definition of this block preconditioner for the Stokes equation can be found in {cite}`kronbichler:etal:2012`.

(parameters:Solver_20parameters/Stokes_20solver_20parameters/Linear_20solver_20tolerance)=
### __Parameter name:__ Linear solver tolerance
**Default value:** 1e-7

**Pattern:** [Double 0...1 (inclusive)]

**Documentation:** A relative tolerance up to which the linear Stokes systems in each time or nonlinear step should be solved. The absolute tolerance will then be $\| M x_0 - F \| \cdot \text{tol}$, where $x_0 = (0,p_0)$ is the initial guess of the pressure, $M$ is the system matrix, $F$ is the right-hand side, and tol is the parameter specified here. We include the initial guess of the pressure to remove the dependency of the tolerance on the static pressure. A given tolerance value of 1 would mean that a zero solution vector is an acceptable solution since in that case the norm of the residual of the linear system equals the norm of the right hand side. A given tolerance of 0 would mean that the linear system has to be solved exactly, since this is the only way to obtain a zero residual.

In practice, you should choose the value of this parameter to be so that if you make it smaller the results of your simulation do not change any more (qualitatively) whereas if you make it larger, they do. For most cases, the default value should be sufficient. In fact, a tolerance of 1e-4 might be accurate enough.

(parameters:Solver_20parameters/Stokes_20solver_20parameters/Maximum_20number_20of_20expensive_20Stokes_20solver_20steps)=
### __Parameter name:__ Maximum number of expensive Stokes solver steps
**Default value:** 1000

**Pattern:** [Integer range 0...2147483647 (inclusive)]

**Documentation:** This sets the maximum number of iterations used in the expensive Stokes solver. If this value is set too low for the size of the problem, the Stokes solver will not converge and return an error message pointing out that the user didn&rsquo;t allow a sufficiently large number of iterations for the iterative solver to converge.

(parameters:Solver_20parameters/Stokes_20solver_20parameters/Number_20of_20cheap_20Stokes_20solver_20steps)=
### __Parameter name:__ Number of cheap Stokes solver steps
**Default value:** 1000

**Pattern:** [Integer range 0...2147483647 (inclusive)]

**Documentation:** As explained in the paper that describes ASPECT (Kronbichler, Heister, and Bangerth, 2012, see {cite}`kronbichler:etal:2012`) we first try to solve the Stokes system in every time step using a GMRES iteration with a poor but cheap preconditioner. By default, we try whether we can converge the GMRES solver in 200 such iterations before deciding that we need a better preconditioner. This is sufficient for simple problems with variable viscosity and we never need the second phase with the more expensive preconditioner. On the other hand, for more complex problems, and in particular for problems with strongly nonlinear viscosity, the 200 cheap iterations don&rsquo;t actually do very much good and one might skip this part right away. In that case, this parameter can be set to zero, i.e., we immediately start with the better but more expensive preconditioner.

(parameters:Solver_20parameters/Stokes_20solver_20parameters/Stokes_20solver_20type)=
### __Parameter name:__ Stokes solver type
**Default value:** default solver

**Pattern:** [Selection default solver|block AMG|direct solver|block GMG ]

**Documentation:** This is the type of solver used on the Stokes system. The block geometric multigrid solver currently has a limited implementation and therefore may trigger Asserts in the code when used. If this is the case, please switch to &rsquo;block AMG&rsquo;. Additionally, the block GMG solver requires using material model averaging. The &rsquo;default solver&rsquo; chooses the geometric multigrid solver if supported, otherwise the AMG solver.

(parameters:Solver_20parameters/Stokes_20solver_20parameters/Use_20direct_20solver_20for_20Stokes_20system)=
### __Parameter name:__ Use direct solver for Stokes system
**Default value:** false

**Pattern:** [Bool]

**Documentation:** If set to true the linear system for the Stokes equation will be solved using Trilinos klu, otherwise an iterative Schur complement solver is used. The direct solver is only efficient for small problems.

(parameters:Solver_20parameters/Stokes_20solver_20parameters/Use_20full_20A_20block_20as_20preconditioner)=
### __Parameter name:__ Use full A block as preconditioner
**Default value:** false

**Pattern:** [Bool]

**Documentation:** This parameter determines whether we use an simplified approximation of the $A$ block as preconditioner for the Stokes solver, or the full $A$ block. The simplified approximation only contains the terms that describe the coupling of identical components (plus boundary conditions) as described in {cite}`kronbichler:etal:2012`. The full block is closer to the description in {cite}`rudi2017weighted`.

There is no clear way to determine which preconditioner performs better. The default value (simplified approximation) requires more outer GMRES iterations, but is faster to apply in each iteration. The full block needs less assembly time (because the block is available anyway), converges in less GMRES iterations, but requires more time per iteration. There are also differences in the amount of memory consumption between the two approaches.

The default value should be good for relatively simple models, but in particular for very strong viscosity contrasts the full $A$ block can be advantageous. This parameter is always set to true when using the GMG solver.

(parameters:Solver_20parameters/Stokes_20solver_20parameters/Use_20weighted_20BFBT_20for_20Schur_20complement)=
### __Parameter name:__ Use weighted BFBT for Schur complement
**Default value:** false

**Pattern:** [Bool]

**Documentation:** If set to true, the Schur complement approximation in the Block preconditioner uses the weighted BFBT preconditioner, otherwise a weighted mass matrix will be used. The BFBT preconditioner is more expensive, but works better for large viscosity variations.
