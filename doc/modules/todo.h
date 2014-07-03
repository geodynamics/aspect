/**
 * @page TODO TODOs -- things that will eventually need to be addressed
 *
 * <h3>Extensions we'd like to implement</h3>
 *
 * <ol>
 * <li>Mark compositional fields as to their meaning. This may include:
 * <ol>
 * <li>concentration: $0 \leq c \leq 1$
 * <li>binary: $c \in \{0,1\}$
 * <li>positive: $c \geq 0$
 * <li>level set: $c \in \mathbb{R}$
 * <li>arbitrary: $c \in \mathbb{R}$
 * </ol>
 *
 * Depending on the meaning of a field, we may want to use different
 * artificial viscosity approaches, compression, etc. We may also want to
 * project back into the feasible set when passing compositional values to
 * material models.
 *
 * <li>Implement "open boundary conditions" as described in "Using open
 * sidewalls for modelling self-consistent lithosphere subduction dynamics",
 * M. V. Chertova, T. Geenen, A. van den Berg, and W. Spakman see discussion
 * on the aspect mailing list "[aspect-devel] strain rate corner values and
 * velocity boundary conditions"
 *
 * <li>Make snapshotting as atomic as possible so that a job that is
 * terminated while creating a snapshot doesn't get us in trouble
 *
 * <li>Provide an option to store not only the last snapshot but many
 *
 * <li>Move the adiabatic heating term (the last one in the temperature
 * equation) to the left hand side and make it implicit
 *
 * <li>The nonlinear solvers are at best lightly tested. Do this more
 * systemtically.
 *
 * <li>The same can be said of the compressible solvers
 *
 * <li>Rewrite the tracer code
 *
 * <li>Self gravity (Ian Rose has something initial)
 *
 * <li>We need a scheme to verify that plugins are compiled against the same
 * version of ASPECT that they are running under.
 *
 * <li>More benchmarks. One possibility is the benchmark by King that has been
 * described (apart from the original paper by Scott) in the thesis at
 * https://www10.informatik.uni-erlangen.de/Research/Projects/terraneo/docs/schlag-thesis-2014.pdf
 * starting at page 22.
 * </ol>
 *
 *
 * <h3>Bugs we know of</h3>
 *
 * <ol>
 * <li>The temperature equation must contain adiabatic terms in the
 * compressible case.
 * <li>Name the "table" material model better than "table"
 * <li>Re-write the geometry section of the manual since information on the
 * boundary components is no longer there
 *
 * </ol>
 *
 *
 * <h3>Features missing tests</h3>
 *
 * <ol>
 * <li>tests for the various refinement criteria
 * <li>write a testcase for friction heating (Poiseulle flow) and adiabatic
 * heating (constant downward velocity)
 * <li>write tests for the depth averaging functions (test adaptive
 * refinement, averaged composition, non-zero average/sinking velocity)
 *
 * </ol>
 *
 *
 * <h3>Other issues to address</h3>
 *
 * <ol>
 * <li>We should collect a list of publications and posters created with the
 * help of Aspect
 *
 * <li>Unify readme.html in aspect/ and from the webpage again
 *
 * <li>Use extrapolated solution for grid refinement
 *
 * <li>Make initial conditions class say whether it implements a temperature
 * or a temperature perturbation
 *
 * </ol>
 *
 *
 */
