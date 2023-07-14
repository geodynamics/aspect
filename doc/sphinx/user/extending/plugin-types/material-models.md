(sec:extending:plugin-types:material-models)=
# Material models

The material model is responsible for describing the various coefficients in
the equations that ASPECT solves. To implement
a new material model, you need to overload the
[aspect::MaterialModel::Interface](https://aspect.geodynamics.org/doc/doxygen/namespaceaspect_1_1MaterialModel.html)
class and use the
`ASPECT_REGISTER_MATERIAL_MODEL` macro to register your new class. The
implementation of the new class should be in namespace
`aspect::MaterialModel`. An example of a material model implemented this way
is given in {ref}`sec:benchmarks:davies_et_al:case2.3`.

Specifically, your new class needs to implement the following interface:

```{code-block} c++
template <int dim>
    class aspect::MaterialModel::Interface
    {
      public:
        // Physical parameters used in the basic equations
        virtual void evaluate(const MaterialModelInputs &in, MaterialModelOutputs &out) const=0;

        virtual bool is_compressible () const = 0;


        // Reference quantities
        virtual double reference_viscosity () const = 0;


        // Functions used in dealing with run-time parameters
        static void
        declare_parameters (ParameterHandler &prm);

        virtual void
        parse_parameters (ParameterHandler &prm);


        // Optional:
        virtual void initialize ();

        virtual void update ();
}
```

The main properties of the material are computed in the function evaluate()
that takes a struct of type MaterialModelInputs and is supposed to fill a
MaterialModelOutputs structure. For performance reasons this function is
handling lookups at an arbitrary number of positions, so for each variable
(for example viscosity), a std::vector is returned. The following members of
MaterialModelOutputs need to be filled:

```{code-block} c++
struct MaterialModelOutputs
{
          std::vector<double> viscosities;
          std::vector<double> densities;
          std::vector<double> thermal_expansion_coefficients;
          std::vector<double> specific_heat;
          std::vector<double> thermal_conductivities;
          std::vector<double> compressibilities;
}
```

The variables refer to the coefficients $\eta,C_p,k,\rho$ in equations
{math:numref}`eq:stokes-1`&ndash;{math:numref}`eq:temperature`, each as a function of
temperature, pressure, position, compositional fields and, in the case of the
viscosity, the strain rate (all handed in by MaterialModelInputs).
Implementations of evaluate() may of course choose to ignore dependencies on
any of these arguments. In writing a new material model, you should consider
coefficient self-consistency
({ref}`sec:coefficient_self_consistency`).

The remaining functions are used in postprocessing as well as handling
run-time parameters. The exact meaning of these member functions is documented
in the [aspect::MaterialModel::Interface](https://aspect.geodynamics.org/doc/doxygen/namespaceaspect_1_1MaterialModel.html)
class documentation. Note that some of the
functions listed above have a default implementation, as discussed on the
documentation page just mentioned.

The function `is_compressible` returns whether we should consider the material
as compressible or not, see {ref}`sec:methods:approximate-equations:ba` on the
Boussinesq model. As discussed there, incompressibility as described by this
function does not necessarily imply that the density is constant; rather, it
may still depend on temperature or pressure. In the current context,
compressibility simply means whether we should solve the continuity equation
as $\nabla \cdot (\rho \mathbf u)=0$ (compressible Stokes) or as
$\nabla \cdot \mathbf{u}=0$ (incompressible Stokes).

The purpose of the parameter handling functions has been discussed in the
general overview of plugins above.

The functions initialize() and update() can be implemented if desired (the
default implementation does nothing) and are useful if the material model has
internal state. The function initialize() is called once during the
initialization of ASPECT and can be used to
allocate memory, initialize state, or read information from an external file.
The function update() is called at the beginning of every time step.

Additionally, every material model has a member variable
"modeldependence," declared in the Interface class, which can be
accessed from the plugin as `this$\rightarrow$modeldependence`.
This structure describes the nonlinear dependence of the various coefficients
on pressure, temperature, composition or strain rate. This information will be
used in future versions of ASPECT to implement
a fully nonlinear solution scheme based on, for example, a Newton iteration.
The initialization of this variable is optional, but only plugins that declare
correct dependencies can benefit from these solver types. All packaged
material models declare their dependencies in the parseparameters() function
and can be used as a starting point for implementations of new material
models.

Older versions of ASPECT used to have
individual functions like `viscosity()` instead of the `evaluate()` function
discussed above. This old interface is no longer supported, restructure your
plugin to implement `evaluate()` instead (even if this function only calls the
old functions).
