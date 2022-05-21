# Heating models

The heating model is responsible for describing the various terms in the
energy equation&nbsp;{math:numref}`eq:temperature`3], using the coefficients provided
by the material model. These can be source terms such as radiogenic heat
production or shear heating, they can be terms on the left-hand side of the
equation, such as part of the latent heating terms, or they can be heating
processes related to reactions. Each of these terms is described by a
"heating model," and a simulation can have none, one, or many
heating models that are active throughout a simulation, with each heating
model usually only implementing the terms for one specific heating process.
One can then decide in the input file which heating processes should be
included in the computation by providing a list of heating models in the input
file.

When the equations are assembled and solved, the heating terms from all
heating models used in the computation are added up.

To implement a new heating model, you need to overload the
[aspect::HeatingModel::Interface][] class and use the
`ASPECT_REGISTER_HEATING_MODEL` macro to register your new class. The
implementation of the new class should be in namespace `aspect::HeatingModel`.

Specifically, your new class needs to implement the following basic interface:

```{code-block} c++
template <int dim>
    class aspect::HeatingModel::Interface
    {
      public:
        // compute heating terms used in the energy equation
        virtual
        void
        evaluate (const MaterialModel::MaterialModelInputs<dim> &material_model_inputs,
                  const MaterialModel::MaterialModelOutputs<dim> &material_model_outputs,
                  HeatingModel::HeatingModelOutputs &heating_model_outputs) const;

        // All the following functions are optional:
        virtual
        void
        initialize ();

        virtual
        void
        update ();

        // Functions used in dealing with run-time parameters
        static
        void
        declare_parameters (ParameterHandler &prm);

        virtual
        void
        parse_parameters (ParameterHandler &prm);

        // Allow the heating model to attach additional material model outputs in case it needs
        // them to compute the heating terms
        virtual
        void
        create_additional_material_model_outputs(MaterialModel::MaterialModelOutputs<dim> &) const;
    };
```

The main properties of the material are computed in the function `evaluate()`
that takes references to `MaterialModelInputs` and `MaterialModelOutputs`
objects and is supposed to fill the `HeatingModelOutputs` structure. As in the
material model, this function is handling lookups at an arbitrary number of
positions, so for each heating term (for example the heating source terms), a
`std::vector` is returned. The following members of `HeatingModelOutputs` need
to be filled:

```{code-block} c++
struct HeatingModelOutputs
{
       std::vector<double> heating_source_terms;
       std::vector<double> lhs_latent_heat_terms;

       // optional:
       std::vector<double> rates_of_temperature_change;
}
```

Heating source terms are terms on the right-hand side of the equations, such
as the adiabatic heating $\alpha T \left( \mathbf u \cdot \nabla p \right)$ in
equation {math:numref}`eq:temperature`3]. An example for a left-hand side heating term
is the temperature-derivative term
$\rho T \Delta S \frac{\partial X}{\partial T}$ that is part of latent heat
production (see equation {math:numref}`eq:temperature-reformulated`19]).[4] Rates of
temperature change[5] are used when the heating term is related to a reaction
process, happening on a faster time scale than the temperature advection. All
of these terms can depend on any of the material model inputs or outputs.
Implementations of `evaluate()` may of course choose to ignore dependencies on
any of these arguments.

The remaining functions are used in postprocessing as well as handling
run-time parameters. The exact meaning of these member functions is documented
in the [aspect::HeatingModel::Interface class
documentation][aspect::HeatingModel::Interface]. Note that some of the
functions listed above have a default implementation, as discussed on the
documentation page just mentioned.

Just like for material models, the functions `initialize()` and `update()` can
be implemented if desired (the default implementation does nothing) and are
useful if the heating model has an internal state. The function `initialize()`
is called once during the initialization of
ASPECT and can be used to allocate memory for the
heating model, initialize state, or read information from an external file.
The function `update()` is called at the beginning of every time step.
