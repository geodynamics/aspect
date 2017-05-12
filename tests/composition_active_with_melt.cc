#include <aspect/material_model/interface.h>
#include <aspect/melt.h>
#include <deal.II/base/parameter_handler.h>

/**
 * This material model extends the simple material model
 * to include melt transport.
 */

namespace aspect
{
  using namespace dealii;

  template <int dim>
  class SimpleWithMelt:
    public MaterialModel::Interface<dim>
  {
    public:
      virtual bool is_compressible () const
      {
        return false;
      }

      virtual double reference_viscosity () const
      {
        return 0.2;
      }

      /**
       * Declare the parameters this class takes through input files.
       */
      static
      void
      declare_parameters (ParameterHandler &prm);

      /**
       * Read the parameters this class declares from the parameter file.
       */
      virtual
      void
      parse_parameters (ParameterHandler &prm);

      virtual void evaluate(const typename MaterialModel::Interface<dim>::MaterialModelInputs &in,
                            typename MaterialModel::Interface<dim>::MaterialModelOutputs &out) const
      {
        for (unsigned int i=0; i<in.position.size(); ++i)
          {
            out.viscosities[i] = 1.0;
            out.thermal_expansion_coefficients[i] = 0.01;
            out.specific_heat[i] = 1250.0;
            out.thermal_conductivities[i] = 1e-6;
            out.compressibilities[i] = 0.0;
            out.densities[i] = 1.0 * (1 - out.thermal_expansion_coefficients[i] * in.temperature[i])
                               + compositional_delta_rho*std::max(in.composition[i][0],0.0);

            // Pressure derivative of entropy at the given positions.
            out.entropy_derivative_pressure[i] = 0.0;
            // Temperature derivative of entropy at the given positions.
            out.entropy_derivative_temperature[i] = 0.0;
            for (unsigned int c=0; c<in.composition[i].size(); ++c)
              out.reaction_terms[i][c] = 0.0;
          }

        // fill melt outputs if they exist
        aspect::MaterialModel::MeltOutputs<dim> *melt_out = out.template get_additional_output<aspect::MaterialModel::MeltOutputs<dim> >();

        if (melt_out != NULL)
          {

            for (unsigned int i=0; i<in.position.size(); ++i)
              {
                double porosity = in.composition[i][0];
                melt_out->compaction_viscosities[i] = 1.0;
                melt_out->fluid_viscosities[i] = 1.0;
                melt_out->permeabilities[i] = 0.0; //1e-30*porosity * porosity;
                melt_out->fluid_density_gradients[i] = Tensor<1,dim>();
                melt_out->fluid_densities[i] = 1.0;
              }
          }
      }

    private:
      double compositional_delta_rho;

  };

  template <int dim>
  void
  SimpleWithMelt<dim>::declare_parameters (ParameterHandler &prm)
  {
    prm.enter_subsection("Material model");
    {
      prm.enter_subsection("Simple with melt");
      {
        prm.declare_entry ("Density differential for compositional field 1", "100.0",
                           Patterns::Double(),
                           "If compositional fields are used, then one would frequently want "
                           "to make the density depend on these fields. In this simple material "
                           "model, we make the following assumptions: if no compositional fields "
                           "are used in the current simulation, then the density is simply the usual "
                           "one with its linear dependence on the temperature. If there are compositional "
                           "fields, then the density only depends on the first one in such a way that "
                           "the density has an additional term of the kind $+\\Delta \\rho \\; c_1(\\mathbf x)$. "
                           "This parameter describes the value of $\\Delta \\rho$. Units: $kg/m^3/\\textrm{unit "
                           "change in composition}$.");
      }
      prm.leave_subsection();
    }
    prm.leave_subsection();
  }



  template <int dim>
  void
  SimpleWithMelt<dim>::parse_parameters (ParameterHandler &prm)
  {
    prm.enter_subsection("Material model");
    {
      prm.enter_subsection("Simple with melt");
      {
        compositional_delta_rho    = prm.get_double ("Density differential for compositional field 1");
      }
      prm.leave_subsection();
    }
    prm.leave_subsection();
  }
}


// explicit instantiations
namespace aspect
{
  ASPECT_REGISTER_MATERIAL_MODEL(SimpleWithMelt,
                                 "simple with melt",
                                 "A simple material model that is like the "
                                 "'simple' model, but with melt migration.")
}
