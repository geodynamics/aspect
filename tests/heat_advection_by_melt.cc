#include <aspect/fluid_pressure_boundary_conditions/interface.h>
#include <aspect/material_model/melt_global.h>
#include <aspect/simulator_access.h>
#include <aspect/global.h>
#include <aspect/melt.h>

namespace aspect
{

  template <int dim>
  class PressureBdry:
    public FluidPressureBoundaryConditions::Interface<dim>, public ::aspect::SimulatorAccess<dim>
  {
    public:
      virtual
      void fluid_pressure_gradient (
        const types::boundary_id /*boundary_indicator*/,
        const MaterialModel::MaterialModelInputs<dim> &material_model_inputs,
        const MaterialModel::MaterialModelOutputs<dim> &material_model_outputs,
        const std::vector<Tensor<1,dim> > &normal_vectors,
        std::vector<double> &fluid_pressure_gradient_outputs
      ) const
      {
        const MaterialModel::MeltOutputs<dim> *melt_outputs = material_model_outputs.template get_additional_output<MaterialModel::MeltOutputs<dim> >();
        Assert(melt_outputs != NULL, ExcMessage("Need MeltOutputs from the material model for shear heating with melt."));

        for (unsigned int q=0; q<fluid_pressure_gradient_outputs.size(); ++q)
          {
            const double density = material_model_outputs.densities[q];
            const double melt_density = melt_outputs->fluid_densities[q];
            const double porosity = material_model_inputs.composition[q][this->introspection().compositional_index_for_name("porosity")];

            const double bulk_density = (porosity * melt_density + (1.0 - porosity) * density);
            fluid_pressure_gradient_outputs[q] = this->get_gravity_model().gravity_vector(material_model_inputs.position[q]) * bulk_density * normal_vectors[q];
          }
      }

  };

  namespace MaterialModel
  {
    using namespace dealii;

    template <int dim>
    class Advection : public MaterialModel::MeltGlobal<dim>
    {
      public:

        virtual void evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
                              MaterialModel::MaterialModelOutputs<dim> &out) const;
    };


    template <int dim>
    void
    Advection<dim>::
    evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
             MaterialModel::MaterialModelOutputs<dim> &out) const
    {
      MeltGlobal<dim>::evaluate(in, out);

      // fill melt outputs if they exist
      MeltOutputs<dim> *melt_out = out.template get_additional_output<MeltOutputs<dim> >();
      for (unsigned int i=0; i < in.position.size(); ++i)
        {
          out.densities[i] = 2.0;
        }

      if (melt_out != NULL)
        {
          for (unsigned int i=0; i < in.position.size(); ++i)
            {
              melt_out->fluid_densities[i] = 1.0;
              melt_out->permeabilities[i] = 1.0;
            }
        }



    }

  }

}

// explicit instantiations
namespace aspect
{
  ASPECT_REGISTER_FLUID_PRESSURE_BOUNDARY_CONDITIONS(PressureBdry,
                                                     "PressureBdry",
                                                     "A fluid pressure boundary condition that prescribes the "
                                                     "gradient of the fluid pressure at the boundaries as "
                                                     "calculated in the analytical solution. ")

}

namespace aspect
{
  namespace MaterialModel
  {
    ASPECT_REGISTER_MATERIAL_MODEL(Advection,
                                   "melt global advection",
                                   "A simple material model that is like the "
                                   "'melt simple' model, but has a constant reaction term.")
  }
}


