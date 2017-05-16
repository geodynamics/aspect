#include <aspect/material_model/interface.h>
#include <aspect/boundary_velocity/interface.h>
#include <aspect/simulator_access.h>
#include <aspect/global.h>
#include <aspect/melt.h>

#include <deal.II/dofs/dof_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function_lib.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/vector_tools.h>


namespace aspect
{
  template <int dim>
  class MeltMaterial:
    public MaterialModel::MeltInterface<dim>, public ::aspect::SimulatorAccess<dim>
  {
      virtual bool is_compressible () const
      {
        return false;
      }

      virtual double reference_viscosity () const
      {
        return 1.0;
      }

      virtual double reference_darcy_coefficient () const
      {
        const double porosity = 0.01;
        const double permeability = 1.0 * std::pow(porosity,3) * std::pow(1.0-porosity,2);
        return permeability / 0.1;
      }

      virtual void evaluate(const typename MaterialModel::Interface<dim>::MaterialModelInputs &in,
                            typename MaterialModel::Interface<dim>::MaterialModelOutputs &out) const
      {

        for (unsigned int i=0; i<in.position.size(); ++i)
          {
            out.viscosities[i] = 1.0;
            out.densities[i] = 1.0 + (in.composition[i][0]);
            out.thermal_expansion_coefficients[i] = 1.0;
            out.specific_heat[i] = 1.0;
            out.thermal_conductivities[i] = 1.0;
            out.compressibilities[i] = 0.0;
            for (unsigned int c=0; c<in.composition[i].size(); ++c)
              out.reaction_terms[i][c] = 0.0;
          }

        // fill melt outputs if they exist
        aspect::MaterialModel::MeltOutputs<dim> *melt_out = out.template get_additional_output<aspect::MaterialModel::MeltOutputs<dim> >();

        if (melt_out != NULL)
          {
            const unsigned int porosity_idx = this->introspection().compositional_index_for_name("porosity");
            for (unsigned int i=0; i<in.position.size(); ++i)
              {
                const double porosity = in.composition[i][porosity_idx];

                melt_out->compaction_viscosities[i] = 100.0;
                melt_out->fluid_viscosities[i]=0.1;
                melt_out->permeabilities[i]=1.0 * std::pow(porosity,3) * std::pow(1.0-porosity,2);
                melt_out->fluid_densities[i]=.1;
                melt_out->fluid_density_gradients[i] = Tensor<1,dim>();
              }
          }
      }

  };


}



// explicit instantiations
namespace aspect
{

  ASPECT_REGISTER_MATERIAL_MODEL(MeltMaterial,
                                 "MeltMaterial",
                                 "")
}
