#include <aspect/material_model/melt_interface.h>
#include <aspect/velocity_boundary_conditions/interface.h>
#include <aspect/simulator_access.h>
#include <aspect/global.h>
#include <aspect/fluid_pressure_boundary_conditions/interface.h>
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
      virtual bool
      viscosity_depends_on (const MaterialModel::NonlinearDependence::Dependence dependence) const
      {
        return false;
      }

      virtual bool
      density_depends_on (const MaterialModel::NonlinearDependence::Dependence dependence) const
      {
        return false;
      }


      virtual bool
      compressibility_depends_on (const MaterialModel::NonlinearDependence::Dependence dependence) const
      {
        return false;
      }


      virtual bool
      specific_heat_depends_on (const MaterialModel::NonlinearDependence::Dependence dependence) const
      {
        return false;
      }


      virtual bool
      thermal_conductivity_depends_on (const MaterialModel::NonlinearDependence::Dependence dependence) const
      {
        return false;
      }

      virtual bool is_compressible () const
      {
        return false;
      }

      virtual double reference_viscosity () const
      {
        return 1.0;
      }

      virtual double reference_density () const
      {
        return 1.0;
      }
      virtual void evaluate(const typename MaterialModel::Interface<dim>::MaterialModelInputs &in,
                 typename MaterialModel::Interface<dim>::MaterialModelOutputs &out) const
      {

        for (unsigned int i=0;i<in.position.size();++i)
          {
            out.viscosities[i] = 3.0/4.0;
            out.densities[i] = 3.0;
            out.thermal_expansion_coefficients[i] = 1.0;
            out.specific_heat[i] = 1.0;
            out.thermal_conductivities[i] = 1.0;
            out.compressibilities[i] = 0.0;
            for (unsigned int c=0;c<in.composition[i].size();++c)
              out.reaction_terms[i][c] = 0.0;
          }
      }

      virtual void evaluate_with_melt(const typename MaterialModel::MeltInterface<dim>::MaterialModelInputs &in,
          typename MaterialModel::MeltInterface<dim>::MaterialModelOutputs &out) const
      {
        evaluate(in, out);
        for (unsigned int i=0;i<in.position.size();++i)
          {
            out.compaction_viscosities[i] = 1.0;
            out.fluid_viscosities[i]= 2.0;
            out.permeabilities[i]= 1.0 + std::pow(in.position[i][1],2);
            out.fluid_densities[i]= 1.0;
            out.fluid_compressibilities[i] = 0.0;
          }

      }

  };
  
  template <int dim>
  class PressureBdry:

      public FluidPressureBoundaryConditions::Interface<dim>
  {
    public:
      virtual
      void fluid_pressure_gradient (
        const typename MaterialModel::MeltInterface<dim>::MaterialModelInputs &material_model_inputs,
        const typename MaterialModel::MeltInterface<dim>::MaterialModelOutputs &material_model_outputs,
        std::vector<Tensor<1,dim> > & output
      ) const
        {
          for (unsigned int q=0; q<output.size(); ++q)
            {
              Tensor<1,dim> gradient;
              gradient[0] = 0.0;
              gradient[1] = 0.0;
              output[q] = gradient;
            }
        }
  };
  
}



// explicit instantiations
namespace aspect
{

    ASPECT_REGISTER_MATERIAL_MODEL(MeltMaterial,
                                   "MeltMaterial2",
				   "")

    ASPECT_REGISTER_FLUID_PRESSURE_BOUNDARY_CONDITIONS(PressureBdry,
                                                       "PressureBdry",
                                                       "A fluid pressure boundary condition that prescribes the "
                                                       "gradient of the fluid pressure at the boundaries as "
                                                       "calculated in the analytical solution. ")
}
