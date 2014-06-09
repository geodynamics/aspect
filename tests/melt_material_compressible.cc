#include <aspect/material_model/melt_interface.h>
#include <aspect/velocity_boundary_conditions/interface.h>
#include <aspect/simulator_access.h>
#include <aspect/global.h>

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
            out.viscosities[i] = 3.0;
            out.thermal_expansion_coefficients[i] = 1.0;
            out.specific_heat[i] = 1.0;
            out.thermal_conductivities[i] = 1.0;
            out.compressibilities[i] = 1.0;
            out.densities[i] = rho_0 * std::exp(-out.compressibilities[i] * in.position[i][1]);
          }
      }

      double rho_0 = 1.0;
      double rho_f_0 = 1.0;

      double permeability (const Point<dim> &position,
                           double porosity,
                           double shear_viscosity,
                           double compaction_viscosity)
      {
        double gravity = this->get_gravity_model().gravity_vector(position).norm();
        z = position[dim-1];
        x = position[0];

        double u_0    = x * (1.0-x);
        double u_0_x  = 1.0 - 2.0*x;
        double u_0_xx = -2.0;

        double A = u_0_xx * (compaction_viscosity + 1.0/3.0 * shear_viscosity);
        double B = 2.0 * gravity * (1.0 - porosity) * (rho_f_0 - rho_0);
        double C = u_0_x * 7.0/3.0 * shear_viscosity - u_0 * compaction_viscosity;
        double D = -B/2.0;

        return std::exp(- 4.0 * z + ((A/C + 2.0) * std::log(2.0*C*std::exp(2.0*z) - B)));
      }

      virtual void evaluate_with_melt(const typename MaterialModel::MeltInterface<dim>::MaterialModelInputs &in,
          typename MaterialModel::MeltInterface<dim>::MaterialModelOutputs &out) const
      {
        evaluate(in, out);
        const unsigned int porosity_idx = this->introspection().compositional_index_for_name("porosity");

        for (unsigned int i=0;i<in.position.size();++i)
          {
            double porosity = in.composition[i][porosity_idx];
            out.compaction_viscosities[i] = 1.0;
            out.fluid_viscosities[i]= 2.0;
            out.permeabilities[i]= permeability(in.position[i],porosity,out.viscosities[i],out.compaction_viscosities[i]);
            out.fluid_compressibilities[i] = 1.0;
            out.fluid_densities[i] = rho_f_0 * std::exp(-out.fluid_compressibilities[i] * in.position[i][1]);
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
}
