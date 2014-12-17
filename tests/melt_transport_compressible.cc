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
      public:
      virtual bool
      viscosity_depends_on (const MaterialModel::NonlinearDependence::Dependence dependence) const
      {
        return false;
      }

      virtual bool
      density_depends_on (const MaterialModel::NonlinearDependence::Dependence dependence) const
      {
        if ((dependence & MaterialModel::NonlinearDependence::compositional_fields) != MaterialModel::NonlinearDependence::none)
          return true;
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
        return true;
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
        const unsigned int porosity_idx = this->introspection().compositional_index_for_name("porosity");
        for (unsigned int i=0;i<in.position.size();++i)
          {
            double porosity = in.composition[i][porosity_idx];
            out.viscosities[i] = 0.5 * std::exp(2.0 * in.position[i][0]);
            out.thermal_expansion_coefficients[i] = 0.0;
            out.specific_heat[i] = 1.0;
            out.thermal_conductivities[i] = 1.0;
            out.compressibilities[i] = 1.0 / (rho_s_0 * C);
            out.densities[i] = rho_s_0 * std::exp(-in.position[i][1]);
            for (unsigned int c=0;c<in.composition[i].size();++c)
              out.reaction_terms[i][c] = - rho_s_0 * B * D * std::exp(in.position[i][1]);
          }
      }

      virtual void evaluate_with_melt(const typename MaterialModel::MeltInterface<dim>::MaterialModelInputs &in,
          typename MaterialModel::MeltInterface<dim>::MaterialModelOutputs &out) const
      {
        evaluate(in, out);
        const unsigned int porosity_idx = this->introspection().compositional_index_for_name("porosity");

        for (unsigned int i=0;i<in.position.size();++i)
          {
            double porosity = in.composition[i][porosity_idx];
            out.compaction_viscosities[i] = xi_1 * std::exp(-in.position[i][1]) + 2.0/3.0 * std::exp(2.0 * in.position[i][0]) + xi_0;
            out.fluid_viscosities[i] = 1.0;
            out.permeabilities[i] = K_D_0 + 2.0 * B / E - rho_s_0 * B * D / E * (1.0/rho_s_0 - 1.0/rho_f_0) * std::exp(in.position[i][1]);
            out.fluid_compressibilities[i] = 1.0 / (rho_f_0 * C);
            out.fluid_densities[i] = rho_f_0 * std::exp(-in.position[i][1]);
          }

      }

      virtual void initialize ()
      {
          rho_s_0 = 1.2;
          rho_f_0 = 1.0;
          xi_0 = 1.0;
          xi_1 = 1.0;

          // A, B and C are constants from the velocity boundary conditions and gravity model
          // they have to be consistent!
          A = 0.1;
          B = -3.0/4.0 * A;
          C = 1.0;
          D = 0.3;
          E = - 3.0/4.0 * xi_0 * A + C * D *(rho_f_0 - rho_s_0);

          K_D_0 = 2.2;
      }


      private:
        double rho_s_0;
        double rho_f_0;
        double xi_0;
        double xi_1;
        double K_D_0;
        double A;
        double B;
        double C;
        double D;
        double E;


  };
  
  
}

// TODO add postprocessor for error calculation



// explicit instantiations
namespace aspect
{

    ASPECT_REGISTER_MATERIAL_MODEL(MeltMaterial,
                                   "Compressible melt material",
				   "")
}
