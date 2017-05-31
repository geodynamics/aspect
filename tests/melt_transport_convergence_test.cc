#include <aspect/material_model/interface.h>
#include <aspect/boundary_velocity/interface.h>
#include <aspect/boundary_fluid_pressure/interface.h>
#include <aspect/postprocess/interface.h>
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
  class TestMeltMaterial:
    public MaterialModel::MeltInterface<dim>, public ::aspect::SimulatorAccess<dim>
  {
    public:
      virtual bool is_compressible () const
      {
        return false;
      }

      virtual double reference_darcy_coefficient () const
      {
        const double porosity = 0.01;
        const double permeability = porosity * porosity;
        return permeability / 1.0;
      }

      virtual double reference_viscosity () const
      {
        return 1.0;
      }

      virtual void evaluate(const typename MaterialModel::Interface<dim>::MaterialModelInputs &in,
                            typename MaterialModel::Interface<dim>::MaterialModelOutputs &out) const
      {
        const unsigned int porosity_idx = this->introspection().compositional_index_for_name("porosity");
        double c = 1.0;
        for (unsigned int i=0; i<in.position.size(); ++i)
          {
            const double porosity = in.composition[i][porosity_idx];
            const double x = in.position[i](0);
            const double z = in.position[i](1);
            out.viscosities[i] = exp(c*porosity);
            out.thermal_expansion_coefficients[i] = 0.0;
            out.specific_heat[i] = 1.0;
            out.thermal_conductivities[i] = 1.0;
            out.compressibilities[i] = 0.0;
            out.densities[i] = 1.0;
            // This is the RHS we need to use for the manufactured solution.
            // We calculated it by subtracting the first term of the RHS (\Nabla (K_D \rho_f g)) from the LHS
            // we computed using our analytical solution.
            out.reaction_terms[i][porosity_idx] =
              0.1000000000e1 + 0.1000000000e0 * (0.1e0 - 0.5e-1 * atan(pow(x + 0.20e1 * z, 0.2e1))) * (0.5e-1 * (0.2e1 * x + 0.40e1 * z) / (0.1e1 + pow(x + 0.20e1 * z, 0.4e1)) * exp(0.1e0 - 0.5e-1 * atan(pow(x + 0.20e1 * z, 0.2e1))) / (-0.9e0 - 0.5e-1 * atan(pow(x + 0.20e1 * z, 0.2e1))) - 0.5e-1 * exp(0.1e0 - 0.5e-1 * atan(pow(x + 0.20e1 * z, 0.2e1))) * pow(-0.9e0 - 0.5e-1 * atan(pow(x + 0.20e1 * z, 0.2e1)), -0.2e1) * (0.2e1 * x + 0.40e1 * z) / (0.1e1 + pow(x + 0.20e1 * z, 0.4e1))) * (0.2e1 * x + 0.40e1 * z) / (0.1e1 + pow(x + 0.20e1 * z, 0.4e1)) - 0.1000000000e1 * pow(0.1e0 - 0.5e-1 * atan(pow(x + 0.20e1 * z, 0.2e1)), 0.2e1) * (0.10e0 / (0.1e1 + pow(x + 0.20e1 * z, 0.4e1)) * exp(0.1e0 - 0.5e-1 * atan(pow(x + 0.20e1 * z, 0.2e1))) / (-0.9e0 - 0.5e-1 * atan(pow(x + 0.20e1 * z, 0.2e1))) - 0.20e0 * (0.2e1 * x + 0.40e1 * z) * pow(0.1e1 + pow(x + 0.20e1 * z, 0.4e1), -0.2e1) * exp(0.1e0 - 0.5e-1 * atan(pow(x + 0.20e1 * z, 0.2e1))) / (-0.9e0 - 0.5e-1 * atan(pow(x + 0.20e1 * z, 0.2e1))) * pow(x + 0.20e1 * z, 0.3e1) - 0.25e-2 * pow(0.2e1 * x + 0.40e1 * z, 0.2e1) * pow(0.1e1 + pow(x + 0.20e1 * z, 0.4e1), -0.2e1) * exp(0.1e0 - 0.5e-1 * atan(pow(x + 0.20e1 * z, 0.2e1))) / (-0.9e0 - 0.5e-1 * atan(pow(x + 0.20e1 * z, 0.2e1))) + 0.50e-2 * pow(0.2e1 * x + 0.40e1 * z, 0.2e1) * pow(0.1e1 + pow(x + 0.20e1 * z, 0.4e1), -0.2e1) * exp(0.1e0 - 0.5e-1 * atan(pow(x + 0.20e1 * z, 0.2e1))) * pow(-0.9e0 - 0.5e-1 * atan(pow(x + 0.20e1 * z, 0.2e1)), -0.2e1) - 0.50e-2 * exp(0.1e0 - 0.5e-1 * atan(pow(x + 0.20e1 * z, 0.2e1))) * pow(-0.9e0 - 0.5e-1 * atan(pow(x + 0.20e1 * z, 0.2e1)), -0.3e1) * pow(0.2e1 * x + 0.40e1 * z, 0.2e1) * pow(0.1e1 + pow(x + 0.20e1 * z, 0.4e1), -0.2e1) - 0.10e0 * exp(0.1e0 - 0.5e-1 * atan(pow(x + 0.20e1 * z, 0.2e1))) * pow(-0.9e0 - 0.5e-1 * atan(pow(x + 0.20e1 * z, 0.2e1)), -0.2e1) / (0.1e1 + pow(x + 0.20e1 * z, 0.4e1)) + 0.20e0 * exp(0.1e0 - 0.5e-1 * atan(pow(x + 0.20e1 * z, 0.2e1))) * pow(-0.9e0 - 0.5e-1 * atan(pow(x + 0.20e1 * z, 0.2e1)), -0.2e1) * (0.2e1 * x + 0.40e1 * z) * pow(0.1e1 + pow(x + 0.20e1 * z, 0.4e1), -0.2e1) * pow(x + 0.20e1 * z, 0.3e1)) + 0.1000000000e0 * (0.1e0 - 0.5e-1 * atan(pow(x + 0.20e1 * z, 0.2e1))) * (0.5e-1 * (0.40e1 * x + 0.800e1 * z) / (0.1e1 + pow(x + 0.20e1 * z, 0.4e1)) * exp(0.1e0 - 0.5e-1 * atan(pow(x + 0.20e1 * z, 0.2e1))) / (-0.9e0 - 0.5e-1 * atan(pow(x + 0.20e1 * z, 0.2e1))) - 0.5e-1 * exp(0.1e0 - 0.5e-1 * atan(pow(x + 0.20e1 * z, 0.2e1))) * pow(-0.9e0 - 0.5e-1 * atan(pow(x + 0.20e1 * z, 0.2e1)), -0.2e1) * (0.40e1 * x + 0.800e1 * z) / (0.1e1 + pow(x + 0.20e1 * z, 0.4e1)) + 0.5000000000e-1 * exp(0.5000000000e-1 * z)) * (0.40e1 * x + 0.800e1 * z) / (0.1e1 + pow(x + 0.20e1 * z, 0.4e1)) - 0.1000000000e1 * pow(0.1e0 - 0.5e-1 * atan(pow(x + 0.20e1 * z, 0.2e1)), 0.2e1) * (0.4000e0 / (0.1e1 + pow(x + 0.20e1 * z, 0.4e1)) * exp(0.1e0 - 0.5e-1 * atan(pow(x + 0.20e1 * z, 0.2e1))) / (-0.9e0 - 0.5e-1 * atan(pow(x + 0.20e1 * z, 0.2e1))) - 0.400e0 * (0.40e1 * x + 0.800e1 * z) * pow(0.1e1 + pow(x + 0.20e1 * z, 0.4e1), -0.2e1) * exp(0.1e0 - 0.5e-1 * atan(pow(x + 0.20e1 * z, 0.2e1))) / (-0.9e0 - 0.5e-1 * atan(pow(x + 0.20e1 * z, 0.2e1))) * pow(x + 0.20e1 * z, 0.3e1) - 0.25e-2 * pow(0.40e1 * x + 0.800e1 * z, 0.2e1) * pow(0.1e1 + pow(x + 0.20e1 * z, 0.4e1), -0.2e1) * exp(0.1e0 - 0.5e-1 * atan(pow(x + 0.20e1 * z, 0.2e1))) / (-0.9e0 - 0.5e-1 * atan(pow(x + 0.20e1 * z, 0.2e1))) + 0.50e-2 * pow(0.40e1 * x + 0.800e1 * z, 0.2e1) * pow(0.1e1 + pow(x + 0.20e1 * z, 0.4e1), -0.2e1) * exp(0.1e0 - 0.5e-1 * atan(pow(x + 0.20e1 * z, 0.2e1))) * pow(-0.9e0 - 0.5e-1 * atan(pow(x + 0.20e1 * z, 0.2e1)), -0.2e1) - 0.50e-2 * exp(0.1e0 - 0.5e-1 * atan(pow(x + 0.20e1 * z, 0.2e1))) * pow(-0.9e0 - 0.5e-1 * atan(pow(x + 0.20e1 * z, 0.2e1)), -0.3e1) * pow(0.40e1 * x + 0.800e1 * z, 0.2e1) * pow(0.1e1 + pow(x + 0.20e1 * z, 0.4e1), -0.2e1) - 0.4000e0 * exp(0.1e0 - 0.5e-1 * atan(pow(x + 0.20e1 * z, 0.2e1))) * pow(-0.9e0 - 0.5e-1 * atan(pow(x + 0.20e1 * z, 0.2e1)), -0.2e1) / (0.1e1 + pow(x + 0.20e1 * z, 0.4e1)) + 0.400e0 * exp(0.1e0 - 0.5e-1 * atan(pow(x + 0.20e1 * z, 0.2e1))) * pow(-0.9e0 - 0.5e-1 * atan(pow(x + 0.20e1 * z, 0.2e1)), -0.2e1) * (0.40e1 * x + 0.800e1 * z) * pow(0.1e1 + pow(x + 0.20e1 * z, 0.4e1), -0.2e1) * pow(x + 0.20e1 * z, 0.3e1) + 0.2500000000e-2 * exp(0.5000000000e-1 * z)) - 0.5000000000e-1 * (0.1e0 - 0.5e-1 * atan(pow(x + 0.20e1 * z, 0.2e1))) * (0.1166666666e0 * (0.2e1 * x + 0.40e1 * z) / (0.1e1 + pow(x + 0.20e1 * z, 0.4e1)) * exp(0.1e0 - 0.5e-1 * atan(pow(x + 0.20e1 * z, 0.2e1))) + 0.500e-1 * (0.40e1 * x + 0.800e1 * z) / (0.1e1 + pow(x + 0.20e1 * z, 0.4e1)) * exp(0.1e0 - 0.5e-1 * atan(pow(x + 0.20e1 * z, 0.2e1))) * (cos(z) + 0.1e1) + 0.100e1 * exp(0.1e0 - 0.5e-1 * atan(pow(x + 0.20e1 * z, 0.2e1))) * sin(z) + 0.5e-1 * (0.2e1 * x + 0.40e1 * z) / (0.1e1 + pow(x + 0.20e1 * z, 0.4e1)) * exp(0.1e0 - 0.5e-1 * atan(pow(x + 0.20e1 * z, 0.2e1))) / (-0.9e0 - 0.5e-1 * atan(pow(x + 0.20e1 * z, 0.2e1))) - 0.5e-1 * exp(0.1e0 - 0.5e-1 * atan(pow(x + 0.20e1 * z, 0.2e1))) * pow(-0.9e0 - 0.5e-1 * atan(pow(x + 0.20e1 * z, 0.2e1)), -0.2e1) * (0.2e1 * x + 0.40e1 * z) / (0.1e1 + pow(x + 0.20e1 * z, 0.4e1))) / (0.95e0 + 0.25e-1 * atan(pow(x + 0.20e1 * z, 0.2e1))) * (0.2e1 * x + 0.40e1 * z) / (0.1e1 + pow(x + 0.20e1 * z, 0.4e1)) + 0.5000000000e0 * pow(0.1e0 - 0.5e-1 * atan(pow(x + 0.20e1 * z, 0.2e1)), 0.2e1) * (0.2333333332e0 / (0.1e1 + pow(x + 0.20e1 * z, 0.4e1)) * exp(0.1e0 - 0.5e-1 * atan(pow(x + 0.20e1 * z, 0.2e1))) - 0.4666666664e0 * (0.2e1 * x + 0.40e1 * z) * pow(0.1e1 + pow(x + 0.20e1 * z, 0.4e1), -0.2e1) * exp(0.1e0 - 0.5e-1 * atan(pow(x + 0.20e1 * z, 0.2e1))) * pow(x + 0.20e1 * z, 0.3e1) - 0.5833333330e-2 * pow(0.2e1 * x + 0.40e1 * z, 0.2e1) * pow(0.1e1 + pow(x + 0.20e1 * z, 0.4e1), -0.2e1) * exp(0.1e0 - 0.5e-1 * atan(pow(x + 0.20e1 * z, 0.2e1))) + 0.20000e0 / (0.1e1 + pow(x + 0.20e1 * z, 0.4e1)) * exp(0.1e0 - 0.5e-1 * atan(pow(x + 0.20e1 * z, 0.2e1))) * (cos(z) + 0.1e1) - 0.2000e0 * (0.40e1 * x + 0.800e1 * z) * pow(0.1e1 + pow(x + 0.20e1 * z, 0.4e1), -0.2e1) * exp(0.1e0 - 0.5e-1 * atan(pow(x + 0.20e1 * z, 0.2e1))) * (cos(z) + 0.1e1) * pow(x + 0.20e1 * z, 0.3e1) - 0.2500e-2 * (0.40e1 * x + 0.800e1 * z) * pow(0.1e1 + pow(x + 0.20e1 * z, 0.4e1), -0.2e1) * (0.2e1 * x + 0.40e1 * z) * exp(0.1e0 - 0.5e-1 * atan(pow(x + 0.20e1 * z, 0.2e1))) * (cos(z) + 0.1e1) - 0.500e-1 * (0.2e1 * x + 0.40e1 * z) / (0.1e1 + pow(x + 0.20e1 * z, 0.4e1)) * exp(0.1e0 - 0.5e-1 * atan(pow(x + 0.20e1 * z, 0.2e1))) * sin(z) + 0.10e0 / (0.1e1 + pow(x + 0.20e1 * z, 0.4e1)) * exp(0.1e0 - 0.5e-1 * atan(pow(x + 0.20e1 * z, 0.2e1))) / (-0.9e0 - 0.5e-1 * atan(pow(x + 0.20e1 * z, 0.2e1))) - 0.20e0 * (0.2e1 * x + 0.40e1 * z) * pow(0.1e1 + pow(x + 0.20e1 * z, 0.4e1), -0.2e1) * exp(0.1e0 - 0.5e-1 * atan(pow(x + 0.20e1 * z, 0.2e1))) / (-0.9e0 - 0.5e-1 * atan(pow(x + 0.20e1 * z, 0.2e1))) * pow(x + 0.20e1 * z, 0.3e1) - 0.25e-2 * pow(0.2e1 * x + 0.40e1 * z, 0.2e1) * pow(0.1e1 + pow(x + 0.20e1 * z, 0.4e1), -0.2e1) * exp(0.1e0 - 0.5e-1 * atan(pow(x + 0.20e1 * z, 0.2e1))) / (-0.9e0 - 0.5e-1 * atan(pow(x + 0.20e1 * z, 0.2e1))) + 0.50e-2 * pow(0.2e1 * x + 0.40e1 * z, 0.2e1) * pow(0.1e1 + pow(x + 0.20e1 * z, 0.4e1), -0.2e1) * exp(0.1e0 - 0.5e-1 * atan(pow(x + 0.20e1 * z, 0.2e1))) * pow(-0.9e0 - 0.5e-1 * atan(pow(x + 0.20e1 * z, 0.2e1)), -0.2e1) - 0.50e-2 * exp(0.1e0 - 0.5e-1 * atan(pow(x + 0.20e1 * z, 0.2e1))) * pow(-0.9e0 - 0.5e-1 * atan(pow(x + 0.20e1 * z, 0.2e1)), -0.3e1) * pow(0.2e1 * x + 0.40e1 * z, 0.2e1) * pow(0.1e1 + pow(x + 0.20e1 * z, 0.4e1), -0.2e1) - 0.10e0 * exp(0.1e0 - 0.5e-1 * atan(pow(x + 0.20e1 * z, 0.2e1))) * pow(-0.9e0 - 0.5e-1 * atan(pow(x + 0.20e1 * z, 0.2e1)), -0.2e1) / (0.1e1 + pow(x + 0.20e1 * z, 0.4e1)) + 0.20e0 * exp(0.1e0 - 0.5e-1 * atan(pow(x + 0.20e1 * z, 0.2e1))) * pow(-0.9e0 - 0.5e-1 * atan(pow(x + 0.20e1 * z, 0.2e1)), -0.2e1) * (0.2e1 * x + 0.40e1 * z) * pow(0.1e1 + pow(x + 0.20e1 * z, 0.4e1), -0.2e1) * pow(x + 0.20e1 * z, 0.3e1)) / (0.95e0 + 0.25e-1 * atan(pow(x + 0.20e1 * z, 0.2e1))) - 0.1250000000e-1 * pow(0.1e0 - 0.5e-1 * atan(pow(x + 0.20e1 * z, 0.2e1)), 0.2e1) * (0.1166666666e0 * (0.2e1 * x + 0.40e1 * z) / (0.1e1 + pow(x + 0.20e1 * z, 0.4e1)) * exp(0.1e0 - 0.5e-1 * atan(pow(x + 0.20e1 * z, 0.2e1))) + 0.500e-1 * (0.40e1 * x + 0.800e1 * z) / (0.1e1 + pow(x + 0.20e1 * z, 0.4e1)) * exp(0.1e0 - 0.5e-1 * atan(pow(x + 0.20e1 * z, 0.2e1))) * (cos(z) + 0.1e1) + 0.100e1 * exp(0.1e0 - 0.5e-1 * atan(pow(x + 0.20e1 * z, 0.2e1))) * sin(z) + 0.5e-1 * (0.2e1 * x + 0.40e1 * z) / (0.1e1 + pow(x + 0.20e1 * z, 0.4e1)) * exp(0.1e0 - 0.5e-1 * atan(pow(x + 0.20e1 * z, 0.2e1))) / (-0.9e0 - 0.5e-1 * atan(pow(x + 0.20e1 * z, 0.2e1))) - 0.5e-1 * exp(0.1e0 - 0.5e-1 * atan(pow(x + 0.20e1 * z, 0.2e1))) * pow(-0.9e0 - 0.5e-1 * atan(pow(x + 0.20e1 * z, 0.2e1)), -0.2e1) * (0.2e1 * x + 0.40e1 * z) / (0.1e1 + pow(x + 0.20e1 * z, 0.4e1))) * pow(0.95e0 + 0.25e-1 * atan(pow(x + 0.20e1 * z, 0.2e1)), -0.2e1) * (0.2e1 * x + 0.40e1 * z) / (0.1e1 + pow(x + 0.20e1 * z, 0.4e1)) - 0.5000000000e-1 * (0.1e0 - 0.5e-1 * atan(pow(x + 0.20e1 * z, 0.2e1))) * (0.1666666667e-1 * (0.40e1 * x + 0.800e1 * z) / (0.1e1 + pow(x + 0.20e1 * z, 0.4e1)) * exp(0.1e0 - 0.5e-1 * atan(pow(x + 0.20e1 * z, 0.2e1))) + 0.500e-1 * (0.2e1 * x + 0.40e1 * z) / (0.1e1 + pow(x + 0.20e1 * z, 0.4e1)) * exp(0.1e0 - 0.5e-1 * atan(pow(x + 0.20e1 * z, 0.2e1))) * (cos(z) + 0.1e1) + 0.5e-1 * (0.40e1 * x + 0.800e1 * z) / (0.1e1 + pow(x + 0.20e1 * z, 0.4e1)) * exp(0.1e0 - 0.5e-1 * atan(pow(x + 0.20e1 * z, 0.2e1))) / (-0.9e0 - 0.5e-1 * atan(pow(x + 0.20e1 * z, 0.2e1))) - 0.5e-1 * exp(0.1e0 - 0.5e-1 * atan(pow(x + 0.20e1 * z, 0.2e1))) * pow(-0.9e0 - 0.5e-1 * atan(pow(x + 0.20e1 * z, 0.2e1)), -0.2e1) * (0.40e1 * x + 0.800e1 * z) / (0.1e1 + pow(x + 0.20e1 * z, 0.4e1)) + 0.5000000000e-1 * exp(0.5000000000e-1 * z)) / (0.95e0 + 0.25e-1 * atan(pow(x + 0.20e1 * z, 0.2e1))) * (0.40e1 * x + 0.800e1 * z) / (0.1e1 + pow(x + 0.20e1 * z, 0.4e1)) + 0.5000000000e0 * pow(0.1e0 - 0.5e-1 * atan(pow(x + 0.20e1 * z, 0.2e1)), 0.2e1) * (0.1333333334e0 / (0.1e1 + pow(x + 0.20e1 * z, 0.4e1)) * exp(0.1e0 - 0.5e-1 * atan(pow(x + 0.20e1 * z, 0.2e1))) - 0.1333333334e0 * (0.40e1 * x + 0.800e1 * z) * pow(0.1e1 + pow(x + 0.20e1 * z, 0.4e1), -0.2e1) * exp(0.1e0 - 0.5e-1 * atan(pow(x + 0.20e1 * z, 0.2e1))) * pow(x + 0.20e1 * z, 0.3e1) - 0.8333333335e-3 * pow(0.40e1 * x + 0.800e1 * z, 0.2e1) * pow(0.1e1 + pow(x + 0.20e1 * z, 0.4e1), -0.2e1) * exp(0.1e0 - 0.5e-1 * atan(pow(x + 0.20e1 * z, 0.2e1))) + 0.20000e0 / (0.1e1 + pow(x + 0.20e1 * z, 0.4e1)) * exp(0.1e0 - 0.5e-1 * atan(pow(x + 0.20e1 * z, 0.2e1))) * (cos(z) + 0.1e1) - 0.40000e0 * (0.2e1 * x + 0.40e1 * z) * pow(0.1e1 + pow(x + 0.20e1 * z, 0.4e1), -0.2e1) * exp(0.1e0 - 0.5e-1 * atan(pow(x + 0.20e1 * z, 0.2e1))) * (cos(z) + 0.1e1) * pow(x + 0.20e1 * z, 0.3e1) - 0.2500e-2 * (0.40e1 * x + 0.800e1 * z) * pow(0.1e1 + pow(x + 0.20e1 * z, 0.4e1), -0.2e1) * (0.2e1 * x + 0.40e1 * z) * exp(0.1e0 - 0.5e-1 * atan(pow(x + 0.20e1 * z, 0.2e1))) * (cos(z) + 0.1e1) - 0.500e-1 * (0.2e1 * x + 0.40e1 * z) / (0.1e1 + pow(x + 0.20e1 * z, 0.4e1)) * exp(0.1e0 - 0.5e-1 * atan(pow(x + 0.20e1 * z, 0.2e1))) * sin(z) + 0.4000e0 / (0.1e1 + pow(x + 0.20e1 * z, 0.4e1)) * exp(0.1e0 - 0.5e-1 * atan(pow(x + 0.20e1 * z, 0.2e1))) / (-0.9e0 - 0.5e-1 * atan(pow(x + 0.20e1 * z, 0.2e1))) - 0.400e0 * (0.40e1 * x + 0.800e1 * z) * pow(0.1e1 + pow(x + 0.20e1 * z, 0.4e1), -0.2e1) * exp(0.1e0 - 0.5e-1 * atan(pow(x + 0.20e1 * z, 0.2e1))) / (-0.9e0 - 0.5e-1 * atan(pow(x + 0.20e1 * z, 0.2e1))) * pow(x + 0.20e1 * z, 0.3e1) - 0.25e-2 * pow(0.40e1 * x + 0.800e1 * z, 0.2e1) * pow(0.1e1 + pow(x + 0.20e1 * z, 0.4e1), -0.2e1) * exp(0.1e0 - 0.5e-1 * atan(pow(x + 0.20e1 * z, 0.2e1))) / (-0.9e0 - 0.5e-1 * atan(pow(x + 0.20e1 * z, 0.2e1))) + 0.50e-2 * pow(0.40e1 * x + 0.800e1 * z, 0.2e1) * pow(0.1e1 + pow(x + 0.20e1 * z, 0.4e1), -0.2e1) * exp(0.1e0 - 0.5e-1 * atan(pow(x + 0.20e1 * z, 0.2e1))) * pow(-0.9e0 - 0.5e-1 * atan(pow(x + 0.20e1 * z, 0.2e1)), -0.2e1) - 0.50e-2 * exp(0.1e0 - 0.5e-1 * atan(pow(x + 0.20e1 * z, 0.2e1))) * pow(-0.9e0 - 0.5e-1 * atan(pow(x + 0.20e1 * z, 0.2e1)), -0.3e1) * pow(0.40e1 * x + 0.800e1 * z, 0.2e1) * pow(0.1e1 + pow(x + 0.20e1 * z, 0.4e1), -0.2e1) - 0.4000e0 * exp(0.1e0 - 0.5e-1 * atan(pow(x + 0.20e1 * z, 0.2e1))) * pow(-0.9e0 - 0.5e-1 * atan(pow(x + 0.20e1 * z, 0.2e1)), -0.2e1) / (0.1e1 + pow(x + 0.20e1 * z, 0.4e1)) + 0.400e0 * exp(0.1e0 - 0.5e-1 * atan(pow(x + 0.20e1 * z, 0.2e1))) * pow(-0.9e0 - 0.5e-1 * atan(pow(x + 0.20e1 * z, 0.2e1)), -0.2e1) * (0.40e1 * x + 0.800e1 * z) * pow(0.1e1 + pow(x + 0.20e1 * z, 0.4e1), -0.2e1) * pow(x + 0.20e1 * z, 0.3e1) + 0.2500000000e-2 * exp(0.5000000000e-1 * z)) / (0.95e0 + 0.25e-1 * atan(pow(x + 0.20e1 * z, 0.2e1))) - 0.1250000000e-1 * pow(0.1e0 - 0.5e-1 * atan(pow(x + 0.20e1 * z, 0.2e1)), 0.2e1) * (0.1666666667e-1 * (0.40e1 * x + 0.800e1 * z) / (0.1e1 + pow(x + 0.20e1 * z, 0.4e1)) * exp(0.1e0 - 0.5e-1 * atan(pow(x + 0.20e1 * z, 0.2e1))) + 0.500e-1 * (0.2e1 * x + 0.40e1 * z) / (0.1e1 + pow(x + 0.20e1 * z, 0.4e1)) * exp(0.1e0 - 0.5e-1 * atan(pow(x + 0.20e1 * z, 0.2e1))) * (cos(z) + 0.1e1) + 0.5e-1 * (0.40e1 * x + 0.800e1 * z) / (0.1e1 + pow(x + 0.20e1 * z, 0.4e1)) * exp(0.1e0 - 0.5e-1 * atan(pow(x + 0.20e1 * z, 0.2e1))) / (-0.9e0 - 0.5e-1 * atan(pow(x + 0.20e1 * z, 0.2e1))) - 0.5e-1 * exp(0.1e0 - 0.5e-1 * atan(pow(x + 0.20e1 * z, 0.2e1))) * pow(-0.9e0 - 0.5e-1 * atan(pow(x + 0.20e1 * z, 0.2e1)), -0.2e1) * (0.40e1 * x + 0.800e1 * z) / (0.1e1 + pow(x + 0.20e1 * z, 0.4e1)) + 0.5000000000e-1 * exp(0.5000000000e-1 * z)) * pow(0.95e0 + 0.25e-1 * atan(pow(x + 0.20e1 * z, 0.2e1)), -0.2e1) * (0.40e1 * x + 0.800e1 * z) / (0.1e1 + pow(x + 0.20e1 * z, 0.4e1));

          }

        // fill melt outputs if they exist
        aspect::MaterialModel::MeltOutputs<dim> *melt_out = out.template get_additional_output<aspect::MaterialModel::MeltOutputs<dim> >();

        if (melt_out != NULL)
          {
            double c = 1.0;
            const unsigned int porosity_idx = this->introspection().compositional_index_for_name("porosity");

            for (unsigned int i=0; i<in.position.size(); ++i)
              {
                double porosity = in.composition[i][porosity_idx];
                melt_out->compaction_viscosities[i] = exp(c*porosity);
                melt_out->fluid_viscosities[i] = 1.0;
                melt_out->permeabilities[i] = porosity * porosity;
                melt_out->fluid_density_gradients[i] = Tensor<1,dim>();
                melt_out->fluid_densities[i] = 0.5;
              }
          }
      }

  };




  template <int dim>
  class RefFunction : public Function<dim>
  {
    public:
      RefFunction () : Function<dim>(2*dim+5) {}
      virtual void vector_value (const Point< dim >   &p,
                                 Vector< double >   &values) const
      {
        double x = p(0);
        double z = p(1);

        values[0] = x + sin(z);
        values[1] = x;
        values[2] = -exp(0.1e0 - 0.5e-1 * atan(pow(x + 0.20e1 * z, 0.2e1))) / (-0.9e0 - 0.5e-1 * atan(pow(x + 0.20e1 * z, 0.2e1))) + exp(0.5000000000e-1 * z);
        values[3] = -exp(0.1e0 - 0.5e-1 * atan(pow(x + 0.20e1 * z, 0.2e1)));
        values[4] = x + sin(z) - (0.1e0 - 0.5e-1 * atan(pow(x + 0.20e1 * z, 0.2e1))) * (0.5e-1 * (0.2e1 * x + 0.40e1 * z) / (0.1e1 + pow(x + 0.20e1 * z, 0.4e1)) * exp(0.1e0 - 0.5e-1 * atan(pow(x + 0.20e1 * z, 0.2e1))) / (-0.9e0 - 0.5e-1 * atan(pow(x + 0.20e1 * z, 0.2e1))) - 0.5e-1 * exp(0.1e0 - 0.5e-1 * atan(pow(x + 0.20e1 * z, 0.2e1))) * pow(-0.9e0 - 0.5e-1 * atan(pow(x + 0.20e1 * z, 0.2e1)), -0.2e1) * (0.2e1 * x + 0.40e1 * z) / (0.1e1 + pow(x + 0.20e1 * z, 0.4e1)) - 0.5e0 * (0.1166666666e0 * (0.2e1 * x + 0.40e1 * z) / (0.1e1 + pow(x + 0.20e1 * z, 0.4e1)) * exp(0.1e0 - 0.5e-1 * atan(pow(x + 0.20e1 * z, 0.2e1))) + 0.500e-1 * (0.40e1 * x + 0.800e1 * z) / (0.1e1 + pow(x + 0.20e1 * z, 0.4e1)) * exp(0.1e0 - 0.5e-1 * atan(pow(x + 0.20e1 * z, 0.2e1))) * (cos(z) + 0.1e1) + 0.100e1 * exp(0.1e0 - 0.5e-1 * atan(pow(x + 0.20e1 * z, 0.2e1))) * sin(z) + 0.5e-1 * (0.2e1 * x + 0.40e1 * z) / (0.1e1 + pow(x + 0.20e1 * z, 0.4e1)) * exp(0.1e0 - 0.5e-1 * atan(pow(x + 0.20e1 * z, 0.2e1))) / (-0.9e0 - 0.5e-1 * atan(pow(x + 0.20e1 * z, 0.2e1))) - 0.5e-1 * exp(0.1e0 - 0.5e-1 * atan(pow(x + 0.20e1 * z, 0.2e1))) * pow(-0.9e0 - 0.5e-1 * atan(pow(x + 0.20e1 * z, 0.2e1)), -0.2e1) * (0.2e1 * x + 0.40e1 * z) / (0.1e1 + pow(x + 0.20e1 * z, 0.4e1))) / (0.95e0 + 0.25e-1 * atan(pow(x + 0.20e1 * z, 0.2e1))));
        values[5] = x - (0.1e0 - 0.5e-1 * atan(pow(x + 0.20e1 * z, 0.2e1))) * (0.5e-1 * (0.40e1 * x + 0.800e1 * z) / (0.1e1 + pow(x + 0.20e1 * z, 0.4e1)) * exp(0.1e0 - 0.5e-1 * atan(pow(x + 0.20e1 * z, 0.2e1))) / (-0.9e0 - 0.5e-1 * atan(pow(x + 0.20e1 * z, 0.2e1))) - 0.5e-1 * exp(0.1e0 - 0.5e-1 * atan(pow(x + 0.20e1 * z, 0.2e1))) * pow(-0.9e0 - 0.5e-1 * atan(pow(x + 0.20e1 * z, 0.2e1)), -0.2e1) * (0.40e1 * x + 0.800e1 * z) / (0.1e1 + pow(x + 0.20e1 * z, 0.4e1)) + 0.5000000000e-1 * exp(0.5000000000e-1 * z) - 0.5e0 * (0.1666666667e-1 * (0.40e1 * x + 0.800e1 * z) / (0.1e1 + pow(x + 0.20e1 * z, 0.4e1)) * exp(0.1e0 - 0.5e-1 * atan(pow(x + 0.20e1 * z, 0.2e1))) + 0.500e-1 * (0.2e1 * x + 0.40e1 * z) / (0.1e1 + pow(x + 0.20e1 * z, 0.4e1)) * exp(0.1e0 - 0.5e-1 * atan(pow(x + 0.20e1 * z, 0.2e1))) * (cos(z) + 0.1e1) + 0.5e-1 * (0.40e1 * x + 0.800e1 * z) / (0.1e1 + pow(x + 0.20e1 * z, 0.4e1)) * exp(0.1e0 - 0.5e-1 * atan(pow(x + 0.20e1 * z, 0.2e1))) / (-0.9e0 - 0.5e-1 * atan(pow(x + 0.20e1 * z, 0.2e1))) - 0.5e-1 * exp(0.1e0 - 0.5e-1 * atan(pow(x + 0.20e1 * z, 0.2e1))) * pow(-0.9e0 - 0.5e-1 * atan(pow(x + 0.20e1 * z, 0.2e1)), -0.2e1) * (0.40e1 * x + 0.800e1 * z) / (0.1e1 + pow(x + 0.20e1 * z, 0.4e1)) + 0.5000000000e-1 * exp(0.5000000000e-1 * z)) / (0.95e0 + 0.25e-1 * atan(pow(x + 0.20e1 * z, 0.2e1))));
        values[6] = exp(0.5000000000e-1 * z);
        values[7] = 0;
        values[8] = 0.1e0 - 0.5e-1 * atan(pow(x + 0.20e1 * z, 0.2e1));

      }
  };



  /**
    * A postprocessor that evaluates the accuracy of the solution
    * by using the L2 norm.
    */
  template <int dim>
  class ConvergenceMeltPostprocessor : public Postprocess::Interface<dim>, public ::aspect::SimulatorAccess<dim>
  {
    public:
      /**
       * Generate graphical output from the current solution.
       */
      virtual
      std::pair<std::string,std::string>
      execute (TableHandler &statistics);

  };

  template <int dim>
  std::pair<std::string,std::string>
  ConvergenceMeltPostprocessor<dim>::execute (TableHandler &statistics)
  {
    RefFunction<dim> ref_func;
    const QGauss<dim> quadrature_formula (this->introspection().polynomial_degree.velocities +2);
    const unsigned int n_total_comp = this->introspection().n_components;

    Vector<float> cellwise_errors_u (this->get_triangulation().n_active_cells());
    Vector<float> cellwise_errors_p_f (this->get_triangulation().n_active_cells());
    Vector<float> cellwise_errors_p_c (this->get_triangulation().n_active_cells());
    Vector<float> cellwise_errors_u_f (this->get_triangulation().n_active_cells());
    Vector<float> cellwise_errors_p (this->get_triangulation().n_active_cells());
    Vector<float> cellwise_errors_porosity (this->get_triangulation().n_active_cells());

    ComponentSelectFunction<dim> comp_u(std::pair<unsigned int, unsigned int>(0,dim),
                                        n_total_comp);
    ComponentSelectFunction<dim> comp_p_f(dim, n_total_comp);
    ComponentSelectFunction<dim> comp_p_c(dim+1, n_total_comp);
    ComponentSelectFunction<dim> comp_u_f(std::pair<unsigned int, unsigned int>(dim+2,dim+2+
                                                                                dim),
                                          n_total_comp);
    ComponentSelectFunction<dim> comp_p(dim+2+dim, n_total_comp);
    ComponentSelectFunction<dim> comp_porosity(dim+2+dim+2, n_total_comp);

    VectorTools::integrate_difference (this->get_mapping(),this->get_dof_handler(),
                                       this->get_solution(),
                                       ref_func,
                                       cellwise_errors_u,
                                       quadrature_formula,
                                       VectorTools::L2_norm,
                                       &comp_u);
    VectorTools::integrate_difference (this->get_mapping(),this->get_dof_handler(),
                                       this->get_solution(),
                                       ref_func,
                                       cellwise_errors_p_f,
                                       quadrature_formula,
                                       VectorTools::L2_norm,
                                       &comp_p_f);
    VectorTools::integrate_difference (this->get_mapping(),this->get_dof_handler(),
                                       this->get_solution(),
                                       ref_func,
                                       cellwise_errors_p,
                                       quadrature_formula,
                                       VectorTools::L2_norm,
                                       &comp_p);
    VectorTools::integrate_difference (this->get_mapping(),this->get_dof_handler(),
                                       this->get_solution(),
                                       ref_func,
                                       cellwise_errors_p_c,
                                       quadrature_formula,
                                       VectorTools::L2_norm,
                                       &comp_p_c);
    VectorTools::integrate_difference (this->get_mapping(),this->get_dof_handler(),
                                       this->get_solution(),
                                       ref_func,
                                       cellwise_errors_porosity,
                                       quadrature_formula,
                                       VectorTools::L2_norm,
                                       &comp_porosity);
    VectorTools::integrate_difference (this->get_mapping(),this->get_dof_handler(),
                                       this->get_solution(),
                                       ref_func,
                                       cellwise_errors_u_f,
                                       quadrature_formula,
                                       VectorTools::L2_norm,
                                       &comp_u_f);

    const double u_l2 = std::sqrt(Utilities::MPI::sum(cellwise_errors_u.norm_sqr(),this->get_mpi_communicator()));
    const double p_l2 = std::sqrt(Utilities::MPI::sum(cellwise_errors_p.norm_sqr(),this->get_mpi_communicator()));
    const double p_f_l2 = std::sqrt(Utilities::MPI::sum(cellwise_errors_p_f.norm_sqr(),this->get_mpi_communicator()));
    const double p_c_l2 = std::sqrt(Utilities::MPI::sum(cellwise_errors_p_c.norm_sqr(),this->get_mpi_communicator()));
    const double poro_l2 = std::sqrt(Utilities::MPI::sum(cellwise_errors_porosity.norm_sqr(),this->get_mpi_communicator()));
    const double u_f_l2 = std::sqrt(Utilities::MPI::sum(cellwise_errors_u_f.norm_sqr(),this->get_mpi_communicator()));


    std::ostringstream os;
    os << std::scientific << u_l2
       << ", " << p_l2
       << ", " << p_f_l2
       << ", " << p_c_l2
       << ", " << poro_l2
       << ", " << u_f_l2;

    return std::make_pair("Errors u_L2, p_L2, p_f_L2, p_c_L2, porosity_L2, u_f_L2:", os.str());
  }


  template <int dim>
  class PressureBdry:

    public BoundaryFluidPressure::Interface<dim>
  {
    public:
      virtual
      void fluid_pressure_gradient (
        const types::boundary_id boundary_indicator,
        const MaterialModel::MaterialModelInputs<dim> &material_model_inputs,
        const MaterialModel::MaterialModelOutputs<dim> &material_model_outputs,
        const std::vector<Tensor<1,dim> > &normal_vectors,
        std::vector<double> &output
      ) const
      {
        for (unsigned int q=0; q<output.size(); ++q)
          {
            const double x = material_model_inputs.position[q][0];
            const double z = material_model_inputs.position[q][1];
            Tensor<1,dim> gradient;
            gradient[0] = 0.5e-1 * (0.2e1 * x + 0.40e1 * z) / (0.1e1 + pow(x + 0.20e1 * z, 0.4e1)) * exp(0.1e0 - 0.5e-1 * atan(pow(x + 0.20e1 * z, 0.2e1))) / (-0.9e0 - 0.5e-1 * atan(pow(x + 0.20e1 * z, 0.2e1))) - 0.5e-1 * exp(0.1e0 - 0.5e-1 * atan(pow(x + 0.20e1 * z, 0.2e1))) * pow(-0.9e0 - 0.5e-1 * atan(pow(x + 0.20e1 * z, 0.2e1)), -0.2e1) * (0.2e1 * x + 0.40e1 * z) / (0.1e1 + pow(x + 0.20e1 * z, 0.4e1));
            gradient[1] = 0.5e-1 * (0.40e1 * x + 0.800e1 * z) / (0.1e1 + pow(x + 0.20e1 * z, 0.4e1)) * exp(0.1e0 - 0.5e-1 * atan(pow(x + 0.20e1 * z, 0.2e1))) / (-0.9e0 - 0.5e-1 * atan(pow(x + 0.20e1 * z, 0.2e1))) - 0.5e-1 * exp(0.1e0 - 0.5e-1 * atan(pow(x + 0.20e1 * z, 0.2e1))) * pow(-0.9e0 - 0.5e-1 * atan(pow(x + 0.20e1 * z, 0.2e1)), -0.2e1) * (0.40e1 * x + 0.800e1 * z) / (0.1e1 + pow(x + 0.20e1 * z, 0.4e1)) + 0.5000000000e-1 * exp(0.5000000000e-1 * z);
            output[q] = gradient * normal_vectors[q];
          }
      }



  };

}

// explicit instantiations
namespace aspect
{

  ASPECT_REGISTER_MATERIAL_MODEL(TestMeltMaterial,
                                 "test melt material",
                                 "")


  ASPECT_REGISTER_POSTPROCESSOR(ConvergenceMeltPostprocessor,
                                "melt error calculation",
                                "A postprocessor that compares the numerical solution to the analytical "
                                "solution derived for incompressible melt transport in a 2D box as described "
                                "in the manuscript and reports the error.")

  ASPECT_REGISTER_BOUNDARY_FLUID_PRESSURE_MODEL(PressureBdry,
                                                "PressureBdry",
                                                "A fluid pressure boundary condition that prescribes the "
                                                "gradient of the fluid pressure at the boundaries as "
                                                "calculated in the analytical solution. ")

}
