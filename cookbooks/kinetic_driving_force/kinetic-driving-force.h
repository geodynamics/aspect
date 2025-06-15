#ifndef COOKBOOKS_KINETIC_DRIVING_FORCE_KINETIC_DRIVING_FORCE_H_
#define COOKBOOKS_KINETIC_DRIVING_FORCE_KINETIC_DRIVING_FORCE_H_

#include <deal.II/base/patterns.h>
#include <deal.II/base/types.h>

#include <aspect/adiabatic_conditions/interface.h>
#include <aspect/geometry_model/interface.h>
#include <aspect/material_model/equation_of_state/interface.h>
#include <aspect/material_model/interface.h>
#include <aspect/simulator_access.h>
#include <aspect/utilities.h>

#include <vector>

namespace aspect
{
  namespace MaterialModel
  {
    template <int dim>
    class KineticDrivingForce : public MaterialModel::Interface<dim>,
      public ::aspect::SimulatorAccess<dim>
    {
      public:
        KineticDrivingForce();

        void
        initialize() override;

        bool
        is_compressible() const override;

        void
        evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
                 MaterialModel::MaterialModelOutputs<dim> &out) const override;

        static void
        declare_parameters(ParameterHandler &prm);

        void
        parse_parameters(ParameterHandler &prm) override;

        void
        create_additional_named_outputs(
          MaterialModel::MaterialModelOutputs<dim> &out) const override;

      private:
        // Data reader
        Utilities::AsciiDataProfile<dim> profile;

        // Cached column indices
        unsigned int rho_a_idx;
        unsigned int rho_b_idx;
        unsigned int alpha_a_idx;
        unsigned int alpha_b_idx;
        unsigned int beta_a_idx;
        unsigned int beta_b_idx;
        unsigned int cp_a_idx;
        unsigned int cp_b_idx;
        unsigned int dG_idx;
        unsigned int dS_idx;
        unsigned int dV_idx;

        // Rheological parameters
        double              viscosity;
        double              thermal_viscosity_exponent;
        std::vector<double> transition_depths;
        std::vector<double> viscosity_prefactors;
        double              minimum_viscosity;
        double              maximum_viscosity;

        // Material properties
        double k;
        bool   use_adiabatic_pressure_for_density;

        // Kinetic parameters
        double Q_kinetic_prefactor;
        double reaction_cutoff_temperature;
        bool   use_temperature_cutoff;
    };
  } // namespace MaterialModel
} // namespace aspect

#endif // COOKBOOKS_KINETIC_DRIVING_FORCE_KINETIC_DRIVING_FORCE_H_
