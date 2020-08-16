/*
  Copyright (C) 2011 - 2020 by the authors of the ASPECT code.

  This file is part of ASPECT.

  ASPECT is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2, or (at your option)
  any later version.

  ASPECT is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with ASPECT; see the file LICENSE.  If not see
  <http://www.gnu.org/licenses/>.
*/

#include <aspect/material_model/interface.h>
#include <aspect/material_model/rheology/constant_viscosity.h>
#include <aspect/material_model/rheology/constant_viscosity_prefactors.h>
#include <aspect/simulator_access.h>

#include <iostream>
#include <cmath>


namespace aspect
{
  namespace MaterialModel
  {
    using namespace dealii;

    /**
     * A material model where all output parameters are hard-coded, but where
     * AdditionalNamedMaterialOutputs is filled with rheological properties
     * from several different rheological models.
     * Useful for checking the outputs of the various rheology plugins.
     *
     * @ingroup MaterialModels
     */
    template <int dim>
    class RheologyChecker : public MaterialModel::Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:

        virtual bool is_compressible () const;

        virtual double reference_viscosity () const;

        virtual void evaluate(const typename Interface<dim>::MaterialModelInputs &in,
                              typename Interface<dim>::MaterialModelOutputs &out) const;

        /**
         * @name Functions used in dealing with run-time parameters
         * @{
         */
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

        /**
         * Create additional named outputs
         */
        void
        create_additional_named_outputs (MaterialModel::MaterialModelOutputs<dim> &out) const override;

        /**
         * Make vector of the names of the additional outputs
         */
        std::vector<std::string>
        make_rheology_output_names() const;

        /**
         * Fill additional named outputs
         */
        void
        fill_rheology_outputs(const typename Interface<dim>::MaterialModelInputs &in,
                              NamedAdditionalMaterialOutputs<dim> *rheology_out) const;

        /**
         * @}
         */


      private:
        double ref_strain_rate;
        double ref_stress;

        Rheology::ConstantViscosity constant_viscosity;
        Rheology::ConstantViscosityPrefactors<dim> constant_viscosity_prefactors;
    };



    template <int dim>
    bool
    RheologyChecker<dim>::
    is_compressible () const
    {
      return false;
    }



    template <int dim>
    double
    RheologyChecker<dim>::
    reference_viscosity () const
    {
      return 5e24;
    }



    template <int dim>
    std::vector<std::string>
    RheologyChecker<dim>::make_rheology_output_names() const
    {
      const std::vector<std::string> names = {"constant viscosity strain rate",
                                              "constant viscosity strain rate derivative"
                                             };
      return names;
    }



    template <int dim>
    void
    RheologyChecker<dim>::fill_rheology_outputs(const typename Interface<dim>::MaterialModelInputs &in,
                                                NamedAdditionalMaterialOutputs<dim> *rheology_out) const
    {
      std::vector<std::vector<double>> outputs(rheology_out->get_names().size(), std::vector<double>(in.n_evaluation_points(), 0.));

      for (unsigned int i=0; i<in.n_evaluation_points(); ++i)
        {
          const std::pair strain_and_deriv_constant_viscosity = constant_viscosity_prefactors.compute_strain_rate_and_derivative (ref_stress, constant_viscosity.compute_viscosity(), 0);
          outputs[0][i] = strain_and_deriv_constant_viscosity.first;
          outputs[1][i] = strain_and_deriv_constant_viscosity.second;
        }

      rheology_out->output_values = outputs;
    }



    template <int dim>
    void
    RheologyChecker<dim>::
    evaluate(const typename Interface<dim>::MaterialModelInputs &in, typename Interface<dim>::MaterialModelOutputs &out ) const
    {
      for (unsigned int i=0; i<in.n_evaluation_points(); ++i)
        {
          out.viscosities[i] = 5e24;
          out.densities[i] = 3300. * (1.0 - 2.e-5 * (in.temperature[i] - 293));
          out.thermal_expansion_coefficients[i] = 2.e-5;
          out.specific_heat[i] = 1250.;
          out.thermal_conductivities[i] = 4.7;
          out.compressibilities[i] = 0.0;
        }

      if (MaterialModel::NamedAdditionalMaterialOutputs<dim> *rheology_out = out.template get_additional_output<NamedAdditionalMaterialOutputs<dim>>())
        fill_rheology_outputs(in, rheology_out);;
    }



    template <int dim>
    void
    RheologyChecker<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Rheology checker model");
        {
          prm.declare_entry ("Reference strain rate","1.0e-15",Patterns::Double (0.),
                             "Reference strain rate for calculation of viscosities. Units: \\si{\\per\\second}.");
          prm.declare_entry ("Reference stress","1.0e6",Patterns::Double (0.),
                             "Reference stress for calculation of strain rate and derivatives. Units: \\si{\\pascal}.");

          // Constant viscosity parameters
          Rheology::ConstantViscosity::declare_parameters(prm);

          // Constant viscosity prefactor parameters
          Rheology::ConstantViscosityPrefactors<dim>::declare_parameters(prm);
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    RheologyChecker<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Rheology checker model");
        {
          ref_strain_rate = prm.get_double("Reference strain rate");
          ref_stress = prm.get_double("Reference stress");

          const unsigned int n_phase_transitions_for_each_composition = 1;


          // Constant viscosity parameters
          constant_viscosity.parse_parameters(prm);

          // Constant viscosity prefactor parameters
          constant_viscosity_prefactors.initialize_simulator (this->get_simulator());
          constant_viscosity_prefactors.parse_parameters(prm);

        }
        prm.leave_subsection();
      }
      prm.leave_subsection();

      // Declare dependencies on solution variables
      this->model_dependence.viscosity = NonlinearDependence::none;
      this->model_dependence.density = NonlinearDependence::none;
      this->model_dependence.compressibility = NonlinearDependence::none;
      this->model_dependence.specific_heat = NonlinearDependence::none;
      this->model_dependence.thermal_conductivity = NonlinearDependence::none;
    }



    template <int dim>
    void
    RheologyChecker<dim>::create_additional_named_outputs (MaterialModel::MaterialModelOutputs<dim> &out) const
    {
      if (out.template get_additional_output<NamedAdditionalMaterialOutputs<dim> >() == nullptr)
        {
          const unsigned int n_points = out.n_evaluation_points();
          out.additional_outputs.push_back(
            std_cxx14::make_unique<MaterialModel::NamedAdditionalMaterialOutputs<dim>> (make_rheology_output_names(), n_points));
        }
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
    ASPECT_REGISTER_MATERIAL_MODEL(RheologyChecker,
                                   "rheology checker",
                                   "A material model that is like the ``simpler'' model but "
                                   "with additional outputs.")
  }
}
