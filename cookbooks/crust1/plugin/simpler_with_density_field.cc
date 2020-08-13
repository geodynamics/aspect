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
#include <aspect/simulator_access.h>
#include <deal.II/base/parameter_handler.h>
#include <aspect/utilities.h>


namespace aspect
{
  namespace MaterialModel
  {
    using namespace dealii;

    /**
     * A material model similar to the "simpler" material model, but where the
     * viscosity has two different values dependent on whether we are above or
     * below a line at a certain z-value, i.e., representing a crustal layer.
     *
     * @ingroup MaterialModels
     */
    template <int dim>
    class SimplerWithDensityField : public Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        /**
         * Return whether the model is compressible or not. This model is
         * incompressible.
         */
        virtual bool is_compressible () const;

        /**
         * Return a reference value typical of the viscosities that appear in
         * this model. This value is not actually used in the material
         * description itself, but is used in scaling variables to the same
         * numerical order of magnitude when solving linear systems.
         * Specifically, the reference viscosity appears in the factor scaling
         * the pressure against the velocity. It is also used in computing
         * dimension-less quantities. You may want to take a look at the
         * Kronbichler, Heister, Bangerth 2012 paper that describes the
         * design of ASPECT for a description of this pressure scaling.
         */
        virtual double reference_viscosity () const;

        /**
         * Function to compute the material properties in @p out given the
         * inputs in @p in.
         */
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
         * @}
         */

      private:
        double eta;
        double thermal_expansion_coefficient;
        double reference_specific_heat;
        double thermal_conductivity;
        unsigned int density_idx;
    };



    template <int dim>
    bool
    SimplerWithDensityField<dim>::
    is_compressible () const
    {
      return false;
    }



    template <int dim>
    double
    SimplerWithDensityField<dim>::
    reference_viscosity () const
    {
      return eta;
    }



    template <int dim>
    void
    SimplerWithDensityField<dim>::
    evaluate(const typename Interface<dim>::MaterialModelInputs &in,
             typename Interface<dim>::MaterialModelOutputs &out ) const
    {
      for (unsigned int i=0; i<in.n_evaluation_points(); ++i)
        {
          out.viscosities[i] = eta;
          out.densities[i] = std::max(in.composition[i][density_idx], 1.);
          out.thermal_expansion_coefficients[i] = thermal_expansion_coefficient;
          out.specific_heat[i] = reference_specific_heat;
          out.thermal_conductivities[i] = thermal_conductivity;
          out.compressibilities[i] = 0.0;
          out.entropy_derivative_pressure[i] = 0.0;
          out.entropy_derivative_temperature[i] = 0.0;
          for (unsigned int c=0; c<in.composition[i].size(); ++c)
            out.reaction_terms[i][c] = 0.0;
        }
    }



    template <int dim>
    void
    SimplerWithDensityField<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Simpler with density field model");
        {
          prm.declare_entry ("Viscosity", "1e23",
                             Patterns::Double (0),
                             "The value of the viscosity $\\eta$. "
                             "Units: \\si{\\pascal\\second}.");
          prm.declare_entry ("Thermal conductivity", "4.7",
                             Patterns::Double (0),
                             "The value of the thermal conductivity $k$. "
                             "Units: \\si{\\watt\\per\\meter\\per\\kelvin}.");
          prm.declare_entry ("Reference specific heat", "1250",
                             Patterns::Double (0),
                             "The value of the specific heat capacity $c_p$. "
                             "Units: \\si{\\joule\\per\\kelvin\\per\\kilogram}.");
          prm.declare_entry ("Thermal expansion coefficient", "2e-5",
                             Patterns::Double (0),
                             "The value of the thermal expansion coefficient $\\alpha$. "
                             "Units: \\si{\\per\\kelvin}.");
          prm.declare_entry ("Compositional field index for density", "1",
                             Patterns::Integer (0),
                             "The index of the compositional field corresponding to "
                             "the density. Indexing starts from 1."
                             "Units: [].");

        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    SimplerWithDensityField<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Simpler with density field model");
        {
          eta                      = prm.get_double ("Viscosity");
          thermal_conductivity       = prm.get_double ("Thermal conductivity");
          reference_specific_heat    = prm.get_double ("Reference specific heat");
          thermal_expansion_coefficient = prm.get_double ("Thermal expansion coefficient");
          density_idx = prm.get_integer ("Compositional field index for density") - 1;
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
  }
}

// explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
    ASPECT_REGISTER_MATERIAL_MODEL(SimplerWithDensityField,
                                   "simpler with density field",
                                   "A material model that is like the ``simpler'' model but "
                                   "uses a compositional field to assign densities.")
  }
}
