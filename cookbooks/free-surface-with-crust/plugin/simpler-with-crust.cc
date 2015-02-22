/*
  Copyright (C) 2011 - 2015 by the authors of the ASPECT code.

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
  along with ASPECT; see the file doc/COPYING.  If not see
  <http://www.gnu.org/licenses/>.
*/


#include <aspect/material_model/interface.h>
#include <aspect/simulator_access.h>

#include <deal.II/base/parameter_handler.h>

#include <iostream>
#include <cmath>


namespace aspect
{
  namespace MaterialModel
  {
    using namespace dealii;

    /**
     * A material model similar to the "simpler" material model, but where the
     * viscosity has two different values dependent on whether we are above or
     * below a line at a certain z-value, i.e., with a
     * crust.
     *
     * @ingroup MaterialModels
     */
    template <int dim>
    class SimplerWithCrust : public Interface<dim>
    {
      public:

        virtual bool
        viscosity_depends_on (const NonlinearDependence::Dependence dependence) const;

        virtual bool
        density_depends_on (const NonlinearDependence::Dependence dependence) const;

        virtual bool
        compressibility_depends_on (const NonlinearDependence::Dependence dependence) const;

        virtual bool
        specific_heat_depends_on (const NonlinearDependence::Dependence dependence) const;

        virtual bool
        thermal_conductivity_depends_on (const NonlinearDependence::Dependence dependence) const;

        virtual bool is_compressible () const;

        virtual double reference_viscosity () const;

        virtual double reference_density () const;

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
        double reference_rho;
        double reference_T;
        double eta_L;
        double eta_U;
        double jump_height;
        double thermal_alpha;
        double reference_specific_heat;
        double k_value;
    };


    template <int dim>
    bool
    SimplerWithCrust<dim>::
    viscosity_depends_on (const NonlinearDependence::Dependence dependence) const
    {
      return false;
    }


    template <int dim>
    bool
    SimplerWithCrust<dim>::
    density_depends_on (const NonlinearDependence::Dependence dependence) const
    {
      return false;
    }

    template <int dim>
    bool
    SimplerWithCrust<dim>::
    compressibility_depends_on (const NonlinearDependence::Dependence) const
    {
      return false;
    }

    template <int dim>
    bool
    SimplerWithCrust<dim>::
    specific_heat_depends_on (const NonlinearDependence::Dependence) const
    {
      return false;
    }

    template <int dim>
    bool
    SimplerWithCrust<dim>::
    thermal_conductivity_depends_on (const NonlinearDependence::Dependence dependence) const
    {
      return false;
    }


    template <int dim>
    bool
    SimplerWithCrust<dim>::
    is_compressible () const
    {
      return false;
    }

    template <int dim>
    double
    SimplerWithCrust<dim>::
    reference_viscosity () const
    {
      return eta_L;
    }



    template <int dim>
    double
    SimplerWithCrust<dim>::
    reference_density () const
    {
      return reference_rho;
    }

    template <int dim>
    void
    SimplerWithCrust<dim>::
    evaluate(const typename Interface<dim>::MaterialModelInputs &in, typename Interface<dim>::MaterialModelOutputs &out ) const
    {
      for (unsigned int i=0; i<in.position.size(); ++i)
        {
          const double z = in.position[i][1];
          if (z>jump_height)
            out.viscosities[i] = eta_U;
          else
            out.viscosities[i] = eta_L;

          out.densities[i] = reference_rho * (1.0 - thermal_alpha * (in.temperature[i] - reference_T));
          out.thermal_expansion_coefficients[i] = thermal_alpha;
          out.specific_heat[i] = reference_specific_heat;
          out.thermal_conductivities[i] = k_value;
          out.compressibilities[i] = 0.0;
        }
    }


    template <int dim>
    void
    SimplerWithCrust<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Simpler with crust model");
        {
          prm.declare_entry ("Reference density", "3300",
                             Patterns::Double (0),
                             "Reference density $\\rho_0$. Units: $kg/m^3$.");
          prm.declare_entry ("Reference temperature", "293",
                             Patterns::Double (0),
                             "The reference temperature $T_0$. The reference temperature is used "
                             "in the density formula. Units: $K$.");
          prm.declare_entry ("Lower viscosity", "5e24",
                             Patterns::Double (0),
                             "The value of the viscosity $\\eta$L. Units: $kg/m/s$.");
          prm.declare_entry ("Upper viscosity", "5e24",
                             Patterns::Double (0),
                             "The value of the viscosity in the top section $\\eta$U. Units: $kg/m/s$.");
          prm.declare_entry ("Jump height", "100000",
                             Patterns::Double (0),
                             "The height at which the viscosity changes. Units: m.");
          prm.declare_entry ("Thermal conductivity", "4.7",
                             Patterns::Double (0),
                             "The value of the thermal conductivity $k$. "
                             "Units: $W/m/K$.");
          prm.declare_entry ("Reference specific heat", "1250",
                             Patterns::Double (0),
                             "The value of the specific heat $cp$. "
                             "Units: $J/kg/K$.");
          prm.declare_entry ("Thermal expansion coefficient", "2e-5",
                             Patterns::Double (0),
                             "The value of the thermal expansion coefficient $\\beta$. "
                             "Units: $1/K$.");

        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    SimplerWithCrust<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Simpler with crust model");
        {
          reference_rho              = prm.get_double ("Reference density");
          reference_T                = prm.get_double ("Reference temperature");
          eta_L                      = prm.get_double ("Lower viscosity");
          eta_U                      = prm.get_double ("Upper viscosity");
          jump_height                = prm.get_double ("Jump height");
          k_value                    = prm.get_double ("Thermal conductivity");
          reference_specific_heat    = prm.get_double ("Reference specific heat");
          thermal_alpha              = prm.get_double ("Thermal expansion coefficient");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
    ASPECT_REGISTER_MATERIAL_MODEL(SimplerWithCrust,
                                   "simpler with crust",
                                   "A material model that is like the ``simpler'' model but "
                                   "has a jump in the viscosity.")
  }
}
