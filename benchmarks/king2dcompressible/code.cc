/*
  Copyright (C) 2011 - 2022 by the authors of the ASPECT code.

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

#include <aspect/postprocess/interface.h>
#include <aspect/material_model/simple.h>
#include <aspect/boundary_temperature/interface.h>
#include <aspect/adiabatic_conditions/interface.h>
#include <aspect/simulator_access.h>
#include <aspect/global.h>
#include <aspect/geometry_model/box.h>

#include <deal.II/dofs/dof_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function_lib.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/vector_tools.h>



namespace aspect
{
  using namespace dealii;

  /**
   * This benchmark is from the article
   * @code
   * @article{KLKLZTTK10,
   *   title={A community benchmark for 2-{D} {C}artesian compressible
   *          convection in the {E}arth's mantle},
   *   author={King, Scott D and Lee, Changyeol and Van Keken, Peter E
   *           and Leng, Wei and Zhong, Shijie and Tan, Eh and Tosi, Nicola
   *           and Kameyama, Masanori C},
   *   journal={Geophysical Journal International},
   *   volume={180},
   *   number={1},
   *   pages={73--87},
   *   year={2010},
   *   publisher={Oxford University Press}
   * }
   * @endcode
   */


  namespace MaterialModel
  {

    template <int dim>
    class Material : public MaterialModel::Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:

        /**
        * Evaluate material properties.
        */
        virtual void evaluate(const MaterialModelInputs<dim> &in,
                              MaterialModelOutputs<dim> &out) const
        {

          for (unsigned int i=0; i < in.n_evaluation_points(); ++i)
            {
              const Point<dim> position = in.position[i];
              const double temperature = in.temperature[i];
              const double pressure = in.pressure[i];

              const double depth = 1.0-position(dim-1);

              out.viscosities[i] = ((Di==0.0)?1.0:Di)/Ra*exp(-b*(temperature- this->get_adiabatic_conditions().temperature(position))+c*depth);

              out.specific_heat[i] = reference_specific_heat;
              out.thermal_conductivities[i] = 1.0;
              out.thermal_expansion_coefficients[i] = (Di==0.0)?1.0:Di;

              double rho = reference_rho * exp(depth * Di/gamma);
              rho *= 1.0 - out.thermal_expansion_coefficients[i] * (temperature - this->get_adiabatic_conditions().temperature(position))
                     + (tala?0.0:1.0)*Di*gamma
                     *  (pressure - this->get_adiabatic_conditions().pressure(position));

              out.densities[i] = rho;

              out.compressibilities[i] = 0.0;
              out.entropy_derivative_pressure[i] = 0.0;
              out.entropy_derivative_temperature[i] = 0.0;
              // Change in composition due to chemical reactions at the
              // given positions. The term reaction_terms[i][c] is the
              // change in compositional field c at point i.
              for (unsigned int c=0; c<in.composition[i].size(); ++c)
                out.reaction_terms[i][c] = 0.0;
            }
        }

        virtual bool is_compressible () const
        {
          return true;
        }


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

        bool tala;

        /**
         * Parameters describing temperature and depth-dependence
         * of viscosity.
         */
        double b, c;

        /**
         * The surface density.
         */
        double reference_rho;

        /**
         * The nondimensional numbers (Dissipation number,
         * Rayleigh number, grueneisen parameter).
         */
        double Di, Ra, gamma;

        /**
         * The constant specific heat
         */
        double reference_specific_heat;
    };



    template <int dim>
    void
    Material<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("King model");
        {
          prm.declare_entry ("Reference density", "3300",
                             Patterns::Double (0),
                             "Reference density $\\rho_0$. "
                             "Units: \\si{\\kilogram\\per\\meter\\cubed}.");
          prm.declare_entry ("Ra", "1e4",
                             Patterns::Double (0),
                             "");
          prm.declare_entry ("Di", "0.0",
                             Patterns::Double (0),
                             "");
          prm.declare_entry ("gamma", "1.0",
                             Patterns::Double (0),
                             "");
          prm.declare_entry ("Reference specific heat", "1250",
                             Patterns::Double (0),
                             "The value of the specific heat $cp$. "
                             "Units: \\si{\\joule\\per\\kelvin\\per\\kilogram}.");
          prm.declare_entry ("b", "6.907755279",
                             Patterns::Double (0),
                             "");
          prm.declare_entry ("c", "0",
                             Patterns::Double (0),
                             "");
          prm.declare_entry ("Use TALA", "false",
                             Patterns::Bool (),
                             "Whether to use the TALA instead of the ALA "
                             "approximation.");

        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    Material<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("King model");
        {
          reference_rho   = prm.get_double ("Reference density");
          b               = prm.get_double ("b");
          c               = prm.get_double ("c");
          Di              = prm.get_double ("Di");
          Ra              = prm.get_double ("Ra");
          gamma           = prm.get_double ("gamma");

          tala            = prm.get_bool ("Use TALA");

          reference_specific_heat = prm.get_double ("Reference specific heat");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();

      // Declare dependencies on solution variables
      this->model_dependence.viscosity = NonlinearDependence::temperature;
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
    ASPECT_REGISTER_MATERIAL_MODEL(Material,
                                   "king material",
                                   "A simple compressible material model based on a benchmark "
                                   "from the paper of King et al. (2010). It uses the "
                                   "nondimensional numbers Di, Ra and gamma to define material "
                                   "properties.")
  }

}
