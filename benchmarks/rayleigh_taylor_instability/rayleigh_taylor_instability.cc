/*
  Copyright (C) 2011 - 2023 by the authors of the ASPECT code.

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

#include <aspect/simulator.h>
#include <aspect/material_model/simple.h>
#include <aspect/boundary_velocity/interface.h>
#include <aspect/simulator_access.h>
#include <aspect/global.h>
#include <aspect/utilities.h>

#include <deal.II/dofs/dof_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function_lib.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/vector_tools.h>

namespace aspect
{
  namespace RTinstabilityBenchmark
  {
    /**
     * @note This benchmark only talks about the flow field, not about a
     * temperature field. All quantities related to the temperature are
     * therefore set to zero in the implementation of this class.
     *
     * @ingroup MaterialModels
     */
    template <int dim>
    class RTinstabilityMaterial : public MaterialModel::Interface<dim>
    {
      public:
        /**
         * @name Physical parameters used in the basic equations
         * @{
         */
        void evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
                      MaterialModel::MaterialModelOutputs<dim> &out) const override
        {
          for (unsigned int i=0; i < in.n_evaluation_points(); ++i)
            {
              const Point<dim> &pos = in.position[i];

              if ( pos[1]< 256e3 + amplitude*std::cos(pos[0]/lambda*2.0*numbers::PI))
                {
                  out.viscosities[i] = eta2;
                  out.densities[i] = rho2;
                }
              else
                {
                  out.viscosities[i] = eta1;
                  out.densities[i] = rho1;
                }
              out.compressibilities[i] = 0;
              out.specific_heat[i] = 0;
              out.thermal_expansion_coefficients[i] = 0;
              out.thermal_conductivities[i] = 0.0;
            }
        }
        /**
         * @}
         */

        /**
         * @name Qualitative properties one can ask a material model
         * @{
         */

        /**
         * Return whether the model is compressible or not.
         * Incompressibility does not necessarily imply that the density is
         * constant; rather, it may still depend on temperature or pressure.
         * In the current context, compressibility means whether we should
         * solve the continuity equation as $\nabla \cdot (\rho \mathbf u)=0$
         * (compressible Stokes) or as $\nabla \cdot \mathbf{u}=0$
         * (incompressible Stokes).
         */
        bool is_compressible () const override;
        /**
         * @}
         */
        /**
         * Declare the parameters this class takes through input files.
         */
        static
        void
        declare_parameters (ParameterHandler &prm)
        {
          prm.enter_subsection("Material model");
          {
            prm.enter_subsection("RTinstability");
            {
              prm.declare_entry ("eta1", "1e21",
                                 Patterns::Double (0),
                                 "Viscosity of the mantle.");
              prm.declare_entry ("eta2", "1e23",
                                 Patterns::Double (0),
                                 "Viscosity in the inclusion.");
              prm.declare_entry ("rho1", "3200",
                                 Patterns::Double (0),
                                 "density of the mantle.");
              prm.declare_entry ("rho2", "3232",
                                 Patterns::Double (0),
                                 "density in the inclusion.");
              prm.declare_entry ("amplitude", "1000",
                                 Patterns::Double (0),
                                 "perturbation amplitude.");
              prm.declare_entry ("lambda", "256e3",
                                 Patterns::Double (0),
                                 "perturbation wavelength.");
            }
            prm.leave_subsection();
          }
          prm.leave_subsection();
        }

        /**
         * Read the parameters this class declares from the parameter file.
         */
        void
        parse_parameters (ParameterHandler &prm) override
        {
          prm.enter_subsection("Material model");
          {
            prm.enter_subsection("RTinstability");
            {
              eta1 = prm.get_double ("eta1");
              eta2 = prm.get_double ("eta2");
              rho1 = prm.get_double ("rho1");
              rho2 = prm.get_double ("rho2");
              amplitude = prm.get_double ("amplitude");
              lambda = prm.get_double ("lambda");
            }
            prm.leave_subsection();
          }
          prm.leave_subsection();

          // Declare dependencies on solution variables
          this->model_dependence.viscosity = MaterialModel::NonlinearDependence::none;
          this->model_dependence.density = MaterialModel::NonlinearDependence::none;
          this->model_dependence.compressibility = MaterialModel::NonlinearDependence::none;
          this->model_dependence.specific_heat = MaterialModel::NonlinearDependence::none;
          this->model_dependence.thermal_conductivity = MaterialModel::NonlinearDependence::none;
        }

      private:
        double eta1;
        double eta2;
        double rho1;
        double rho2;
        double amplitude;
        double lambda;

    };


    template <int dim>
    bool
    RTinstabilityMaterial<dim>::
    is_compressible () const
    {
      return false;
    }

  }
}



// explicit instantiations
namespace aspect
{
  namespace RTinstabilityBenchmark
  {
    ASPECT_REGISTER_MATERIAL_MODEL(RTinstabilityMaterial,
                                   "RTinstabilityMaterial",
                                   "A material model that corresponds to the `RTinstability' benchmark. "
                                   "See the manual for more information.")
  }
}
