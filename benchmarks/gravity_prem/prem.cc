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
#include <aspect/utilities.h>

namespace aspect
{
  namespace GravityPremBenchmark
  {
    using namespace dealii;

    /**
     * @note This benchmark only talks about the flow field, not about a
     * temperature field. All quantities related to the temperature are
     * therefore set to zero in the implementation of this class.
     *
     * @ingroup MaterialModels
     */
    template <int dim>
    class PremDensity : public MaterialModel::Interface<dim>
    {
      public:
        /**
         * @name Physical parameters used in the basic equations
         * @{
         */
        virtual void evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
                              MaterialModel::MaterialModelOutputs<dim> &out) const
        {
          for (unsigned int i=0; i < in.n_evaluation_points(); ++i)
            {
              const std::array<double,dim> position_point_sph = Utilities::Coordinates::cartesian_to_spherical_coordinates<dim>(in.position[i]);
              const double normalized_r = position_point_sph[0] / 6371e3;

              if (position_point_sph[0] > 6371e3)
                out.densities[i]=1;
              else if (position_point_sph[0] <= 1221.5e3)
                out.densities[i]=(13.0885-8.8381*std::pow(normalized_r,2))*1e3;
              else if ((position_point_sph[0] > 1221.5e3) && (position_point_sph[0] < 3480.e3))
                out.densities[i]=(12.5815-1.2638*normalized_r-3.6426*std::pow(normalized_r,2)-5.5281*std::pow(normalized_r,3))*1e3;
              else if ((position_point_sph[0] >= 3480.e3) && (position_point_sph[0] <= 5701.e3))
                out.densities[i]=(7.9565-6.4761*normalized_r+5.5283*std::pow(normalized_r,2)-3.0807*std::pow(normalized_r,3))*1e3;
              else if ((position_point_sph[0] > 5701.e3) && (position_point_sph[0] <= 5771.e3))
                out.densities[i]=(5.3197-1.4836*normalized_r)*1e3;
              else if ((position_point_sph[0] > 5771.e3) && (position_point_sph[0] <= 5971.e3))
                out.densities[i]=(11.2494-8.0298*normalized_r)*1e3;
              else if ((position_point_sph[0] > 5971.e3) && (position_point_sph[0] <= 6151.e3))
                out.densities[i]=(7.1089-3.8045*normalized_r)*1e3;
              else if ((position_point_sph[0] > 6151.e3) && (position_point_sph[0] <= 6346.e3))
                out.densities[i]=(2.6910+0.6924*normalized_r)*1e3;
              else if ((position_point_sph[0] > 6346.e3) && (position_point_sph[0] <= 6356.e3))
                out.densities[i]=2.900*1e3;
              else if ((position_point_sph[0] > 6356.e3) && (position_point_sph[0] <= 6368.e3))
                out.densities[i]=2.600*1e3;
              else
                out.densities[i]=1.020*1e3;

              out.compressibilities[i] = 0;
              out.specific_heat[i] = 0;
              out.thermal_expansion_coefficients[i] = 0;
              out.thermal_conductivities[i] = 0.0;
              out.viscosities[i] = 0;
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
        virtual bool is_compressible () const
        {
          return false;
        }
        /**
         * @}
         */

        /**
         * Read the parameters this class declares from the parameter file.
         */
        virtual
        void
        parse_parameters (ParameterHandler &)
        {
          // Declare dependencies on solution variables
          this->model_dependence.viscosity = MaterialModel::NonlinearDependence::none;
          this->model_dependence.density = MaterialModel::NonlinearDependence::none;
          this->model_dependence.compressibility = MaterialModel::NonlinearDependence::none;
          this->model_dependence.specific_heat = MaterialModel::NonlinearDependence::none;
          this->model_dependence.thermal_conductivity = MaterialModel::NonlinearDependence::none;
        }

        /**
         * @name Reference quantities
         * @{
         */
        virtual double reference_viscosity () const
        {
          return 1.e21;
        }
        /**
         * @}
         */

    };

  }
}



// explicit instantiations
namespace aspect
{
  namespace GravityPremBenchmark
  {
    ASPECT_REGISTER_MATERIAL_MODEL(PremDensity,
                                   "premDensity",
                                   "A material model that corresponds to the `gravity_prem' benchmark. "
                                   "See the manual for more information.")
  }
}

