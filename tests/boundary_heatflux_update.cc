/*
  Copyright (C) 2022 by the authors of the ASPECT code.

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

#include <aspect/boundary_heat_flux/interface.h>
#include <aspect/simulator_access.h>
#include <aspect/utilities.h>

#include <deal.II/base/parsed_function.h>
#include <aspect/global.h>
#include <deal.II/base/signaling_nan.h>

namespace aspect
{
  namespace BoundaryHeatFlux
  {
    /**
     * A class that implements heat flux boundary conditions based on a
     * functional description provided in the input file.
     *
     * @ingroup BoundaryHeatFlux
     */
    template <int dim>
    class Update : public Interface<dim>, public SimulatorAccess<dim>
    {
      public:
        /**
         * Return the boundary heat flux as a function of position.
         */
        virtual
        std::vector<Tensor<1,dim>>
        heat_flux (const types::boundary_id /*boundary_indicator*/,
                   const MaterialModel::MaterialModelInputs<dim> &material_model_inputs,
                   const MaterialModel::MaterialModelOutputs<dim> &/*material_model_outputs*/,
                   const std::vector<Tensor<1,dim>> &normal_vectors) const
        {
          std::vector<Tensor<1,dim>> heat_flux(normal_vectors);
          const unsigned int n_evaluation_points = material_model_inputs.position.size();
          for (unsigned int i=0; i<n_evaluation_points; ++i)
            heat_flux[i] *= timestep * (-1.e-4);
          return heat_flux;
        }

        /**
         * Initialization function. This function is called once at the
         * beginning of the program after parse_parameters is run and after the
         * SimulatorAccess (if applicable) is initialized.
         */
        virtual
        void
        initialize ()
        {
          std::cout << "Initializing boundary heat flux!" << std::endl;
          timestep = 0;
        };

        /**
         * A function that is called at the beginning of each time step to
         * indicate what the model time is for which the boundary values will
         * next be evaluated. For the current class, the function passes to
         * the parsed function what the current time is.
         */
        virtual
        void
        update ()
        {
          timestep = this->get_timestep_number();
          std::cout << "Updating boundary heat flux to time step "
                    << timestep
                    << '!'
                    << std::endl;
        }

      private:
        unsigned int timestep;
    };

  }
}

// explicit instantiations
namespace aspect
{
  namespace BoundaryHeatFlux
  {
    ASPECT_REGISTER_BOUNDARY_HEAT_FLUX_MODEL(Update,
                                             "heatflux update",
                                             "Heat flux boundary to test the update and initialize functions.")
  }
}
