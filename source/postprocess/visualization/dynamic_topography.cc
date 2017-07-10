/*
  Copyright (C) 2011 - 2016 by the authors of the ASPECT code.

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
#include <aspect/postprocess/visualization/dynamic_topography.h>
#include <aspect/postprocess/dynamic_topography.h>

namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      template <int dim>
      std::pair<std::string, Vector<float> *>
      DynamicTopography<dim>::execute() const
      {
        Postprocess::DynamicTopography<dim> *dynamic_topography =
          this->template find_postprocessor<Postprocess::DynamicTopography<dim> >();
        AssertThrow(dynamic_topography != NULL,
                    ExcMessage("Could not find the DynamicTopography postprocessor."));
        std::pair<std::string, Vector<float> *>
        return_value ("dynamic_topography",
                      new Vector<float>(dynamic_topography->cellwise_topography()));

        return return_value;
      }

      /**
       * Register the other postprocessor that we need: DynamicTopography
       */
      template <int dim>
      std::list<std::string>
      DynamicTopography<dim>::required_other_postprocessors() const
      {
        return std::list<std::string> (1, "dynamic topography");
      }
    }
  }
}


// explicit instantiations
namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      ASPECT_REGISTER_VISUALIZATION_POSTPROCESSOR(DynamicTopography,
                                                  "dynamic topography",
                                                  "A visualization output object that generates output "
                                                  "for the dynamic topography at the top and bottom of the model space. The approach to determine the "
                                                  "dynamic topography requires us to compute the stress tensor and "
                                                  "evaluate the component of it in the direction in which "
                                                  "gravity acts. In other words, we compute "
                                                  "$\\sigma_{rr}={\\hat g}^T(2 \\eta \\varepsilon(\\mathbf u)-\\frac 13 (\\textrm{div}\\;\\mathbf u)I)\\hat g - p_d$ "
                                                  "where $\\hat g = \\mathbf g/\\|\\mathbf g\\|$ is the direction of "
                                                  "the gravity vector $\\mathbf g$ and $p_d=p-p_a$ is the dynamic "
                                                  "pressure computed by subtracting the adiabatic pressure $p_a$ "
                                                  "from the total pressure $p$ computed as part of the Stokes "
                                                  "solve. From this, the dynamic "
                                                  "topography is computed using the formula "
                                                  "$h=\\frac{\\sigma_{rr}}{(\\mathbf g \\cdot \\mathbf n)  \\rho}$ where $\\rho$ "
                                                  "is the density at the cell center. For the bottom surface we chose the convection "
                                                  "that positive values are up (out) and negative values are in (down), analogous to "
                                                  "the deformation of the upper surface. "
                                                  "Note that this implementation takes "
                                                  "the direction of gravity into account, which means that reversing the flow "
                                                  "in backward advection calculations will not reverse the instantaneous topography "
                                                  "because the reverse flow will be divided by the reverse surface gravity."
                                                  "\n\n"
                                                  "Strictly speaking, the dynamic topography is of course a "
                                                  "quantity that is only of interest at the surface. However, "
                                                  "we compute it everywhere to make things fit into the framework "
                                                  "within which we produce data for visualization. You probably "
                                                  "only want to visualize whatever data this postprocessor generates "
                                                  "at the surface of your domain and simply ignore the rest of the "
                                                  "data generated.")
    }
  }
}
