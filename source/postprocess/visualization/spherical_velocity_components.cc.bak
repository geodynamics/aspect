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


#include <aspect/postprocess/visualization/spherical_velocity_components.h>
#include <aspect/utilities.h>


namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      template <int dim>
      SphericalVelocityComponents<dim>::
      SphericalVelocityComponents ()
        :
        DataPostprocessorVector<dim> ("spherical_velocity_components",
                                      update_values | update_quadrature_points),
        Interface<dim>()    // unknown units at construction time, will be filled by a separate function
      {}



      template <int dim>
      std::string
      SphericalVelocityComponents<dim>::
      get_physical_units () const
      {
        if (this->convert_output_to_years())
          return "m/year";
        else
          return "m/s";
      }



      template <int dim>
      void
      SphericalVelocityComponents<dim>::
      evaluate_vector_field(const DataPostprocessorInputs::Vector<dim> &input_data,
                            std::vector<Vector<double>> &computed_quantities) const
      {
        const unsigned int n_quadrature_points = input_data.evaluation_points.size();
        Assert (computed_quantities.size() == n_quadrature_points,    ExcInternalError());
        Assert (computed_quantities[0].size() == dim,    ExcInternalError());

        const double velocity_scaling_factor =
          this->convert_output_to_years() ? year_in_seconds : 1.0;

        Tensor<1,dim> velocity;
        for (unsigned int q=0; q<n_quadrature_points; ++q)
          {
            for (unsigned int d = 0; d < dim; ++d)
              {
                velocity[d] = input_data.solution_values[q][this->introspection().component_indices.velocities[d]] * velocity_scaling_factor;
              }

            std::array<double,dim> scoord = aspect::Utilities::Coordinates::cartesian_to_spherical_coordinates(input_data.evaluation_points[q]);

            if (dim==2)
              {
                computed_quantities[q](0) =  std::cos(scoord[1])*velocity[0]+std::sin(scoord[1])*velocity[1]; // v_r
                computed_quantities[q](1) = -std::sin(scoord[1])*velocity[0]+std::cos(scoord[1])*velocity[1]; // v_scoord[1]
              }

            if (dim==3)
              {
                computed_quantities[q](0) = velocity[0]*(std::sin(scoord[2])*std::cos(scoord[1]))
                                            + velocity[1]*(std::sin(scoord[2])*std::sin(scoord[1]))
                                            + velocity[2]*(std::cos(scoord[2]));                // v_r
                computed_quantities[q](1) = velocity[0]*(-std::sin(scoord[1]))
                                            + velocity[1]*(std::cos(scoord[1]));                // v_phi
                computed_quantities[q](2) = velocity[0]*(std::cos(scoord[2])*std::cos(scoord[1]))
                                            + velocity[1]*(std::cos(scoord[2])*std::sin(scoord[1]))
                                            + velocity[2]*(-std::sin(scoord[2]));               // v_theta
              }
          }
      }


      template <int dim>
      std::vector<std::string>
      SphericalVelocityComponents<dim>::get_names () const
      {
        std::vector<std::string> names;
        names.emplace_back("v_r");
        names.emplace_back("v_phi");
        if (dim==3)
          names.emplace_back("v_theta");
        return names;
      }


      template <int dim>
      std::vector<DataComponentInterpretation::DataComponentInterpretation>
      SphericalVelocityComponents<dim>::get_data_component_interpretation () const
      {
        return
          std::vector<DataComponentInterpretation::DataComponentInterpretation>
          (dim,DataComponentInterpretation::component_is_scalar);
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
      ASPECT_REGISTER_VISUALIZATION_POSTPROCESSOR(SphericalVelocityComponents,
                                                  "spherical velocity components",
                                                  "A visualization output object that outputs the polar coordinates "
                                                  "components $v_r$ and $v_\\phi$ of the velocity field in 2d and the "
                                                  "spherical coordinates components $v_r$, $v_{\\phi}$ and $v_{\\theta}$ "
                                                  "of the velocity field in 3d."
                                                  "\n\n"
                                                  "Physical units: $\\frac{\\text{m}}{\\text{s}}$ or "
                                                  "$\\frac{\\text{m}}{\\text{year}}$, depending on settings in the input file.")
    }
  }
}
